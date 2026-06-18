use crate::config::*;
use crate::syncmer::*;
use crate::index::*;
use crate::fastq::*;
use crate::chain::*;
use std::io::{BufWriter, Write};
use std::fs::File;
use std::time::Instant;
use std::sync::atomic::Ordering;
use rayon::prelude::*;

pub(crate) struct MapResult {
    pub(crate) chr:    u8,
    pub(crate) pos:    i64,   // inferred reference start (may be negative near contig ends)
    pub(crate) strand: bool,  // true = forward
    pub(crate) votes:  u32,
    pub(crate) mapq:   u8,    // 0 = ambiguous, 60 = uniquely placed
}

/// Inner mapping pass: collect anchors on both strands and vote.
///
/// MAPQ correctness note: a read from (say) chr1:P forward strand also produces
/// a near-identical cluster on the RC strand at the same genomic locus. Counting
/// that as a competing "second-best" would collapse MAPQ to ~0 for every
/// uniquely-mapped read, so we only promote a locus to second_score when it is
/// genuinely different (different chromosome, or offset > CLUSTER_WINDOW).
pub(crate) fn map_read_with_occ(fwd: &[u8], rc: &[u8], index: &GIndex,
                     k: usize, s: usize, t: usize, mode: SeedMode,
                     max_occ: usize, max_occ_l1: usize)
    -> (Option<MapResult>, u32)
{
    let mut best: Option<MapResult> = None;
    let mut second_score: u32 = 0;

    for (seq, strand) in [(&fwd[..], true), (&rc[..], false)] {
        let mut anchors = collect_anchors(seq, index, k, s, t, mode,
                                          max_occ, max_occ_l1);
        let (top, second_here) = timed(&PROF_SELECT, || vote_locus(&mut anchors));
        second_score = second_score.max(second_here);

        if let Some((chr, pos, score)) = top {
            let candidate = MapResult { chr, pos, strand, votes: score, mapq: 0 };
            match &best {
                None => best = Some(candidate),
                Some(b) if score > b.votes => {
                    let same_locus = chr == b.chr && (pos - b.pos).abs() < CLUSTER_WINDOW;
                    if !same_locus { second_score = second_score.max(b.votes); }
                    best = Some(candidate);
                }
                Some(b) => {
                    let same_locus = chr == b.chr && (pos - b.pos).abs() < CLUSTER_WINDOW;
                    if !same_locus { second_score = second_score.max(score); }
                }
            }
        }
    }
    (best, second_score)
}

/// Map a single read by diagonal voting on variable-span anchors.
/// MAPQ = (best − second) × 60 / best — the fraction of the best score that is
/// uncontested. Naturally low when two genomic loci score equally well.
pub(crate) fn map_read(fwd: &[u8], index: &GIndex, k: usize, s: usize, t: usize, mode: SeedMode,
            max_occ: usize, max_occ_l1: usize) -> Option<MapResult> {
    let rc = timed(&PROF_RC, || revcomp(fwd));
    let (mut best, second_score) =
        map_read_with_occ(fwd, &rc, index, k, s, t, mode, max_occ, max_occ_l1);
    if let Some(ref mut b) = best {
        if b.votes > second_score {
            b.mapq = (b.votes.saturating_sub(second_score)
                .saturating_mul(60) / b.votes.max(1)).min(60) as u8;
            return best;
        }
    }
    None
}

#[derive(Default)]
pub(crate) struct MapStats {
    pub(crate) total:    u64,
    pub(crate) mapped:   u64,
    pub(crate) bases:    u64,
}

pub(crate) fn write_paf_line(
    paf:       &mut Option<BufWriter<File>>,
    name:      &str,
    len:       u64,
    result:    &Option<MapResult>,
    chr_names: &[String],
) {
    let Some(w) = paf else { return };
    if let Some(mr) = result {
        let chr_name  = chr_names.get(mr.chr as usize).map(|s| s.as_str()).unwrap_or("*");
        let strand    = if mr.strand { '+' } else { '-' };
        let ref_start = mr.pos.max(0) as u64;
        let ref_end   = ref_start + len;
        writeln!(w, "{name}\t{len}\t0\t{len}\t{strand}\t{chr_name}\t0\t{ref_start}\t{ref_end}\t{len}\t{len}\t{}", mr.mapq).unwrap();
    } else {
        writeln!(w, "{name}\t{len}\t0\t{len}\t*\t*\t0\t0\t0\t0\t0\t0").unwrap();
    }
}

pub(crate) fn run_mapping(reads_path: &str, genome_or_idx: &str, paf_out: Option<&str>,
               k: usize, s: usize, t: usize, max_occ: usize, max_occ_l1: usize,
               mode: SeedMode)
{
    // ── Load or build index ───────────────────────────────────────────────────
    // When loading a prebuilt .idx file the k/s/t/mode stored inside are
    // authoritative; CLI values are ignored so the user doesn't have to
    // repeat them on the command line.
    let (index, k, s, t, mode, idx_elapsed) = if genome_or_idx.ends_with(".idx") {
        print!("  Loading index from {} ... ", genome_or_idx);
        std::io::stdout().flush().ok();
        let t0 = Instant::now();
        let (idx, ik, is, it, imode) = load_index(genome_or_idx);
        let elapsed = t0.elapsed();
        println!("done  ({:.2}s)", elapsed.as_secs_f64());
        (idx, ik, is, it, imode, elapsed)
    } else {
        // In-memory index for direct FASTA mapping: build to the same level floor
        // we will query at (min_lvl(), default 3).
        println!("  Building index from {} (levels ≥ {}) ...", genome_or_idx, min_lvl());
        let (idx, elapsed) = build_index(genome_or_idx, k, s, t, MAX_LEVELS, mode);
        (idx, k, s, t, mode, elapsed)
    };
    for (li, lv) in index.levels.iter().enumerate() {
        println!("  L{} anchors indexed : {:>16}", li, commas(lv.len() as u64));
    }
    println!("  Chromosomes        : {:>16}", index.chr_names.len());
    println!("  Index build/load   : {:>13.2}s", idx_elapsed.as_secs_f64());
    println!();

    // ── Optional PAF output writer ────────────────────────────────────────────
    let mut paf_writer: Option<BufWriter<File>> = paf_out.map(|p| {
        BufWriter::with_capacity(1 << 20,
            File::create(p).unwrap_or_else(|e| panic!("Cannot create {p}: {e}")))
    });

    // ── Map reads (parallel — batched rayon) ─────────────────────────────────
    println!("  Mapping reads from {} ...", reads_path);
    if paf_out.is_some() { println!("  PAF output        → {}", paf_out.unwrap()); }
    let t0 = Instant::now();
    let mut stats = MapStats::default();

    // Batch by sequence-byte budget (--batch-mb) so peak memory is independent of
    // read length; 0 = unbounded (single batch). BATCH_MAX caps tiny-read counts.
    let budget = batch_bytes_budget();
    const BATCH_MAX: usize = 200_000;      // read-count safety cap
    let mut batch: Vec<(String, Vec<u8>)> = Vec::new();
    let mut batch_bytes: usize = 0;

    // flush_batch: map all reads in parallel, then write PAF sequentially.
    let flush_batch = |batch: &Vec<(String, Vec<u8>)>,
                           stats: &mut MapStats,
                           paf: &mut Option<BufWriter<File>>| {
        let results: Vec<(String, u64, Option<MapResult>)> =
            batch.par_iter()
            .map(|(name, seq)| {
                let r = map_read(seq, &index, k, s, t, mode, max_occ, max_occ_l1);
                (name.clone(), seq.len() as u64, r)
            })
            .collect();

        for (name, len, result) in results {
            stats.total += 1;
            stats.bases += len;
            if result.is_some() { stats.mapped += 1; }
            write_paf_line(paf, &name, len, &result, &index.chr_names);
        }

        if stats.total % 1_000_000 == 0 {
            print!("\r  mapped {}/{} ({:.1}%)   ",
                stats.mapped, stats.total,
                100.0 * stats.mapped as f64 / stats.total as f64);
            std::io::stdout().flush().ok();
        }
    };

    for read in FastqReader::open(reads_path) {
        batch_bytes += read.1.len();
        batch.push(read);
        if batch_bytes >= budget || batch.len() >= BATCH_MAX {
            flush_batch(&batch, &mut stats, &mut paf_writer);
            batch.clear();
            batch_bytes = 0;
        }
    }
    if !batch.is_empty() {
        flush_batch(&batch, &mut stats, &mut paf_writer);
    }
    if let Some(ref mut w) = paf_writer { w.flush().unwrap(); }

    let map_elapsed = t0.elapsed();

    println!("\r                                                  ");
    let sep = "═".repeat(64);
    println!("\n{sep}");
    println!("  Syncmer LCP — Read Mapping Results");
    let _ = mode;
    println!("  k={k}  s={s}  t={t}  syncmer  max_occ={max_occ}  max_occ_l1={max_occ_l1}");
    println!("{sep}");
    println!("  Total reads        : {:>16}", commas(stats.total));
    println!("  Mapped             : {:>16}  ({:.2}%)",
        commas(stats.mapped),
        100.0 * stats.mapped as f64 / stats.total as f64);
    println!("  Unmapped           : {:>16}  ({:.2}%)",
        commas(stats.total - stats.mapped),
        100.0 * (stats.total - stats.mapped) as f64 / stats.total as f64);
    println!("  Total bases        : {:>16}", commas(stats.bases));
    println!();
    println!("  Index build time   : {:>13.2}s", idx_elapsed.as_secs_f64());
    println!("  Mapping time       : {:>13.2}s", map_elapsed.as_secs_f64());
    let rps = stats.total as f64 / map_elapsed.as_secs_f64();
    let bps = stats.bases as f64 / map_elapsed.as_secs_f64();
    println!("  Throughput         : {:>10.0} reads/s  ({:.2} Mbp/s)",
        rps, bps / 1e6);
    println!();

    // ── Per-stage profile (PROFILE=1) ────────────────────────────────────────
    if profiling() {
        let rc  = PROF_RC.load(Ordering::Relaxed);
        let hi  = PROF_HIER.load(Ordering::Relaxed);
        let em  = PROF_EMIT.load(Ordering::Relaxed);
        let pr  = PROF_PRUNE.load(Ordering::Relaxed);
        let se  = PROF_SELECT.load(Ordering::Relaxed);
        let sum = (rc + hi + em + pr + se).max(1);
        let nthreads = rayon::current_num_threads() as f64;
        let ms = |ns: u64| ns as f64 / 1e6;
        let pct = |ns: u64| 100.0 * ns as f64 / sum as f64;
        let sd = PROF_SEED.load(Ordering::Relaxed);
        let l1 = PROF_L1.load(Ordering::Relaxed);
        let up = PROF_UPPER.load(Ordering::Relaxed);
        println!("  Per-stage CPU time (summed over {} threads):", nthreads);
        println!("    revcomp          : {:>10.0} ms  ({:>5.1}%)", ms(rc), pct(rc));
        println!("    seeding+LCP+hier : {:>10.0} ms  ({:>5.1}%)", ms(hi), pct(hi));
        println!("      ├ syncmer seed : {:>10.0} ms  ({:>5.1}%)", ms(sd), pct(sd));
        println!("      ├ L1 parse     : {:>10.0} ms  ({:>5.1}%)", ms(l1), pct(l1));
        println!("      └ L2..L6 recur : {:>10.0} ms  ({:>5.1}%)", ms(up), pct(up));
        println!("    index lookups    : {:>10.0} ms  ({:>5.1}%)", ms(em), pct(em));
        println!("    anchor prune     : {:>10.0} ms  ({:>5.1}%)", ms(pr), pct(pr));
        println!("    select (vote/dp) : {:>10.0} ms  ({:>5.1}%)", ms(se), pct(se));
        println!("    ─ profiled total : {:>10.0} ms  (≈ {:.2}s wall ÷ {:.0} thr)",
                 ms(sum), ms(sum)/1000.0/nthreads, nthreads);
        println!();
    }
}


// ─────────────────────────────────────────────────────────────────────────────
// 15.  Main
// ─────────────────────────────────────────────────────────────────────────────
