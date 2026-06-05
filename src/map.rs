#![allow(unused_imports, dead_code)]
use crate::config::*;
use crate::hash::*;
use crate::syncmer::*;
use crate::lcp::*;
use crate::index::*;
use crate::fastq::*;
use crate::align::*;
use crate::chain::*;
use std::io::{self, BufRead, BufReader, BufWriter, Read, Write};
use std::fs::File;
use std::time::Instant;
use std::sync::atomic::{AtomicBool, Ordering};
use std::sync::OnceLock;
use std::cell::RefCell;
use std::collections::HashSet;
use flate2::read::MultiGzDecoder;
use rayon::prelude::*;

pub(crate) struct MapResult {
    pub(crate) chr:    u8,
    pub(crate) pos:    i64,   // inferred reference start (may be negative near contig ends)
    pub(crate) strand: bool,  // true = forward
    pub(crate) votes:  u32,
    pub(crate) mapq:   u8,    // 0 = ambiguous, 60 = uniquely placed
}

// ── Anchor collection + chain DP ─────────────────────────────────────────────

/// Inner mapping pass for one set of occ thresholds.
///
/// MAPQ correctness note: a read from (say) chr1:P forward strand will also
/// produce a near-identical chain on the RC strand at the same genomic locus.
/// Counting that as a competing "second-best" chain would collapse MAPQ to ~0
/// for every uniquely-mapped read.  We avoid this by only promoting a chain
/// to second_score when it maps to a *genuinely different* locus
/// (different chromosome, or offset differing by more than CLUSTER_WINDOW).
pub(crate) fn map_read_with_occ(fwd: &[u8], rc: &[u8], index: &GIndex,
                     k: usize, s: usize, t: usize, mode: SeedMode,
                     max_occ: usize, max_occ_l1: usize)
    -> (Option<MapResult>, u32)
{
    let mut best: Option<MapResult> = None;
    let mut second_score: u32 = 0;

    for (seq, strand) in [(&fwd[..], true), (&rc[..], false)] {
        let mut anchors = collect_anchors(seq, index, k, s, t, mode,
                                          max_occ, max_occ_l1,
                                          None, None, None);
        let (top, second_here) = if use_chaining() {
            chain_dp(&mut anchors)
        } else {
            vote_locus(&mut anchors)
        };
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

/// Map a single read using chain DP on variable-span anchors.
/// MAPQ = (best − second) × 60 / best — the fraction of the best chain score
/// that is uncontested.  Naturally low when two genomic loci chain equally well.
pub(crate) fn map_read(fwd: &[u8], index: &GIndex, k: usize, s: usize, t: usize, mode: SeedMode,
            max_occ: usize, max_occ_l1: usize) -> Option<MapResult> {
    let rc = revcomp(fwd);

    let try_pass = |mo: usize, mo1: usize| -> Option<MapResult> {
        let (mut best, second_score) =
            map_read_with_occ(fwd, &rc, index, k, s, t, mode, mo, mo1);
        if let Some(ref mut b) = best {
            if b.votes > second_score {
                b.mapq = (b.votes.saturating_sub(second_score)
                    .saturating_mul(60) / b.votes.max(1)).min(60) as u8;
                return best;
            }
        }
        None
    };

    // Pass 1: default thresholds — the fast path for the cleanly-mapped majority.
    if let Some(r) = try_pass(max_occ, max_occ_l1) { return Some(r); }

    // Pass 2 (--rescue only): reads that fail pass 1 are often in repeat-rich
    // regions where ~half their true-locus L1/L2 blocks were filtered as too
    // frequent.  Retry with relaxed thresholds so those blocks survive — gated
    // strictly to first-pass failures (never disturbs a mapped read).  Recovered
    // reads come out at MAPQ 0 (flagged uncertain).
    if rescue_pass() {
        if let Some(r) = try_pass(max_occ * 4, max_occ_l1 * 4) { return Some(r); }
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
    paf:           &mut Option<BufWriter<File>>,
    name:          &str,
    len:           u64,
    result:        &Option<MapResult>,
    chr_names:     &[String],
    cigar:         Option<&str>,
    actual_ref_end: Option<u64>,
) {
    let Some(w) = paf else { return };
    if let Some(mr) = result {
        let chr_name  = chr_names.get(mr.chr as usize).map(|s| s.as_str()).unwrap_or("*");
        let strand    = if mr.strand { '+' } else { '-' };
        let ref_start = mr.pos.max(0) as u64;
        let ref_end   = actual_ref_end.unwrap_or(ref_start + len);
        let aln_len   = ref_end.saturating_sub(ref_start).max(len);
        // Count = bases from CIGAR; fall back to len for the matches field.
        let matches = cigar.map(|cg| {
            let mut m = 0u64; let mut n = 0u64;
            for b in cg.bytes() {
                if b.is_ascii_digit() { n = n * 10 + (b - b'0') as u64; }
                else { if b == b'=' { m += n; } n = 0; }
            }
            m
        }).unwrap_or(len);
        if let Some(cg) = cigar {
            writeln!(w, "{name}\t{len}\t0\t{len}\t{strand}\t{chr_name}\t0\t{ref_start}\t{ref_end}\t{matches}\t{aln_len}\t{}\tcg:Z:{cg}", mr.mapq).unwrap();
        } else {
            writeln!(w, "{name}\t{len}\t0\t{len}\t{strand}\t{chr_name}\t0\t{ref_start}\t{ref_end}\t{matches}\t{aln_len}\t{}", mr.mapq).unwrap();
        }
    } else {
        writeln!(w, "{name}\t{len}\t0\t{len}\t*\t*\t0\t0\t0\t0\t0\t0").unwrap();
    }
}

pub(crate) fn run_mapping(reads_path: &str, genome_or_idx: &str, paf_out: Option<&str>,
               k: usize, s: usize, t: usize, max_occ: usize, max_occ_l1: usize,
               mode: SeedMode, cigar_band: Option<usize>, ref_path_override: Option<&str>)
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
        let (idx, elapsed) = build_index(genome_or_idx, k, s, t, MAX_LEVELS, min_lvl(), mode);
        (idx, k, s, t, mode, elapsed)
    };
    for (li, lv) in index.levels.iter().enumerate() {
        println!("  L{} anchors indexed : {:>16}", li, commas(lv.len() as u64));
    }
    println!("  Chromosomes        : {:>16}", index.chr_names.len());
    println!("  Index build/load   : {:>13.2}s", idx_elapsed.as_secs_f64());
    println!();

    // ── Optional reference sequences for CIGAR alignment ─────────────────────
    // Reference sequences are loaded only when CIGAR output is requested.
    let ref_seqs: Option<Vec<Vec<u8>>> = if cigar_band.is_some() {
        let genome_path = ref_path_override.unwrap_or_else(|| {
            if genome_or_idx.ends_with(".idx") { "" } else { genome_or_idx }
        });
        if genome_path.is_empty() {
            eprintln!("  [warn] --cigar with .idx requires --ref <genome.fa>; CIGAR disabled.");
            None
        } else {
            print!("  Loading genome for alignment ... ");
            std::io::stdout().flush().ok();
            let chrs = read_fasta_chrs(genome_path);
            println!("{} chromosomes loaded", chrs.len());
            Some(chrs.into_iter().map(|(_, seq)| seq).collect())
        }
    } else {
        None
    };

    // ── Optional PAF output writer ────────────────────────────────────────────
    let mut paf_writer: Option<BufWriter<File>> = paf_out.map(|p| {
        BufWriter::with_capacity(1 << 20,
            File::create(p).unwrap_or_else(|e| panic!("Cannot create {p}: {e}")))
    });

    // ── Map reads (parallel — batched rayon) ─────────────────────────────────
    println!("  Mapping reads from {} ...", reads_path);
    if paf_out.is_some() { println!("  PAF output        → {}", paf_out.unwrap()); }
    if cigar_band.is_some() && ref_seqs.is_some() {
        println!("  CIGAR alignment   : enabled  (half-band = {})", cigar_band.unwrap());
    }
    let t0 = Instant::now();
    let mut stats = MapStats::default();

    const BATCH: usize = 100_000;
    let mut batch: Vec<(String, Vec<u8>)> = Vec::with_capacity(BATCH);

    // flush_batch: map all reads in parallel, collect results + optional CIGAR
    let flush_batch = |batch: &Vec<(String, Vec<u8>)>,
                           stats: &mut MapStats,
                           paf: &mut Option<BufWriter<File>>| {
        // Each result carries: name, read_len, mapping, Option<(cigar, actual_ref_end)>
        let results: Vec<(String, u64, Option<MapResult>, Option<(String, u64)>)> =
            batch.par_iter()
            .map(|(name, seq)| {
                let r = map_read(seq, &index, k, s, t, mode, max_occ, max_occ_l1);
                let cg: Option<(String, u64)> = match (&r, &ref_seqs, cigar_band) {
                    (Some(mr), Some(refs), Some(band)) => {
                        let hb = band.max(50).max(seq.len() / 50);
                        Some(cigar_for_mapping(
                            seq, mr, refs,
                            &index, k, s, t, mode, max_occ, max_occ_l1,
                            hb,
                        ))
                    }
                    _ => None,
                };
                (name.clone(), seq.len() as u64, r, cg)
            })
            .collect();

        for (name, len, result, cigar) in results {
            stats.total += 1;
            stats.bases += len;
            if result.is_some() { stats.mapped += 1; }
            let (cg_str, actual_ref_end) = match &cigar {
                Some((s, e)) => (Some(s.as_str()), Some(*e)),
                None         => (None, None),
            };
            write_paf_line(paf, &name, len, &result, &index.chr_names,
                           cg_str, actual_ref_end);
        }

        if stats.total % 1_000_000 == 0 {
            print!("\r  mapped {}/{} ({:.1}%)   ",
                stats.mapped, stats.total,
                100.0 * stats.mapped as f64 / stats.total as f64);
            std::io::stdout().flush().ok();
        }
    };

    for read in FastqReader::open(reads_path) {
        batch.push(read);
        if batch.len() == BATCH {
            flush_batch(&batch, &mut stats, &mut paf_writer);
            batch.clear();
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
}


// ─────────────────────────────────────────────────────────────────────────────
// 15.  Main
// ─────────────────────────────────────────────────────────────────────────────
