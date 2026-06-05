//! syncmer-hifi — hierarchical syncmer-LCP read mapper for PacBio HiFi.
//!
//! Pipeline (see README for the full method description):
//!   reads → open syncmers → locally-consistent parsing (LCP) into a 6-level
//!   block hierarchy → index lookup → conservation-of-mass weighted anchors →
//!   diagonal voting (or chain DP) → MAPQ → PAF (optional chain-guided CIGAR).
//!
//! Source is split into modules: config, hash, syncmer, lcp, index, fastq,
//! align, chain, map (see README "Source layout").
//!
//! Default preset: k=19, s=10 (syncmer density 1/(k-s+1) = 1/10), k-mer atom.
//!
//!   cargo run --release -- --build-index genome.fa.gz index.idx
//!   cargo run --release -- --map reads.fq.gz index.idx -o out.paf

mod config; mod hash; mod syncmer; mod lcp; mod index;
mod fastq; mod align; mod chain; mod map;

use config::*;
use syncmer::*;
use index::*;
use map::*;
use std::env;
use std::io::Write;
use std::time::Instant;
use std::sync::atomic::Ordering;

const USAGE: &str = "\
syncmer-hifi — hierarchical syncmer-LCP mapper for PacBio HiFi reads

USAGE:
  Build an index from a reference FASTA:
    syncmer-hifi --build-index <genome.fa[.gz]> <out.idx> [options]

  Map HiFi reads against a prebuilt index (or a FASTA):
    syncmer-hifi --map <reads.fq[.gz]> <index.idx|genome.fa[.gz]> -o <out.paf> [options]
    syncmer-hifi --map <reads.fq[.gz]> <index.idx> --ref <genome.fa[.gz]> --cigar -o <out.paf>

OPTIONS:
  --k N            k-mer length          (default 19)
  --s N            syncmer s-mer length  (default 10; density = 1/(k-s+1))
  --threads N      worker threads        (default: all cores)
  --cigar [BAND]   emit base-level CIGAR (cg:Z:) — needs the reference;
                   BAND = half-band in bp (default auto = max(100, len/50))
  --ref <fa>       reference FASTA for --cigar when mapping against a .idx
  --max-occ N      max genomic occurrences per L2+ block (default 500)
  --index-min-level N  (build only) lowest block level to index (default 3).
                   Matches --min-level; skipping L0-L2 shrinks the index and
                   speeds the build ~10×. Build with 0 to support mapping noisy
                   reads at --min-level 0.
  --min-level N    lowest block level allowed to anchor a read (default 3).
                   Default 3 is tuned for HiFi (≤~0.5% error): higher accuracy,
                   near-zero wrong-chromosome, ~50% faster. Use 0 for noisy
                   (>1% error) reads, which need the finer-block fallback.
  --vote           place reads by diagonal voting (heaviest anchor cluster).
                   This is the default selector — slightly more accurate.
  --chaining       place reads by colinear chain-DP instead of voting.
  --rescue         second relaxed-filter pass for reads that fail the first;
                   maps a few % more reads in repeat-rich regions, emitted at
                   MAPQ 0 (flagged uncertain).  Off by default.

Recommended HiFi preset (the default): --k 19 --s 10
";

fn main() {
    let args: Vec<String> = env::args().skip(1).collect();
    if args.is_empty() || args.iter().any(|a| a == "-h" || a == "--help") {
        print!("{USAGE}");
        return;
    }

    // --threads N  (must be parsed before rayon is first used)
    if let Some(pos) = args.iter().position(|a| a == "--threads") {
        let n: usize = args.get(pos + 1)
            .and_then(|s| s.parse().ok())
            .expect("--threads requires a number");
        rayon::ThreadPoolBuilder::new().num_threads(n).build_global().unwrap();
    }

    // HiFi always uses open-syncmer seeds with the full k-mer as the LCP atom
    // (the k-mer atom gives minimap2-level block specificity at HiFi error rates).
    let mode = SeedMode::Syncmer;
    KMER_ATOM.store(true, Ordering::Relaxed);
    if args.iter().any(|a| a == "--rescue") {
        RESCUE.store(true, Ordering::Relaxed);
    }
    // Locus selector: --vote (default, diagonal voting) vs --chaining (chain-DP).
    if args.iter().any(|a| a == "--chaining") { set_chaining(true); }
    else if args.iter().any(|a| a == "--vote") { set_chaining(false); }

    // --min-level N  (default 3): minimum hierarchy level allowed to anchor a read.
    // Blocks below N are ignored, placing reads via the coarser, more-unique
    // high-level blocks. Default 3 is tuned for HiFi (≤~0.5% error): higher
    // accuracy, near-zero wrong-chromosome, ~50% faster. Use --min-level 0 for
    // noisy (>1%) reads, where the finer-block fallback is needed for sensitivity.
    if let Some(v) = args.iter().position(|a| a == "--min-level")
        .and_then(|p| args.get(p + 1)).and_then(|v| v.parse::<usize>().ok())
    {
        set_min_lvl(v);
    }

    let k: usize = args.iter().position(|a| a == "--k")
        .and_then(|p| args.get(p + 1)).and_then(|v| v.parse().ok()).unwrap_or(19);
    let s: usize = args.iter().position(|a| a == "--s")
        .and_then(|p| args.get(p + 1)).and_then(|v| v.parse().ok()).unwrap_or(10);
    assert!(k > s, "--k ({k}) must be greater than --s ({s})");
    let t: usize = (k - s) / 2;  // middle s-mer position

    let max_occ: usize = args.iter().position(|a| a == "--max-occ")
        .and_then(|p| args.get(p + 1)).and_then(|s| s.parse().ok())
        .unwrap_or(MAX_OCC_DEFAULT);
    let max_occ_l1: usize = (max_occ * 2 / 5).max(10);

    // --index-min-level N (default 3): lowest block level to actually index.
    // Default matches the mapping default (--min-level 3): the bulk L0-L2 blocks
    // are never used at HiFi error rates, so skipping them shrinks the index and
    // speeds the build ~10×. Build with 0 if you also map noisy reads at
    // --min-level 0 (those passes need the finer levels present in the index).
    let idx_min_level: usize = args.iter().position(|a| a == "--index-min-level")
        .and_then(|p| args.get(p + 1)).and_then(|v| v.parse().ok())
        .unwrap_or(3);

    // --build-index <genome.fa[.gz]> <output.idx>
    if let Some(bi_pos) = args.iter().position(|a| a == "--build-index") {
        let genome_path = args.get(bi_pos + 1)
            .expect("--build-index requires a genome path");
        let idx_path = args[bi_pos + 2..].iter()
            .find(|a| !a.starts_with('-'))
            .map(|s| s.as_str()).unwrap_or("syncmer-hifi.idx");
        let nthreads = rayon::current_num_threads();
        println!("\n  Building L{MAX_LEVELS} index (levels ≥ {idx_min_level})  k={k}  s={s}  threads={nthreads}");
        println!("  genome : {genome_path}");
        println!("  output : {idx_path}\n");
        let (index, elapsed) = build_index(genome_path, k, s, t, MAX_LEVELS, idx_min_level, mode);
        println!("  Build time : {:.2}s", elapsed.as_secs_f64());
        print!("  Saving to {idx_path} ... ");
        std::io::stdout().flush().ok();
        let t0 = Instant::now();
        save_index(&index, idx_path, k, s, t, mode);
        println!("done  ({:.1}s)", t0.elapsed().as_secs_f64());
        println!("  Index size : {:.0} MB",
            std::fs::metadata(idx_path).map(|m| m.len()).unwrap_or(0) as f64 / 1e6);
        return;
    }

    // --cigar [BAND]  + --ref <genome.fa>
    let cigar_band: Option<usize> = if args.iter().any(|a| a == "--cigar") {
        let explicit = args.iter().position(|a| a == "--cigar")
            .and_then(|p| args.get(p + 1))
            .and_then(|v| v.parse::<usize>().ok());
        Some(explicit.unwrap_or(0))  // 0 = auto (max(100, qlen/50) per read)
    } else {
        None
    };
    let ref_override = args.iter().position(|a| a == "--ref")
        .and_then(|p| args.get(p + 1)).map(|s| s.as_str());

    // --map <reads.fastq[.gz]> <genome.fa[.gz]|index.idx> [-o out.paf]
    if let Some(map_pos) = args.iter().position(|a| a == "--map") {
        let reads_path  = args.get(map_pos + 1)
            .expect("--map requires a reads path argument");
        let genome_path = args.get(map_pos + 2)
            .expect("--map requires a genome/index path argument");
        let paf_out = args.iter().position(|a| a == "-o")
            .and_then(|p| args.get(p + 1)).map(|s| s.as_str());
        run_mapping(reads_path, genome_path, paf_out,
                    k, s, t, max_occ, max_occ_l1, mode, cigar_band, ref_override);
        return;
    }

    print!("{USAGE}");
}
