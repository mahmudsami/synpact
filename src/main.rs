//! synpact — hierarchical syncmer-LCP read mapper for PacBio HiFi.
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
mod fastq; mod chain; mod map;

use config::*;
use syncmer::*;
use index::*;
use map::*;
use std::env;
use std::io::Write;
use std::time::Instant;

const USAGE: &str = "\
synpact — hierarchical syncmer-LCP mapper for PacBio HiFi reads

USAGE:
  Build an index from a reference FASTA:
    synpact --build-index <genome.fa[.gz]> <out.idx> [options]

  Map HiFi reads against a prebuilt index (or a FASTA):
    synpact --map <reads.fq[.gz]> <index.idx|genome.fa[.gz]> -o <out.paf> [options]

OPTIONS:
  --k N            k-mer length          (default 19, must be ≤ 32)
  --s N            syncmer s-mer length  (default 10; density = 1/(k-s+1))
  --threads N      worker threads        (default: all cores)
  --batch-mb N     per-batch read-sequence budget in MiB (default 256); keeps
                   peak memory independent of read length. 0 = unbounded.
  --min-level N    lowest block level used to anchor a read AND the lowest level
                   stored in the index (default 3; one unified value). Default 3
                   is tuned for HiFi: high accuracy, near-zero wrong-chromosome.
  --max-occ N      max genomic occurrences per L2+ block before it is skipped as
                   too-repetitive (default 10000). Higher surfaces more anchors in
                   segmental-dup regions (better accuracy); very large values can
                   slow highly repetitive genomes.

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

    // HiFi always uses open-syncmer seeds with the full k-mer as the LCP atom.
    let mode = SeedMode::Syncmer;

    // --min-level N (default 3): the single floor used for both index build and
    // mapping. Blocks below it are neither stored nor used.
    if let Some(v) = args.iter().position(|a| a == "--min-level")
        .and_then(|p| args.get(p + 1)).and_then(|v| v.parse::<usize>().ok())
    {
        set_min_lvl(v);
    }

    // --batch-mb N (default 256): per-batch read-sequence budget in MiB.
    if let Some(v) = args.iter().position(|a| a == "--batch-mb")
        .and_then(|p| args.get(p + 1)).and_then(|v| v.parse::<usize>().ok())
    {
        set_batch_mb(v);
    }

    let k: usize = args.iter().position(|a| a == "--k")
        .and_then(|p| args.get(p + 1)).and_then(|v| v.parse().ok()).unwrap_or(19);
    let s: usize = args.iter().position(|a| a == "--s")
        .and_then(|p| args.get(p + 1)).and_then(|v| v.parse().ok()).unwrap_or(10);
    assert!(k > s, "--k ({k}) must be greater than --s ({s})");
    assert!(k <= 32, "--k ({k}) must be ≤ 32 (k-mer atom must fit in u64)");
    let t: usize = (k - s) / 2;  // middle s-mer position

    // --max-occ N (default 10000): per-L2+-block occurrence cap before a block is
    // skipped as too-repetitive. The L1 cap scales with it (2/5, min 10).
    let max_occ: usize = args.iter().position(|a| a == "--max-occ")
        .and_then(|p| args.get(p + 1)).and_then(|v| v.parse().ok())
        .unwrap_or(MAX_OCC_DEFAULT);
    let max_occ_l1: usize = (max_occ * 2 / 5).max(10);

    // --build-index <genome.fa[.gz]> <output.idx>
    if let Some(bi_pos) = args.iter().position(|a| a == "--build-index") {
        let genome_path = args.get(bi_pos + 1)
            .expect("--build-index requires a genome path");
        let idx_path = args[bi_pos + 2..].iter()
            .find(|a| !a.starts_with('-'))
            .map(|s| s.as_str()).unwrap_or("synpact.idx");
        let nthreads = rayon::current_num_threads();
        println!("\n  Building L{MAX_LEVELS} index (levels ≥ {})  k={k}  s={s}  threads={nthreads}", min_lvl());
        println!("  genome : {genome_path}");
        println!("  output : {idx_path}\n");
        let (index, elapsed) = build_index(genome_path, k, s, t, MAX_LEVELS, mode);
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

    // --map <reads.fastq[.gz]> <genome.fa[.gz]|index.idx> [-o out.paf]
    if let Some(map_pos) = args.iter().position(|a| a == "--map") {
        let reads_path  = args.get(map_pos + 1)
            .expect("--map requires a reads path argument");
        let genome_path = args.get(map_pos + 2)
            .expect("--map requires a genome/index path argument");
        let paf_out = args.iter().position(|a| a == "-o")
            .and_then(|p| args.get(p + 1)).map(|s| s.as_str());
        run_mapping(reads_path, genome_path, paf_out,
                    k, s, t, max_occ, max_occ_l1, mode);
        return;
    }

    print!("{USAGE}");
}
