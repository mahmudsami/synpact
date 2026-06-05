#![allow(unused_imports, dead_code)]
use crate::hash::*;
use crate::syncmer::*;
use crate::lcp::*;
use crate::index::*;
use crate::fastq::*;
use crate::align::*;
use crate::chain::*;
use crate::map::*;
use std::io::{self, BufRead, BufReader, BufWriter, Read, Write};
use std::fs::File;
use std::time::Instant;
use std::sync::atomic::{AtomicBool, Ordering};
use std::sync::OnceLock;
use std::cell::RefCell;
use std::collections::HashSet;
use flate2::read::MultiGzDecoder;
use rayon::prelude::*;

// ── Per-stage profiling (env PROFILE=1) ─────────────────────────────────────
// Each counter accumulates wall-nanoseconds summed across all worker threads.
// Zero overhead when PROFILE is unset (the Instant calls are skipped).
use std::sync::atomic::AtomicU64;
pub(crate) static PROF_RC:     AtomicU64 = AtomicU64::new(0); // reverse-complement
pub(crate) static PROF_HIER:   AtomicU64 = AtomicU64::new(0); // syncmers+LCP+hierarchy (total)
pub(crate) static PROF_SEED:   AtomicU64 = AtomicU64::new(0); //   ├ syncmer selection+hashing
pub(crate) static PROF_L1:     AtomicU64 = AtomicU64::new(0); //   ├ L1 parse + block build
pub(crate) static PROF_UPPER:  AtomicU64 = AtomicU64::new(0); //   └ L2..L6 recursion
pub(crate) static PROF_EMIT:   AtomicU64 = AtomicU64::new(0); // index lookups / anchors
pub(crate) static PROF_PRUNE:  AtomicU64 = AtomicU64::new(0); // anchor pruning
pub(crate) static PROF_SELECT: AtomicU64 = AtomicU64::new(0); // voting / chain-DP

#[inline]
pub(crate) fn profiling() -> bool {
    static V: OnceLock<bool> = OnceLock::new();
    *V.get_or_init(|| std::env::var("PROFILE").map(|v| v == "1").unwrap_or(false))
}
/// Run `f`, and if profiling is on add its elapsed time to `ctr`. Returns f's value.
#[inline]
pub(crate) fn timed<T>(ctr: &AtomicU64, f: impl FnOnce() -> T) -> T {
    if profiling() {
        let t = Instant::now();
        let r = f();
        ctr.fetch_add(t.elapsed().as_nanos() as u64, Ordering::Relaxed);
        r
    } else {
        f()
    }
}

/// Hierarchy-atom selector.  When true, LCP uses the full k-mer as each leaf
/// value (high specificity, best for low-error HiFi).  When false (default),
/// it uses the middle s-mer (more error-tolerant, best for ONT).
/// Set once at startup from the CLI / loaded index, before any parallel work.
pub(crate) static KMER_ATOM: AtomicBool = AtomicBool::new(false);
#[inline] pub(crate) fn kmer_atom() -> bool { KMER_ATOM.load(Ordering::Relaxed) }

/// Canonical-atom selector.  When true, each atom value is min(encode(atom),
/// encode(revcomp(atom))) so a sequence and its reverse-complement share a key.
/// Halves the atom value space (less specificity) — tested as an alternative.
pub(crate) static CANON_ATOM: AtomicBool = AtomicBool::new(false);
#[inline] pub(crate) fn canon_atom() -> bool { CANON_ATOM.load(Ordering::Relaxed) }

/// Disable the L0 raw-syncmer fallback at query time (env NO_L0=1).
pub(crate) fn no_l0() -> bool {
    use std::sync::OnceLock; static V: OnceLock<bool> = OnceLock::new();
    *V.get_or_init(|| std::env::var("NO_L0").map(|v| v == "1").unwrap_or(false))
}

/// Minimum hierarchy level whose blocks may produce anchors (`--min-level N`,
/// default 3). Blocks below this level are neither emitted nor recursed into — the
/// read is placed only by the coarser, more-unique high-level blocks, which at HiFi
/// error rates removes paralog/segmental-dup ambiguity (higher accuracy + faster).
/// A floor > 0 also disables the L0 raw-syncmer fallback.
/// Set once from the CLI before mapping; falls back to env MIN_LVL then default 3.
pub(crate) static MIN_LVL_CELL: std::sync::OnceLock<usize> = std::sync::OnceLock::new();

pub(crate) fn set_min_lvl(v: usize) { let _ = MIN_LVL_CELL.set(v); }

pub(crate) fn min_lvl() -> usize {
    *MIN_LVL_CELL.get_or_init(|| {
        std::env::var("MIN_LVL").ok().and_then(|v| v.parse().ok()).unwrap_or(3)
    })
}

/// Locus selector. Default is diagonal **voting** (`--vote`): pick the heaviest
/// anchor cluster on a diagonal. `--chaining` (env CHAINING=1) instead uses the
/// colinear gap-penalised chain-DP. Voting is slightly more accurate on contested
/// reads; both run at the same speed (chaining is not the bottleneck).
pub(crate) static USE_CHAINING_CELL: std::sync::OnceLock<bool> = std::sync::OnceLock::new();

pub(crate) fn set_chaining(v: bool) { let _ = USE_CHAINING_CELL.set(v); }

pub(crate) fn use_chaining() -> bool {
    *USE_CHAINING_CELL.get_or_init(|| std::env::var("CHAINING").map(|v| v == "1").unwrap_or(false))
}

/// Second-pass relaxed-filter rescue for first-pass failures (--rescue).
/// Off by default; recovered reads are emitted at MAPQ 0 (flagged uncertain).
pub(crate) static RESCUE: AtomicBool = AtomicBool::new(false);
#[inline] pub(crate) fn rescue_pass() -> bool { RESCUE.load(Ordering::Relaxed) }

pub(crate) fn commas(n: u64) -> String {
    let s = n.to_string();
    let mut r = String::new();
    for (i, c) in s.chars().rev().enumerate() {
        if i > 0 && i % 3 == 0 { r.push(','); }
        r.push(c);
    }
    r.chars().rev().collect()
}

// ─────────────────────────────────────────────────────────────────────────────
// 7.  Segment and chromosome processing
// ─────────────────────────────────────────────────────────────────────────────



// ─────────────────────────────────────────────────────────────────────────────
// 8.  FASTA reader (handles .gz and plain)
// ─────────────────────────────────────────────────────────────────────────────


// ─────────────────────────────────────────────────────────────────────────────
// 9.  Built-in demo (no FASTA file)
// ─────────────────────────────────────────────────────────────────────────────



// ─────────────────────────────────────────────────────────────────────────────
// 10.  L2 block extraction (shared by indexing and read mapping)
// ─────────────────────────────────────────────────────────────────────────────

// ─────────────────────────────────────────────────────────────────────────────
// 10.  Recursive N-level block extraction
// ─────────────────────────────────────────────────────────────────────────────
