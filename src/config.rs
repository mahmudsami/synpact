use std::time::Instant;
use std::sync::atomic::Ordering;
use std::sync::OnceLock;

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

/// Minimum hierarchy level whose blocks may produce anchors, and the lowest level
/// the index stores (`--min-level N`, default 3 — one unified value). Blocks below
/// this level are neither stored nor used; the read is placed only by the coarser,
/// more-unique high-level blocks, which at HiFi error rates removes paralog /
/// segmental-dup ambiguity (higher accuracy + faster). Set once from the CLI.
pub(crate) static MIN_LVL_CELL: std::sync::OnceLock<usize> = std::sync::OnceLock::new();

pub(crate) fn set_min_lvl(v: usize) { let _ = MIN_LVL_CELL.set(v); }

pub(crate) fn min_lvl() -> usize { *MIN_LVL_CELL.get_or_init(|| 3) }

/// Maximum hierarchy level: the highest level the index builds and stores, and the
/// highest level a read is anchored from (`--max-level N`, default `MAX_LEVELS` = 6
/// — one unified value, the ceiling to `--min-level`'s floor). Anchoring starts at
/// the lowest of this, the read's own top level and the index's top level, then
/// falls back down to `--min-level`. Set once from the CLI.
pub(crate) static MAX_LVL_CELL: std::sync::OnceLock<usize> = std::sync::OnceLock::new();

pub(crate) fn set_max_lvl(v: usize) { let _ = MAX_LVL_CELL.set(v); }

pub(crate) fn max_lvl() -> usize { *MAX_LVL_CELL.get_or_init(|| crate::index::MAX_LEVELS) }

/// Near-tie margin, in parts per million of the best locus's votes (`--mapq0-margin F`,
/// default 0 = off). When it is on, a read whose runner-up locus scores at least
/// (1 − F) × the best is placed at its best locus with MAPQ 0, exact ties included
/// (which are otherwise left unmapped). Set once from the CLI.
pub(crate) static MAPQ0_MARGIN_CELL: std::sync::OnceLock<u32> = std::sync::OnceLock::new();

pub(crate) fn set_mapq0_margin(f: f64) { let _ = MAPQ0_MARGIN_CELL.set((f * 1e6).round() as u32); }

pub(crate) fn mapq0_margin_ppm() -> u32 { *MAPQ0_MARGIN_CELL.get_or_init(|| 0) }

/// Near-tie tiebreak (`--tiebreak votes|floor|chain`, default votes): which of the
/// near-tied loci (within `--mapq0-margin`) is reported. `Votes` keeps the heaviest
/// window (the first found on an exact tie); `Floor` re-ranks them by floor-level
/// support (chain::floor_support); `Chain` by colinear chaining score
/// (chain::chain_score). Set once from the CLI.
#[derive(Clone, Copy, PartialEq, Eq, Debug)]
pub(crate) enum Tiebreak { Votes, Floor, Chain }

pub(crate) static TIEBREAK_CELL: std::sync::OnceLock<Tiebreak> = std::sync::OnceLock::new();

pub(crate) fn set_tiebreak(t: Tiebreak) { let _ = TIEBREAK_CELL.set(t); }

pub(crate) fn tiebreak() -> Tiebreak { *TIEBREAK_CELL.get_or_init(|| Tiebreak::Votes) }

/// Read-span switch (`--span-switch F,R`, default off): when the runner-up has at
/// least F × the best's votes, the loci with at least F × the best's votes on both
/// strands are compared by how much of the read their anchors span (first to last
/// anchor); the widest is reported instead of the vote-best if its span is at least
/// R × the vote-best's, at MAPQ 0. Stored as (F in parts per million, R). Set once
/// from the CLI.
pub(crate) static SPAN_SWITCH_CELL: std::sync::OnceLock<Option<(u32, f64)>> = std::sync::OnceLock::new();

pub(crate) fn set_span_switch(f: f64, r: f64) { let _ = SPAN_SWITCH_CELL.set(Some(((f * 1e6).round() as u32, r))); }

pub(crate) fn span_switch() -> Option<(u32, f64)> { *SPAN_SWITCH_CELL.get_or_init(|| None) }

/// Per-batch read-sequence budget in MiB (`--batch-mb N`, default 128). Reads are
/// buffered until this many bytes of sequence accumulate, then mapped in parallel
/// while the next batch is read (so two batches are resident).
/// Batching by bytes (not by read count) keeps peak memory independent of read
/// length. `0` = unbounded (load the whole input as one batch). Set from the CLI
/// before mapping; falls back to env BATCH_MB then default 128.
pub(crate) static BATCH_MB_CELL: std::sync::OnceLock<usize> = std::sync::OnceLock::new();

pub(crate) fn set_batch_mb(mb: usize) { let _ = BATCH_MB_CELL.set(mb); }

/// Resolved per-batch budget in *bytes* (`usize::MAX` when unbounded).
pub(crate) fn batch_bytes_budget() -> usize {
    let mb = *BATCH_MB_CELL.get_or_init(|| {
        std::env::var("BATCH_MB").ok().and_then(|v| v.parse().ok()).unwrap_or(128)
    });
    if mb == 0 { usize::MAX } else { mb << 20 }
}

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
