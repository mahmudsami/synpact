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

/// Per-batch read-sequence budget in MiB (`--batch-mb N`, default 256). Reads are
/// buffered until this many bytes of sequence accumulate, then mapped in parallel.
/// Batching by bytes (not by read count) keeps peak memory independent of read
/// length. `0` = unbounded (load the whole input as one batch). Set from the CLI
/// before mapping; falls back to env BATCH_MB then default 256.
pub(crate) static BATCH_MB_CELL: std::sync::OnceLock<usize> = std::sync::OnceLock::new();

pub(crate) fn set_batch_mb(mb: usize) { let _ = BATCH_MB_CELL.set(mb); }

/// Resolved per-batch budget in *bytes* (`usize::MAX` when unbounded).
pub(crate) fn batch_bytes_budget() -> usize {
    let mb = *BATCH_MB_CELL.get_or_init(|| {
        std::env::var("BATCH_MB").ok().and_then(|v| v.parse().ok()).unwrap_or(256)
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
