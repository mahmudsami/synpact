//! syncmer-hifi — hierarchical syncmer-LCP read mapper for PacBio HiFi.
//!
//! Pipeline (see README for the full method description):
//!   reads → open syncmers → locally-consistent parsing (LCP) into a 6-level
//!   block hierarchy → index lookup → conservation-of-mass weighted anchors →
//!   colinear chain DP → MAPQ → PAF (optional chain-guided affine-gap CIGAR).
//!
//! Default preset: k=19, s=10 (syncmer density 1/(k-s+1) = 1/10), k-mer atom.
//!
//!   cargo run --release -- --build-index genome.fa.gz index.idx
//!   cargo run --release -- --map reads.fq.gz index.idx -o out.paf

use std::env;
use std::fs::File;
use std::io::{BufRead, BufReader, BufWriter, Read, Write};
use std::time::Instant;
use std::sync::atomic::{AtomicBool, Ordering};
use flate2::read::MultiGzDecoder;
use rayon::prelude::*;

/// Hierarchy-atom selector.  When true, LCP uses the full k-mer as each leaf
/// value (high specificity, best for low-error HiFi).  When false (default),
/// it uses the middle s-mer (more error-tolerant, best for ONT).
/// Set once at startup from the CLI / loaded index, before any parallel work.
static KMER_ATOM: AtomicBool = AtomicBool::new(false);
#[inline] fn kmer_atom() -> bool { KMER_ATOM.load(Ordering::Relaxed) }

/// Canonical-atom selector.  When true, each atom value is min(encode(atom),
/// encode(revcomp(atom))) so a sequence and its reverse-complement share a key.
/// Halves the atom value space (less specificity) — tested as an alternative.
static CANON_ATOM: AtomicBool = AtomicBool::new(false);
#[inline] fn canon_atom() -> bool { CANON_ATOM.load(Ordering::Relaxed) }

/// Disable the L0 raw-syncmer fallback at query time (env NO_L0=1).
fn no_l0() -> bool {
    use std::sync::OnceLock; static V: OnceLock<bool> = OnceLock::new();
    *V.get_or_init(|| std::env::var("NO_L0").map(|v| v == "1").unwrap_or(false))
}

/// Minimum hierarchy level whose blocks may produce anchors (`--min-level N`,
/// default 3). Blocks below this level are neither emitted nor recursed into — the
/// read is placed only by the coarser, more-unique high-level blocks, which at HiFi
/// error rates removes paralog/segmental-dup ambiguity (higher accuracy + faster).
/// A floor > 0 also disables the L0 raw-syncmer fallback.
/// Set once from the CLI before mapping; falls back to env MIN_LVL then default 3.
static MIN_LVL_CELL: std::sync::OnceLock<usize> = std::sync::OnceLock::new();
fn set_min_lvl(v: usize) { let _ = MIN_LVL_CELL.set(v); }
fn min_lvl() -> usize {
    *MIN_LVL_CELL.get_or_init(|| {
        std::env::var("MIN_LVL").ok().and_then(|v| v.parse().ok()).unwrap_or(3)
    })
}

/// Second-pass relaxed-filter rescue for first-pass failures (--rescue).
/// Off by default; recovered reads are emitted at MAPQ 0 (flagged uncertain).
static RESCUE: AtomicBool = AtomicBool::new(false);
#[inline] fn rescue_pass() -> bool { RESCUE.load(Ordering::Relaxed) }

/// Base-4 encode the reverse-complement of `atom` without allocating.
#[inline]
fn encode_revcomp(atom: &[u8]) -> u64 {
    let mut h = 0u64;
    for &b in atom.iter().rev() {
        let comp = match b.to_ascii_uppercase() {
            b'A' => 3, b'C' => 2, b'G' => 1, b'T' => 0, _ => 0, // complement, base-4
        };
        h = h * 4 + comp;
    }
    h
}

/// The atom value for a window, honoring the canonical flag.
#[inline]
fn atom_value(atom: &[u8]) -> u64 {
    let f = encode_smer(atom);
    if canon_atom() { f.min(encode_revcomp(atom)) } else { f }
}

// ─────────────────────────────────────────────────────────────────────────────
// 1.  Forward-only DNA rolling hash  (replaces the nthash crate)
//
//  Scheme — rotation-XOR, same structure as NT-hash's forward component but
//  with custom seeds and no reverse-complement component:
//
//    H(i) = seed[seq[i]]   rotated (w-1)
//         ^ seed[seq[i+1]] rotated (w-2)
//         ^ …
//         ^ seed[seq[i+w-1]]   (no rotation)
//
//  Rolling O(1) update:
//    H(i+1) = rotate_left(H(i), 1)
//            ^ rotate_left(seed[seq[i]], w)   ← removes outgoing base
//            ^ seed[seq[i+w]]                 ← adds incoming base
//
//  Seeds: the four SplitMix64/φ constants used throughout the codebase.
//  They have excellent 64-bit avalanche (≈50% bit-flip per input bit) and
//  are well-separated in Hamming distance from each other.
// ─────────────────────────────────────────────────────────────────────────────

/// Per-base seeds for the rolling hash. Index by raw byte value.
/// Non-ACGT bytes (including 'N') map to 0.
const DNA_SEED: [u64; 256] = {
    let mut t = [0u64; 256];
    t[b'A' as usize] = 0x9e3779b97f4a7c15;  // floor(2^64 / φ)
    t[b'a' as usize] = 0x9e3779b97f4a7c15;
    t[b'C' as usize] = 0x6c62272e07bb0142;  // FNV-1a basis
    t[b'c' as usize] = 0x6c62272e07bb0142;
    t[b'G' as usize] = 0xbf58476d1ce4e5b9;  // SplitMix64 multiplier 1
    t[b'g' as usize] = 0xbf58476d1ce4e5b9;
    t[b'T' as usize] = 0x94d049bb133111eb;  // SplitMix64 multiplier 2
    t[b't' as usize] = 0x94d049bb133111eb;
    t
};

/// Rolling forward-only DNA hash over a window of `w` consecutive bases.
/// Yields one `u64` per starting position.  An N inside the current window
/// makes that window's hash `u64::MAX` (so it is never selected as a minimum).
struct DnaHashFwd<'a> {
    seq:  &'a [u8],
    w:    usize,
    h:    u64,
    n_in_window: u32,  // count of N-like bases currently inside the window
    pos:  usize,
}

impl<'a> DnaHashFwd<'a> {
    /// Returns `None` when `seq.len() < w`.
    fn new(seq: &'a [u8], w: usize) -> Option<Self> {
        if seq.len() < w { return None; }
        let mut h = 0u64;
        let mut n_count = 0u32;
        for &b in &seq[..w] {
            if DNA_SEED[b as usize] == 0 { n_count += 1; }
            h = h.rotate_left(1) ^ DNA_SEED[b as usize];
        }
        Some(DnaHashFwd { seq, w, h, n_in_window: n_count, pos: 0 })
    }
}

impl<'a> Iterator for DnaHashFwd<'a> {
    type Item = u64;
    #[inline]
    fn next(&mut self) -> Option<u64> {
        if self.pos + self.w > self.seq.len() { return None; }
        let h = if self.n_in_window > 0 { u64::MAX } else { self.h };
        // Roll forward one position (prepare hash for pos+1)
        if self.pos + self.w < self.seq.len() {
            let out = self.seq[self.pos];
            let in_ = self.seq[self.pos + self.w];
            if DNA_SEED[out as usize] == 0 { self.n_in_window -= 1; }
            if DNA_SEED[in_  as usize] == 0 { self.n_in_window += 1; }
            self.h = self.h.rotate_left(1)
                ^ DNA_SEED[out as usize].rotate_left(self.w as u32)
                ^ DNA_SEED[in_  as usize];
        }
        self.pos += 1;
        Some(h)
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// 2.  s-mer direct encoding  (A=0, C=1, G=2, T=3, base-4 positional)
// ─────────────────────────────────────────────────────────────────────────────

#[inline]
fn base4(b: u8) -> u64 {
    match b.to_ascii_uppercase() {
        b'A' => 0, b'C' => 1, b'G' => 2, b'T' => 3, _ => 0,
    }
}

#[inline]
fn encode_smer(smer: &[u8]) -> u64 {
    smer.iter().fold(0u64, |acc, &b| acc * 4 + base4(b))
}

// ─────────────────────────────────────────────────────────────────────────────
// 3.  Level-1 block hash
//
// A block is represented by its ordered sequence of base-4 s-mer values.
// Two blocks with the same s-mer sequence (anywhere in the genome) get the
// same hash, making the hash a canonical block identifier.
//
// Design to minimise collisions:
//   • Length is folded in first → blocks of different sizes can't collide
//     unless the hash itself collides.
//   • Each value is mixed via a full SplitMix64 step → ~50% bit-flip on any
//     single-bit change (avalanche).
//   • Position-sensitivity: the rotation before adding each value means the
//     same values in different order produce different hashes.
//
// 64-bit output.  Expected collisions for N=10^8 distinct blocks ≈ N²/2⁶⁴ ≈ 0.5%.
// ─────────────────────────────────────────────────────────────────────────────

// Domain constants keep every level's hash space disjoint.
// Index 0 → L1, index 1 → L2, … up to L8; beyond that we derive more.
const LEVEL_DOMAINS: [u64; 8] = [
    0x9e3779b97f4a7c15, // L1  (φ-based)
    0x6c62272e07bb0142, // L2  (FNV prime)
    0xd2a98b26625eee7b, // L3
    0xa3b195354a2b7623, // L4
    0x1b03738712fad5c9, // L5
    0xc4ceb9fe1a85ec53, // L6
    0x517cc1b727220a95, // L7
    0x9be1aa58ba6a9f81, // L8
];

#[inline]
fn level_domain(level_0idx: usize) -> u64 {
    LEVEL_DOMAINS.get(level_0idx).copied().unwrap_or_else(||
        LEVEL_DOMAINS[7].wrapping_add(
            (level_0idx as u64).wrapping_mul(0x6c62272e07bb0142)))
}

/// Hash a block's values at the given 0-indexed level.
#[inline]
fn block_hash_for_level(values: &[u64], level_0idx: usize) -> u64 {
    block_hash_with_domain(values, level_domain(level_0idx))
}

// Backward-compatible aliases used by the demo/stats display code.

fn block_hash_with_domain(values: &[u64], domain: u64) -> u64 {
    // Seed mixes in both length and domain so that:
    //   • blocks of different lengths can't collide (even with same prefix)
    //   • L1 and L2 occupy disjoint hash spaces
    let mut h: u64 = (values.len() as u64)
        .wrapping_mul(0xbf58476d1ce4e5b9)
        ^ domain;
    for &v in values {
        h = h.rotate_left(13).wrapping_add(v);
        // SplitMix64 avalanche step
        h ^= h >> 30;
        h = h.wrapping_mul(0xbf58476d1ce4e5b9);
        h ^= h >> 27;
        h = h.wrapping_mul(0x94d049bb133111eb);
        h ^= h >> 31;
    }
    h
}


// ─────────────────────────────────────────────────────────────────────────────
// 4.  Syncmer types
// ─────────────────────────────────────────────────────────────────────────────


/// Lightweight record for genome-scale processing.
#[derive(Clone)]
struct SyncmerLight {
    pos:   u32,
    value: u64,   // s-mer base-4 encoding — used for L1+ block hashing
}


fn select_syncmers_light(seq: &[u8], k: usize, s: usize, t: usize) -> Vec<SyncmerLight> {
    if seq.len() < k { return vec![]; }
    // Precompute all s-mer NT-hashes in one rolling pass.
    let smer_hashes: Vec<u64> = match DnaHashFwd::new(seq, s) {
        Some(iter) => iter.collect(),
        None => return vec![],
    };
    let mut out = Vec::new();
    for i in 0..=(seq.len() - k) {
        let kmer = &seq[i..i + k];
        if kmer.iter().any(|&b| matches!(b.to_ascii_uppercase(), b'N')) { continue; }
        let (_, min_j) = (0..=(k - s)).map(|j| (smer_hashes[i + j], j)).min().unwrap();
        if min_j == t {
            // Atom: full k-mer (high specificity, HiFi) or middle s-mer
            // (error-tolerant, ONT default).  k≤32 keeps the base-4 fit in u64.
            let atom: &[u8] = if kmer_atom() && k <= 32 { kmer } else { &kmer[t..t + s] };
            out.push(SyncmerLight { pos: i as u32, value: atom_value(atom) });
        }
    }
    out
}


/// Seed selection mode.  Kept as a single-variant enum so the on-disk index
/// format (which stores a seed-mode byte) stays stable and self-describing.
#[derive(Clone, Copy, PartialEq, Eq, Debug)]
enum SeedMode {
    /// Open syncmer: min s-mer at position t within the k-mer.
    Syncmer,
}

/// Seed extraction — open syncmers (s = s-mer size, t = position).
#[inline]
fn select_seeds_light(seq: &[u8], k: usize, s: usize, t: usize, _mode: SeedMode)
    -> Vec<SyncmerLight>
{
    select_syncmers_light(seq, k, s, t)
}

// ─────────────────────────────────────────────────────────────────────────────
// 5.  Locally consistent parsing — level-1 blocks
// ─────────────────────────────────────────────────────────────────────────────

#[derive(Debug)]
struct Block {
    indices: Vec<usize>,
    rule:    &'static str,
}

/// Cole–Vishkin deterministic-coin-tossing tag between two consecutive values.
/// = 2 * (index of lowest bit where they differ) + (that bit's value in `cur`).
/// Tags vary non-monotonically even over a sorted run, so their local maxima
/// give content-determined (locally consistent) split points.
#[inline]
fn dct_tag(prev: u64, cur: u64) -> u32 {
    let d = prev ^ cur;
    if d == 0 { return 0; }
    let bit = d.trailing_zeros();
    2 * bit + ((cur >> bit) & 1) as u32
}

/// Split a monotone run of indices (whose `values` are strictly monotone) into
/// consecutive blocks of size ≤3, deterministically, by Cole–Vishkin DCT.
///
/// Boundaries = run start, run end, and every interior position whose DCT tag
/// is a local maximum (tag[i] > tag[i-1] && tag[i] >= tag[i+1]).  Consecutive
/// boundaries are emitted as blocks INCLUSIVE of both ends (adjacent blocks
/// share their boundary index — same intentional-intersection style as the
/// other rules).  Any gap that still exceeds 2 is chunked left-to-right in
/// steps of 2 (measured from the local boundary, so still locally determined),
/// guaranteeing every emitted block spans ≤3 indices.
fn split_monotone_run(run: &[usize], values: &[u64]) -> Vec<Vec<usize>> {
    let m = run.len();
    if m <= 3 { return vec![run.to_vec()]; }

    // DCT tag of each run position relative to its predecessor within the run.
    let tag: Vec<u32> = (0..m)
        .map(|i| if i == 0 { 0 } else { dct_tag(values[run[i-1]], values[run[i]]) })
        .collect();

    // Boundary positions (indices into `run`): start, local-max tags, end.
    let mut bounds = vec![0usize];
    for i in 1..m - 1 {
        if tag[i] > tag[i-1] && tag[i] >= tag[i+1] { bounds.push(i); }
    }
    if *bounds.last().unwrap() != m - 1 { bounds.push(m - 1); }

    // Emit inclusive blocks between consecutive boundaries; chunk any gap >2.
    let mut out: Vec<Vec<usize>> = Vec::new();
    for w in bounds.windows(2) {
        let (mut lo, hi) = (w[0], w[1]);
        while hi - lo > 2 {
            out.push(run[lo..=lo + 2].to_vec());  // size 3, shares end with next
            lo += 2;
        }
        out.push(run[lo..=hi].to_vec());           // size ≤3
    }
    out
}

#[inline]
fn is_local_min(v: &[u64], i: usize) -> bool {
    let n = v.len();
    (i == 0 || v[i] < v[i - 1]) && (i == n - 1 || v[i] < v[i + 1])
}
#[inline]
fn is_local_max(v: &[u64], i: usize) -> bool {
    let n = v.len();
    (i == 0 || v[i] > v[i - 1]) && (i == n - 1 || v[i] > v[i + 1])
}

fn locally_consistent_parsing(values: &[u64]) -> Vec<Block> {
    let n = values.len();
    if n == 0 { return vec![]; }
    if n == 1 { return vec![Block { indices: vec![0], rule: "local_min" }]; }

    // Vec<bool> instead of HashSet — better cache behaviour at genome scale
    let is_min: Vec<bool> = (0..n).map(|i| is_local_min(values, i)).collect();
    let is_max: Vec<bool> = (0..n).map(|i| is_local_max(values, i)).collect();

    let mut assigned = vec![false; n];
    let mut blocks: Vec<Block> = Vec::new();

    // Inclusive commit: all positions including already-assigned (all rules)
    macro_rules! commit_incl {
        ($idxs:expr, $rule:expr) => {{
            let all: Vec<usize> = $idxs.iter().copied().filter(|&i| i < n).collect();
            if all.iter().any(|&i| !assigned[i]) {
                for &i in &all { assigned[i] = true; }
                blocks.push(Block { indices: all, rule: $rule });
            }
        }};
    }

    // Rule 1 — local-minimum block: always full triplet {i-1, i, i+1}
    // Uses commit_incl so adjacent minima each get their full triplet
    // (shared boundary position appears in both — intentional intersection).
    for i in 0..n {
        if is_min[i] {
            let lo = i.saturating_sub(1);
            let hi = (i + 1).min(n - 1);
            commit_incl!((lo..=hi).collect::<Vec<_>>(), "local_min");
        }
    }

    // Rule 2 — local-maximum block (no adjacent local minimum)
    // Also commit_incl so {i-1, i, i+1} always present.
    for i in 0..n {
        if is_max[i] {
            let l_is_min = i > 0     && is_min[i - 1];
            let r_is_min = i + 1 < n && is_min[i + 1];
            if !l_is_min && !r_is_min {
                let lo = i.saturating_sub(1);
                let hi = (i + 1).min(n - 1);
                commit_incl!((lo..=hi).collect::<Vec<_>>(), "local_max");
            }
        }
    }

    // Rule 3 — repetition run ≥ 2 plus immediate prev/next neighbours
    {
        let mut i = 0;
        while i < n {
            let mut j = i + 1;
            while j < n && values[j] == values[i] { j += 1; }
            if j - i >= 2 {
                let lo = i.saturating_sub(1);
                let hi = j.min(n - 1); // j is one past the run; j is the next syncmer
                commit_incl!((lo..=hi).collect::<Vec<_>>(), "repetition");
            }
            i = j;
        }
    }

    // Rule 4 — full monotone run [a..=b], anchors included (intentional intersection)
    for a in 0..n.saturating_sub(2) {
        if !assigned[a] { continue; }
        for &sign in &[1i128, -1i128] {
            let mut b = a + 1;
            while b < n {
                let d = values[b] as i128 - values[b - 1] as i128;
                if d * sign > 0 { b += 1; } else { break; }
            }
            b -= 1;
            if b > a + 1 && assigned[b] && (a + 1..b).any(|j| !assigned[j]) {
                let run: Vec<usize> = (a..=b).collect();
                let rule = if sign == 1 { "mono_inc" } else { "mono_dec" };
                // Split long monotone runs into ≤3-unit blocks via DCT so they
                // aren't one big featureless block (deterministic & locally
                // consistent → genome and read split identically).
                for sub in split_monotone_run(&run, values) {
                    commit_incl!(sub, rule);
                }
            }
        }
    }

    // Panic on unassigned positions — indicates a bug in the algorithm
    let bad: Vec<usize> = (0..n).filter(|&i| !assigned[i]).collect();
    if !bad.is_empty() {
        panic!(
            "LCP bug: {} position(s) unassigned (first few: {:?})\nvalues={:?}",
            bad.len(), &bad[..bad.len().min(5)], values
        );
    }

    blocks.sort_by_key(|b| b.indices[0]);
    blocks
}

// ─────────────────────────────────────────────────────────────────────────────
// 6.  Statistics
// ─────────────────────────────────────────────────────────────────────────────



fn commas(n: u64) -> String {
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

/// Extract all hierarchy levels from one N-free segment.
/// Returns result[0] = L1 anchors, result[1] = L2, … up to max_levels.
/// Used during index construction.
fn extract_all_levels(seq: &[u8], k: usize, s: usize, t: usize, max_levels: usize,
                      mode: SeedMode) -> Vec<Vec<(u64, u32)>>
{
    let syncmers = select_seeds_light(seq, k, s, t, mode);
    if syncmers.is_empty() { return vec![]; }

    let smer_vals: Vec<u64> = syncmers.iter().map(|sm| sm.value).collect();
    let l1_raw = locally_consistent_parsing(&smer_vals);

    let mut cur_hashes: Vec<u64> = Vec::with_capacity(l1_raw.len());
    let mut cur_pos:    Vec<u32> = Vec::with_capacity(l1_raw.len());
    let mut first_out: Vec<(u64, u32)> = Vec::with_capacity(l1_raw.len());

    for blk in &l1_raw {
        let bvals: Vec<u64> = blk.indices.iter().map(|&i| smer_vals[i]).collect();
        let h   = block_hash_for_level(&bvals, 0);
        let pos = syncmers[blk.indices[0]].pos;
        cur_hashes.push(h);
        cur_pos.push(pos);
        first_out.push((h, pos));
    }
    // all[0] = L1 blocks, all[1] = L2 blocks, …
    let mut all: Vec<Vec<(u64, u32)>> = vec![first_out];

    for level_1idx in 2..=max_levels {
        if cur_hashes.len() < 2 { break; }
        let next_raw = locally_consistent_parsing(&cur_hashes);
        if next_raw.is_empty() { break; }

        let mut next_hashes = Vec::with_capacity(next_raw.len());
        let mut next_pos    = Vec::with_capacity(next_raw.len());
        let mut level_out   = Vec::with_capacity(next_raw.len());

        for blk in &next_raw {
            let bvals: Vec<u64> = blk.indices.iter().map(|&i| cur_hashes[i]).collect();
            let h   = block_hash_for_level(&bvals, level_1idx - 1);
            let pos = cur_pos[blk.indices[0]];
            next_hashes.push(h);
            next_pos.push(pos);
            level_out.push((h, pos));
        }

        all.push(level_out);
        cur_hashes = next_hashes;
        cur_pos    = next_pos;
    }

    all
}

/// A block at any hierarchy level, with its children at level-1.
/// level is 1-indexed: 1 = L1 (leaf), 2 = L2, …
/// Children are empty for L1 leaves.
/// `pos`..`end` is the half-open query interval covered by this block (in bp).
#[derive(Clone)]
struct HierBlock {
    level:    usize,
    hash:     u64,
    pos:      u32,   // query start (inclusive)
    end:      u32,   // query end   (exclusive); span = end - pos
    #[allow(dead_code)]
    rule:     &'static str,  // LCP rule that formed this block (provenance)
    /// Conservation-of-mass weight: every syncmer carries mass 1; a block's mass
    /// is the sum of its units' masses, with a unit shared by `m` blocks
    /// contributing 1/m to each.  Total mass is conserved across all levels.
    mass:     f32,
    children: Vec<HierBlock>,
}

/// Build the recursive HierBlock tree for one read strand.
/// Returns top-level blocks (highest level reachable from this read).
/// num_levels is taken from the loaded index so we don't exceed what was indexed.
fn extract_hier_blocks_n(seq: &[u8], k: usize, s: usize, t: usize, num_levels: usize,
                         mode: SeedMode) -> Vec<HierBlock>
{
    let syncmers = select_seeds_light(seq, k, s, t, mode);
    if syncmers.is_empty() { return vec![]; }

    let smer_vals: Vec<u64> = syncmers.iter().map(|sm| sm.value).collect();
    let l1_raw = locally_consistent_parsing(&smer_vals);

    // ── L1 masses ────────────────────────────────────────────────────────────
    // Each syncmer carries mass 1.  A syncmer shared by `m` L1 blocks splits its
    // mass 1/m across them.  A block's mass = sum of its syncmers' shares.
    let mut sync_membership = vec![0u32; syncmers.len()];
    for blk in &l1_raw { for &i in &blk.indices { sync_membership[i] += 1; } }

    let mut cur_blocks: Vec<HierBlock> = l1_raw.iter().map(|blk| {
        let bvals: Vec<u64> = blk.indices.iter().map(|&i| smer_vals[i]).collect();
        let h   = block_hash_for_level(&bvals, 0);
        let pos = syncmers[blk.indices[0]].pos;
        let end = syncmers[*blk.indices.last().unwrap()].pos + k as u32;
        let mass: f32 = blk.indices.iter()
            .map(|&i| 1.0 / sync_membership[i] as f32).sum();
        HierBlock { level: 1, hash: h, pos, end, rule: blk.rule,
                    mass, children: vec![] }
    }).collect();

    if cur_blocks.is_empty() { return vec![]; }

    // Iteratively build L2, L3, … up to num_levels
    for level_1idx in 2..=num_levels {
        if cur_blocks.len() < 2 { break; }

        // Snapshot hashes+positions+ends+masses before moving ownership
        let cur_hashes: Vec<u64> = cur_blocks.iter().map(|b| b.hash).collect();
        let cur_pos:    Vec<u32> = cur_blocks.iter().map(|b| b.pos).collect();
        let cur_end:    Vec<u32> = cur_blocks.iter().map(|b| b.end).collect();
        let cur_mass:   Vec<f32> = cur_blocks.iter().map(|b| b.mass).collect();
        let next_raw = locally_consistent_parsing(&cur_hashes);
        if next_raw.is_empty() { break; }

        // Membership: how many parent blocks each child block belongs to.
        let mut child_membership = vec![0u32; cur_blocks.len()];
        for blk in &next_raw { for &i in &blk.indices { child_membership[i] += 1; } }

        let prev_blocks = std::mem::take(&mut cur_blocks);
        cur_blocks = next_raw.iter().map(|blk| {
            let bvals: Vec<u64> = blk.indices.iter().map(|&i| cur_hashes[i]).collect();
            let h        = block_hash_for_level(&bvals, level_1idx - 1);
            let pos      = cur_pos[blk.indices[0]];
            let end      = cur_end[*blk.indices.last().unwrap()];
            // Parent mass = sum of each child's mass / (#parents sharing that child).
            let mass: f32 = blk.indices.iter()
                .map(|&i| cur_mass[i] / child_membership[i] as f32).sum();
            let children = blk.indices.iter().map(|&i| prev_blocks[i].clone()).collect();
            HierBlock { level: level_1idx, hash: h, pos, end, rule: blk.rule,
                        mass, children }
        }).collect();
    }

    cur_blocks  // top-level blocks for this read
}

// ─────────────────────────────────────────────────────────────────────────────
// 11.  Genome index (sorted flat array, binary-searched at query time)
// ─────────────────────────────────────────────────────────────────────────────

struct GIndex {
    /// levels[0] = L0 raw-syncmer entries (k-mer NT-hash, built sequentially).
    /// levels[1] = L1 block entries, levels[2] = L2, …
    /// Each level is stored struct-of-arrays (SoA): three parallel arrays sorted
    /// jointly by hash.  This is 13 bytes/entry (vs 16 for the AoS `(u64,u8,u32)`
    /// which pads to 16), an 18.75% memory saving, and — more importantly — the
    /// binary search touches only the contiguous `hashes` array (8 B/entry
    /// scanned instead of 16), roughly halving cache pressure during lookup.
    levels:    Vec<Level>,
    chr_names: Vec<String>,
}

/// One hierarchy level, struct-of-arrays.  `hashes` is sorted ascending; the
/// other two arrays are permuted in lockstep so index i refers to one entry.
#[derive(Default)]
struct Level {
    hashes: Vec<u64>,
    chrs:   Vec<u8>,
    poss:   Vec<u32>,
}

impl Level {
    #[inline] fn len(&self) -> usize { self.hashes.len() }
    #[inline] fn is_empty(&self) -> bool { self.hashes.is_empty() }

    /// Build a SoA level from an AoS entry vector (used during index build).
    fn from_aos(entries: Vec<(u64, u8, u32)>) -> Self {
        let n = entries.len();
        let mut hashes = Vec::with_capacity(n);
        let mut chrs   = Vec::with_capacity(n);
        let mut poss   = Vec::with_capacity(n);
        for (h, c, p) in entries {
            hashes.push(h); chrs.push(c); poss.push(p);
        }
        Level { hashes, chrs, poss }
    }
}

/// A view over the matching (chr, pos) entries for a looked-up hash.
struct Hits<'a> {
    chrs: &'a [u8],
    poss: &'a [u32],
}

impl<'a> Hits<'a> {
    #[inline] fn is_empty(&self) -> bool { self.chrs.is_empty() }
    /// Iterate (chr_id, genome_pos) pairs.
    #[inline]
    fn iter(&self) -> impl Iterator<Item = (u8, u32)> + '_ {
        self.chrs.iter().copied().zip(self.poss.iter().copied())
    }
}

impl GIndex {
    fn num_levels(&self) -> usize { self.levels.len() }

    /// Look up hash h. Returns (hits, too_repetitive):
    ///   - hits non-empty when found and occ ≤ max_occ
    ///   - too_repetitive=true means hash exists but occ > max_occ
    ///     (children are likely also repetitive — caller should NOT fall back)
    /// Out-of-range levels (e.g. a top-level block above the indexed depth)
    /// return empty hits, matching the original AoS `.get()` behaviour.
    /// The binary search touches only the contiguous `hashes` array.
    fn lookup_with_status(&self, level_0idx: usize, h: u64, max_occ: usize)
        -> (Hits<'_>, bool)
    {
        let empty = Hits { chrs: &[], poss: &[] };
        let Some(lvl) = self.levels.get(level_0idx) else { return (empty, false) };
        let lo = lvl.hashes.partition_point(|&x| x < h);
        let hi = lvl.hashes.partition_point(|&x| x <= h);
        if hi - lo > max_occ {
            (empty, true)
        } else {
            (Hits { chrs: &lvl.chrs[lo..hi], poss: &lvl.poss[lo..hi] }, false)
        }
    }

    /// Convenience wrapper — returns empty hits on not-found OR too-repetitive.
    fn lookup(&self, level_0idx: usize, h: u64, max_occ: usize) -> Hits<'_> {
        self.lookup_with_status(level_0idx, h, max_occ).0
    }
}

/// Read all chromosomes from a FASTA file (plain or .gz) into memory.
fn read_fasta_chrs(path: &str) -> Vec<(String, Vec<u8>)> {
    let file = File::open(path).unwrap_or_else(|e| panic!("Cannot open {path}: {e}"));
    let reader: Box<dyn BufRead> = if path.ends_with(".gz") {
        Box::new(BufReader::with_capacity(1 << 20, MultiGzDecoder::new(file)))
    } else {
        Box::new(BufReader::with_capacity(1 << 20, file))
    };
    let mut chrs: Vec<(String, Vec<u8>)> = Vec::new();
    let mut name = String::new();
    let mut seq: Vec<u8> = Vec::with_capacity(1 << 28);
    for line in reader.lines() {
        let line = line.expect("read error");
        let bytes = line.trim_end().as_bytes();
        if bytes.first() == Some(&b'>') {
            if !seq.is_empty() {
                chrs.push((name.clone(), std::mem::take(&mut seq)));
            }
            name = std::str::from_utf8(&bytes[1..])
                .unwrap_or("?").split_whitespace().next().unwrap_or("?")
                .to_string();
        } else {
            seq.extend(bytes.iter().map(|b| b.to_ascii_uppercase()));
        }
    }
    if !seq.is_empty() { chrs.push((name, seq)); }
    chrs
}

/// Extract entries for all N levels from one chromosome (N-split on 'N' runs).
/// Returns a Vec of N vecs: result[0] = L1 entries, result[1] = L2, …
fn chr_to_entries_n(chr_id: u8, seq: &[u8], k: usize, s: usize, t: usize, max_levels: usize,
                    mode: SeedMode) -> Vec<Vec<(u64, u8, u32)>>
{
    // slots 0..max_levels-1 = L1-L{max_levels} blocks
    let mut per_level: Vec<Vec<(u64, u8, u32)>> = vec![Vec::new(); max_levels];
    let mut seg_start: Option<usize> = None;

    let flush = |s0: usize, end: usize, per_level: &mut Vec<Vec<(u64, u8, u32)>>| {
        if end - s0 < k { return; }
        let levels = extract_all_levels(&seq[s0..end], k, s, t, max_levels, mode);
        for (li, entries) in levels.into_iter().enumerate() {
            if li < per_level.len() {
                for (h, pos) in entries {
                    per_level[li].push((h, chr_id, s0 as u32 + pos));
                }
            }
        }
    };

    for (i, &b) in seq.iter().enumerate() {
        match (seg_start, matches!(b, b'N')) {
            (None, false) => seg_start = Some(i),
            (Some(s0), true) => { flush(s0, i, &mut per_level); seg_start = None; }
            _ => {}
        }
    }
    if let Some(s0) = seg_start { flush(s0, seq.len(), &mut per_level); }

    per_level
}

/// Build the L0 (raw-syncmer) index level sequentially in fixed-size chunks.
///
/// Processes one chromosome at a time, each chromosome in CHUNK_SIZE-bp slices.
/// Peak memory per chunk ≈ 2 × CHUNK_SIZE × 8 bytes (smer + kmer hash arrays) = 160 MB.
/// The growing `all_l0` accumulator is the dominant memory user (~700M × 13 bytes ≈ 9 GB
/// for the human genome before filtering); it is safe on a 32 GB machine.
///
/// Returns a sorted, occ-filtered slice ready for `GIndex.levels[0]`.
fn extract_l0_sequential(
    chrs: &[(String, Vec<u8>)],
    k: usize, s: usize, t: usize,
    mode: SeedMode,
    max_occ: usize,
) -> Vec<(u64, u8, u32)>
{
    const CHUNK_SIZE: usize = 10_000_000; // 10 Mbp → ~80 MB per hash array
    let overlap = k - 1; // left extension so syncmers near chunk boundaries aren't missed

    let mut all_l0: Vec<(u64, u8, u32)> = Vec::new();

    for (chr_id, (_name, seq)) in chrs.iter().enumerate() {
        if seq.len() < k { continue; }
        let chr_id = chr_id as u8;

        let mut owned_start = 0usize;
        while owned_start < seq.len() {
            // Extend chunk leftward so the full k-mer window is available for
            // syncmers whose left edge is just at `owned_start`.
            let chunk_start = owned_start.saturating_sub(overlap);
            let chunk_end   = (owned_start + CHUNK_SIZE).min(seq.len());

            let chunk = &seq[chunk_start..chunk_end];
            if chunk.len() < k { owned_start += CHUNK_SIZE; continue; }

            // Syncmer positions within this chunk (N-containing k-mers already skipped).
            let syncmers = select_seeds_light(chunk, k, s, t, mode);

            // k-mer NT-hashes for the chunk — the lookup keys stored in L0.
            let kmer_hashes: Vec<u64> = match DnaHashFwd::new(chunk, k) {
                Some(it) => it.collect(),
                None     => { owned_start += CHUNK_SIZE; continue; }
            };

            for sm in &syncmers {
                let abs_pos = chunk_start + sm.pos as usize;
                // Assign each syncmer to exactly one chunk (the one that "owns" its abs pos).
                if abs_pos >= owned_start && abs_pos < chunk_end {
                    let local = sm.pos as usize;
                    if local < kmer_hashes.len() {
                        all_l0.push((kmer_hashes[local], chr_id, abs_pos as u32));
                    }
                }
            }

            owned_start += CHUNK_SIZE;
        }
    }

    // Sort by hash for binary-search lookup, then filter repetitive k-mers.
    all_l0.sort_unstable_by_key(|e| e.0);
    let before = all_l0.len();
    let mut filtered: Vec<(u64, u8, u32)> = Vec::with_capacity(before);
    let mut i = 0;
    while i < all_l0.len() {
        let h  = all_l0[i].0;
        let mut j = i + 1;
        while j < all_l0.len() && all_l0[j].0 == h { j += 1 }
        if j - i <= max_occ { filtered.extend_from_slice(&all_l0[i..j]); }
        i = j;
    }
    eprintln!("      L0 raw syncmers: {} total → {} after occ≤{} filter",
              commas(before as u64), commas(filtered.len() as u64), max_occ);
    filtered
}

fn build_index(genome_path: &str, k: usize, s: usize, t: usize, max_levels: usize,
               mode: SeedMode) -> (GIndex, std::time::Duration)
{
    let t0 = Instant::now();

    print!("    reading sequences ... ");
    std::io::stdout().flush().ok();
    let chrs = read_fasta_chrs(genome_path);
    println!("{} chromosomes", chrs.len());

    let chr_names: Vec<String> = chrs.iter().map(|(n, _)| n.clone()).collect();

    print!("    extracting {max_levels} levels of LCP blocks (parallel) ... ");
    std::io::stdout().flush().ok();

    // Parallel per-chromosome extraction
    let all: Vec<Vec<Vec<(u64, u8, u32)>>> = chrs.par_iter().enumerate()
        .map(|(chr_id, (_, seq))| chr_to_entries_n(chr_id as u8, seq, k, s, t, max_levels, mode))
        .collect();

    // Merge across chromosomes
    let mut merged: Vec<Vec<(u64, u8, u32)>> = vec![Vec::new(); max_levels];
    for chr_levels in all {
        for (li, entries) in chr_levels.into_iter().enumerate() {
            if li < merged.len() { merged[li].extend(entries); }
        }
    }
    // Drop trailing empty levels (can happen for short genomes / large k)
    while merged.last().map_or(false, |v| v.is_empty()) { merged.pop(); }

    print!("    sorting L1-L{max_levels} ... ");
    std::io::stdout().flush().ok();
    for level in &mut merged {
        level.sort_unstable_by_key(|e| e.0);
    }
    println!("done");

    // Build L0 (raw syncmers) sequentially in 10 Mbp chunks to stay within
    // memory budget.  L1-L{max_levels} are already sorted in `merged`.
    print!("    building L0 raw-syncmer level (sequential, 10 Mbp chunks) ... ");
    std::io::stdout().flush().ok();
    let l0 = extract_l0_sequential(&chrs, k, s, t, mode, MAX_OCC_L1_DEFAULT);
    println!("done  ({} entries, occ≤{})", commas(l0.len() as u64), MAX_OCC_L1_DEFAULT);

    // Prepend L0 so GIndex.levels[0]=L0, levels[1]=L1, …
    merged.insert(0, l0);

    let level_counts: Vec<String> = merged.iter().enumerate()
        .map(|(li, v)| format!("L{}: {}", li, commas(v.len() as u64)))
        .collect();
    println!("    {}", level_counts.join("  "));

    // Convert AoS build buffers → SoA levels.
    let levels: Vec<Level> = merged.into_iter().map(Level::from_aos).collect();

    let elapsed = t0.elapsed();
    (GIndex { levels, chr_names }, elapsed)
}

// ── Index serialisation ───────────────────────────────────────────────────────
// Format v7  ("SYNCL2\x07\x00"):
//   8  bytes  magic
//   4  bytes  k  (u32 LE)
//   4  bytes  s  (u32 LE)
//   4  bytes  t  (u32 LE)
//   1  byte   seed_mode  (0 = Syncmer)
//   1  byte   atom flags (bit0 = k-mer atom, bit1 = canonical)
//   4  bytes  num_chrs  (u32 LE)
//   for each chr: 4-byte len + UTF-8 bytes
//   4  bytes  num_levels  (u32 LE)
//   for level 0 … num_levels-1:
//       8  bytes  num_entries  (u64 LE)
//       num_entries × 13 bytes  (8-byte hash + 1-byte chr_id + 4-byte pos)
//   Level 0 = L0 raw-syncmer fallback, level 1 = L1 blocks, …
// Older formats (v4–v6) still load (atom flags default off).

/// Write one level in the packed 13-bytes/entry on-disk format (unchanged from
/// the AoS era, so existing v6 index files remain compatible).
fn write_entries(w: &mut impl Write, lvl: &Level) {
    w.write_all(&(lvl.len() as u64).to_le_bytes()).unwrap();
    for i in 0..lvl.len() {
        w.write_all(&lvl.hashes[i].to_le_bytes()).unwrap();
        w.write_all(&[lvl.chrs[i]]).unwrap();
        w.write_all(&lvl.poss[i].to_le_bytes()).unwrap();
    }
}

fn save_index(idx: &GIndex, path: &str, k: usize, s: usize, t: usize, mode: SeedMode) {
    let file = File::create(path)
        .unwrap_or_else(|e| panic!("Cannot create {path}: {e}"));
    let mut w = BufWriter::with_capacity(1 << 23, file);
    // v7: like v6 but adds a 1-byte atom flag (1 = k-mer atom, 0 = s-mer atom)
    // immediately after the seed_mode byte.  v6 loads with atom defaulted to 0.
    w.write_all(b"SYNCL2\x07\x00").unwrap();
    w.write_all(&(k as u32).to_le_bytes()).unwrap();
    w.write_all(&(s as u32).to_le_bytes()).unwrap();
    w.write_all(&(t as u32).to_le_bytes()).unwrap();
    w.write_all(&[mode as u8]).unwrap();
    // v7 flag byte: bit0 = k-mer atom, bit1 = canonical atom.
    let flags = (kmer_atom() as u8) | ((canon_atom() as u8) << 1);
    w.write_all(&[flags]).unwrap();
    w.write_all(&(idx.chr_names.len() as u32).to_le_bytes()).unwrap();
    for name in &idx.chr_names {
        let b = name.as_bytes();
        w.write_all(&(b.len() as u32).to_le_bytes()).unwrap();
        w.write_all(b).unwrap();
    }
    w.write_all(&(idx.levels.len() as u32).to_le_bytes()).unwrap();
    for lvl in &idx.levels {
        write_entries(&mut w, lvl);
    }
    w.flush().unwrap();
}

#[inline]
fn read_u32le(r: &mut impl Read) -> u32 {
    let mut b = [0u8; 4]; r.read_exact(&mut b).unwrap(); u32::from_le_bytes(b)
}
fn load_index(path: &str) -> (GIndex, usize, usize, usize, SeedMode) {
    let file = File::open(path)
        .unwrap_or_else(|e| panic!("Cannot open {path}: {e}"));
    let mut r = BufReader::with_capacity(1 << 23, file);

    let mut magic = [0u8; 8];
    r.read_exact(&mut magic).unwrap();
    assert!(&magic[..6] == b"SYNCL2", "Not a valid SYNCL2 index file");

    let k = read_u32le(&mut r) as usize;
    let s = read_u32le(&mut r) as usize;
    let t = read_u32le(&mut r) as usize;

    let version = magic[6];
    let mode = if version >= 5 {
        let mut b = [0u8; 1];
        r.read_exact(&mut b).unwrap();   // seed-mode byte (always syncmer)
        SeedMode::Syncmer
    } else {
        SeedMode::Syncmer
    };
    // v7+: flag byte (bit0 = k-mer atom, bit1 = canonical).  Older = s-mer, fwd.
    if version >= 7 {
        let mut b = [0u8; 1];
        r.read_exact(&mut b).unwrap();
        KMER_ATOM.store(b[0] & 1 != 0, Ordering::Relaxed);
        CANON_ATOM.store(b[0] & 2 != 0, Ordering::Relaxed);
    } else {
        KMER_ATOM.store(false, Ordering::Relaxed);
        CANON_ATOM.store(false, Ordering::Relaxed);
    }

    let num_chrs = read_u32le(&mut r) as usize;
    let mut chr_names = Vec::with_capacity(num_chrs);
    for _ in 0..num_chrs {
        let len = read_u32le(&mut r) as usize;
        let mut buf = vec![0u8; len];
        r.read_exact(&mut buf).unwrap();
        chr_names.push(String::from_utf8(buf).unwrap());
    }

    // Read one level directly into SoA arrays from the packed 13-byte format.
    fn read_level(r: &mut impl Read) -> Level {
        let mut b8 = [0u8; 8];
        r.read_exact(&mut b8).unwrap();
        let num = u64::from_le_bytes(b8) as usize;
        let mut hashes = Vec::with_capacity(num);
        let mut chrs   = Vec::with_capacity(num);
        let mut poss   = Vec::with_capacity(num);
        const CHUNK: usize = 65536;
        let mut buf = vec![0u8; CHUNK * 13];
        let mut remaining = num;
        while remaining > 0 {
            let n = remaining.min(CHUNK);
            r.read_exact(&mut buf[..n * 13]).unwrap();
            for i in 0..n {
                let base = i * 13;
                hashes.push(u64::from_le_bytes(buf[base..base+8].try_into().unwrap()));
                chrs.push(buf[base+8]);
                poss.push(u32::from_le_bytes(buf[base+9..base+13].try_into().unwrap()));
            }
            remaining -= n;
        }
        Level { hashes, chrs, poss }
    }

    let levels: Vec<Level> = if version >= 6 {
        // v6: explicit num_levels + entries in level order (L0, L1, L2, …)
        let num_levels = read_u32le(&mut r) as usize;
        (0..num_levels).map(|_| read_level(&mut r)).collect()
    } else if version >= 4 {
        // v4/v5: no L0 level — prepend an empty L0 slot so existing level
        // numbering (L1 at GIndex.levels[1], etc.) remains consistent.
        let num_levels = read_u32le(&mut r) as usize;
        std::iter::once(Level::default())
            .chain((0..num_levels).map(|_| read_level(&mut r)))
            .collect()
    } else if version == 3 {
        // v3 legacy: L2, L1, L3 — prepend empty L0
        let l2 = read_level(&mut r);
        let l1 = read_level(&mut r);
        let l3 = read_level(&mut r);
        vec![Level::default(), l1, l2, l3]
    } else if version == 2 {
        // v2 legacy: L2, L1 — prepend empty L0
        let l2 = read_level(&mut r);
        let l1 = read_level(&mut r);
        vec![Level::default(), l1, l2]
    } else {
        // v1 legacy: L2 only — prepend empty L0 and empty L1
        let l2 = read_level(&mut r);
        vec![Level::default(), Level::default(), l2]
    };

    (GIndex { levels, chr_names }, k, s, t, mode)
}

// ─────────────────────────────────────────────────────────────────────────────
// 12.  FASTQ reader + reverse complement
// ─────────────────────────────────────────────────────────────────────────────

fn revcomp(seq: &[u8]) -> Vec<u8> {
    seq.iter().rev().map(|&b| match b.to_ascii_uppercase() {
        b'A' => b'T', b'T' => b'A', b'C' => b'G', b'G' => b'C', x => x,
    }).collect()
}

/// Yields raw sequence bytes for every FASTQ record in `path` (gzip or plain).
struct FastqReader {
    inner: Box<dyn BufRead>,
    buf:   String,
}

impl FastqReader {
    fn open(path: &str) -> Self {
        let file = File::open(path)
            .unwrap_or_else(|e| panic!("Cannot open {path}: {e}"));
        let inner: Box<dyn BufRead> = if path.ends_with(".gz") {
            Box::new(BufReader::with_capacity(1 << 20, MultiGzDecoder::new(file)))
        } else {
            Box::new(BufReader::with_capacity(1 << 20, file))
        };
        FastqReader { inner, buf: String::new() }
    }
}

impl Iterator for FastqReader {
    type Item = (String, Vec<u8>);   // (read_name, sequence)
    fn next(&mut self) -> Option<(String, Vec<u8>)> {
        // line 1: @header
        self.buf.clear();
        if self.inner.read_line(&mut self.buf).unwrap_or(0) == 0 { return None; }
        if !self.buf.starts_with('@') { return None; }
        // read name = first whitespace-delimited token after '@'
        let name = self.buf[1..].trim_end()
            .split_whitespace().next().unwrap_or("").to_string();
        // line 2: sequence
        self.buf.clear();
        self.inner.read_line(&mut self.buf).unwrap_or(0);
        let seq: Vec<u8> = self.buf.trim_end().as_bytes().iter().map(|b| b.to_ascii_uppercase()).collect();
        // line 3: + separator
        self.buf.clear();
        self.inner.read_line(&mut self.buf).unwrap_or(0);
        // line 4: quality
        self.buf.clear();
        self.inner.read_line(&mut self.buf).unwrap_or(0);
        Some((name, seq))
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// 13.  Read mapping
// ─────────────────────────────────────────────────────────────────────────────

const MAX_OCC_DEFAULT:    usize = 500;
const MAX_OCC_L1_DEFAULT: usize = 200;

/// Maximum number of hierarchy levels to build (L1 … L_MAX_LEVELS).
/// MEASURED block spans (k=19,s=10): each level grows ~2.27× (≈avg children
/// per block), NOT 3×.  Actual mean spans:
///   L1≈40 bp → L2≈91 → L3≈207 → L4≈478 → L5≈1.1 kb → L6≈2.5 kb.
/// (The old "3×, up to 10 kb" estimate was a stale k=15 guess — overstated ~4×.)
/// Tested 7 levels + smaller k: same accuracy but 2× slower (denser index +
/// the L7 level is only ~1% usable, doesn't short-circuit) → kept at 6.
const MAX_LEVELS: usize = 6;

/// Vote weight for a match at hierarchy level L (0-indexed: L0=1, L1=1, L2=3, L3=9, …).
/// L0 (raw syncmer fallback) shares the same weight as L1 blocks.
/// Each level above L1 is 3× the one below it.
#[inline]
#[allow(dead_code)]
fn vote_weight(level: usize) -> u32 {
    if level == 0 { 1 } else { 3u32.pow((level - 1) as u32) }
}

/// Maximum occurrences allowed at each level.
/// Level is 0-indexed: 0 = L0 raw syncmers, 1 = L1 blocks, 2 = L2, …
///   L0/L1: base_l1   L2/L3: base   L4: base/5   L5: base/20   L6+: base/100
#[inline]
fn max_occ_for_level(level: usize, base: usize, base_l1: usize) -> usize {
    match level {
        0 | 1 => base_l1,
        2 | 3 => base,
        4 => (base / 5).max(10),
        5 => (base / 20).max(5),
        _ => (base / 100).max(2),
    }
}

const CLUSTER_WINDOW: i64  = 500;
const MIN_VOTES:      u32  = 1;

struct MapResult {
    chr:    u8,
    pos:    i64,   // inferred reference start (may be negative near contig ends)
    strand: bool,  // true = forward
    votes:  u32,
    mapq:   u8,    // 0 = ambiguous, 60 = uniquely placed
}

// ── Anchor collection + chain DP ─────────────────────────────────────────────

/// A seed hit: a HierBlock from the query matched at a specific reference position.
/// q_pos..q_pos+span  covers the query interval.
/// r_pos..r_pos+span  is the corresponding reference interval (if colinear).
#[derive(Clone, Copy)]
struct Anchor {
    chr:    u8,
    q_pos:  u32,
    r_pos:  u32,
    span:   u32,
    weight: u32,
}

/// Recursively emit anchors for one HierBlock.
/// Tries the block's own level first (highest weight); falls back to children.
/// `filter` is (chr, lo, hi) in reference-start-offset coords (r_pos - q_pos).
fn emit_anchors(
    block: &HierBlock,
    index: &GIndex,
    max_occ: usize,
    max_occ_l1: usize,
    filter: Option<(u8, i64, i64)>,
    anchors: &mut Vec<Anchor>,
    visited: &mut std::collections::HashSet<(u8, u32, u64)>,
) {
    // Level floor: ignore blocks below MIN_LVL (don't emit, don't recurse deeper).
    if block.level < min_lvl() { return; }

    // LCP rules intentionally share boundary sub-blocks between adjacent parents.
    // A shared child is cloned into both parents, so when both fall through to
    // their children it would otherwise be visited (and its mass emitted) TWICE.
    // Deduplicate by block identity (level, query-pos, hash) → emit once.
    if !visited.insert((block.level as u8, block.pos, block.hash)) { return; }

    // block.level is now 0-indexed: 0 = L0 raw syncmer, 1 = L1 block, …
    let mocc = max_occ_for_level(block.level, max_occ, max_occ_l1);
    let hits = index.lookup(block.level, block.hash, mocc);

    if !hits.is_empty() {
        // Conservation-of-mass weight: a block's vote = its accumulated syncmer
        // mass (every syncmer contributes 1, split across the blocks sharing it).
        // Scaled ×MASS_SCALE so fractional masses keep integer-DP resolution.
        const MASS_SCALE: f32 = 16.0;
        let w = (block.mass * MASS_SCALE).round().max(1.0) as u32;
        let span = block.end - block.pos;
        for (chr_id, genome_pos) in hits.iter() {
            if let Some((fc, flo, fhi)) = filter {
                let vote_pos = genome_pos as i64 - block.pos as i64;
                if chr_id != fc || vote_pos < flo || vote_pos > fhi { continue; }
            }
            anchors.push(Anchor {
                chr:    chr_id,
                q_pos:  block.pos,
                r_pos:  genome_pos,
                span,
                weight: w,
            });
        }
        return;
    }
    // Block not found or too repetitive — try finer-grained children.
    for child in &block.children {
        emit_anchors(child, index, max_occ, max_occ_l1, filter, anchors, visited);
    }
}




/// Minimum distinct anchor positions from the L1+ hierarchy below which
/// we fall back to individual L0 k-mer anchors.
/// Minimum distinct anchor count from L1+ hierarchy below which the L0 raw-syncmer
/// fallback is activated.  Reads with enough L1+ anchors skip the L0 pass entirely.
const L0_FALLBACK_THRESHOLD: usize = 1;

/// Collect anchors for one strand of a read.
///
/// Strategy:
///  1. Hierarchical pass (L1-L6): fast, high-weight anchors from LCP blocks.
///  2. L0 fallback pass: if the hierarchical pass yields fewer than
///     L0_FALLBACK_THRESHOLD distinct genome positions, query individual
///     syncmer k-mer hashes against GIndex.levels[0].  This rescues reads
///     that fall in regions where all LCP blocks are too repetitive or where
///     the read is too short to form blocks.
///     L0 is skipped if the index has no L0 entries (e.g. legacy v4/v5 index).
fn collect_anchors(
    seq: &[u8], index: &GIndex, k: usize, s: usize, t: usize, mode: SeedMode,
    max_occ: usize, max_occ_l1: usize,
    filter_chr: Option<u8>, filter_lo: Option<i64>, filter_hi: Option<i64>,
) -> Vec<Anchor> {
    let filter = match (filter_chr, filter_lo, filter_hi) {
        (Some(c), Some(lo), Some(hi)) => Some((c, lo, hi)),
        _ => None,
    };
    let hier = extract_hier_blocks_n(seq, k, s, t, index.num_levels(), mode);
    let mut anchors = Vec::new();
    let mut visited: std::collections::HashSet<(u8, u32, u64)> = std::collections::HashSet::new();
    for block in &hier {
        emit_anchors(block, index, max_occ, max_occ_l1, filter, &mut anchors, &mut visited);
    }

    // ── L0 fallback pass ─────────────────────────────────────────────────────
    // Only runs when the L1+ hierarchy gives very few anchors AND the index
    // actually has an L0 level (levels[0] non-empty).
    let l0_populated = index.levels.first().map_or(false, |l| !l.is_empty());
    if l0_populated && !no_l0() && min_lvl() == 0 && anchors.len() < L0_FALLBACK_THRESHOLD {
        // Compute k-mer hashes for the read (small — reads are ~15 kbp).
        if let Some(kh_iter) = DnaHashFwd::new(seq, k) {
            let kmer_hashes: Vec<u64> = kh_iter.collect();
            let syncmers = select_seeds_light(seq, k, s, t, mode);
            for sm in &syncmers {
                let pos = sm.pos as usize;
                let Some(&h) = kmer_hashes.get(pos) else { continue };
                if h == u64::MAX { continue }  // N-containing k-mer
                let hits = index.lookup(0, h, max_occ_l1);
                let span = k as u32;
                for (chr_id, genome_pos) in hits.iter() {
                    if let Some((fc, flo, fhi)) = filter {
                        let vote_pos = genome_pos as i64 - sm.pos as i64;
                        if chr_id != fc || vote_pos < flo || vote_pos > fhi { continue; }
                    }
                    anchors.push(Anchor {
                        chr:   chr_id,
                        q_pos: sm.pos,
                        r_pos: genome_pos,
                        span,
                        weight: 16,   // one syncmer = mass 1 × MASS_SCALE
                    });
                }
            }
        }
    }

    prune_anchors(&mut anchors);
    anchors
}

/// Prune anchors for chain DP so O(n²) stays fast.
///
/// Strategy:
///   1. Walk the cascade [L6=243, L5=81, L4=27, L3=9, L2=3, L1=1], scoring
///      (chr, diag-bucket) positions by accumulated anchor weight at each level.
///      Stop at the first level that identifies ≥1 candidate position.
///      For L2/L1 on large anchor sets, sample uniformly (~5000 anchors) to
///      keep per-level cost O(SAMPLE_TARGET) instead of O(n).
///   2. Keep only the top MAX_BUCKETS positions (by total weight).
///   3. Discard any anchor whose ref-offset bucket was not in the top set.
///   4. Within each surviving bucket, keep at most CAP_PER_BUCKET anchors
///      (highest-weight first) to bound total chain DP input.
fn prune_anchors(anchors: &mut Vec<Anchor>) {
    const CAP_PER_BUCKET: usize = 100;
    const MAX_BUCKETS:    usize = 20;
    const SAMPLE_TARGET:  usize = 5_000;

    if anchors.is_empty() { return; }

    // ── Single sort: (chr, diagonal_bucket, weight desc) ─────────────────────
    // Shared by both the cascade scan and the per-bucket cap filter.
    // Sorting once here eliminates the HashMap entirely: buckets are now
    // contiguous runs that a linear sweep can aggregate in O(n).
    anchors.sort_unstable_by_key(|a| {
        let diag = (a.r_pos as i64 - a.q_pos as i64).div_euclid(CLUSTER_WINDOW);
        (a.chr, diag, u32::MAX - a.weight)   // weight desc within each bucket
    });

    // ── Cascade: find top candidate positions via linear sweep ────────────────
    // For each mass threshold (high → low) sum anchor weights per diagonal
    // bucket. Stop at the first threshold that produces ≥1 candidate.
    // Weights are mass×16: a single syncmer ≈ 16; large multi-syncmer blocks
    // reach several hundred.  Thresholds descend through that range so the
    // densest (highest-mass) candidate positions are found first.
    // For low thresholds on large sets, sample every n/SAMPLE_TARGET-th anchor.
    let mut top: Vec<((u8, i64), u32)> = Vec::new();   // (bucket, score)

    'cascade: for &min_w in &[256u32, 128, 64, 32, 16, 1] {
        top.clear();
        let step = if min_w < 32 && anchors.len() > SAMPLE_TARGET {
            anchors.len() / SAMPLE_TARGET
        } else { 1 };

        let mut cur:   (u8, i64) = (255, i64::MAX);
        let mut score: u32       = 0;

        for a in anchors.iter().step_by(step) {
            if a.weight < min_w { continue; }
            let diag   = (a.r_pos as i64 - a.q_pos as i64).div_euclid(CLUSTER_WINDOW);
            let bucket = (a.chr, diag);
            if bucket != cur {
                if score > 0 { top.push((cur, score)); }
                cur   = bucket;
                score = a.weight;
            } else {
                score = score.saturating_add(a.weight);
            }
        }
        if score > 0 { top.push((cur, score)); }
        if !top.is_empty() { break 'cascade; }
    }

    if top.is_empty() {
        // No signal: hard-cap by weight and return.
        anchors.sort_unstable_by(|a, b| b.weight.cmp(&a.weight));
        anchors.truncate(CAP_PER_BUCKET * MAX_BUCKETS);
        return;
    }

    // Keep top MAX_BUCKETS by accumulated score, then sort by bucket key for
    // O(log MAX_BUCKETS) membership test via binary search.
    if top.len() > MAX_BUCKETS {
        top.select_nth_unstable_by(MAX_BUCKETS, |a, b| b.1.cmp(&a.1));
        top.truncate(MAX_BUCKETS);
    }
    top.sort_unstable_by_key(|x| x.0);   // sort by bucket key for bsearch

    // ── Filter + per-bucket cap ───────────────────────────────────────────────
    // Anchors are already sorted by (chr, diag, weight desc) — one pass suffices.
    let mut result: Vec<Anchor> = Vec::with_capacity(MAX_BUCKETS * CAP_PER_BUCKET);
    let mut cur:   (u8, i64) = (255, i64::MAX);
    let mut count: usize     = 0;

    for &a in anchors.iter() {
        let diag   = (a.r_pos as i64 - a.q_pos as i64).div_euclid(CLUSTER_WINDOW);
        let bucket = (a.chr, diag);
        if top.binary_search_by_key(&bucket, |x| x.0).is_err() { continue; }
        if bucket != cur { cur = bucket; count = 0; }
        if count < CAP_PER_BUCKET { result.push(a); count += 1; }
    }

    *anchors = result;
}

// ── Banded alignment (CIGAR generation) ──────────────────────────────────────

use std::cell::RefCell;

struct AlignBufs { dp: Vec<i32>, bt: Vec<u8> }
thread_local! {
    static ALIGN: RefCell<AlignBufs> = RefCell::new(AlignBufs { dp: Vec::new(), bt: Vec::new() });
}

/// Banded semi-global alignment: aligns `query` fully against `target`.
/// `half_band` = maximum allowed diagonal deviation in bases.
/// Scoring: match +2, mismatch −4, linear gap −2.
/// Thread-local DP+BT buffers are reused across calls — no per-read allocation.
/// Returns extended CIGAR (= X I D). Falls back to `{n}M` if band is too narrow.
fn banded_align(query: &[u8], target: &[u8], half_band: usize) -> String {
    let qn = query.len();
    let tn = target.len();
    if qn == 0 { return format!("{}D", tn); }
    if tn == 0 { return format!("{}I", qn); }

    const SMATCH: i32 =  2;
    const SMISS:  i32 = -4;
    // Affine gap (Gotoh): a gap of length L costs GAP_OPEN + L*GAP_EXT.  This
    // models real indels — one long indel is cheap to extend, so the aligner
    // emits a single clean I/D run instead of fragmenting it, and a wrong locus
    // that needs several separate gaps pays a fresh OPEN each time.
    const GAP_OPEN: i32 = -4;
    const GAP_EXT:  i32 = -2;
    const NEG:      i32 = i32::MIN / 4;   // /4 so OPEN+EXT additions never wrap

    let w     = 2 * half_band + 1;
    let cells = (qn + 1) * w;

    // col-offset helper: for row i, column j → index k in [0,w)
    #[inline]
    fn jk(i: usize, j: usize, hb: usize, w: usize) -> Option<usize> {
        let k = j as i64 - i as i64 + hb as i64;
        if k >= 0 && (k as usize) < w { Some(k as usize) } else { None }
    }

    // Three banded matrices: M (match/mismatch end), I (gap in target / query
    // insertion), D (gap in query / deletion).  Back-trace stores, per cell and
    // per matrix, which matrix we came from, packed into one byte per M-cell.
    let ops = ALIGN.with(|cell| {
        let mut bufs = cell.borrow_mut();
        // layout: dp holds 3 interleaved matrices → 3*cells; bt holds M-traceback.
        bufs.dp.resize(3 * cells, NEG);
        bufs.bt.resize(cells, 0u8);
        for v in bufs.dp.iter_mut() { *v = NEG; }
        for v in bufs.bt.iter_mut() { *v = 0;   }
        // index helpers into the three sub-matrices
        let m_off = 0usize; let i_off = cells; let d_off = 2 * cells;

        // ── initialise (0,0) and row 0 (all deletions) ───────────────────────
        bufs.dp[m_off + 0 * w + half_band] = 0;
        for j in 1..=half_band.min(tn) {
            if let Some(k) = jk(0, j, half_band, w) {
                // D run of length j
                bufs.dp[d_off + k] = GAP_OPEN + j as i32 * GAP_EXT;
                bufs.dp[m_off + k] = bufs.dp[d_off + k];
                bufs.bt[k] = 2; // came from D
            }
        }

        for i in 1..=qn {
            let jlo = (i as i64 - half_band as i64).max(0) as usize;
            let jhi = (i + half_band).min(tn);

            if jlo == 0 {
                if let Some(k) = jk(i, 0, half_band, w) {
                    bufs.dp[i_off + i*w + k] = GAP_OPEN + i as i32 * GAP_EXT;
                    bufs.dp[m_off + i*w + k] = bufs.dp[i_off + i*w + k];
                    bufs.bt[i*w + k] = 1; // came from I
                }
            }

            for j in jlo.max(1)..=jhi {
                // I[i][j] = best( M[i-1][j]+OPEN+EXT , I[i-1][j]+EXT )   (query insertion)
                let mut iv = NEG;
                if let Some(k) = jk(i-1, j, half_band, w) {
                    let m = bufs.dp[m_off+(i-1)*w+k];
                    if m != NEG { iv = iv.max(m + GAP_OPEN + GAP_EXT); }
                    let ii = bufs.dp[i_off+(i-1)*w+k];
                    if ii != NEG { iv = iv.max(ii + GAP_EXT); }
                }
                // D[i][j] = best( M[i][j-1]+OPEN+EXT , D[i][j-1]+EXT )   (deletion)
                let mut dv = NEG;
                if let Some(k) = jk(i, j-1, half_band, w) {
                    let m = bufs.dp[m_off+i*w+k];
                    if m != NEG { dv = dv.max(m + GAP_OPEN + GAP_EXT); }
                    let dd = bufs.dp[d_off+i*w+k];
                    if dd != NEG { dv = dv.max(dd + GAP_EXT); }
                }
                // M[i][j] = best( M/I/D[i-1][j-1] ) + match/mismatch
                let mut mv = NEG; let mut from = 0u8;
                if let Some(k) = jk(i-1, j-1, half_band, w) {
                    let sc = if query[i-1].to_ascii_uppercase()
                                == target[j-1].to_ascii_uppercase() { SMATCH } else { SMISS };
                    let pm = bufs.dp[m_off+(i-1)*w+k];
                    let pi = bufs.dp[i_off+(i-1)*w+k];
                    let pd = bufs.dp[d_off+(i-1)*w+k];
                    if pm != NEG && pm+sc > mv { mv = pm+sc; from = 0; }
                    if pi != NEG && pi+sc > mv { mv = pi+sc; from = 1; }
                    if pd != NEG && pd+sc > mv { mv = pd+sc; from = 2; }
                }
                // The M cell can also END on a gap (so traceback can leave via I/D).
                if iv > mv { mv = iv; from = 1; }
                if dv > mv { mv = dv; from = 2; }

                if let Some(k) = jk(i, j, half_band, w) {
                    bufs.dp[i_off+i*w+k] = iv;
                    bufs.dp[d_off+i*w+k] = dv;
                    bufs.dp[m_off+i*w+k] = mv;
                    bufs.bt[i*w+k] = from;
                }
            }
        }

        // ── traceback ─────────────────────────────────────────────────────────
        // Semi-global on reference end: best-scoring M cell in the last query row.
        let jlo_end = (qn as i64 - half_band as i64).max(0) as usize;
        let jhi_end = (qn + half_band).min(tn);
        let j_best = (jlo_end..=jhi_end)
            .filter_map(|j| jk(qn, j, half_band, w).map(|k| (j, bufs.dp[m_off + qn*w + k])))
            .filter(|(_, s)| *s != NEG)
            .max_by_key(|(_, s)| *s)
            .map(|(j, _)| j)
            .unwrap_or(tn.min(qn));

        if jk(qn, j_best, half_band, w).map_or(true, |k| bufs.dp[m_off + qn*w + k] == NEG) {
            return vec![(qn as u32, 4u8)];
        }
        let tn = j_best;

        let mut ops: Vec<(u32, u8)> = Vec::with_capacity(qn / 4);
        let (mut i, mut j) = (qn, tn);
        let mut state = 0u8;  // 0=M, 1=I, 2=D — which matrix we're tracing in
        while i > 0 || j > 0 {
            if i == 0 {            // only deletions can reach the origin
                let n = j as u32; ops.push((n, 3)); break;
            }
            if j == 0 {            // only insertions
                let n = i as u32; ops.push((n, 2)); break;
            }
            let k = match jk(i, j, half_band, w) { Some(k) => k, None => break };
            match state {
                0 => {  // in M: emit =/X, move diagonally, switch to predecessor matrix
                    let from = bufs.bt[i*w + k];
                    if from == 0 {
                        let e = if query[i-1].to_ascii_uppercase()
                                    == target[j-1].to_ascii_uppercase() { 0u8 } else { 1u8 };
                        // emit diagonal then step
                        if ops.last().map_or(false,|x:&(u32,u8)|x.1==e){ops.last_mut().unwrap().0+=1;}
                        else { ops.push((1,e)); }
                        i -= 1; j -= 1;
                    } else {
                        // M ended on a gap → switch to that gap matrix without moving
                        state = from;
                    }
                }
                1 => {  // in I (query insertion): consume query base
                    if ops.last().map_or(false,|x:&(u32,u8)|x.1==2){ops.last_mut().unwrap().0+=1;}
                    else { ops.push((1,2)); }
                    // came from M or I at (i-1, j)?
                    let kk = jk(i-1, j, half_band, w);
                    let stay = kk.map_or(false, |k2| {
                        let ii = bufs.dp[i_off+(i-1)*w+k2];
                        ii != NEG && ii + GAP_EXT == bufs.dp[i_off+i*w+k]
                    });
                    i -= 1;
                    state = if stay { 1 } else { 0 };
                }
                _ => {  // in D (deletion): consume target base
                    if ops.last().map_or(false,|x:&(u32,u8)|x.1==3){ops.last_mut().unwrap().0+=1;}
                    else { ops.push((1,3)); }
                    let kk = jk(i, j-1, half_band, w);
                    let stay = kk.map_or(false, |k2| {
                        let dd = bufs.dp[d_off+i*w+k2];
                        dd != NEG && dd + GAP_EXT == bufs.dp[d_off+i*w+k]
                    });
                    j -= 1;
                    state = if stay { 2 } else { 0 };
                }
            }
        }
        ops.reverse();
        ops
    });

    // encode ops as CIGAR string
    const OP_CHARS: [char; 5] = ['=', 'X', 'I', 'D', 'M']; // M is fallback sentinel
    ops.iter().map(|(n, c)| format!("{}{}", n, OP_CHARS[*c as usize])).collect()
}

/// Given a mapped result and the reference chromosome sequences, extract the
/// reference slice and compute the CIGAR string via banded alignment.
/// `half_band` controls the diagonal band width (good default: max(100, qlen/50)).
/// Reverse CIGAR operation order (used when converting RC-query CIGAR to
/// forward-query CIGAR for PAF strand='-' output).
fn reverse_cigar_ops(cg: &str) -> String {
    let mut ops: Vec<(u32, u8)> = Vec::new();
    let mut n = 0u32;
    for b in cg.bytes() {
        if b.is_ascii_digit() { n = n * 10 + (b - b'0') as u32; }
        else { ops.push((n, b)); n = 0; }
    }
    ops.iter().rev()
        .map(|(n, c)| format!("{}{}", n, *c as char))
        .collect()
}

/// Returns `(cigar_string, actual_ref_end)`.
///
/// Chain-guided path (fast):
///   Re-collects anchors at the mapped position, backtracks the best chain,
///   and aligns only the inter-anchor gaps.  Anchor spans are emitted as `=`.
///   Per-gap band is adaptive (≥ 20, scaled to the gap's indel estimate).
///
/// Fallback (full-read banded DP):
///   Used when the chain produces no anchors (e.g. unmapped region retry,
///   or very short reads with no surviving anchors after re-collection).
fn cigar_for_mapping(
    query_fwd: &[u8],
    result:    &MapResult,
    ref_seqs:  &[Vec<u8>],
    index:     &GIndex,
    k: usize, s: usize, t: usize,
    mode:      SeedMode,
    max_occ: usize, max_occ_l1: usize,
    half_band: usize,            // only used for the full-read fallback
) -> (String, u64) {
    let default_end = result.pos.max(0) as u64 + query_fwd.len() as u64;
    let chr = result.chr as usize;
    let Some(chr_seq) = ref_seqs.get(chr) else {
        return (format!("{}M", query_fwd.len()), default_end);
    };

    // Strand-adjusted query: RC(query) for reverse-strand mappings so that
    // the chain anchors are expressed in the forward-reference frame.
    let seq: Vec<u8> = if result.strand {
        query_fwd.to_vec()
    } else {
        revcomp(query_fwd)
    };

    let ref_start = result.pos.max(0) as usize;

    // ── Chain-guided CIGAR ────────────────────────────────────────────────────
    // Re-collect anchors restricted to a ±500 bp window around the mapped
    // position so the chain DP is cheap and produces clean anchors.
    let slack = 500i64;
    let mut anchors = collect_anchors(
        &seq, index, k, s, t, mode, max_occ, max_occ_l1,
        Some(result.chr),
        Some(result.pos - slack),
        Some(result.pos + slack),
    );
    anchors.sort_unstable_by_key(|a| (a.chr, a.q_pos));

    if let Some((_, _, _, chain)) = single_chain_with_trace(&anchors) {
        if !chain.is_empty() {
            let (cg_raw, ref_end) = cigar_from_chain_anchors(&seq, &chain, chr_seq, ref_start);
            // For reverse strand: reverse the CIGAR ops to convert from
            // "RC(query) vs forward ref" to "query vs RC(ref)" (PAF convention).
            let cg = if result.strand { cg_raw } else { reverse_cigar_ops(&cg_raw) };
            return (cg, ref_end);
        }
    }

    // ── Fallback: full-read banded DP ─────────────────────────────────────────
    let ref_end_s = (ref_start + seq.len() + half_band).min(chr_seq.len());
    let ref_slice = &chr_seq[ref_start..ref_end_s];
    let cg_raw = banded_align(&seq, ref_slice, half_band);
    let ref_consumed = ref_bases_in_cigar(&cg_raw) as u64;
    let cg = if result.strand { cg_raw } else { reverse_cigar_ops(&cg_raw) };
    (cg, ref_start as u64 + ref_consumed)
}

// ── Chain DP constants ────────────────────────────────────────────────────────

/// Max query gap between consecutive anchors in a chain (bases).
const CHAIN_MAX_GAP: u32 = 5_000;
/// Allowed absolute difference |dq - dr| before a chain extension is rejected.
/// Also admits 10% of max(dq,dr) for proportional tolerance.
const CHAIN_GAP_TOL: u32 = 150;
/// Penalty per base of gap inconsistency |dq - dr|.
const CHAIN_GAP_SCALE: u32 = 1;

/// Single forward chain-DP pass on a slice of anchors pre-sorted by (chr, q_pos).
/// Returns (chr, ref_start_offset, best_score) or None if nothing exceeds MIN_VOTES.
/// ref_start_offset = r_pos - q_pos at the chain's first anchor.
fn single_chain(anchors: &[Anchor]) -> Option<(u8, i64, u32)> {
    let n = anchors.len();
    if n == 0 { return None; }

    let mut dp   = vec![0u32; n];
    let mut from = vec![usize::MAX; n];

    for j in 0..n {
        let aj = anchors[j];
        dp[j] = aj.weight;

        // Walk backwards; once we exit the MAX_CHAIN_GAP window, stop
        let mut i = j;
        while i > 0 {
            i -= 1;
            let ai = anchors[i];
            if ai.chr != aj.chr { break; }

            // q-gap from end of ai to start of aj
            let q_end_i = ai.q_pos.saturating_add(ai.span);
            if aj.q_pos < q_end_i { continue; }         // overlapping in query
            let dq = aj.q_pos - q_end_i;
            if dq > CHAIN_MAX_GAP { break; }             // too far back

            // r-gap from end of ai to start of aj
            if aj.r_pos < ai.r_pos.saturating_add(ai.span) { continue; } // ref goes backward
            let dr = aj.r_pos - ai.r_pos.saturating_add(ai.span);
            if dr > CHAIN_MAX_GAP + CHAIN_GAP_TOL { continue; }

            // Colinearity: |dq - dr| must be small
            let gap_diff = if dq > dr { dq - dr } else { dr - dq };
            let tol = CHAIN_GAP_TOL + dq / 10;                  // 10% proportional slack
            if gap_diff > tol { continue; }

            let penalty = gap_diff.saturating_mul(CHAIN_GAP_SCALE);
            let score   = dp[i].saturating_add(aj.weight).saturating_sub(penalty);
            if score > dp[j] {
                dp[j]   = score;
                from[j] = i;
            }
        }
    }

    // Best endpoint
    let best_j = (0..n).max_by_key(|&i| dp[i]).unwrap();
    if dp[best_j] < MIN_VOTES { return None; }

    // Backtrack to chain start for accurate ref_start
    let mut j = best_j;
    while from[j] != usize::MAX { j = from[j]; }
    let first = anchors[j];
    let ref_offset = first.r_pos as i64 - first.q_pos as i64;
    Some((anchors[best_j].chr, ref_offset, dp[best_j]))
}

/// Like single_chain but also backtracks and returns the chain anchors in query order.
/// Used by chain-guided CIGAR generation.
fn single_chain_with_trace(anchors: &[Anchor]) -> Option<(u8, i64, u32, Vec<Anchor>)> {
    let n = anchors.len();
    if n == 0 { return None; }

    let mut dp   = vec![0u32; n];
    let mut from = vec![usize::MAX; n];

    for j in 0..n {
        let aj = anchors[j];
        dp[j] = aj.weight;
        let mut i = j;
        while i > 0 {
            i -= 1;
            let ai = anchors[i];
            if ai.chr != aj.chr { break; }
            let q_end_i = ai.q_pos.saturating_add(ai.span);
            if aj.q_pos < q_end_i { continue; }
            let dq = aj.q_pos - q_end_i;
            if dq > CHAIN_MAX_GAP { break; }
            if aj.r_pos < ai.r_pos.saturating_add(ai.span) { continue; }
            let dr = aj.r_pos - ai.r_pos.saturating_add(ai.span);
            if dr > CHAIN_MAX_GAP + CHAIN_GAP_TOL { continue; }
            let gap_diff = if dq > dr { dq - dr } else { dr - dq };
            let tol = CHAIN_GAP_TOL + dq / 10;
            if gap_diff > tol { continue; }
            let penalty = gap_diff.saturating_mul(CHAIN_GAP_SCALE);
            let score   = dp[i].saturating_add(aj.weight).saturating_sub(penalty);
            if score > dp[j] { dp[j] = score; from[j] = i; }
        }
    }

    let best_j = (0..n).max_by_key(|&i| dp[i]).unwrap();
    if dp[best_j] < MIN_VOTES { return None; }

    // Backtrack full chain
    let mut chain: Vec<Anchor> = Vec::new();
    let mut j = best_j;
    loop {
        chain.push(anchors[j]);
        if from[j] == usize::MAX { break; }
        j = from[j];
    }
    chain.reverse();  // now sorted by q_pos ascending

    let first = chain[0];
    let ref_offset = first.r_pos as i64 - first.q_pos as i64;
    Some((anchors[best_j].chr, ref_offset, dp[best_j], chain))
}

/// Count reference bases consumed by a CIGAR string (=, X, D, M operators).
#[inline]
fn ref_bases_in_cigar(cg: &str) -> usize {
    let mut total = 0usize;
    let mut n = 0usize;
    for b in cg.bytes() {
        if b.is_ascii_digit() { n = n * 10 + (b - b'0') as usize; }
        else { if matches!(b, b'=' | b'X' | b'D' | b'M') { total += n; } n = 0; }
    }
    total
}

/// Build a CIGAR string using chain anchors as alignment seeds.
///
/// Anchor spans are emitted as exact-match ops (`=`): LCP block hash matches
/// guarantee the same s-mer sequences, making the anchor region near-identical.
/// Only the inter-anchor gaps (and the leading / trailing unanchored regions)
/// are aligned with banded DP using an adaptive per-gap band.
///
/// `seq`       — strand-adjusted query (RC(query) if reverse-strand mapping)
/// `chain`     — chain anchors sorted by q_pos (output of single_chain_with_trace)
/// `chr_seq`   — full chromosome sequence
/// `ref_start` — absolute ref position for seq[0], i.e. result.pos.max(0)
fn cigar_from_chain_anchors(
    seq:       &[u8],
    chain:     &[Anchor],
    chr_seq:   &[u8],
    ref_start: usize,
) -> (String, u64) {
    const MIN_BAND: usize = 20;
    let qlen = seq.len();
    let mut cg = String::with_capacity(qlen / 3);

    let mut q_cur: usize = 0;
    let mut r_cur: usize = ref_start;

    /// Align a small query/reference slice with an adaptive band.
    fn align_gap(q: &[u8], r: &[u8]) -> String {
        if q.is_empty() && r.is_empty() { return String::new(); }
        if q.is_empty() { return format!("{}D", r.len()); }
        if r.is_empty() { return format!("{}I", q.len()); }
        let net = (q.len() as i64 - r.len() as i64).unsigned_abs() as usize;
        banded_align(q, r, (20usize).max(net + 10))
    }

    for anchor in chain {
        let aq = anchor.q_pos as usize;
        let ar = anchor.r_pos as usize;
        let sp = anchor.span as usize;

        // ── gap before this anchor ────────────────────────────────────────────
        if aq > q_cur || ar > r_cur {
            let q_gap = &seq[q_cur..aq.min(qlen)];
            let r_gap_end = ar.min(chr_seq.len());
            let r_gap = if r_cur < r_gap_end { &chr_seq[r_cur..r_gap_end] } else { &[] };
            let gc = align_gap(q_gap, r_gap);
            // r_cur is re-anchored to the anchor's absolute ref coord below, so
            // the gap's ref-base advance does not need to be tracked here.
            cg.push_str(&gc);
        }

        // ── anchor span: emit as exact match ─────────────────────────────────
        // (q_cur/r_cur are re-anchored to absolute anchor coords below.)
        let q_end = (aq + sp).min(qlen);
        let r_end = (ar + sp).min(chr_seq.len());
        let actual_sp = q_end.saturating_sub(aq).min(r_end.saturating_sub(ar));
        if actual_sp > 0 {
            cg.push_str(&format!("{}=", actual_sp));
        }
        q_cur = aq + actual_sp;
        r_cur = ar + actual_sp;
    }

    // ── trailing region ───────────────────────────────────────────────────────
    if q_cur < qlen {
        let q_tail = &seq[q_cur..qlen];
        let r_tail_end = (r_cur + q_tail.len() + MIN_BAND * 2).min(chr_seq.len());
        let r_tail = if r_cur < r_tail_end { &chr_seq[r_cur..r_tail_end] } else { &[] };
        let gc = align_gap(q_tail, r_tail);
        r_cur += ref_bases_in_cigar(&gc);
        cg.push_str(&gc);
    }

    (cg, r_cur as u64)
}

fn chain_dp(anchors: &mut Vec<Anchor>) -> (Option<(u8, i64, u32)>, u32) {
    if anchors.is_empty() { return (None, 0); }
    anchors.sort_unstable_by_key(|a| (a.chr, a.q_pos));

    let best = single_chain(anchors);
    let Some((bc, bo, _)) = best else { return (None, 0); };

    // Second chain: exclude anchors whose ref_start_offset is near the best chain's
    let rest: Vec<Anchor> = anchors.iter().copied().filter(|a| {
        a.chr != bc || (a.r_pos as i64 - a.q_pos as i64 - bo).abs() > CLUSTER_WINDOW
    }).collect();
    let second_score = single_chain(&rest).map_or(0, |c| c.2);

    (best, second_score)
}

/// Inner mapping pass for one set of occ thresholds.
///
/// MAPQ correctness note: a read from (say) chr1:P forward strand will also
/// produce a near-identical chain on the RC strand at the same genomic locus.
/// Counting that as a competing "second-best" chain would collapse MAPQ to ~0
/// for every uniquely-mapped read.  We avoid this by only promoting a chain
/// to second_score when it maps to a *genuinely different* locus
/// (different chromosome, or offset differing by more than CLUSTER_WINDOW).
fn map_read_with_occ(fwd: &[u8], rc: &[u8], index: &GIndex,
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
        let (top, second_here) = chain_dp(&mut anchors);
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
fn map_read(fwd: &[u8], index: &GIndex, k: usize, s: usize, t: usize, mode: SeedMode,
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
struct MapStats {
    total:    u64,
    mapped:   u64,
    bases:    u64,
}

fn write_paf_line(
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

fn run_mapping(reads_path: &str, genome_or_idx: &str, paf_out: Option<&str>,
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
        println!("  Building index from {} ...", genome_or_idx);
        let (idx, elapsed) = build_index(genome_or_idx, k, s, t, MAX_LEVELS, mode);
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
  --min-level N    lowest block level allowed to anchor a read (default 3).
                   Default 3 is tuned for HiFi (≤~0.5% error): higher accuracy,
                   near-zero wrong-chromosome, ~50% faster. Use 0 for noisy
                   (>1% error) reads, which need the finer-block fallback.
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

    // --build-index <genome.fa[.gz]> <output.idx>
    if let Some(bi_pos) = args.iter().position(|a| a == "--build-index") {
        let genome_path = args.get(bi_pos + 1)
            .expect("--build-index requires a genome path");
        let idx_path = args[bi_pos + 2..].iter()
            .find(|a| !a.starts_with('-'))
            .map(|s| s.as_str()).unwrap_or("syncmer-hifi.idx");
        let nthreads = rayon::current_num_threads();
        println!("\n  Building L{MAX_LEVELS} index  k={k}  s={s}  threads={nthreads}");
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
