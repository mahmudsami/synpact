#![allow(unused_imports, dead_code)]
use crate::config::*;
use crate::hash::*;
use crate::syncmer::*;
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

#[derive(Debug)]
pub(crate) struct Block {
    pub(crate) indices: Vec<usize>,
    pub(crate) rule:    &'static str,
}

/// Cole–Vishkin deterministic-coin-tossing tag between two consecutive values.
/// = 2 * (index of lowest bit where they differ) + (that bit's value in `cur`).
/// Tags vary non-monotonically even over a sorted run, so their local maxima
/// give content-determined (locally consistent) split points.
#[inline]
pub(crate) fn dct_tag(prev: u64, cur: u64) -> u32 {
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
pub(crate) fn split_monotone_run(run: &[usize], values: &[u64]) -> Vec<Vec<usize>> {
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
pub(crate) fn is_local_min(v: &[u64], i: usize) -> bool {
    let n = v.len();
    (i == 0 || v[i] < v[i - 1]) && (i == n - 1 || v[i] < v[i + 1])
}

#[inline]
pub(crate) fn is_local_max(v: &[u64], i: usize) -> bool {
    let n = v.len();
    (i == 0 || v[i] > v[i - 1]) && (i == n - 1 || v[i] > v[i + 1])
}

pub(crate) fn locally_consistent_parsing(values: &[u64]) -> Vec<Block> {
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

/// Extract all hierarchy levels from one N-free segment.
/// Returns result[0] = L1 anchors, result[1] = L2, … up to max_levels.
/// Used during index construction.
pub(crate) fn extract_all_levels(seq: &[u8], k: usize, s: usize, t: usize, max_levels: usize,
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

/// A block at any hierarchy level. `level` is 1-indexed: 1 = L1 (leaf), 2 = L2, …
/// `pos`..`end` is the half-open query interval covered by this block (in bp).
///
/// `children` holds **indices into the level below** (`forest.levels[level-1]`),
/// not owned subtrees — LCP shares boundary sub-blocks between adjacent parents,
/// so an owned-children tree had to deep-clone shared subtrees at every level
/// (~32 % of map time). Index references store each block exactly once.
pub(crate) struct HierNode {
    pub(crate) level:    usize,
    pub(crate) hash:     u64,
    pub(crate) pos:      u32,   // query start (inclusive)
    pub(crate) end:      u32,   // query end   (exclusive); span = end - pos
    /// Conservation-of-mass weight: every syncmer carries mass 1; a block's mass
    /// is the sum of its units' masses, with a unit shared by `m` blocks
    /// contributing 1/m to each.  Total mass is conserved across all levels.
    pub(crate) mass:     f32,
    pub(crate) children: Vec<u32>,  // indices into the level below
}

/// The full per-read block hierarchy, stored flat: `levels[L]` is every block at
/// level L (1-indexed; `levels[0]` is an unused placeholder). The top-level blocks
/// to anchor from are `levels[top_level]`.
pub(crate) struct HierForest {
    pub(crate) levels:    Vec<Vec<HierNode>>,
    pub(crate) top_level: usize,
}

/// Build the per-read block hierarchy as a flat `HierForest`.
/// num_levels is taken from the loaded index so we don't exceed what was indexed.
pub(crate) fn extract_hier_blocks_n(seq: &[u8], k: usize, s: usize, t: usize, num_levels: usize,
                         mode: SeedMode) -> HierForest
{
    let empty = || HierForest { levels: vec![Vec::new()], top_level: 0 };

    let syncmers = timed(&PROF_SEED, || select_seeds_light(seq, k, s, t, mode));
    if syncmers.is_empty() { return empty(); }

    let smer_vals: Vec<u64> = syncmers.iter().map(|sm| sm.value).collect();
    let l1_raw = timed(&PROF_L1, || locally_consistent_parsing(&smer_vals));

    // ── L1 masses ────────────────────────────────────────────────────────────
    // Each syncmer carries mass 1.  A syncmer shared by `m` L1 blocks splits its
    // mass 1/m across them.  A block's mass = sum of its syncmers' shares.
    let mut sync_membership = vec![0u32; syncmers.len()];
    for blk in &l1_raw { for &i in &blk.indices { sync_membership[i] += 1; } }

    let l1: Vec<HierNode> = l1_raw.iter().map(|blk| {
        let bvals: Vec<u64> = blk.indices.iter().map(|&i| smer_vals[i]).collect();
        let h   = block_hash_for_level(&bvals, 0);
        let pos = syncmers[blk.indices[0]].pos;
        let end = syncmers[*blk.indices.last().unwrap()].pos + k as u32;
        let mass: f32 = blk.indices.iter()
            .map(|&i| 1.0 / sync_membership[i] as f32).sum();
        HierNode { level: 1, hash: h, pos, end, mass, children: Vec::new() }
    }).collect();

    if l1.is_empty() { return empty(); }

    // levels[0] is an unused placeholder so level numbers index directly.
    let mut levels: Vec<Vec<HierNode>> = vec![Vec::new(), l1];
    let mut top_level = 1usize;

    // Iteratively build L2, L3, … up to num_levels. Each new node stores its
    // children as indices into the level below — no subtree cloning.
    timed(&PROF_UPPER, || {
        for level_1idx in 2..=num_levels {
            let cur = &levels[level_1idx - 1];
            if cur.len() < 2 { break; }

            // Snapshot the fields we need so we can push to `levels` afterwards.
            let cur_hashes: Vec<u64> = cur.iter().map(|b| b.hash).collect();
            let cur_pos:    Vec<u32> = cur.iter().map(|b| b.pos).collect();
            let cur_end:    Vec<u32> = cur.iter().map(|b| b.end).collect();
            let cur_mass:   Vec<f32> = cur.iter().map(|b| b.mass).collect();
            let next_raw = locally_consistent_parsing(&cur_hashes);
            if next_raw.is_empty() { break; }

            // Membership: how many parent blocks each child block belongs to.
            let mut child_membership = vec![0u32; cur_hashes.len()];
            for blk in &next_raw { for &i in &blk.indices { child_membership[i] += 1; } }

            let new_level: Vec<HierNode> = next_raw.iter().map(|blk| {
                let bvals: Vec<u64> = blk.indices.iter().map(|&i| cur_hashes[i]).collect();
                let h        = block_hash_for_level(&bvals, level_1idx - 1);
                let pos      = cur_pos[blk.indices[0]];
                let end      = cur_end[*blk.indices.last().unwrap()];
                // Parent mass = sum of each child's mass / (#parents sharing that child).
                let mass: f32 = blk.indices.iter()
                    .map(|&i| cur_mass[i] / child_membership[i] as f32).sum();
                let children: Vec<u32> = blk.indices.iter().map(|&i| i as u32).collect();
                HierNode { level: level_1idx, hash: h, pos, end, mass, children }
            }).collect();

            levels.push(new_level);
            top_level = level_1idx;
        }
    });

    HierForest { levels, top_level }
}

// ─────────────────────────────────────────────────────────────────────────────
// 11.  Genome index (sorted flat array, binary-searched at query time)
// ─────────────────────────────────────────────────────────────────────────────
