use crate::config::*;
use crate::hash::*;
use crate::syncmer::*;

/// One LCP block: the contiguous half-open range of unit indices `start..end`.
///
/// Every rule in `locally_consistent_parsing` emits a contiguous `lo..=hi` range,
/// so a block never needs an explicit index list — storing it as two `u32`s
/// instead of a `Vec<usize>` removes millions of tiny heap allocations at genome
/// scale (the build's dominant per-thread scratch) and is also faster.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) struct Block {
    pub(crate) start: u32,
    pub(crate) end:   u32,   // exclusive
}

impl Block {
    #[inline] pub(crate) fn range(&self) -> std::ops::Range<usize> {
        self.start as usize .. self.end as usize
    }
    #[inline] pub(crate) fn first(&self) -> usize { self.start as usize }
    #[inline] pub(crate) fn last(&self)  -> usize { self.end as usize - 1 }
}

/// Cole–Vishkin deterministic coin-tossing: reduce a *proper* colouring (adjacent
/// entries distinct) of `vals` to a proper **3-colouring** over {0,1,2}, following
/// "Pearls of Algorithm Engineering", §4.3.2 (List Ranking → Deterministic
/// Coin-Tossing).  A strictly-monotone run is automatically a proper colouring.
///
/// Two phases:
///   1. *Get six colours.*  Repeat, until every colour `< 6`:
///        coin'(i) = 2·π(i) + z(i),
///      where π(i) is the lowest bit position at which coin(i) and its successor
///      coin(i+1) differ, and z(i) is that bit of coin(i).  Each round keeps the
///      colouring proper and shrinks the alphabet logarithmically (u64 → <6 in a
///      handful of rounds).  The last element has no successor, so it is given any
///      colour distinct from its (new) predecessor.
///   2. *Get three colours.*  For v = 3,4,5 in turn, recolour every item of
///      colour v to the smallest colour in {0,1,2} differing from both neighbours.
///      No two adjacent items share colour v (the colouring is proper), so this is
///      well-defined and keeps the colouring proper.
///
/// Each output colour depends only on a bounded window of `vals` (≈ log* rounds,
/// successor-directed), so the colouring — and the split derived from it — is
/// locally consistent: a single changed value perturbs only nearby boundaries.
fn dct_three_coloring(vals: &[u64]) -> Vec<u8> {
    let m = vals.len();
    debug_assert!(m >= 2);
    let mut coin: Vec<u64> = vals.to_vec();

    // ── Phase 1: reduce to ≤ 6 colours ───────────────────────────────────────
    while coin.iter().copied().max().unwrap_or(0) >= 6 {
        let mut next = vec![0u64; m];
        for i in 0..m - 1 {
            let p = (coin[i] ^ coin[i + 1]).trailing_zeros() as u64; // first differing bit
            next[i] = 2 * p + ((coin[i] >> p) & 1);
        }
        // Tail (no successor): any colour ≠ its new predecessor keeps it proper.
        next[m - 1] = if next[m - 2] == 0 { 1 } else { 0 };
        coin = next;
    }

    // ── Phase 2: collapse {3,4,5} → {0,1,2}, preserving adjacent-distinctness ──
    let mut col: Vec<u8> = coin.iter().map(|&c| c as u8).collect();
    for v in [3u8, 4, 5] {
        for i in 0..m {
            if col[i] != v { continue; }
            let l = if i > 0     { col[i - 1] } else { u8::MAX };
            let r = if i + 1 < m { col[i + 1] } else { u8::MAX };
            col[i] = (0u8..=2).find(|&c| c != l && c != r).unwrap();
        }
    }
    col
}

/// Split a strictly-monotone run of indices (whose `values` are strictly
/// monotone, hence repetition-free) into blocks of size ≤ 3, deterministically
/// and locally-consistently.
///
/// The run is reduced to a proper 3-colouring by deterministic coin-tossing
/// (`dct_three_coloring`), then parsed by the **same local-minimum and
/// local-maximum rules used by the rest of the LCP** (rules 1 and 2): a triplet
/// `{i-1, i, i+1}` around every local minimum, and around every local maximum
/// whose neighbours are not minima.  This keeps the monotone-run blocks identical
/// in form to all other blocks (centred triplets with shared boundaries / the
/// same intentional-intersection and conservation-of-mass semantics).
///
/// Completeness: over the alphabet {0,1,2} no two slope points are adjacent (four
/// strictly-monotone colours cannot fit in three), so every slope point neighbours
/// a local maximum and every position is covered by some min/max triplet — the
/// repetition (rule 3) and monotone (rule 4) cases are unreachable on a 3-colour
/// proper sequence.  Every emitted block therefore spans ≤ 3 indices; no chunking.
/// Returns contiguous half-open index ranges `start..end` (in the same index space
/// as `run`); the run itself is contiguous, so every triplet is too.
pub(crate) fn split_monotone_run(run: &[usize], values: &[u64]) -> Vec<(u32, u32)> {
    let m = run.len();
    if m <= 3 { return vec![(run[0] as u32, run[m - 1] as u32 + 1)]; }

    // Colour the run by its *content* (the values) so read and reference split
    // identically wherever the same run occurs.
    let vals: Vec<u64> = run.iter().map(|&i| values[i]).collect();
    let col: Vec<u64> = dct_three_coloring(&vals).iter().map(|&c| c as u64).collect();

    let triplet = |i: usize| -> (u32, u32) {
        let lo = i.saturating_sub(1);
        let hi = (i + 1).min(m - 1);
        (run[lo] as u32, run[hi] as u32 + 1)
    };

    let mut out: Vec<(u32, u32)> = Vec::new();
    // Rule 1 — local-minimum triplet.
    for i in 0..m {
        if is_local_min(&col, i) { out.push(triplet(i)); }
    }
    // Rule 2 — local-maximum triplet, when no adjacent local minimum (that min's
    // own triplet already covers the maximum).
    for i in 0..m {
        if is_local_max(&col, i) {
            let l_min = i > 0     && is_local_min(&col, i - 1);
            let r_min = i + 1 < m && is_local_min(&col, i + 1);
            if !l_min && !r_min { out.push(triplet(i)); }
        }
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
    if n == 1 { return vec![Block { start: 0, end: 1 }]; }

    // Vec<bool> instead of HashSet — better cache behaviour at genome scale
    let is_min: Vec<bool> = (0..n).map(|i| is_local_min(values, i)).collect();
    let is_max: Vec<bool> = (0..n).map(|i| is_local_max(values, i)).collect();

    let mut assigned = vec![false; n];
    let mut blocks: Vec<Block> = Vec::new();

    // Commit the contiguous inclusive range `lo..=hi` (clamped to the sequence) as
    // one block, but only if it contributes at least one not-yet-assigned position
    // (shared boundaries between adjacent blocks are intentional intersections).
    macro_rules! commit_range {
        ($lo:expr, $hi:expr) => {{
            let lo = $lo;
            let hi = ($hi).min(n - 1);
            if (lo..=hi).any(|i| !assigned[i]) {
                for i in lo..=hi { assigned[i] = true; }
                blocks.push(Block { start: lo as u32, end: hi as u32 + 1 });
            }
        }};
    }

    // Rule 1 — local-minimum block: always full triplet {i-1, i, i+1}.
    // Adjacent minima each get their full triplet (shared boundary appears in
    // both — intentional intersection).
    for i in 0..n {
        if is_min[i] {
            commit_range!(i.saturating_sub(1), i + 1);
        }
    }

    // Rule 2 — local-maximum block (no adjacent local minimum), {i-1, i, i+1}.
    for i in 0..n {
        if is_max[i] {
            let l_is_min = i > 0     && is_min[i - 1];
            let r_is_min = i + 1 < n && is_min[i + 1];
            if !l_is_min && !r_is_min {
                commit_range!(i.saturating_sub(1), i + 1);
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
                // j is one past the run; the macro clamps it to the next syncmer.
                commit_range!(i.saturating_sub(1), j);
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
                // Split long monotone runs into ≤3-unit blocks via DCT so they
                // aren't one big featureless block (deterministic & locally
                // consistent → genome and read split identically).
                for (s, e) in split_monotone_run(&run, values) {
                    commit_range!(s as usize, e as usize - 1);
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

    blocks.sort_by_key(|b| b.start);
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
    // Only levels ≥ min_level are stored in the index. Lower levels must still be
    // *computed* (each level is parsed from the one below), but their (hash, pos)
    // entry tuples are never emitted — building them only to discard them wastes
    // the build's largest entry vec (L1 has the most blocks).
    let min_level = min_lvl();

    let syncmers = select_seeds_light(seq, k, s, t, mode);
    if syncmers.is_empty() { return vec![]; }

    let smer_vals: Vec<u64> = syncmers.iter().map(|sm| sm.value).collect();
    let l1_raw = locally_consistent_parsing(&smer_vals);

    let mut cur_hashes: Vec<u64> = Vec::with_capacity(l1_raw.len());
    let mut cur_pos:    Vec<u32> = Vec::with_capacity(l1_raw.len());
    let store_l1 = 1 >= min_level;
    let mut first_out: Vec<(u64, u32)> =
        if store_l1 { Vec::with_capacity(l1_raw.len()) } else { Vec::new() };

    for blk in &l1_raw {
        let h   = block_hash_for_level(&smer_vals[blk.range()], 0);
        let pos = syncmers[blk.first()].pos;
        cur_hashes.push(h);
        cur_pos.push(pos);
        if store_l1 { first_out.push((h, pos)); }
    }
    // The syncmer atoms, their s-mer values and the L1 block index lists are dead
    // once the L1 block hashes/positions exist. These are the heaviest per-segment
    // allocations (millions of tiny `Vec<usize>` in `l1_raw`); free them before
    // building L2…L6 so they don't stay resident through the upper levels and, more
    // importantly, don't stack up across worker threads into the build's peak.
    drop(l1_raw);
    drop(smer_vals);
    drop(syncmers);

    // all[0] = L1 blocks, all[1] = L2 blocks, … (empty for levels below min_level).
    let mut all: Vec<Vec<(u64, u32)>> = vec![first_out];

    for level_1idx in 2..=max_levels {
        if cur_hashes.len() < 2 { break; }
        let next_raw = locally_consistent_parsing(&cur_hashes);
        if next_raw.is_empty() { break; }

        let mut next_hashes = Vec::with_capacity(next_raw.len());
        let mut next_pos    = Vec::with_capacity(next_raw.len());
        let store = level_1idx >= min_level;
        let mut level_out: Vec<(u64, u32)> =
            if store { Vec::with_capacity(next_raw.len()) } else { Vec::new() };

        for blk in &next_raw {
            let h   = block_hash_for_level(&cur_hashes[blk.range()], level_1idx - 1);
            let pos = cur_pos[blk.first()];
            next_hashes.push(h);
            next_pos.push(pos);
            if store { level_out.push((h, pos)); }
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
    for blk in &l1_raw { for i in blk.range() { sync_membership[i] += 1; } }

    let l1: Vec<HierNode> = l1_raw.iter().map(|blk| {
        let h   = block_hash_for_level(&smer_vals[blk.range()], 0);
        let pos = syncmers[blk.first()].pos;
        let end = syncmers[blk.last()].pos + k as u32;
        let mass: f32 = blk.range()
            .map(|i| 1.0 / sync_membership[i] as f32).sum();
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
            for blk in &next_raw { for i in blk.range() { child_membership[i] += 1; } }

            let new_level: Vec<HierNode> = next_raw.iter().map(|blk| {
                let h        = block_hash_for_level(&cur_hashes[blk.range()], level_1idx - 1);
                let pos      = cur_pos[blk.first()];
                let end      = cur_end[blk.last()];
                // Parent mass = sum of each child's mass / (#parents sharing that child).
                let mass: f32 = blk.range()
                    .map(|i| cur_mass[i] / child_membership[i] as f32).sum();
                let children: Vec<u32> = (blk.start..blk.end).collect();
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

#[cfg(test)]
mod dct_tests {
    use super::{dct_three_coloring, split_monotone_run};

    // Tiny deterministic xorshift RNG (no external dep).
    struct Rng(u64);
    impl Rng {
        fn next(&mut self) -> u64 {
            let mut x = self.0;
            x ^= x << 13; x ^= x >> 7; x ^= x << 17;
            self.0 = x; x
        }
    }

    /// dct_three_coloring must return a proper 3-colouring of any
    /// adjacent-distinct sequence.
    #[test]
    fn three_coloring_is_proper() {
        let mut rng = Rng(0x9E3779B97F4A7C15);
        for _ in 0..20_000 {
            let m = 2 + (rng.next() % 180) as usize;
            // Build an adjacent-distinct sequence (proper colouring).
            let mut vals = Vec::with_capacity(m);
            let mut last = u64::MAX;
            for _ in 0..m {
                let mut v = rng.next();
                if v == last { v ^= 1; }        // force distinct from neighbour
                vals.push(v);
                last = v;
            }
            let col = dct_three_coloring(&vals);
            assert_eq!(col.len(), m);
            assert!(col.iter().all(|&c| c <= 2), "colour out of {{0,1,2}}");
            for i in 0..m - 1 {
                assert_ne!(col[i], col[i + 1], "colouring not proper at {i}: {col:?}");
            }
        }
    }

    /// split_monotone_run: every block ≤ 3 and consecutive, the union of blocks
    /// covers the whole run (rule-1/2 triplets, so blocks may overlap rather than
    /// stitch), and it is deterministic — for strictly increasing and decreasing
    /// runs.
    #[test]
    fn monotone_split_blocks_le_3_and_cover() {
        let mut rng = Rng(0x1234_5678_9ABC_DEF0);
        for _ in 0..20_000 {
            let m = 1 + (rng.next() % 200) as usize;
            let decreasing = rng.next() & 1 == 0;

            // Strictly monotone values; run indices are 0..m (identity).
            let run: Vec<usize> = (0..m).collect();
            let mut values = Vec::with_capacity(m);
            let mut acc: u64 = (rng.next() % 1000) + 1;
            for _ in 0..m {
                values.push(acc);
                acc += 1 + (rng.next() % 50); // strictly increasing steps
            }
            if decreasing { values.reverse(); } // strictly decreasing

            let blocks = split_monotone_run(&run, &values);
            assert!(!blocks.is_empty());

            // Determinism.
            assert_eq!(blocks, split_monotone_run(&run, &values));

            // Size bound: each block is a contiguous range of 1..=3 indices.
            for &(s, e) in &blocks {
                let len = e - s;
                assert!((1..=3).contains(&len), "block size {len} (m={m}): {s}..{e}");
            }

            // Full coverage: the union of all blocks must hit every run index.
            let mut covered = vec![false; m];
            for &(s, e) in &blocks { for i in s..e { covered[i as usize] = true; } }
            assert!(covered.iter().all(|&c| c),
                    "uncovered index (m={m}, decreasing={decreasing})");
        }
    }
}
