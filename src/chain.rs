// Anchor collection from the LCP block hierarchy + diagonal-voting locus
// selection. (Open-syncmer seeds → LCP blocks → index lookups → weighted
// anchors → heaviest diagonal window = locus.)
use crate::config::*;
use crate::syncmer::*;
use crate::lcp::*;
use crate::index::*;

pub(crate) const CLUSTER_WINDOW: i64 = 500;
pub(crate) const MIN_VOTES:      u32 = 1;

/// A seed hit: a query block matched a reference position.
///
/// Voting only needs a hit's chromosome and diagonal (r_pos − q_pos), and sorts
/// by the pair, so both are packed at construction into one integer `key` whose
/// order is the (chr, diagonal) order: one u64 compare per sort step, and a
/// radix sort for the large anchor sets of repetitive reads (see sort_anchors).
#[derive(Clone, Copy, Default)]
pub(crate) struct Anchor {
    key:               u64,
    pub(crate) q_pos:  u32,
    pub(crate) weight: u32,
}

/// Diagonals lie in (−2³², 2³²); biased by 2³² they fill 33 bits and are ≥ 1.
const DIAG_BIAS: i64 = 1 << 32;
/// The chromosome sits above a one-bit gap over the biased diagonal, so keys on
/// different chromosomes are always > 2³³ (≫ VOTE_BAND) apart: "same chromosome
/// and diagonals ≤ VOTE_BAND apart" is just "keys ≤ VOTE_BAND apart".
const CHR_SHIFT: u32 = 34;

impl Anchor {
    #[inline]
    pub(crate) fn new(chr: u8, q_pos: u32, r_pos: u32, weight: u32) -> Anchor {
        let biased_diag = (r_pos as i64 - q_pos as i64 + DIAG_BIAS) as u64;
        Anchor { key: (chr as u64) << CHR_SHIFT | biased_diag, q_pos, weight }
    }
    #[inline]
    fn chr(&self) -> u8 { (self.key >> CHR_SHIFT) as u8 }
    /// Diagonal = r_pos − q_pos.
    #[inline]
    fn diag(&self) -> i64 { (self.key & ((1 << CHR_SHIFT) - 1)) as i64 - DIAG_BIAS }
}

/// Collect anchors for one strand of a read from its LCP block hierarchy.
///
/// Levels are processed top-down. Each block first tries its own level (highest
/// weight) and falls back to its finer children on a miss / too-repetitive.
///
/// `share[i]` is the fraction of block i's mass still to be emitted. LCP shares
/// boundary sub-blocks between adjacent parents, and a parent's mass holds only
/// its share of such a child (child mass / #parents containing it). A missed
/// parent therefore passes each child `share / #parents`, so the weight emitted
/// under it equals its own mass — whether its neighbour hit or missed. (Emitting
/// the child's full mass instead counts the part a neighbour that hit already
/// emitted a second time, so a locus that matches the read worse — more misses,
/// more fallbacks — gains weight.) Shares from all missed parents are summed
/// before the child is looked up, so every block is looked up at most once.
pub(crate) fn collect_anchors(
    syncmers: &[SyncmerLight], index: &GIndex, k: usize,
    max_occ: usize, max_occ_l1: usize,
) -> Vec<Anchor> {
    // Build the read's levels only up to --max-level (and no higher than the index
    // holds), so anchoring starts at min(read's top level, --max-level, index top).
    let top = index.max_level().min(max_lvl());
    let forest = timed(&PROF_HIER, || extract_hier_blocks_n(syncmers, k, top));
    let mut anchors = Vec::new();
    timed(&PROF_EMIT, || {
        let floor = min_lvl().max(1);
        let mut level = forest.top_level;
        if level < floor { return; }
        // Two-stage prefetch over every block that may be looked up, before any
        // is: first each block's bucket-table entry, then (once those have had
        // time to arrive) the entries of the buckets whose filter admits it. The
        // cache misses of all levels thus overlap instead of stalling the
        // lookups one by one. Blocks under a parent that hits are prefetched in
        // vain, but most parents miss.
        let looked_up = || forest.levels[floor..=level].iter().flatten();
        for block in looked_up() { index.prefetch_bucket(block.level, block.hash); }
        for block in looked_up() { index.prefetch_entries(block.level, block.hash); }
        let mut share = vec![1.0f32; forest.levels[level].len()];
        // memb[c] = number of blocks at this level that contain child c;
        // child_share[c] accumulates the share passed down to child c. Both are
        // reused across levels.
        let mut memb: Vec<u32> = Vec::new();
        let mut child_share: Vec<f32> = Vec::new();
        loop {
            let below = if level > floor { forest.levels[level - 1].len() } else { 0 };
            memb.clear(); memb.resize(below, 0);
            child_share.clear(); child_share.resize(below, 0.0);
            if below > 0 {
                for p in &forest.levels[level] { for c in p.children.range() { memb[c] += 1; } }
            }
            for (block, &sh) in forest.levels[level].iter().zip(&share) {
                if sh == 0.0 { continue; }
                let mocc = max_occ_for_level(block.level, max_occ, max_occ_l1);
                let hits = index.lookup(block.level, block.hash, mocc);
                if !hits.is_empty() {
                    // Vote = the block's syncmer mass times the share this path
                    // emits, scaled ×16 so fractional masses keep integer resolution.
                    let w = (block.mass * sh * 16.0).round().max(1.0) as u32;
                    anchors.extend(hits.iter().map(|(chr, r_pos)|
                        Anchor::new(chr, block.pos, r_pos, w)));
                } else if below > 0 {
                    for c in block.children.range() {
                        child_share[c] += sh / memb[c] as f32;
                    }
                }
            }
            if below == 0 { break; }
            std::mem::swap(&mut share, &mut child_share);
            level -= 1;
        }
    });
    anchors
}

const VOTE_BAND: i64 = CLUSTER_WINDOW;

/// Anchor sets at least this large are radix sorted; below it a comparison sort
/// is cheaper than the radix sort's fixed histogram work.
const RADIX_MIN: usize = 256;
/// Radix digit width: 2¹¹ counters per pass (8 KB) stay in L1.
const RADIX_BITS: u32 = 11;

/// Sort anchors by key, i.e. by (chr, diagonal).
///
/// The order of anchors with equal keys is left unspecified: vote_locus's
/// result does not depend on it when no window sum saturates (see there). A few reads hit repeats and carry
/// tens of thousands of anchors — they dominate the sorting time — so large
/// sets get an LSD radix sort over only the bits in which their keys differ.
/// `scratch` is the radix sort's second buffer; reusing it across calls saves
/// an allocation (the two buffers may trade places).
fn sort_anchors(anchors: &mut Vec<Anchor>, scratch: &mut Vec<Anchor>) {
    let n = anchors.len();
    if n < RADIX_MIN {
        anchors.sort_unstable_by_key(|a| a.key);
        return;
    }
    let (mut min, mut max) = (u64::MAX, 0u64);
    for a in anchors.iter() { min = min.min(a.key); max = max.max(a.key); }
    // Keys are < 2⁴², so at most 4 digits of key − min.
    let passes = (u64::BITS - (max - min).leading_zeros()).div_ceil(RADIX_BITS) as usize;
    const MASK: u64 = (1 << RADIX_BITS) - 1;
    let mut counts = [[0u32; 1 << RADIX_BITS]; 4];
    for a in anchors.iter() {
        let k = a.key - min;
        for (p, c) in counts.iter_mut().enumerate() {
            c[(k >> (p as u32 * RADIX_BITS) & MASK) as usize] += 1;
        }
    }
    scratch.clear();
    scratch.resize(n, Anchor::default());
    for (p, c) in counts.iter_mut().enumerate().take(passes) {
        let shift = p as u32 * RADIX_BITS;
        // A digit shared by every key leaves the order unchanged.
        if c[((anchors[0].key - min) >> shift & MASK) as usize] as usize == n { continue; }
        let mut start = 0u32;
        for slot in c.iter_mut() { let cnt = *slot; *slot = start; start += cnt; }
        for a in anchors.iter() {
            let slot = &mut c[((a.key - min) >> shift & MASK) as usize];
            scratch[*slot as usize] = *a;
            *slot += 1;
        }
        std::mem::swap(anchors, scratch);
    }
}

/// Heaviest window of `anchors` (sorted by key) spanning at most VOTE_BAND
/// diagonals on one chromosome: (weight, first, last) of the leftmost.
/// `record(sum)` is called with the weight of each window [lo, hi], in hi order.
///
/// Keys of one chromosome differ exactly by their diagonals, and keys of
/// different chromosomes by more than VOTE_BAND (see CHR_SHIFT), so a single
/// key difference decides whether `lo` is still in the window of `hi`.
#[inline]
fn heaviest_window(anchors: &[Anchor], mut record: impl FnMut(u32)) -> (u32, usize, usize) {
    let (mut best_sum, mut best_lo, mut best_hi) = (0u32, 0usize, 0usize);
    let mut lo = 0usize;
    let mut sum = 0u32;
    for hi in 0..anchors.len() {
        sum = sum.saturating_add(anchors[hi].weight);
        while anchors[hi].key - anchors[lo].key > VOTE_BAND as u64 {
            sum = sum.saturating_sub(anchors[lo].weight);
            lo += 1;
        }
        record(sum);
        if sum > best_sum { best_sum = sum; best_lo = lo; best_hi = hi; }
    }
    (best_sum, best_lo, best_hi)
}

/// Buffers for vote_locus, created once per read and shared by both strands.
#[derive(Default)]
pub(crate) struct VoteBuf {
    /// The radix sort's second buffer (see sort_anchors).
    scratch: Vec<Anchor>,
    /// sums[hi] = weight of pass 1's window ending at anchor hi.
    sums:    Vec<u32>,
}

/// Locus selection by diagonal voting: anchors are sorted by (chr, diagonal =
/// r_pos − q_pos) and a sliding window of width VOTE_BAND accumulates anchor
/// weight; the heaviest window is the locus. Returns (best (chr, ref_offset,
/// score), second_score) — second_score is the heaviest window at a genuinely
/// different locus (used for MAPQ). O(n log n).
///
/// Unless window sums saturate, the result does not depend on how the sort
/// orders anchors with equal keys: a window's start only moves on key
/// comparisons, so it always sits at the first anchor of a key group; and as
/// every weight is ≥ 1 the window sum grows within a group, so the first window
/// to reach the maximum ends at the last anchor of a group. The best window is
/// thus a union of whole key groups, its sum fixed, and the min-q_pos
/// representative's diagonal is that of its group. A total weight that fits a
/// u32 rules saturation out; otherwise (never seen in practice) the comparison
/// sort and two full scans give the same result as always.
pub(crate) fn vote_locus(anchors: &mut Vec<Anchor>, buf: &mut VoteBuf)
    -> (Option<(u8, i64, u32)>, u32)
{
    if anchors.is_empty() { return (None, 0); }
    let total: u64 = anchors.iter().map(|x| x.weight as u64).sum();
    let unsaturated = total <= u32::MAX as u64;
    if unsaturated {
        sort_anchors(anchors, &mut buf.scratch);
    } else {
        anchors.sort_unstable_by_key(|a| a.key);
    }

    // Pass 1: heaviest diagonal window.
    let sums = &mut buf.sums;
    sums.clear();
    let (best_sum, best_lo, best_hi) = heaviest_window(anchors, |s| sums.push(s));
    if best_sum < MIN_VOTES { return (None, 0); }

    // Representative locus: offset of the smallest-q_pos anchor in the window.
    let best_chr = anchors[best_hi].chr();
    let mut off = anchors[best_lo].diag();
    let mut min_q = u32::MAX;
    for a in &anchors[best_lo..=best_hi] {
        if a.q_pos < min_q { min_q = a.q_pos; off = a.diag(); }
    }

    // Pass 2: heaviest window at a different locus → second_score for MAPQ. The
    // anchors within VOTE_BAND of the window's diagonal are one contiguous run
    // [a, b); no window spans it (its two sides are > 2·VOTE_BAND apart), so
    // each side is scanned on its own. Every key is ≥ (chr << CHR_SHIFT) + 1,
    // and the previous chromosome's keys end 2³³ below that, so the run's bounds
    // (VOTE_BAND below / above the window's key) never reach another chromosome.
    let win_key = anchors[best_hi].key;
    let a = anchors.partition_point(|x| x.key < win_key.saturating_sub(VOTE_BAND as u64));
    let b = anchors.partition_point(|x| x.key <= win_key + VOTE_BAND as u64);
    let second = if unsaturated {
        // No sum saturates, so most of pass 1's window sums can be reused rather
        // than rescanning both sides. Windows ending before a lie wholly before
        // a. Of those ending at or after b, only the ones still reaching back
        // across b — ending before h, the first anchor more than VOTE_BAND past
        // anchor b − 1 — are cut short there, and only they are rescanned.
        let h = match b {
            0 => 0,
            _ => b + anchors[b..].partition_point(|x| x.key - anchors[b - 1].key <= VOTE_BAND as u64),
        };
        let max = |s: &[u32]| s.iter().copied().max().unwrap_or(0);
        max(&sums[..a]).max(max(&sums[h..])).max(heaviest_window(&anchors[b..h], |_| {}).0)
    } else {
        heaviest_window(&anchors[..a], |_| {}).0.max(heaviest_window(&anchors[b..], |_| {}).0)
    };

    (Some((best_chr, off, best_sum)), second)
}

#[cfg(test)]
mod tests {
    use super::*;

    struct Rng(u64);
    impl Rng {
        fn next(&mut self) -> u64 {
            let mut x = self.0;
            x ^= x << 13; x ^= x >> 7; x ^= x << 17;
            self.0 = x; x
        }
    }

    /// The voting as written before anchors carried a packed key: comparison
    /// sort on (chr, diagonal) and two full window scans.
    fn vote_reference(mut anchors: Vec<(u8, u32, u32, u32)>) -> (Option<(u8, i64, u32)>, u32) {
        let diag = |a: &(u8, u32, u32, u32)| a.2 as i64 - a.1 as i64;
        let window = |v: &[(u8, u32, u32, u32)]| {
            let (mut best_sum, mut best_lo, mut best_hi) = (0u32, 0usize, 0usize);
            let (mut lo, mut sum) = (0usize, 0u32);
            for hi in 0..v.len() {
                sum = sum.saturating_add(v[hi].3);
                while v[lo].0 != v[hi].0 || diag(&v[hi]) - diag(&v[lo]) > VOTE_BAND {
                    sum = sum.saturating_sub(v[lo].3);
                    lo += 1;
                }
                if sum > best_sum { best_sum = sum; best_lo = lo; best_hi = hi; }
            }
            (best_sum, best_lo, best_hi)
        };
        if anchors.is_empty() { return (None, 0); }
        anchors.sort_unstable_by_key(|a| (a.0, diag(a)));
        let (best_sum, best_lo, best_hi) = window(&anchors);
        if best_sum < MIN_VOTES { return (None, 0); }
        let best_chr = anchors[best_hi].0;
        let mut off = diag(&anchors[best_lo]);
        let mut min_q = u32::MAX;
        for a in &anchors[best_lo..=best_hi] {
            if a.1 < min_q { min_q = a.1; off = diag(a); }
        }
        let win_diag = diag(&anchors[best_hi]);
        let a = anchors.partition_point(|x| (x.0, diag(x)) < (best_chr, win_diag - VOTE_BAND));
        let b = anchors.partition_point(|x| (x.0, diag(x)) <= (best_chr, win_diag + VOTE_BAND));
        let second = window(&anchors[..a]).0.max(window(&anchors[b..]).0);
        (Some((best_chr, off, best_sum)), second)
    }

    /// vote_locus (packed keys, radix sort, reused window sums) must agree with
    /// the reference exactly: clustered diagonals with many equal keys, several
    /// chromosomes (incl. 0 and 255), extreme positions, both sort paths, and
    /// weights large enough to saturate the window sums.
    #[test]
    fn vote_matches_reference() {
        let mut rng = Rng(0x0DDB1A5E5BAD5EED);
        let mut buf = VoteBuf::default();
        for iter in 0..3000 {
            let n = match iter % 4 { 0 => rng.next() % 40, 1 => rng.next() % 600, _ => rng.next() % 4000 } as usize;
            let chrs = [0u8, 1, 2, 7, 254, 255];
            let n_chr = 1 + (rng.next() % 6) as usize;
            let n_centres = 1 + (rng.next() % 12) as usize;
            let centres: Vec<(u8, i64)> = (0..n_centres).map(|_| {
                let c = chrs[(rng.next() % n_chr as u64) as usize];
                let d = match rng.next() % 4 {
                    0 => -(u32::MAX as i64) + (rng.next() % 2000) as i64,   // r_pos 0, q_pos huge
                    1 => u32::MAX as i64 - (rng.next() % 2000) as i64,      // r_pos huge, q_pos 0
                    _ => (rng.next() % 200_000) as i64 - 50_000,
                };
                (c, d)
            }).collect();
            let saturate = iter % 10 == 0;
            let mut raw: Vec<(u8, u32, u32, u32)> = Vec::with_capacity(n);
            for _ in 0..n {
                let (c, d) = centres[(rng.next() % n_centres as u64) as usize];
                let d = d + (rng.next() % 1500) as i64 - 750;   // within ±3 bands
                // q_pos so that r_pos = q_pos + d stays in u32.
                let q_lo = (-d).max(0);
                let q_hi = (u32::MAX as i64 - d).min(u32::MAX as i64);
                if q_lo > q_hi { continue; }
                let (q_lo, q_hi) = (q_lo as u64, q_hi as u64);
                let q = q_lo + rng.next() % (q_hi - q_lo + 1).min(30_000);
                let r = (q as i64 + d) as u32;
                let w = if saturate { 1 + (rng.next() % (u32::MAX as u64 / 8)) as u32 }
                        else { 1 + (rng.next() % 3000) as u32 };
                raw.push((c, q as u32, r, w));
                if rng.next() % 3 == 0 {   // an equal key with another q_pos / weight
                    let q2 = (q + rng.next() % 50).min(q_hi);
                    raw.push((c, q2 as u32, (q2 as i64 + d) as u32, 1 + (rng.next() % 3000) as u32));
                }
            }
            let mut anchors: Vec<Anchor> = raw.iter().map(|&(c, q, r, w)| Anchor::new(c, q, r, w)).collect();
            assert_eq!(vote_locus(&mut anchors, &mut buf), vote_reference(raw), "iter {iter}");
        }
    }
}
