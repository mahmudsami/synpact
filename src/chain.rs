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
#[derive(Clone, Copy)]
pub(crate) struct Anchor {
    pub(crate) chr:    u8,
    pub(crate) q_pos:  u32,
    pub(crate) r_pos:  u32,
    pub(crate) weight: u32,
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
    seq: &[u8], index: &GIndex, k: usize, s: usize, t: usize, mode: SeedMode,
    max_occ: usize, max_occ_l1: usize,
) -> Vec<Anchor> {
    let forest = timed(&PROF_HIER, || extract_hier_blocks_n(seq, k, s, t, index.max_level(), mode));
    let mut anchors = Vec::new();
    timed(&PROF_EMIT, || {
        let floor = min_lvl().max(1);
        let mut level = forest.top_level;
        if level < floor { return; }
        let mut share = vec![1.0f32; forest.levels[level].len()];
        loop {
            let below = if level > floor { forest.levels[level - 1].len() } else { 0 };
            // memb[c] = number of blocks at this level that contain child c.
            let mut memb = vec![0u32; below];
            let mut child_share = vec![0.0f32; below];
            if below > 0 {
                for p in &forest.levels[level] { for c in p.children.range() { memb[c] += 1; } }
            }
            // Touch every block's bucket before looking any up, so that their
            // cache misses overlap rather than stall the lookups one by one.
            for block in &forest.levels[level] { index.touch(block.level, block.hash); }
            for (block, &sh) in forest.levels[level].iter().zip(&share) {
                if sh == 0.0 { continue; }
                let mocc = max_occ_for_level(block.level, max_occ, max_occ_l1);
                let hits = index.lookup(block.level, block.hash, mocc);
                if !hits.is_empty() {
                    // Vote = the block's syncmer mass times the share this path
                    // emits, scaled ×16 so fractional masses keep integer resolution.
                    let w = (block.mass * sh * 16.0).round().max(1.0) as u32;
                    for (chr_id, genome_pos) in hits.iter() {
                        anchors.push(Anchor { chr: chr_id, q_pos: block.pos, r_pos: genome_pos, weight: w });
                    }
                } else if below > 0 {
                    for c in block.children.range() {
                        child_share[c] += sh / memb[c] as f32;
                    }
                }
            }
            if below == 0 { break; }
            share = child_share;
            level -= 1;
        }
    });
    anchors
}

/// Locus selection by diagonal voting: anchors are sorted by (chr, diagonal =
/// r_pos − q_pos) and a sliding window of width VOTE_BAND accumulates anchor
/// weight; the heaviest window is the locus. Returns (best (chr, ref_offset,
/// score), second_score) — second_score is the heaviest window at a genuinely
/// different locus (used for MAPQ). O(n log n).
pub(crate) fn vote_locus(anchors: &mut Vec<Anchor>) -> (Option<(u8, i64, u32)>, u32) {
    let n = anchors.len();
    if n == 0 { return (None, 0); }
    const VOTE_BAND: i64 = CLUSTER_WINDOW;
    let diag = |a: &Anchor| a.r_pos as i64 - a.q_pos as i64;
    anchors.sort_unstable_by_key(|a| (a.chr, diag(a)));

    // Pass 1: heaviest diagonal window.
    let mut best_sum = 0u32;
    let (mut best_lo, mut best_hi) = (0usize, 0usize);
    let mut lo = 0usize;
    let mut sum = 0u32;
    for hi in 0..n {
        sum = sum.saturating_add(anchors[hi].weight);
        while anchors[lo].chr != anchors[hi].chr
            || diag(&anchors[hi]) - diag(&anchors[lo]) > VOTE_BAND
        {
            sum = sum.saturating_sub(anchors[lo].weight);
            lo += 1;
        }
        if sum > best_sum { best_sum = sum; best_lo = lo; best_hi = hi; }
    }
    if best_sum < MIN_VOTES { return (None, 0); }

    // Representative locus: offset of the smallest-q_pos anchor in the window.
    let best_chr = anchors[best_hi].chr;
    let mut off = diag(&anchors[best_lo]);
    let mut min_q = u32::MAX;
    for a in &anchors[best_lo..=best_hi] {
        if a.q_pos < min_q { min_q = a.q_pos; off = diag(a); }
    }

    // Pass 2: heaviest window at a different locus → second_score for MAPQ.
    let win_diag = diag(&anchors[best_hi]);
    let rest: Vec<Anchor> = anchors.iter().copied()
        .filter(|a| a.chr != best_chr || (diag(a) - win_diag).abs() > VOTE_BAND)
        .collect();
    let mut second = 0u32;
    if !rest.is_empty() {
        let mut lo = 0usize; let mut sum = 0u32;
        for hi in 0..rest.len() {
            sum = sum.saturating_add(rest[hi].weight);
            while rest[lo].chr != rest[hi].chr
                || diag(&rest[hi]) - diag(&rest[lo]) > VOTE_BAND
            {
                sum = sum.saturating_sub(rest[lo].weight);
                lo += 1;
            }
            if sum > second { second = sum; }
        }
    }

    (Some((best_chr, off, best_sum)), second)
}
