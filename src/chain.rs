// Anchor collection from the LCP block hierarchy + diagonal-voting locus
// selection. (Open-syncmer seeds → LCP blocks → index lookups → weighted
// anchors → heaviest diagonal window = locus.)
use crate::config::*;
use crate::syncmer::*;
use crate::lcp::*;
use crate::index::*;
use std::collections::HashSet;

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

/// Recursively emit anchors for one block: try its own level first (highest
/// weight), falling back to its finer children on a miss/too-repetitive.
/// `filter` is (chr, lo, hi) in reference-start-offset coords (r_pos - q_pos).
pub(crate) fn emit_anchors(
    forest: &HierForest,
    level: usize,
    idx: usize,
    index: &GIndex,
    max_occ: usize,
    max_occ_l1: usize,
    filter: Option<(u8, i64, i64)>,
    anchors: &mut Vec<Anchor>,
    visited: &mut HashSet<(u8, u32, u64)>,
) {
    // Level floor: ignore blocks below --min-level (don't emit, don't recurse).
    if level < min_lvl() { return; }

    let block = &forest.levels[level][idx];

    // LCP rules share boundary sub-blocks between adjacent parents; deduplicate
    // by block identity (level, query-pos, hash) so its mass is emitted once.
    if !visited.insert((block.level as u8, block.pos, block.hash)) { return; }

    let mocc = max_occ_for_level(block.level, max_occ, max_occ_l1);
    let hits = index.lookup(block.level, block.hash, mocc);

    if !hits.is_empty() {
        // Conservation-of-mass weight: a block's vote = its accumulated syncmer
        // mass (each syncmer contributes 1, split across the blocks sharing it),
        // scaled ×MASS_SCALE so fractional masses keep integer resolution.
        const MASS_SCALE: f32 = 16.0;
        let w = (block.mass * MASS_SCALE).round().max(1.0) as u32;
        for (chr_id, genome_pos) in hits.iter() {
            if let Some((fc, flo, fhi)) = filter {
                let vote_pos = genome_pos as i64 - block.pos as i64;
                if chr_id != fc || vote_pos < flo || vote_pos > fhi { continue; }
            }
            anchors.push(Anchor { chr: chr_id, q_pos: block.pos, r_pos: genome_pos, weight: w });
        }
        return;
    }
    // Block not found or too repetitive — recurse into finer children (level-1).
    for ci in 0..block.children.len() {
        let c = forest.levels[level][idx].children[ci] as usize;
        emit_anchors(forest, level - 1, c, index,
                     max_occ, max_occ_l1, filter, anchors, visited);
    }
}

/// Collect anchors for one strand of a read from its LCP block hierarchy.
pub(crate) fn collect_anchors(
    seq: &[u8], index: &GIndex, k: usize, s: usize, t: usize, mode: SeedMode,
    max_occ: usize, max_occ_l1: usize,
    filter_chr: Option<u8>, filter_lo: Option<i64>, filter_hi: Option<i64>,
) -> Vec<Anchor> {
    let filter = match (filter_chr, filter_lo, filter_hi) {
        (Some(c), Some(lo), Some(hi)) => Some((c, lo, hi)),
        _ => None,
    };
    let forest = timed(&PROF_HIER, || extract_hier_blocks_n(seq, k, s, t, index.num_levels(), mode));
    let mut anchors = Vec::new();
    let mut visited: HashSet<(u8, u32, u64)> = HashSet::new();
    timed(&PROF_EMIT, || {
        let tl = forest.top_level;
        if tl >= 1 {
            for idx in 0..forest.levels[tl].len() {
                emit_anchors(&forest, tl, idx, index, max_occ, max_occ_l1,
                             filter, &mut anchors, &mut visited);
            }
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
