#![allow(unused_imports, dead_code)]
use crate::config::*;
use crate::hash::*;
use crate::syncmer::*;
use crate::lcp::*;
use crate::index::*;
use crate::fastq::*;
use crate::align::*;
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

pub(crate) const CLUSTER_WINDOW: i64  = 500;

pub(crate) const MIN_VOTES:      u32  = 1;

/// A seed hit: a HierBlock from the query matched at a specific reference position.
/// q_pos..q_pos+span  covers the query interval.
/// r_pos..r_pos+span  is the corresponding reference interval (if colinear).
#[derive(Clone, Copy)]
pub(crate) struct Anchor {
    pub(crate) chr:    u8,
    pub(crate) q_pos:  u32,
    pub(crate) r_pos:  u32,
    pub(crate) span:   u32,
    pub(crate) weight: u32,
}

/// Recursively emit anchors for one HierBlock.
/// Tries the block's own level first (highest weight); falls back to children.
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
    visited: &mut std::collections::HashSet<(u8, u32, u64)>,
) {
    // Level floor: ignore blocks below MIN_LVL (don't emit, don't recurse deeper).
    if level < min_lvl() { return; }

    let block = &forest.levels[level][idx];

    // LCP rules intentionally share boundary sub-blocks between adjacent parents,
    // so a shared child is reachable from two parents that both fell through to
    // their children. Deduplicate by block identity (level, query-pos, hash) so
    // its mass is emitted exactly once.
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
    // Block not found or too repetitive — try finer-grained children (level-1).
    // `forest` is only ever shared-borrowed, so the recursion needs no clone.
    for ci in 0..block.children.len() {
        let c = forest.levels[level][idx].children[ci] as usize;
        emit_anchors(forest, level - 1, c, index,
                     max_occ, max_occ_l1, filter, anchors, visited);
    }
}

/// Minimum distinct anchor positions from the L1+ hierarchy below which
/// we fall back to individual L0 k-mer anchors.
/// Minimum distinct anchor count from L1+ hierarchy below which the L0 raw-syncmer
/// fallback is activated.  Reads with enough L1+ anchors skip the L0 pass entirely.
pub(crate) const L0_FALLBACK_THRESHOLD: usize = 1;

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
    let mut visited: std::collections::HashSet<(u8, u32, u64)> = std::collections::HashSet::new();
    timed(&PROF_EMIT, || {
        let tl = forest.top_level;
        if tl >= 1 {
            for idx in 0..forest.levels[tl].len() {
                emit_anchors(&forest, tl, idx, index, max_occ, max_occ_l1,
                             filter, &mut anchors, &mut visited);
            }
        }
    });

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

    timed(&PROF_PRUNE, || prune_anchors(&mut anchors));
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
pub(crate) fn prune_anchors(anchors: &mut Vec<Anchor>) {
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

/// Max query gap between consecutive anchors in a chain (bases).
pub(crate) const CHAIN_MAX_GAP: u32 = 5_000;

/// Allowed absolute difference |dq - dr| before a chain extension is rejected.
/// Also admits 10% of max(dq,dr) for proportional tolerance.
pub(crate) const CHAIN_GAP_TOL: u32 = 150;

/// Penalty per base of gap inconsistency |dq - dr|.
pub(crate) const CHAIN_GAP_SCALE: u32 = 1;

/// Single forward chain-DP pass on a slice of anchors pre-sorted by (chr, q_pos).
/// Returns (chr, ref_start_offset, best_score) or None if nothing exceeds MIN_VOTES.
/// ref_start_offset = r_pos - q_pos at the chain's first anchor.
pub(crate) fn single_chain(anchors: &[Anchor]) -> Option<(u8, i64, u32)> {
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
pub(crate) fn single_chain_with_trace(anchors: &[Anchor]) -> Option<(u8, i64, u32, Vec<Anchor>)> {
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

pub(crate) fn chain_dp(anchors: &mut Vec<Anchor>) -> (Option<(u8, i64, u32)>, u32) {
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

/// Fast locus selection by diagonal voting (no chaining).
///
/// Anchors are sorted by (chr, diagonal = r_pos − q_pos) and a sliding window of
/// width VOTE_BAND accumulates anchor weight; the heaviest window is the locus.
/// Returns (best (chr, ref_offset, score), second_score) like `chain_dp`, where
/// second_score is the heaviest window at a *genuinely different* locus (used for
/// MAPQ). O(n log n) — no O(n²) gap-penalised DP. Default selector (`--vote`).
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

    // Representative locus: offset of the smallest-q_pos anchor in the window
    // (mirrors chain_dp's "ref start = first anchor's offset").
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
