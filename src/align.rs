#![allow(unused_imports, dead_code)]
use crate::config::*;
use crate::hash::*;
use crate::syncmer::*;
use crate::lcp::*;
use crate::index::*;
use crate::fastq::*;
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

pub(crate) struct AlignBufs { pub(crate) dp: Vec<i32>, pub(crate) bt: Vec<u8> }
thread_local! {
    static ALIGN: RefCell<AlignBufs> = RefCell::new(AlignBufs { dp: Vec::new(), bt: Vec::new() });
}

/// Banded semi-global alignment: aligns `query` fully against `target`.
/// `half_band` = maximum allowed diagonal deviation in bases.
/// Scoring: match +2, mismatch −4, linear gap −2.
/// Thread-local DP+BT buffers are reused across calls — no per-read allocation.
/// Returns extended CIGAR (= X I D). Falls back to `{n}M` if band is too narrow.
pub(crate) fn banded_align(query: &[u8], target: &[u8], half_band: usize) -> String {
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
pub(crate) fn reverse_cigar_ops(cg: &str) -> String {
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
pub(crate) fn cigar_for_mapping(
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

/// Count reference bases consumed by a CIGAR string (=, X, D, M operators).
#[inline]
pub(crate) fn ref_bases_in_cigar(cg: &str) -> usize {
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
pub(crate) fn cigar_from_chain_anchors(
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
