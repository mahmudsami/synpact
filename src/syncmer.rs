use crate::hash::*;

/// Lightweight record for genome-scale processing.
#[derive(Clone)]
pub(crate) struct SyncmerLight {
    pub(crate) pos:   u32,
    pub(crate) value: u64,   // s-mer base-4 encoding — used for L1+ block hashing
}

/// 2-bit code per byte (A=0, C=1, G=2, T=3, anything else 0) — same mapping as
/// `base4`, as a table so the rolling k-mer code needs no branches.
const BASE4: [u8; 256] = {
    let mut t = [0u8; 256];
    t[b'C' as usize] = 1; t[b'c' as usize] = 1;
    t[b'G' as usize] = 2; t[b'g' as usize] = 2;
    t[b'T' as usize] = 3; t[b't' as usize] = 3;
    t
};

/// Ring buffer of s-mer hashes. The window w = k − s + 1 ≤ 32 because k ≤ 32,
/// so 64 slots always hold the whole window.
const RING: usize = 64;
const RING_MASK: usize = RING - 1;

/// Open syncmers: k-mers whose minimal s-mer (by `DnaHashFwd` hash, leftmost on
/// ties) sits at offset `t`.
///
/// One pass over the bases computes, in the same step, the rolling s-mer hash,
/// the rolling 2-bit k-mer code (the syncmer's `value`), and the positions of the
/// most recent non-ACGT / 'N' bases. The window minimum is tracked as a single
/// (index, value) pair and only rescanned when the minimum leaves the window
/// (≈ 2/(w+1) of positions for random hashes). Output is identical to the
/// earlier monotone-deque version (same hash, tie-break and N handling).
pub(crate) fn select_syncmers_light(seq: &[u8], k: usize, s: usize, t: usize) -> Vec<SyncmerLight> {
    let n = seq.len();
    if n < k || s == 0 || k > 32 || s > k { return vec![]; }
    let w = k - s + 1;
    debug_assert!(w <= RING);
    let kmask: u64 = if k == 32 { u64::MAX } else { (1u64 << (2 * k)) - 1 };
    let s_rot = s as u32;

    let mut out = Vec::with_capacity((n - k + 1) / w + 1);
    let mut ring = [0u64; RING];

    let mut h: u64 = 0;           // rolling s-mer hash (raw, before N masking)
    let mut code: u64 = 0;        // rolling 2-bit code of the last k bases
    // One past the position of the most recent byte that is non-ACGT (masks the
    // s-mer hash to u64::MAX) / that is 'N' or 'n' (rejects the k-mer). 0 = none.
    let mut bad_end: usize = 0;
    let mut n_end: usize = 0;
    // Leftmost minimum of the current window of s-mer hashes.
    let mut min_idx: usize = 0;
    let mut min_val: u64 = u64::MAX;

    for p in 0..n {
        let b = seq[p];
        let seed_in = DNA_SEED[b as usize];
        if seed_in == 0 { bad_end = p + 1; }
        if b == b'N' || b == b'n' { n_end = p + 1; }
        code = ((code << 2) | BASE4[b as usize] as u64) & kmask;
        h = h.rotate_left(1) ^ seed_in;
        if p >= s {
            h ^= DNA_SEED[seq[p - s] as usize].rotate_left(s_rot);
        }
        if p + 1 < s { continue; }

        // s-mer j = seq[j..=p]
        let j = p + 1 - s;
        let hj = if bad_end > j { u64::MAX } else { h };
        ring[j & RING_MASK] = hj;

        // Window of s-mers [i, j]; the k-mer is seq[i..i+k] = seq[i..=p].
        if j + 1 >= w {
            let i = j + 1 - w;
            if min_idx < i {
                // The minimum fell out: rescan the window (leftmost on ties).
                min_idx = i;
                min_val = ring[i & RING_MASK];
                for q in i + 1..=j {
                    let v = ring[q & RING_MASK];
                    if v < min_val { min_val = v; min_idx = q; }
                }
            } else if hj < min_val {
                min_val = hj; min_idx = j;
            }
            if min_idx == i + t && n_end <= i {
                out.push(SyncmerLight { pos: i as u32, value: code });
            }
        } else if hj < min_val {
            // Still filling the first window.
            min_val = hj; min_idx = j;
        }
    }
    out
}

/// Seed selection mode.  Kept as a single-variant enum so the on-disk index
/// format (which stores a seed-mode byte) stays stable and self-describing.
#[derive(Clone, Copy, PartialEq, Eq, Debug)]
pub(crate) enum SeedMode {
    /// Open syncmer: min s-mer at position t within the k-mer.
    Syncmer,
}

/// Seed extraction — open syncmers (s = s-mer size, t = position).
#[inline]
pub(crate) fn select_seeds_light(seq: &[u8], k: usize, s: usize, t: usize, _mode: SeedMode)
    -> Vec<SyncmerLight>
{
    select_syncmers_light(seq, k, s, t)
}
