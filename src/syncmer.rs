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
/// The default window (w = k − s + 1 = 10, t = 4, i.e. k=19 s=10) runs a
/// branchless kernel specialised at compile time; other parameters use the
/// general tracked-minimum version. Both give identical output.
pub(crate) fn select_syncmers_light(seq: &[u8], k: usize, s: usize, t: usize) -> Vec<SyncmerLight> {
    let n = seq.len();
    if n < k || s == 0 || k > 32 || s > k { return vec![]; }
    match (k - s + 1, t) {
        (10, 4) => select_syncmers_window::<10, 4>(seq, k, s),
        _ => select_syncmers_tracked(seq, k, s, t),
    }
}

/// Branchless open-syncmer selection for a compile-time window of `W` s-mers
/// with the syncmer s-mer at offset `T`.
///
/// The s-mer at offset `T` is the window's leftmost minimum iff it is strictly
/// below every s-mer to its left and ≤ every s-mer to its right, so each window
/// is tested with two fixed-length (fully unrolled) min-reductions over the ring
/// instead of an incrementally tracked minimum, whose "minimum fell out / new
/// minimum / is a syncmer" branches are data-dependent and mispredict often on
/// random hashes. Every k-mer is written to the output unconditionally and the
/// write cursor advances only for syncmers, so the loop has no data-dependent
/// branch. ~1.4× faster than `select_syncmers_tracked` for W=10, T=4.
///
/// Requires n ≥ k (checked by the caller) and k − s + 1 == W.
fn select_syncmers_window<const W: usize, const T: usize>(seq: &[u8], k: usize, s: usize)
    -> Vec<SyncmerLight>
{
    debug_assert!(W <= RING && T < W && k - s + 1 == W && seq.len() >= k);
    let n = seq.len();
    let kmask: u64 = if k == 32 { u64::MAX } else { (1u64 << (2 * k)) - 1 };
    let s_rot = s as u32;
    let nk = n - k + 1;

    // One slot per k-mer; `cnt` only advances past slots holding a syncmer.
    let mut out: Vec<SyncmerLight> = Vec::with_capacity(nk);
    let slots = out.spare_capacity_mut();
    let mut cnt: usize = 0;
    let mut ring = [0u64; RING];

    let mut h: u64 = 0;
    let mut code: u64 = 0;
    let mut bad_end: usize = 0;
    let mut n_end: usize = 0;

    for p in 0..n {
        let b = seq[p];
        let seed_in = DNA_SEED[b as usize];
        bad_end = if seed_in == 0 { p + 1 } else { bad_end };
        n_end = if b == b'N' || b == b'n' { p + 1 } else { n_end };
        code = ((code << 2) | BASE4[b as usize] as u64) & kmask;
        h = h.rotate_left(1) ^ seed_in;
        if p >= s {
            h ^= DNA_SEED[seq[p - s] as usize].rotate_left(s_rot);
        }
        if p + 1 < s { continue; }

        // s-mer j = seq[j..=p]
        let j = p + 1 - s;
        ring[j & RING_MASK] = if bad_end > j { u64::MAX } else { h };
        if j + 1 < W { continue; }

        // Window of s-mers [i, j]; the k-mer is seq[i..i+k] = seq[i..=p].
        let i = j + 1 - W;
        let c = ring[(i + T) & RING_MASK];
        let mut lmin = u64::MAX;
        for q in 0..T { lmin = lmin.min(ring[(i + q) & RING_MASK]); }
        let mut rmin = u64::MAX;
        for q in T + 1..W { rmin = rmin.min(ring[(i + q) & RING_MASK]); }
        let is_sync = (T == 0 || c < lmin) & (c <= rmin) & (n_end <= i);
        // i < nk always; the check only lets the compiler drop the bounds panic.
        if cnt < nk { slots[cnt].write(SyncmerLight { pos: i as u32, value: code }); }
        cnt += is_sync as usize;
    }
    // SAFETY: slots[..cnt] were all initialised above (cnt ≤ number of k-mers = nk).
    unsafe { out.set_len(cnt); }
    out
}

/// General open-syncmer selection for any k/s/t.
///
/// One pass over the bases computes, in the same step, the rolling s-mer hash,
/// the rolling 2-bit k-mer code (the syncmer's `value`), and the positions of the
/// most recent non-ACGT / 'N' bases. The window minimum is tracked as a single
/// (index, value) pair and only rescanned when the minimum leaves the window
/// (≈ 2/(w+1) of positions for random hashes). Output is identical to the
/// earlier monotone-deque version (same hash, tie-break and N handling).
fn select_syncmers_tracked(seq: &[u8], k: usize, s: usize, t: usize) -> Vec<SyncmerLight> {
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

    /// The branchless window kernel must match the tracked-minimum version
    /// exactly, including N runs, other non-ACGT bytes, lowercase, low-complexity
    /// repeats (equal hashes → tie-breaking) and every k that yields w = 10.
    #[test]
    fn branchless_matches_reference() {
        let mut rng = Rng(0xDEADBEEF12345678);
        let alphabet = b"ACGTACGTACGTacgtNRY";
        for iter in 0..4000 {
            let n = (rng.next() % 600) as usize;
            let mut seq: Vec<u8> = (0..n).map(|_| {
                let r = rng.next();
                if r % 50 == 0 { alphabet[(r >> 8) as usize % alphabet.len()] }
                else { b"ACGT"[(r >> 8) as usize % 4] }
            }).collect();
            if iter % 7 == 0 && n > 40 {             // homopolymer / dinucleotide run
                let a = (rng.next() as usize) % (n - 30);
                for q in a..a + 30 { seq[q] = if iter % 2 == 0 { b'A' } else { b"AC"[q % 2] }; }
            }
            // Every k with s = k − 9 has the specialised window (w=10); t=4 selects it.
            let k = 10 + (rng.next() % 23) as usize;  // 10..=32
            for (k, s) in [(k, k - 9), (19, 10)] {
                if n < k { continue; }   // the dispatcher returns early for n < k
                let a = select_syncmers_window::<10, 4>(&seq, k, s);
                let b = select_syncmers_tracked(&seq, k, s, 4);
                assert_eq!(a.len(), b.len(), "count k={k} s={s} n={n}");
                for (x, y) in a.iter().zip(&b) {
                    assert_eq!((x.pos, x.value), (y.pos, y.value), "k={k} s={s}");
                }
            }
        }
    }
}
