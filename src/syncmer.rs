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

/// 2-bit code of each byte's complement (3 − `BASE4` for ACGT in either case,
/// anything else 0) — the code the byte gets in the reverse complement, since
/// `revcomp` maps every non-ACGT byte to a non-ACGT byte.
const BASE4_COMP: [u8; 256] = {
    let mut t = [0u8; 256];
    t[b'A' as usize] = 3; t[b'a' as usize] = 3;
    t[b'C' as usize] = 2; t[b'c' as usize] = 2;
    t[b'G' as usize] = 1; t[b'g' as usize] = 1;
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
/// block-wise kernel specialised at compile time; other parameters use the
/// general tracked-minimum version. Both give identical output.
pub(crate) fn select_syncmers_light(seq: &[u8], k: usize, s: usize, t: usize) -> Vec<SyncmerLight> {
    let n = seq.len();
    if n < k || s == 0 || k > 32 || s > k { return vec![]; }
    match (k - s + 1, t) {
        (10, 4) => select_syncmers_chunked::<10, 4, false>(seq, k, s).0,
        _ => select_syncmers_tracked(seq, k, s, t),
    }
}

/// Open syncmers of `seq` and of its reverse complement, from one pass over `seq`:
/// exactly `(select_syncmers_light(seq, ..), select_syncmers_light(&revcomp(seq), ..))`,
/// the second in ascending reverse-complement coordinates — without building the
/// reverse complement, and sharing the per-base work of both strands.
#[allow(dead_code)]   // the read mapper can use it in place of revcomp + two calls
pub(crate) fn select_syncmers_both(seq: &[u8], k: usize, s: usize, t: usize)
    -> (Vec<SyncmerLight>, Vec<SyncmerLight>)
{
    let n = seq.len();
    if n < k || s == 0 || k > 32 || s > k { return (vec![], vec![]); }
    match (k - s + 1, t) {
        (10, 4) => select_syncmers_chunked::<10, 4, true>(seq, k, s),
        _ => (select_syncmers_tracked(seq, k, s, t),
              select_syncmers_tracked(&crate::fastq::revcomp(seq), k, s, t)),
    }
}

/// k-mers handled per block by `select_syncmers_chunked`. Its buffers (~50 B per
/// k-mer with both strands) stay in L1.
const CHUNK: usize = 256;

/// True iff every byte is A/C/G/T in either case — i.e. no byte masks an s-mer
/// hash (`DNA_SEED` = 0) or rejects a k-mer ('N'). A branch-free fold, so it
/// vectorises (≈ 0.05 ns/base).
#[inline]
fn all_acgt(seq: &[u8]) -> bool {
    seq.iter().fold(0u8, |acc, &b| {
        let u = b & 0xDF;   // ASCII uppercase for letters; no other byte maps to ACGT
        acc | !((u == b'A') | (u == b'C') | (u == b'G') | (u == b'T')) as u8
    }) == 0
}

/// Largest power of two ≤ x (x ≥ 1).
const fn pow2_floor(x: usize) -> usize { 1 << (usize::BITS - 1 - x.leading_zeros()) }

/// Minimum of the `len` values starting at `a`, given `mp[y]` = minimum of the `p`
/// values starting at y (p ≤ len, or len = 0). Runs may overlap: min is idempotent.
#[inline(always)]
fn run_min(mp: &[u64], a: usize, len: usize, p: usize) -> u64 {
    let mut m = u64::MAX;
    if len == 0 { return m; }
    let mut o = 0;
    while o + p < len { m = m.min(mp[a + o]); o += p; }
    m.min(mp[a + len - p])
}

/// `mp[y]` = minimum of `h[y .. y + p]` for y in 0..=len−p (p a power of two), by
/// doubling. Each pass is a plain elementwise min, so it vectorises.
#[inline(always)]
fn sliding_min(h: &[u64], mp: &mut [u64], len: usize, p: usize) {
    if p == 1 { mp[..len].copy_from_slice(&h[..len]); return; }
    for y in 0..len - 1 { mp[y] = h[y].min(h[y + 1]); }
    let mut d = 2;
    while d < p {
        for y in 0..len - 2 * d + 1 { mp[y] = mp[y].min(mp[y + d]); }
        d *= 2;
    }
}

/// Visits, in ascending order, every x whose flag byte (0/1) is 1. `flags.len()`
/// must be a multiple of 8. 64 flags are gathered into a bit mask and only its
/// set bits (~1 in w) are visited, so a non-syncmer k-mer costs a fraction of an
/// instruction.
#[inline(always)]
fn for_each_flagged(flags: &[u8], mut f: impl FnMut(usize)) {
    for (g, f64) in flags.chunks(64).enumerate() {
        let mut m: u64 = 0;
        for (q, f8) in f64.chunks_exact(8).enumerate() {
            let v = u64::from_le_bytes(f8.try_into().unwrap());
            // Bytes 0/1 → bits 0..7: the multiply moves byte i's bit to bit 56 + i
            // (all partial products land on distinct bits, so nothing carries).
            m |= (v.wrapping_mul(0x0102_0408_1020_4080) >> 56) << (8 * q);
        }
        while m != 0 {
            f(64 * g + m.trailing_zeros() as usize);
            m &= m - 1;
        }
    }
}

/// Per-call seed tables (zero for non-ACGT, like `DNA_SEED`), pre-rotated so each
/// roll is one rotate and two xors: outgoing forward base rotated left by s,
/// incoming RC base (complement) left by s−1, outgoing RC base right by 1.
struct SeedTables { f_out: [u64; 256], r_in: [u64; 256], r_out: [u64; 256] }

impl SeedTables {
    fn new<const BOTH: bool>(s: usize) -> Self {
        let mut t = SeedTables { f_out: [0; 256], r_in: [0; 256], r_out: [0; 256] };
        for (b, c) in [(b'A', b'T'), (b'C', b'G'), (b'G', b'C'), (b'T', b'A')] {
            for b in [b, b.to_ascii_lowercase()] {
                t.f_out[b as usize] = DNA_SEED[b as usize].rotate_left(s as u32);
                if BOTH {
                    t.r_in[b as usize] = DNA_SEED[c as usize].rotate_left(s as u32 - 1);
                    t.r_out[b as usize] = DNA_SEED[c as usize].rotate_right(1);
                }
            }
        }
        t
    }
}

/// Rolls the hashes and codes over a block of plain A/C/G/T bases: base x enters
/// as `ins[x]` while `outs[x]` leaves the s-mer; writes the s-mer hash and k-mer
/// code ending there (and, with `BOTH`, the RC ones). State = (h, rh, code, rcode).
///
/// Two bases per step, each state advanced straight from the one two bases back,
/// so the serial chains take one link per two bases; the in-between state hangs
/// off the chain. The two-base hash step is written rotl(h ⊕ rotr(Δ, 2), 2), not
/// the equal rotl(h, 2) ⊕ Δ: LLVM flattens a xor tree and puts the loop-carried h
/// first, which would chain h through every term of Δ. Kept out of line so its
/// ~20 live values get registers of their own (inlined, the RC state spilled).
#[inline(never)]
#[allow(clippy::too_many_arguments)]
fn roll_clean<const BOTH: bool>(
    (mut h, mut rh, mut code, mut rcode): (u64, u64, u64, u64), tabs: &SeedTables, top: u32,
    ins: &[u8], outs: &[u8], hs: &mut [u64], cs: &mut [u64], rhs: &mut [u64], rcs: &mut [u64],
) -> (u64, u64, u64, u64) {
    let c = ins.len();
    let (f_out, r_in, r_out) = (&tabs.f_out, &tabs.r_in, &tabs.r_out);
    let (hs, cs) = (&mut hs[..c], &mut cs[..c]);
    let (rhs, rcs) = if BOTH { (&mut rhs[..c], &mut rcs[..c]) } else { (rhs, rcs) };
    for (x, (i2, o2)) in ins.chunks_exact(2).zip(outs[..c].chunks_exact(2)).enumerate() {
        let (b0, b1) = (i2[0] as usize, i2[1] as usize);
        let (o0, o1) = (o2[0] as usize, o2[1] as usize);
        let d0 = DNA_SEED[b0] ^ f_out[o0];
        let d1 = DNA_SEED[b1] ^ f_out[o1];
        let (e0, e1) = (BASE4[b0] as u64, BASE4[b1] as u64);
        hs[2 * x] = h.rotate_left(1) ^ d0;
        cs[2 * x] = (code << 2) | e0;
        h = (h ^ (d0.rotate_right(1) ^ d1.rotate_right(2))).rotate_left(2);
        code = (code << 4) | ((e0 << 2) | e1);
        hs[2 * x + 1] = h;
        cs[2 * x + 1] = code;
        if BOTH {
            let r0 = r_in[b0] ^ r_out[o0];
            let r1 = r_in[b1] ^ r_out[o1];
            let (g0, g1) = (BASE4_COMP[b0] as u64, BASE4_COMP[b1] as u64);
            rhs[2 * x] = rh.rotate_right(1) ^ r0;
            rcs[2 * x] = (rcode >> 2) | g0 << top;
            rh = (rh ^ (r0.rotate_left(1) ^ r1.rotate_left(2))).rotate_right(2);
            rcode = (rcode >> 4) | g0 << (top - 2) | g1 << top;
            rhs[2 * x + 1] = rh;
            rcs[2 * x + 1] = rcode;
        }
    }
    if c % 2 == 1 {
        let (b, o) = (ins[c - 1] as usize, outs[c - 1] as usize);
        code = (code << 2) | BASE4[b] as u64;
        h = h.rotate_left(1) ^ (DNA_SEED[b] ^ f_out[o]);
        hs[c - 1] = h;
        cs[c - 1] = code;
        if BOTH {
            rcode = (rcode >> 2) | (BASE4_COMP[b] as u64) << top;
            rh = rh.rotate_right(1) ^ (r_in[b] ^ r_out[o]);
            rhs[c - 1] = rh;
            rcs[c - 1] = rcode;
        }
    }
    (h, rh, code, rcode)
}

/// Open-syncmer selection for a compile-time window of `W` s-mers with the
/// syncmer s-mer at offset `T`, block by block (`CHUNK` k-mers); with `BOTH` it
/// also returns the reverse complement's syncmers (else an empty `Vec`).
///
/// Per block, three simple loops replace one loop doing everything per k-mer:
/// 1. Roll the s-mer hashes and k-mer codes into flat buffers. Blocks of plain
///    A/C/G/T (checked up front; nearly all of them) skip the non-ACGT/'N'
///    tracking, and advance two bases per step so the serial rotate-xor and
///    shift-or chains take one link per two bases.
/// 2. Window test from a sliding minimum over P = 4 s-mers (for W=10, T=4: the
///    4 s-mers left of offset T are one run, the 5 to its right two overlapping
///    runs) — elementwise loops over flat buffers, which vectorise.
/// 3. Write out only the flagged k-mers (`for_each_flagged`).
///
/// Reverse complement: the RC s-mer over forward bases j..j+s hashes to
/// ⊕_r rotl(seed[comp(seq[j+r])], r), which rolls the other way (rotate right,
/// new base rotated by s−1); the RC k-mer code rolls as (code >> 2) | comp << 2(k−1).
/// RC window offset q is forward offset W−1−q, so leftmost-minimum-at-T in RC
/// order is: s-mer W−1−T strictly below the T s-mers to its right and ≤ the
/// W−1−T to its left. RC k-mer at forward i sits at RC position n−k−i; those are
/// found in descending RC order, so the RC output is reversed at the end.
/// Non-ACGT bytes: `revcomp` maps them to non-ACGT bytes ('n' to 'N'), so the
/// same s-mers are masked and the same k-mers rejected on both strands.
///
/// Requires n ≥ k (checked by the caller) and k − s + 1 == W.
fn select_syncmers_chunked<const W: usize, const T: usize, const BOTH: bool>(
    seq: &[u8], k: usize, s: usize,
) -> (Vec<SyncmerLight>, Vec<SyncmerLight>) {
    debug_assert!(W >= 2 && W <= RING && T < W && k - s + 1 == W && seq.len() >= k);
    let n = seq.len();
    let nk = n - k + 1;
    let kmask: u64 = if k == 32 { u64::MAX } else { (1u64 << (2 * k)) - 1 };
    let top = 2 * (k as u32 - 1);   // bit offset of the RC code's newest base
    // Left run (offsets 0..T) and right run (T+1..W) of the forward test; the RC
    // test has them swapped. Both are covered by one sliding minimum of P s-mers.
    let (lt, rt) = (T, W - 1 - T);
    let p_run = pow2_floor(match (lt, rt) { (0, r) => r.max(1), (l, 0) => l, (l, r) => l.min(r) });

    // Two syncmers are more than min(T, W−1−T) k-mers apart (the later one's
    // left run would hold the earlier one's s-mer and vice versa, needing
    // a < b ≤ a), which bounds the output size.
    let gap = lt.min(rt) + 1;
    let cap = nk.div_ceil(gap);
    let mut fwd: Vec<SyncmerLight> = Vec::with_capacity(cap);
    let mut rev: Vec<SyncmerLight> = Vec::with_capacity(if BOTH { cap } else { 0 });
    let fslots = &mut fwd.spare_capacity_mut()[..cap];
    let rslots = &mut rev.spare_capacity_mut()[..if BOTH { cap } else { 0 }];
    let (mut nf, mut nr) = (0usize, 0usize);

    let tabs = SeedTables::new::<BOTH>(s);
    let (f_out, r_in, r_out) = (&tabs.f_out, &tabs.r_in, &tabs.r_out);

    // hb[x] = (masked) hash of s-mer i0 + x: the W − 1 + c s-mers of the block's
    // k-mers i0..i0+c. cb[x] = 2-bit code of k-mer i0 + x (unmasked: bits above
    // 2k are dropped on output). r-prefixed = the same on the reverse complement
    // (RC of the same forward span). okb[x] = k-mer i0 + x holds no 'N'.
    let mut hb = [0u64; CHUNK + RING];
    let mut mb = [0u64; CHUNK + RING];
    let mut cb = [0u64; CHUNK];
    let mut fl = [0u8; CHUNK];
    let mut rhb = [0u64; CHUNK + RING];
    let mut rmb = [0u64; CHUNK + RING];
    let mut rcb = [0u64; CHUNK];
    let mut rfl = [0u8; CHUNK];
    let mut okb = [false; CHUNK];

    let (mut h, mut rh, mut code, mut rcode) = (0u64, 0u64, 0u64, 0u64);
    // One past the most recent byte that is non-ACGT (masks an s-mer hash) /
    // 'N' or 'n' (rejects a k-mer); 0 = none. Tracked only in blocks with such bytes.
    let mut bad_end: usize = 0;
    let mut n_end: usize = 0;

    // One base, with full non-ACGT tracking (warm-up and blocks with non-ACGT).
    macro_rules! roll_tracked { ($p:expr) => {{
        let p: usize = $p;
        let b = seq[p];
        let bo = if p >= s { seq[p - s] } else { b'N' };   // 'N': seed 0 = nothing out
        if DNA_SEED[b as usize] == 0 { bad_end = p + 1; }
        if b == b'N' || b == b'n' { n_end = p + 1; }
        code = (code << 2) | BASE4[b as usize] as u64;
        h = h.rotate_left(1) ^ (DNA_SEED[b as usize] ^ f_out[bo as usize]);
        if BOTH {
            rcode = (rcode >> 2) | (BASE4_COMP[b as usize] as u64) << top;
            rh = rh.rotate_right(1) ^ (r_in[b as usize] ^ r_out[bo as usize]);
        }
    }}}

    // Warm-up: bases 0..k-1 complete s-mers 0..W-1 (the first block's carry).
    for p in 0..k - 1 {
        roll_tracked!(p);
        if p + 1 >= s {
            let j = p + 1 - s;
            let bad = bad_end > j;
            hb[j] = if bad { u64::MAX } else { h };
            if BOTH { rhb[j] = if bad { u64::MAX } else { rh }; }
        }
    }

    let mut i0 = 0;
    while i0 < nk {
        let c = CHUNK.min(nk - i0);
        // The block's k-mers span seq[i0 .. i0 + c + k − 1]; k-mer (and its last
        // s-mer) x ends at base p0 + x.
        let p0 = i0 + k - 1;
        let clean = all_acgt(&seq[i0..p0 + c]);
        let ins = &seq[p0..p0 + c];
        let outs = &seq[p0 - s..p0 - s + c];

        // ── 1. hashes and codes ──
        if clean {
            let st = (h, rh, code, rcode);
            (h, rh, code, rcode) = roll_clean::<BOTH>(st, &tabs, top, ins, outs,
                &mut hb[W - 1..W - 1 + c], &mut cb[..c], &mut rhb[W - 1..W - 1 + c], &mut rcb[..c]);
        } else {
            for x in 0..c {
                let p = p0 + x;
                roll_tracked!(p);
                let bad = bad_end > p + 1 - s;   // s-mer p+1−s holds a non-ACGT byte
                hb[W - 1 + x] = if bad { u64::MAX } else { h };
                cb[x] = code;
                if BOTH {
                    rhb[W - 1 + x] = if bad { u64::MAX } else { rh };
                    rcb[x] = rcode;
                }
                okb[x] = n_end <= i0 + x;
            }
        }

        // ── 2. window tests ──
        let len = c + W - 1;
        sliding_min(&hb, &mut mb, len, p_run);
        for x in 0..c {
            let cv = hb[x + T];
            let lmin = run_min(&mb, x, lt, p_run);
            let rmin = run_min(&mb, x + T + 1, rt, p_run);
            fl[x] = ((T == 0 || cv < lmin) & (cv <= rmin)) as u8;
        }
        if !clean { for x in 0..c { fl[x] &= okb[x] as u8; } }
        if BOTH {
            sliding_min(&rhb, &mut rmb, len, p_run);
            for x in 0..c {
                let cv = rhb[x + rt];
                let lmin = run_min(&rmb, x, rt, p_run);
                let rmin = run_min(&rmb, x + rt + 1, lt, p_run);
                rfl[x] = ((cv <= lmin) & (T == 0 || cv < rmin)) as u8;
            }
            if !clean { for x in 0..c { rfl[x] &= okb[x] as u8; } }
        }

        // ── 3. output ── (flags past c, up to a multiple of 8, cleared)
        let c8 = c.next_multiple_of(8);
        fl[c..c8].fill(0);
        for_each_flagged(&fl[..c8], |x| {
            fslots[nf].write(SyncmerLight { pos: (i0 + x) as u32, value: cb[x] & kmask });
            nf += 1;
        });
        if BOTH {
            rfl[c..c8].fill(0);
            for_each_flagged(&rfl[..c8], |x| {
                rslots[nr].write(SyncmerLight { pos: (nk - 1 - (i0 + x)) as u32, value: rcb[x] });
                nr += 1;
            });
        }

        hb.copy_within(c..c + W - 1, 0);
        if BOTH { rhb.copy_within(c..c + W - 1, 0); }
        i0 += c;
    }
    // SAFETY: slots[..nf] / [..nr] were all initialised above.
    unsafe {
        fwd.set_len(nf);
        rev.set_len(nr);
    }
    rev.reverse();
    (fwd, rev)
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

    /// Random sequence of length < `max_n`: mostly ACGT, with a per-sequence
    /// density (none / sparse / dense) of lowercase, 'N'/'n', IUPAC and other
    /// bytes, plus homopolymer / dinucleotide runs (equal hashes → tie-breaking)
    /// and N runs.
    fn random_seq(rng: &mut Rng, iter: usize, max_n: u64) -> Vec<u8> {
        let alphabet = b"ACGTACGTacgtacgtNnRYKMSWryx-*";
        let n = (rng.next() % max_n) as usize;
        let odd_rate = [0u64, 0, 2000, 300, 50, 8][iter % 6];   // 1 in odd_rate
        let mut seq: Vec<u8> = (0..n).map(|_| {
            let r = rng.next();
            if odd_rate > 0 && r % odd_rate == 0 { alphabet[(r >> 8) as usize % alphabet.len()] }
            else { b"ACGT"[(r >> 8) as usize % 4] }
        }).collect();
        if iter % 5 == 0 && n > 40 {
            let a = (rng.next() as usize) % (n - 30);
            for q in a..a + 30 { seq[q] = if iter % 2 == 0 { b'A' } else { b"AC"[q % 2] }; }
        }
        if iter % 11 == 0 && n > 40 {
            let a = (rng.next() as usize) % (n - 30);
            for q in a..a + (rng.next() as usize % 30) { seq[q] = b'N'; }
        }
        if iter % 13 == 0 { seq.make_ascii_lowercase(); }
        seq
    }

    fn assert_same(a: &[SyncmerLight], b: &[SyncmerLight], what: &str) {
        assert_eq!(a.len(), b.len(), "count: {what}");
        for (x, y) in a.iter().zip(b) {
            assert_eq!((x.pos, x.value), (y.pos, y.value), "{what}");
        }
    }

    /// The block-wise kernel (both strands) must match the tracked-minimum
    /// version on the sequence and on its reverse complement, for window W / offset T.
    fn check_chunked<const W: usize, const T: usize>(seq: &[u8], k: usize) {
        let s = k + 1 - W;
        if seq.len() < k || s == 0 { return; }   // the dispatchers return early for n < k
        let rc = crate::fastq::revcomp(seq);
        let what = format!("W={W} T={T} k={k} n={}", seq.len());
        let (f, r) = select_syncmers_chunked::<W, T, true>(seq, k, s);
        let (f1, r1) = select_syncmers_chunked::<W, T, false>(seq, k, s);
        assert!(r1.is_empty());
        assert_same(&f, &select_syncmers_tracked(seq, k, s, T), &what);
        assert_same(&f1, &f, &what);
        assert_same(&r, &select_syncmers_tracked(&rc, k, s, T), &format!("RC {what}"));
    }

    /// Exhaustive randomized equality of the block-wise kernel with the reference,
    /// including N runs, other non-ACGT bytes, lowercase, low-complexity repeats,
    /// block boundaries, every k that yields W = 10, and other windows/offsets.
    #[test]
    fn chunked_matches_reference() {
        let mut rng = Rng(0xDEADBEEF12345678);
        for iter in 0..3000 {
            let seq = random_seq(&mut rng, iter, if iter % 4 == 0 { 2000 } else { 600 });
            let k = 10 + (rng.next() % 23) as usize;  // 10..=32
            check_chunked::<10, 4>(&seq, k);
            check_chunked::<10, 4>(&seq, 19);
            let k2 = 12 + (rng.next() % 21) as usize; // 12..=32
            match iter % 7 {
                0 => check_chunked::<10, 0>(&seq, k),
                1 => check_chunked::<10, 9>(&seq, k),
                2 => check_chunked::<2, 0>(&seq, k),
                3 => check_chunked::<3, 1>(&seq, k),
                4 => check_chunked::<6, 2>(&seq, k),
                5 => check_chunked::<12, 5>(&seq, k2),
                _ => check_chunked::<12, 3>(&seq, k2),
            }
        }
    }

    /// `select_syncmers_both` returns exactly the syncmers of the sequence and of
    /// its reverse complement, for the specialised and the general parameters.
    #[test]
    fn both_strands_match_revcomp() {
        let mut rng = Rng(0x0123456789ABCDEF);
        for iter in 0..2000 {
            let seq = random_seq(&mut rng, iter, 1200);
            let rc = crate::fastq::revcomp(&seq);
            for (k, s, t) in [(19, 10, 4), (15, 8, 4), (19, 10, 2), (5, 3, 1), (31, 22, 4)] {
                let (f, r) = select_syncmers_both(&seq, k, s, t);
                let what = format!("k={k} s={s} t={t} n={}", seq.len());
                assert_same(&f, &select_syncmers_light(&seq, k, s, t), &what);
                assert_same(&r, &select_syncmers_light(&rc, k, s, t), &format!("RC {what}"));
            }
        }
    }
}
