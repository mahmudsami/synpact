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
#[cfg(test)]
fn dct_three_coloring(vals: &[u64]) -> Vec<u64> {
    let mut col = Vec::new();
    dct_three_coloring_into(vals, &mut col);
    col
}

/// `dct_three_coloring` writing into a caller-owned buffer, so parsing many
/// runs reuses one allocation.
fn dct_three_coloring_into(vals: &[u64], col: &mut Vec<u64>) {
    let m = vals.len();
    debug_assert!(m >= 2);
    col.clear();

    // ── Phase 1: reduce to ≤ 6 colours ───────────────────────────────────────
    // Entry i is rewritten from itself and its not-yet-rewritten successor. The
    // first round reads `vals` directly instead of a copy of it, and the "any
    // colour ≥ 6" test that decides on another round is folded into the rewrite
    // (the tail colour is < 2, so it never needs testing) instead of rescanning.
    let round = |a: u64, b: u64| {
        let p = (a ^ b).trailing_zeros() as u64; // first differing bit
        2 * p + ((a >> p) & 1)
    };
    // Tail (no successor): any colour ≠ its new predecessor keeps it proper.
    let tail = |prev: u64| if prev == 0 { 1 } else { 0 };
    if vals.iter().any(|&c| c >= 6) {
        let mut wide = false;
        col.extend(vals.windows(2).map(|w| {
            let c = round(w[0], w[1]);
            wide |= c >= 6;
            c
        }));
        col.push(tail(col[m - 2]));
        while wide {
            wide = false;
            for i in 0..m - 1 {
                col[i] = round(col[i], col[i + 1]);
                wide |= col[i] >= 6;
            }
            col[m - 1] = tail(col[m - 2]);
        }
    } else {
        col.extend_from_slice(vals);
    }
    let col = &mut col[..];

    // ── Phase 2: collapse {3,4,5} → {0,1,2}, preserving adjacent-distinctness ──
    // Recolouring only produces colours < 3, so a colour absent now stays absent
    // and its pass can be skipped. Each pass is branch-free: every entry is
    // rewritten, with itself unless it has colour v. (Its neighbours never have
    // colour v then, so reading them mid-pass sees the same values as before.)
    let present = col.iter().fold(0u32, |a, &c| a | 1 << c);   // all colours < 6
    for v in [3, 4, 5] {
        if present & 1 << v == 0 { continue; }
        for i in 0..m {
            let l = if i > 0     { col[i - 1] } else { u64::MAX };
            let r = if i + 1 < m { col[i + 1] } else { u64::MAX };
            // Smallest colour in {0,1,2} differing from both neighbours.
            let c = if l != 0 && r != 0 { 0 } else if l != 1 && r != 1 { 1 } else { 2 };
            col[i] = if col[i] == v { c } else { col[i] };
        }
    }
}

/// Split a strictly-monotone (hence repetition-free) run of values `vals` into
/// blocks of size ≤ 3, deterministically and locally-consistently.
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
/// Each block is passed to `emit` as a half-open index range `start..end` into
/// `vals`.
#[cfg(test)]
pub(crate) fn split_monotone_run(vals: &[u64], emit: impl FnMut(usize, usize)) {
    split_monotone_run_buf(vals, &mut Vec::new(), emit)
}

/// `split_monotone_run` with a caller-owned colouring buffer.
fn split_monotone_run_buf(vals: &[u64], col: &mut Vec<u64>, mut emit: impl FnMut(usize, usize)) {
    let m = vals.len();
    if m <= 3 { return emit(0, m); }

    // Colour the run by its *content* (the values) so read and reference split
    // identically wherever the same run occurs.
    dct_three_coloring_into(vals, col);
    let col = &col[..];

    let mut triplet = |i: usize| emit(i.saturating_sub(1), (i + 1).min(m - 1) + 1);
    // Rule 1 — local-minimum triplet.
    for i in 0..m {
        if is_local_min(col, i) { triplet(i); }
    }
    // Rule 2 — local-maximum triplet, when no adjacent local minimum (that min's
    // own triplet already covers the maximum).
    for i in 0..m {
        if is_local_max(col, i) {
            let l_min = i > 0     && is_local_min(col, i - 1);
            let r_min = i + 1 < m && is_local_min(col, i + 1);
            if !l_min && !r_min { triplet(i); }
        }
    }
}

/// `dct_three_coloring` followed by `triplet_centres` for a run of `4 ≤ m ≤ N`
/// values (`N = 8` covers > 99.9% of the runs rule 4 splits on HiFi reads; most
/// have 4 or 5 units), on a fixed array of `N` lanes so that the rounds unroll into
/// independent register operations instead of a store→load chain through a
/// buffer. Lanes `≥ m` hold padding: phase 1 only ever reads a lane's successor,
/// so padding never flows into a live lane (lane `m-1` is overwritten by the tail
/// rule), and phase 2 and the centre test treat lane `m` as absent.
///
/// Phase 2 reads the colours from before each pass: an entry of colour v has no
/// neighbour of colour v, so its neighbours are not rewritten in that pass, and
/// every entry is independent of the others' rewrites — the same result as the
/// in-place sweep.
fn split_small<const N: usize>(vals: &[u64]) -> (u64, u64) {
    let m = vals.len();
    debug_assert!((4..=N).contains(&m) && N <= 64);
    let mut c = [0u64; N];
    for (i, x) in c.iter_mut().enumerate() { *x = vals[i.min(m - 1)]; }
    let live = (1u64 << (m - 1)) - 1;   // lanes 0..m-2: those with a successor

    let mut wide = vals.iter().any(|&x| x >= 6);
    while wide {
        let mut w = 0u64;
        for i in 0..N - 1 {
            let p = (c[i] ^ c[i + 1]).trailing_zeros();
            c[i] = 2 * p as u64 + ((c[i] >> (p & 63)) & 1);
            w |= ((c[i] >= 6) as u64) << i;
        }
        let tail = if c[m - 2] == 0 { 1 } else { 0 };
        for (i, x) in c.iter_mut().enumerate() { if i == m - 1 { *x = tail; } }
        wide = w & live != 0;
    }

    for v in [3, 4, 5] {
        let old = c;
        for i in 0..N {
            let l = if i > 0     { old[i - 1] } else { u64::MAX };
            let r = if i + 1 < m { old[i + 1] } else { u64::MAX };
            let nc = if l != 0 && r != 0 { 0 } else if l != 1 && r != 1 { 1 } else { 2 };
            if old[i] == v { c[i] = nc; }
        }
    }

    // `triplet_centres` on lanes 0..m.
    let (mut up, mut dn) = (0u64, 0u64);
    for j in 0..N - 1 {
        up |= ((c[j] < c[j + 1]) as u64) << j;
        dn |= ((c[j] > c[j + 1]) as u64) << j;
    }
    let (up, dn) = ((up & live) | 1 << (m - 1), (dn & live) | 1 << (m - 1));
    let mins = ((dn << 1) | 1) & up;
    let maxs = ((up << 1) | 1) & dn;
    (mins, maxs & !(mins << 1) & !(mins >> 1))
}

/// Rules 1–2 of `split_monotone_run_buf` on a colouring of `2 ≤ m ≤ 64` units,
/// as bitmasks: (local minima, local maxima with no adjacent local minimum).
/// The same slope-bit formulation as `locally_consistent_parsing`'s rules 1–2,
/// with the ends always beating their missing neighbour.
fn triplet_centres(col: &[u64]) -> (u64, u64) {
    let m = col.len();
    debug_assert!((2..=64).contains(&m));
    let (mut up, mut dn) = (0u64, 0u64);
    for j in 0..m - 1 {
        up |= ((col[j] < col[j + 1]) as u64) << j;
        dn |= ((col[j] > col[j + 1]) as u64) << j;
    }
    up |= 1 << (m - 1);
    dn |= 1 << (m - 1);
    let mins = ((dn << 1) | 1) & up;
    let maxs = ((up << 1) | 1) & dn;
    (mins, maxs & !(mins << 1) & !(mins >> 1))
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

/// A bitset over unit positions `0..n` (one bit per unit, 64 per word), used for
/// the parser's "assigned" set and its "equal to successor" set.  Range tests and
/// fills touch one or two words, and the rule-4 search for the next unassigned
/// unit skips 64 assigned units per step instead of testing them one by one.
struct Bits<'a>(&'a mut [u64]);

impl Bits<'_> {
    /// Mask of bits `lo..=hi` of one word (`lo ≤ hi < 64`).
    #[inline]
    fn mask(lo: usize, hi: usize) -> u64 { (u64::MAX << lo) & (u64::MAX >> (63 - hi)) }

    /// Whether any bit in `lo..=hi` is clear.
    #[inline]
    fn any_clear(&self, lo: usize, hi: usize) -> bool {
        (lo >> 6..=hi >> 6).any(|w| {
            let m = Self::mask(if w == lo >> 6 { lo & 63 } else { 0 },
                               if w == hi >> 6 { hi & 63 } else { 63 });
            self.0[w] & m != m
        })
    }

    /// Set every bit in `lo..=hi`.
    #[inline]
    fn set_range(&mut self, lo: usize, hi: usize) {
        for w in lo >> 6..=hi >> 6 {
            self.0[w] |= Self::mask(if w == lo >> 6 { lo & 63 } else { 0 },
                                    if w == hi >> 6 { hi & 63 } else { 63 });
        }
    }

    /// Bits `lo..lo + len` (`1 ≤ len ≤ 64`) as the low bits of a word.
    #[inline]
    fn get(&self, lo: usize, len: usize) -> u64 {
        let (w, o) = (lo >> 6, lo & 63);
        let mut x = self.0[w] >> o;
        if o + len > 64 { x |= self.0[w + 1] << (64 - o); }
        x & (u64::MAX >> (64 - len))
    }

    /// OR the low `len` bits of `x` (`1 ≤ len ≤ 64`, higher bits clear) into bits
    /// `lo..lo + len`.
    #[inline]
    fn or(&mut self, lo: usize, len: usize, x: u64) {
        let (w, o) = (lo >> 6, lo & 63);
        self.0[w] |= x << o;
        if o + len > 64 { self.0[w + 1] |= x >> (64 - o); }
    }

    /// Index of the last bit `< before` equal to `bit`, if any.
    #[inline]
    fn prev(&self, before: usize, bit: bool) -> Option<usize> {
        let flip = if bit { 0 } else { u64::MAX };
        if before == 0 { return None; }
        let last = before - 1;
        let mut w = last >> 6;
        let mut word = (self.0[w] ^ flip) & (u64::MAX >> (63 - (last & 63)));
        loop {
            if word != 0 { return Some(w * 64 + 63 - word.leading_zeros() as usize); }
            if w == 0 { return None; }
            w -= 1;
            word = self.0[w] ^ flip;
        }
    }

    /// Index of the first bit `≥ from` equal to `bit`, if any.
    #[inline]
    fn next(&self, from: usize, bit: bool) -> Option<usize> {
        let flip = if bit { 0 } else { u64::MAX };
        let mut w = from >> 6;
        let mut word = (*self.0.get(w)? ^ flip) & (u64::MAX << (from & 63));
        loop {
            if word != 0 { return Some(w * 64 + word.trailing_zeros() as usize); }
            w += 1;
            word = *self.0.get(w)? ^ flip;
        }
    }
}

/// Locally-consistent parsing of one level: partition-with-overlaps of `values`
/// into blocks, by four rules applied in order (a block of a later rule is kept
/// only if it covers a unit no earlier block covers):
///   1. a triplet `{i-1, i, i+1}` around every local minimum;
///   2. a triplet around every local maximum with no adjacent local minimum;
///   3. every run of ≥ 2 equal values plus its two neighbours;
///   4. every maximal strictly-monotone run still holding an uncovered unit,
///      split into ≤ 3-unit blocks by `split_monotone_run_buf`.
/// The blocks are returned in order of their first unit (stable: rule order on
/// ties, rule-4 blocks of one run sorted by start).
///
/// Rules 1–2 cover all but a few percent of units, so they are evaluated 64 units
/// at a time on bitsets rather than unit by unit with data-dependent branches.
/// With the slope bits `up[i] = v[i] < v[i+1]` and `dn[i] = v[i] > v[i+1]`, a
/// local minimum is `dn[i-1] ∧ up[i]` and a local maximum `up[i-1] ∧ dn[i]`; since
/// a maximum's neighbours already slope away from it, "maximum with no adjacent
/// minimum" reduces to `¬dn[i-2] ∧ ¬up[i+1]`, so the triplet centres are
///   centre(i) = (dn[i-1] ∧ up[i]) ∨ (up[i-1] ∧ dn[i] ∧ ¬dn[i-2] ∧ ¬up[i+1]),
/// with sentinels at the ends: `up[-1] = dn[-1] = true` and `dn[-2] = false`
/// (unit 0 has no left neighbour to beat, and there is no unit −1 to be a minimum),
/// `up[n-1] = dn[n-1] = true` and `up[n] = false` likewise on the right.  It is the
/// same predicate as `is_local_min(i) ∨ (is_local_max(i) ∧ ¬is_local_min(i±1))`
/// (checked against the original implementation by the randomized test below).
/// No min/max triplet covers another triplet's centre, so every centre's triplet
/// is kept, and the units within one of a centre form the assigned set that rules
/// 3 and 4 then consult.
///
/// Every kept block covers at least one new unit, so there are at most `n` blocks
/// in all. The min/max triplets are emitted last, straight from the centre bits,
/// interleaved with the (few) rule 3–4 blocks.
pub(crate) fn locally_consistent_parsing(values: &[u64]) -> Vec<Block> {
    let n = values.len();
    if n == 0 { return vec![]; }
    if n == 1 { return vec![Block { start: 0, end: 1 }]; }
    let v = values;

    // ── Rules 1–2 on 64-unit words ───────────────────────────────────────────
    // Slope bits first (bits at positions ≥ n stay 0), then the centre bits and
    // the equal-successor bits (`¬up ∧ ¬dn`). The slope bits are kept: rule 4 finds
    // its monotone runs' ends by scanning them.
    let words = n.div_ceil(64);
    let mut bit_buf = vec![0u64; 5 * words];
    let (up, rest) = bit_buf.split_at_mut(words);
    let (dn, rest) = rest.split_at_mut(words);
    let (ctr, rest) = rest.split_at_mut(words);
    let (eq, asg) = rest.split_at_mut(words);
    let tail = (n - 1) / 64;   // the word holding slope n-1 (the sentinel)
    for w in 0..tail {
        // All 64 slopes of this word are real (64w + 64 ≤ n - 1).
        let c = &v[64 * w..64 * w + 65];
        let (mut u, mut d) = (0u64, 0u64);
        for j in 0..64 {
            u |= ((c[j] < c[j + 1]) as u64) << j;
            d |= ((c[j] > c[j + 1]) as u64) << j;
        }
        (up[w], dn[w]) = (u, d);
    }
    let (mut u, mut d) = (0u64, 0u64);
    for i in 64 * tail..n - 1 {
        u |= ((v[i] < v[i + 1]) as u64) << (i & 63);
        d |= ((v[i] > v[i + 1]) as u64) << (i & 63);
    }
    let last = 1u64 << ((n - 1) & 63);
    (up[tail], dn[tail]) = (u | last, d | last);

    // Centres. `*prev` hold the previous word's slope bits; the initial values put
    // the left sentinels (up[-1] = dn[-1] = 1, dn[-2] = 0) in bits 63 / 62.
    let (mut uprev, mut dprev) = (1u64 << 63, 1u64 << 63);
    for w in 0..words {
        let (u, d) = (up[w], dn[w]);
        let unext = if w + 1 < words { up[w + 1] } else { 0 };   // up[n] = 0
        let dn1 = (d << 1) | (dprev >> 63);   // dn[i-1]
        let dn2 = (d << 2) | (dprev >> 62);   // dn[i-2]
        let up1 = (u << 1) | (uprev >> 63);   // up[i-1]
        let upn = (u >> 1) | (unext << 63);   // up[i+1]
        ctr[w] = (dn1 & u) | (up1 & d & !dn2 & !upn);
        let valid = if w == tail { u64::MAX >> (63 - ((n - 1) & 63)) } else { u64::MAX };
        eq[w] = !(u | d) & valid;
        (uprev, dprev) = (u, d);
    }

    // Assigned = every unit within one of a centre (no min/max triplet covers
    // another's centre, so every centre's triplet is kept). Bits past n are set, so
    // the rule-4 search never runs off the end.
    let mut carry = 0u64;   // top centre bit of the previous word
    for w in 0..words {
        let c = ctr[w];
        let next = if w + 1 < words { ctr[w + 1] << 63 } else { 0 };
        asg[w] = c | (c << 1) | carry | (c >> 1) | next;
        carry = c >> 63;
    }
    if n & 63 != 0 { asg[words - 1] |= u64::MAX << (n & 63); }
    let mut assigned = Bits(asg);
    let eq = Bits(eq);

    // Commit the contiguous inclusive range `lo..=hi` (clamped to the sequence) as
    // one block of `$list`, but only if it contributes at least one not-yet-assigned
    // position (shared boundaries between adjacent blocks are intentional
    // intersections).
    macro_rules! commit_range {
        ($list:ident, $lo:expr, $hi:expr) => {{
            let lo = $lo;
            let hi = ($hi).min(n - 1);
            if assigned.any_clear(lo, hi) {
                assigned.set_range(lo, hi);
                $list.push(Block { start: lo as u32, end: hi as u32 + 1 });
            }
        }};
    }

    // ── Rule 3 — repetition run ≥ 2 plus immediate prev/next neighbours ──────
    // A maximal run of set equal-successor bits st..=e is the value run st..=e+1.
    let mut reps: Vec<Block> = Vec::new();
    let mut i = 0;
    while let Some(st) = eq.next(i, true) {
        // eq[n-1] is clear (no successor), so the run always ends.
        let j = eq.next(st, false).unwrap() + 1;   // one past the value run
        // The block runs to j, the unit after the run; the macro clamps it to n-1.
        commit_range!(reps, st.saturating_sub(1), j);
        i = j;
    }

    // ── Rule 4 — full monotone run [s..=b], anchors included (intentional
    // intersection). A position still unassigned is neither an extremum nor next
    // to a repeat, so it lies inside a maximal strictly monotone run whose ends
    // (extrema, repeat neighbours or sequence ends) are assigned. Split such runs
    // into ≤3-unit blocks via DCT so they aren't one big featureless block
    // (deterministic & locally consistent → genome and read split identically).
    //
    // A run of ≤ 64 units (all but pathological ones) is parsed on a local copy of
    // its assigned bits: a block is kept iff its unit mask is not already a subset,
    // the copy is OR-ed back into the bitset afterwards, and the kept blocks are
    // emitted by centre, i.e. by start. Longer runs take the general path.
    let (up, dn) = (Bits(up), Bits(dn));
    let mut mono_blocks: Vec<Block> = Vec::with_capacity(n / 16 + 4);
    let mut col_buf: Vec<u64> = Vec::new();
    let mut i = 0;
    while let Some(u) = assigned.next(i, false) {
        // The run extends over the slopes in u's direction: `s` is just past the
        // last slope before u that is not, `b` the first slope from u on that is
        // not (the sentinel slope n-1 is set in both directions, hence the cap).
        let dir = if v[u] < v[u + 1] { &up } else { &dn };
        let s = dir.prev(u, false).map_or(0, |j| j + 1);
        let b = dir.next(u, false).map_or(n - 1, |j| j.min(n - 1));
        let m = b - s + 1;
        if m <= 64 {
            let full = u64::MAX >> (64 - m);
            let mut a = assigned.get(s, m);
            if m <= 3 {
                // One block, the whole run (as `split_monotone_run_buf`).
                if a != full { mono_blocks.push(Block { start: s as u32, end: b as u32 + 1 }); }
                a = full;
            } else {
                let (mins, lone_maxs) = if m <= 8 {
                    split_small::<8>(&v[s..=b])
                } else {
                    dct_three_coloring_into(&v[s..=b], &mut col_buf);
                    triplet_centres(&col_buf)
                };
                // Commit order as in `split_monotone_run_buf`: all minima, then
                // the lone maxima.
                let mut kept = 0u64;
                for centres in [mins, lone_maxs] {
                    let mut x = centres;
                    while x != 0 {
                        let c = x.trailing_zeros();
                        let r = ((0b111u128 << c >> 1) as u64) & full;   // c-1..=c+1
                        kept |= ((a & r != r) as u64) << c;
                        a |= r;
                        x &= x - 1;
                    }
                }
                while kept != 0 {
                    let c = kept.trailing_zeros() as usize;
                    mono_blocks.push(Block { start: (s + c.saturating_sub(1)) as u32,
                                             end:   (s + (c + 1).min(m - 1) + 1) as u32 });
                    kept &= kept - 1;
                }
            }
            assigned.or(s, m, a);
        } else {
            let first = mono_blocks.len();
            split_monotone_run_buf(&v[s..=b], &mut col_buf,
                                   |lo, hi| commit_range!(mono_blocks, s + lo, s + hi - 1));
            // A run emits its min triplets before its max triplets; order its
            // blocks by start (`sort_by_key` is stable, as the original global
            // sort was).
            mono_blocks[first..].sort_by_key(|b| b.start);
        }
        // Runs are visited left to right and a later run starts no earlier than
        // where the previous one ended, past every block start of the previous
        // run, so `mono_blocks` stays in start order — exactly as a stable sort by
        // start of the blocks in commit order.
        i = b;
    }
    debug_assert!(assigned.next(0, false).is_none(), "LCP bug: unassigned position\nvalues={values:?}");

    // ── Start-ordered result ─────────────────────────────────────────────────
    // Equivalent to a stable sort by `start` of rules-1/2 ++ reps ++ mono_blocks,
    // all three already in start order (each rule scans upwards; see above for
    // `mono_blocks`). The (few) rule 3–4 blocks are merged stably — rule order on
    // equal starts — and the min/max triplets are then emitted from the centre
    // bits with each rule 3–4 block slotted in after the triplets starting no
    // later than it.
    debug_assert!(mono_blocks.windows(2).all(|w| w[0].start <= w[1].start));
    let merged;
    let extra: &[Block] = match (reps.is_empty(), mono_blocks.is_empty()) {
        (true, false) => &mono_blocks,
        (false, true) => &reps,
        _ => { merged = merge_by_start(&reps, &mono_blocks); &merged }
    };
    let mut out: Vec<Block> = Vec::with_capacity(n);
    // Centre bits still to emit: word `w`, remaining bits `x`.
    let (mut w, mut x) = (0usize, ctr[0]);
    // Emit the triplets of all centres `< limit`.
    let mut emit_until = |out: &mut Vec<Block>, limit: usize| loop {
        while x != 0 {
            let i = 64 * w + x.trailing_zeros() as usize;
            if i >= limit { return; }
            out.push(Block { start: i.saturating_sub(1) as u32, end: (i + 2).min(n) as u32 });
            x &= x - 1;
        }
        w += 1;
        if w >= words { return; }
        x = ctr[w];
    };
    for e in extra {
        // Triplet start max(c-1, 0) ≤ e.start  ⇔  c ≤ e.start + 1.
        emit_until(&mut out, e.start as usize + 2);
        out.push(*e);
    }
    emit_until(&mut out, n);
    debug_assert!(out.len() <= n);
    out
}

/// Stable merge of two start-sorted block lists; on equal starts `a` precedes `b`.
fn merge_by_start(a: &[Block], b: &[Block]) -> Vec<Block> {
    let mut out = Vec::with_capacity(a.len() + b.len());
    let (mut i, mut j) = (0, 0);
    while i < a.len() && j < b.len() {
        if a[i].start <= b[j].start { out.push(a[i]); i += 1; }
        else                        { out.push(b[j]); j += 1; }
    }
    out.extend_from_slice(&a[i..]);
    out.extend_from_slice(&b[j..]);
    out
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
/// `children` is the contiguous index range of this block's children in the
/// level below (`forest.levels[level-1]`); LCP blocks always cover a contiguous
/// run, and adjacent parents may share a boundary child. Empty for L1.
pub(crate) struct HierNode {
    pub(crate) level:    usize,
    pub(crate) hash:     u64,
    pub(crate) pos:      u32,   // query start (inclusive)
    pub(crate) end:      u32,   // query end   (exclusive); span = end - pos
    /// Conservation-of-mass weight: every syncmer carries mass 1; a block's mass
    /// is the sum of its units' masses, with a unit shared by `m` blocks
    /// contributing 1/m to each.  Total mass is conserved across all levels.
    pub(crate) mass:     f32,
    pub(crate) children: Block,  // index range into the level below
}

/// The full per-read block hierarchy, stored flat: `levels[L]` is every block at
/// level L (1-indexed; `levels[0]` is an unused placeholder). The top-level blocks
/// to anchor from are `levels[top_level]`.
pub(crate) struct HierForest {
    pub(crate) levels:    Vec<Vec<HierNode>>,
    pub(crate) top_level: usize,
}

/// Build the per-read block hierarchy as a flat `HierForest`.
/// num_levels is the highest level to anchor from: the lower of `--max-level` and
/// the loaded index's top level (`GIndex::max_level`). No higher level is built,
/// since its blocks would never be looked up.
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
        HierNode { level: 1, hash: h, pos, end, mass, children: Block { start: 0, end: 0 } }
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

            let cur_hashes: Vec<u64> = cur.iter().map(|b| b.hash).collect();
            let next_raw = locally_consistent_parsing(&cur_hashes);
            if next_raw.is_empty() { break; }

            // Membership: how many parent blocks each child block belongs to.
            let mut child_membership = vec![0u32; cur_hashes.len()];
            for blk in &next_raw { for i in blk.range() { child_membership[i] += 1; } }

            let new_level: Vec<HierNode> = next_raw.iter().map(|blk| {
                let h        = block_hash_for_level(&cur_hashes[blk.range()], level_1idx - 1);
                let pos      = cur[blk.first()].pos;
                let end      = cur[blk.last()].end;
                // Parent mass = sum of each child's mass / (#parents sharing that child).
                let mass: f32 = blk.range()
                    .map(|i| cur[i].mass / child_membership[i] as f32).sum();
                HierNode { level: level_1idx, hash: h, pos, end, mass, children: *blk }
            }).collect();

            levels.push(new_level);
            top_level = level_1idx;
        }
    });

    HierForest { levels, top_level }
}

// ─────────────────────────────────────────────────────────────────────────────
// 11.  Genome index (sorted flat array + bucket table, see index.rs)
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

            // Strictly monotone values.
            let mut values = Vec::with_capacity(m);
            let mut acc: u64 = (rng.next() % 1000) + 1;
            for _ in 0..m {
                values.push(acc);
                acc += 1 + (rng.next() % 50); // strictly increasing steps
            }
            if decreasing { values.reverse(); } // strictly decreasing

            let split = || { let mut v = vec![]; split_monotone_run(&values, |s, e| v.push((s, e))); v };
            let blocks = split();
            assert!(!blocks.is_empty());

            // Determinism.
            assert_eq!(blocks, split());

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

/// Verbatim copy of the pre-optimisation parser (the bool-array / repeated
/// `is_local_min` implementation), kept as the oracle that the randomized
/// equivalence test checks `locally_consistent_parsing` against block-for-block.
#[cfg(test)]
mod reference {
    use super::Block;

    fn dct_three_coloring_into(vals: &[u64], col: &mut Vec<u64>) {
        let m = vals.len();
        debug_assert!(m >= 2);
        col.clear();
        col.extend_from_slice(vals);

        // ── Phase 1: reduce to ≤ 6 colours ───────────────────────────────────────
        // In place: entry i is rewritten from itself and its not-yet-rewritten successor.
        while col.iter().any(|&c| c >= 6) {
            for i in 0..m - 1 {
                let p = (col[i] ^ col[i + 1]).trailing_zeros() as u64; // first differing bit
                col[i] = 2 * p + ((col[i] >> p) & 1);
            }
            // Tail (no successor): any colour ≠ its new predecessor keeps it proper.
            col[m - 1] = if col[m - 2] == 0 { 1 } else { 0 };
        }

        // ── Phase 2: collapse {3,4,5} → {0,1,2}, preserving adjacent-distinctness ──
        for v in [3, 4, 5] {
            for i in 0..m {
                if col[i] != v { continue; }
                let l = if i > 0     { col[i - 1] } else { u64::MAX };
                let r = if i + 1 < m { col[i + 1] } else { u64::MAX };
                col[i] = (0..=2).find(|&c| c != l && c != r).unwrap();
            }
        }
    }

    fn split_monotone_run_buf(vals: &[u64], col: &mut Vec<u64>, mut emit: impl FnMut(usize, usize)) {
        let m = vals.len();
        if m <= 3 { return emit(0, m); }

        // Colour the run by its *content* (the values) so read and reference split
        // identically wherever the same run occurs.
        dct_three_coloring_into(vals, col);
        let col = &col[..];

        let mut triplet = |i: usize| emit(i.saturating_sub(1), (i + 1).min(m - 1) + 1);
        // Rule 1 — local-minimum triplet.
        for i in 0..m {
            if is_local_min(col, i) { triplet(i); }
        }
        // Rule 2 — local-maximum triplet, when no adjacent local minimum (that min's
        // own triplet already covers the maximum).
        for i in 0..m {
            if is_local_max(col, i) {
                let l_min = i > 0     && is_local_min(col, i - 1);
                let r_min = i + 1 < m && is_local_min(col, i + 1);
                if !l_min && !r_min { triplet(i); }
            }
        }
    }

    #[inline]
    pub(super) fn is_local_min(v: &[u64], i: usize) -> bool {
        let n = v.len();
        (i == 0 || v[i] < v[i - 1]) && (i == n - 1 || v[i] < v[i + 1])
    }

    #[inline]
    pub(super) fn is_local_max(v: &[u64], i: usize) -> bool {
        let n = v.len();
        (i == 0 || v[i] > v[i - 1]) && (i == n - 1 || v[i] > v[i + 1])
    }

    pub(super) fn locally_consistent_parsing(values: &[u64]) -> Vec<Block> {
        let n = values.len();
        if n == 0 { return vec![]; }
        if n == 1 { return vec![Block { start: 0, end: 1 }]; }

        let mut assigned = vec![false; n];
        // Blocks of rules 1–2, rule 3 and rule 4 are collected separately: the first
        // two lists come out in start order, so the final start-ordered list is a
        // merge rather than a sort of everything (see the end of this function).
        let mut minmax: Vec<Block> = Vec::with_capacity(n / 2 + 1);
        let mut reps:   Vec<Block> = Vec::new();
        let mut mono_blocks: Vec<Block> = Vec::new();

        // Commit the contiguous inclusive range `lo..=hi` (clamped to the sequence) as
        // one block of `$list`, but only if it contributes at least one not-yet-assigned
        // position (shared boundaries between adjacent blocks are intentional
        // intersections).
        macro_rules! commit_range {
            ($list:ident, $lo:expr, $hi:expr) => {{
                let lo = $lo;
                let hi = ($hi).min(n - 1);
                if assigned[lo..=hi].contains(&false) {
                    assigned[lo..=hi].fill(true);
                    $list.push(Block { start: lo as u32, end: hi as u32 + 1 });
                }
            }};
        }

        // Rule 1 — local-minimum block: always full triplet {i-1, i, i+1}.
        // Adjacent minima each get their full triplet (shared boundary appears in
        // both — intentional intersection).
        // Rule 2 — local-maximum block (no adjacent local minimum), {i-1, i, i+1}.
        // No other min/max triplet covers a triplet's centre, so each always commits.
        let is_min = |i: usize| is_local_min(values, i);
        for i in 0..n {
            let lone_max = is_local_max(values, i)
                && !(i > 0 && is_min(i - 1)) && !(i + 1 < n && is_min(i + 1));
            if is_min(i) || lone_max {
                commit_range!(minmax, i.saturating_sub(1), i + 1);
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
                    commit_range!(reps, i.saturating_sub(1), j);
                }
                i = j;
            }
        }

        // Rule 4 — full monotone run [s..=b], anchors included (intentional intersection).
        // A position still unassigned is neither an extremum nor next to a repeat, so it
        // lies inside a maximal strictly monotone run whose ends (extrema, repeat
        // neighbours or sequence ends) are assigned. Split such runs into ≤3-unit blocks
        // via DCT so they aren't one big featureless block (deterministic & locally
        // consistent → genome and read split identically).
        let mut col_buf: Vec<u64> = Vec::new();
        let mut i = 0;
        while let Some(u) = assigned[i..].iter().position(|&x| !x).map(|o| o + i) {
            let up = values[u] < values[u + 1];
            let mono = |j: usize| values[j] != values[j + 1] && (values[j] < values[j + 1]) == up;
            let mut s = u;
            while s > 0 && mono(s - 1) { s -= 1; }
            let mut b = u;
            while b + 1 < n && mono(b) { b += 1; }
            split_monotone_run_buf(&values[s..=b], &mut col_buf,
                                   |lo, hi| commit_range!(mono_blocks, s + lo, s + hi - 1));
            i = b;
        }
        debug_assert!(!assigned.contains(&false), "LCP bug: unassigned position\nvalues={values:?}");

        // Equivalent to a stable sort by `start` of minmax ++ reps ++ mono_blocks:
        // `minmax` and `reps` are already in start order (each rule scans i upwards);
        // `mono_blocks` is not (a run emits its min triplets before its max
        // triplets), so it is stably sorted on its own. A stable 3-way merge preferring
        // the earlier list on equal starts then reproduces that sort exactly.
        debug_assert!(minmax.windows(2).all(|w| w[0].start <= w[1].start));
        debug_assert!(reps.windows(2).all(|w| w[0].start <= w[1].start));
        mono_blocks.sort_by_key(|b| b.start);
        merge3_by_start(&minmax, &reps, &mono_blocks)
    }

    fn merge3_by_start(a: &[Block], b: &[Block], c: &[Block]) -> Vec<Block> {
        let mut out = Vec::with_capacity(a.len() + b.len() + c.len());
        let (mut i, mut j, mut l) = (0, 0, 0);
        loop {
            let ka = a.get(i).map_or(u32::MAX, |x| x.start);
            let kb = b.get(j).map_or(u32::MAX, |x| x.start);
            let kc = c.get(l).map_or(u32::MAX, |x| x.start);
            if i < a.len() && ka <= kb && ka <= kc { out.push(a[i]); i += 1; }
            else if j < b.len() && kb <= kc       { out.push(b[j]); j += 1; }
            else if l < c.len()                   { out.push(c[l]); l += 1; }
            else { break; }
        }
        out
    }
}

#[cfg(test)]
mod lcp_equiv_tests {
    use super::{locally_consistent_parsing, reference};

    struct Rng(u64);
    impl Rng {
        fn next(&mut self) -> u64 {
            let mut x = self.0;
            x ^= x << 13; x ^= x >> 7; x ^= x << 17;
            self.0 = x; x
        }
        fn below(&mut self, m: u64) -> u64 { self.next() % m }
    }

    /// A value array of one of several shapes: uniform random, tiny alphabets
    /// (dense repeats and plateaus), long monotone runs (rule 4) with occasional
    /// repeats, and spliced mixtures of all of these.
    fn random_values(rng: &mut Rng) -> Vec<u64> {
        let n = match rng.below(10) {
            0..=3 => rng.below(12) as usize,               // tiny, incl. 0 and 1
            4..=7 => rng.below(200) as usize,
            _     => 60 + rng.below(700) as usize,         // several bitset words
        };
        let mut v = Vec::with_capacity(n);
        while v.len() < n {
            let len = 1 + rng.below(40) as usize;
            match rng.below(6) {
                0 => for _ in 0..len { v.push(rng.next()); },
                5 => {
                    // Long strictly monotone run: the > 8 and > 64 rule-4 paths.
                    let up = rng.next() & 1 == 0;
                    let mut x = rng.below(1 << 20) + (1 << 30);
                    for _ in 0..len * rng.below(6) as usize {
                        v.push(x);
                        let step = 1 + rng.below(1 << 12);
                        x = if up { x + step } else { x - step };
                    }
                }
                1 => { let a = 1 + rng.below(5); for _ in 0..len { v.push(rng.below(a)); } }
                2 | 3 => {
                    // Strictly monotone run, sometimes with a repeat inside.
                    let up = rng.next() & 1 == 0;
                    let mut x = rng.below(1 << 20) + (1 << 30);
                    for _ in 0..len {
                        v.push(x);
                        let step = if rng.below(12) == 0 { 0 } else { 1 + rng.below(3) };
                        x = if up { x + step } else { x - step };
                    }
                }
                _ => { let x = rng.below(4); for _ in 0..len { v.push(x); } }
            }
        }
        v.truncate(n);
        v
    }

    /// The optimised parser must emit exactly the blocks, in exactly the order, of
    /// the original implementation.
    #[test]
    fn lcp_matches_reference() {
        let mut rng = Rng(0xD1B5_4A32_D192_ED03);
        for it in 0..200_000 {
            let v = random_values(&mut rng);
            let got = locally_consistent_parsing(&v);
            let want = reference::locally_consistent_parsing(&v);
            assert_eq!(got, want, "iteration {it}, values={v:?}");
        }
    }
}
