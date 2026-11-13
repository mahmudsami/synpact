// ─────────────────────────────────────────────────────────────────────────────
// 1.  Forward-only DNA rolling hash  (replaces the nthash crate)
//
//  Scheme — rotation-XOR, same structure as NT-hash's forward component but
//  with custom seeds and no reverse-complement component:
//
//    H(i) = seed[seq[i]]   rotated (w-1)
//         ^ seed[seq[i+1]] rotated (w-2)
//         ^ …
//         ^ seed[seq[i+w-1]]   (no rotation)
//
//  Rolling O(1) update:
//    H(i+1) = rotate_left(H(i), 1)
//            ^ rotate_left(seed[seq[i]], w)   ← removes outgoing base
//            ^ seed[seq[i+w]]                 ← adds incoming base
//
//  Seeds: the four SplitMix64/φ constants used throughout the codebase.
//  They have excellent 64-bit avalanche (≈50% bit-flip per input bit) and
//  are well-separated in Hamming distance from each other.
// ─────────────────────────────────────────────────────────────────────────────

/// Per-base seeds for the rolling hash. Index by raw byte value.
/// Non-ACGT bytes (including 'N') map to 0.
pub(crate) const DNA_SEED: [u64; 256] = {
    let mut t = [0u64; 256];
    t[b'A' as usize] = 0x9e3779b97f4a7c15;  // floor(2^64 / φ)
    t[b'a' as usize] = 0x9e3779b97f4a7c15;
    t[b'C' as usize] = 0x6c62272e07bb0142;  // FNV-1a basis
    t[b'c' as usize] = 0x6c62272e07bb0142;
    t[b'G' as usize] = 0xbf58476d1ce4e5b9;  // SplitMix64 multiplier 1
    t[b'g' as usize] = 0xbf58476d1ce4e5b9;
    t[b'T' as usize] = 0x94d049bb133111eb;  // SplitMix64 multiplier 2
    t[b't' as usize] = 0x94d049bb133111eb;
    t
};

// ─────────────────────────────────────────────────────────────────────────────
// 3.  Level-1 block hash
//
// A block is represented by its ordered sequence of base-4 s-mer values.
// Two blocks with the same s-mer sequence (anywhere in the genome) get the
// same hash, making the hash a canonical block identifier.
//
// Design to minimise collisions:
//   • Length is folded in first → blocks of different sizes can't collide
//     unless the hash itself collides.
//   • Each value is mixed via a full SplitMix64 step → ~50% bit-flip on any
//     single-bit change (avalanche).
//   • Position-sensitivity: the rotation before adding each value means the
//     same values in different order produce different hashes.
//
// 64-bit output.  Expected collisions for N=10^8 distinct blocks ≈ N²/2⁶⁴ ≈ 0.5%.
// ─────────────────────────────────────────────────────────────────────────────

// Domain constants keep every level's hash space disjoint.
// Index 0 → L1, index 1 → L2, … up to L8; beyond that we derive more.
pub(crate) const LEVEL_DOMAINS: [u64; 8] = [
    0x9e3779b97f4a7c15, // L1  (φ-based)
    0x6c62272e07bb0142, // L2  (FNV prime)
    0xd2a98b26625eee7b, // L3
    0xa3b195354a2b7623, // L4
    0x1b03738712fad5c9, // L5
    0xc4ceb9fe1a85ec53, // L6
    0x517cc1b727220a95, // L7
    0x9be1aa58ba6a9f81, // L8
];

#[inline]
pub(crate) fn level_domain(level_0idx: usize) -> u64 {
    LEVEL_DOMAINS.get(level_0idx).copied().unwrap_or_else(||
        LEVEL_DOMAINS[7].wrapping_add(
            (level_0idx as u64).wrapping_mul(0x6c62272e07bb0142)))
}

/// Hash `len` block values (yielded by `vals`) into the given domain. Taking an
/// iterator lets callers hash directly from an index list without materialising
/// a temporary `Vec` of the gathered values.
#[inline]
pub(crate) fn block_hash_iter_domain(len: usize, vals: impl Iterator<Item = u64>, domain: u64) -> u64 {
    // Seed mixes in both length and domain so that:
    //   • blocks of different lengths can't collide (even with same prefix)
    //   • L1 and L2 occupy disjoint hash spaces
    let mut h = block_hash_seed(len, domain);
    for v in vals { h = block_hash_mix(h, v); }
    h
}

/// Initial state of `block_hash_iter_domain` for a block of `len` values.
#[inline(always)]
pub(crate) fn block_hash_seed(len: usize, domain: u64) -> u64 {
    (len as u64).wrapping_mul(0xbf58476d1ce4e5b9) ^ domain
}

/// Fold one block value into the hash state (`block_hash_iter_domain`'s step).
#[inline(always)]
pub(crate) fn block_hash_mix(mut h: u64, v: u64) -> u64 {
    h = h.rotate_left(13).wrapping_add(v);
    // SplitMix64 avalanche step
    h ^= h >> 30;
    h = h.wrapping_mul(0xbf58476d1ce4e5b9);
    h ^= h >> 27;
    h = h.wrapping_mul(0x94d049bb133111eb);
    h ^= h >> 31;
    h
}

/// Hash a block at the given 0-indexed level. A block always covers a *contiguous*
/// run of units (the LCP rules only ever emit `lo..=hi` ranges), so its values are
/// exactly the slice `vals` — no gather needed.
#[inline]
pub(crate) fn block_hash_for_level(vals: &[u64], level_0idx: usize) -> u64 {
    block_hash_iter_domain(vals.len(), vals.iter().copied(), level_domain(level_0idx))
}


// ─────────────────────────────────────────────────────────────────────────────
// 4.  Syncmer types
// ─────────────────────────────────────────────────────────────────────────────
