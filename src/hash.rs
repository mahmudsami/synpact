#![allow(unused_imports, dead_code)]
use crate::config::*;
use crate::syncmer::*;
use crate::lcp::*;
use crate::index::*;
use crate::fastq::*;
use crate::align::*;
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

/// Base-4 encode the reverse-complement of `atom` without allocating.
#[inline]
pub(crate) fn encode_revcomp(atom: &[u8]) -> u64 {
    let mut h = 0u64;
    for &b in atom.iter().rev() {
        let comp = match b.to_ascii_uppercase() {
            b'A' => 3, b'C' => 2, b'G' => 1, b'T' => 0, _ => 0, // complement, base-4
        };
        h = h * 4 + comp;
    }
    h
}

/// The atom value for a window, honoring the canonical flag.
#[inline]
pub(crate) fn atom_value(atom: &[u8]) -> u64 {
    let f = encode_smer(atom);
    if canon_atom() { f.min(encode_revcomp(atom)) } else { f }
}

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

/// Rolling forward-only DNA hash over a window of `w` consecutive bases.
/// Yields one `u64` per starting position.  An N inside the current window
/// makes that window's hash `u64::MAX` (so it is never selected as a minimum).
pub(crate) struct DnaHashFwd<'a> {
    pub(crate) seq:  &'a [u8],
    pub(crate) w:    usize,
    pub(crate) h:    u64,
    pub(crate) n_in_window: u32,  // count of N-like bases currently inside the window
    pub(crate) pos:  usize,
}

impl<'a> DnaHashFwd<'a> {
    /// Returns `None` when `seq.len() < w`.
    pub(crate) fn new(seq: &'a [u8], w: usize) -> Option<Self> {
        if seq.len() < w { return None; }
        let mut h = 0u64;
        let mut n_count = 0u32;
        for &b in &seq[..w] {
            if DNA_SEED[b as usize] == 0 { n_count += 1; }
            h = h.rotate_left(1) ^ DNA_SEED[b as usize];
        }
        Some(DnaHashFwd { seq, w, h, n_in_window: n_count, pos: 0 })
    }
}

impl<'a> Iterator for DnaHashFwd<'a> {
    type Item = u64;
    #[inline]
    fn next(&mut self) -> Option<u64> {
        if self.pos + self.w > self.seq.len() { return None; }
        let h = if self.n_in_window > 0 { u64::MAX } else { self.h };
        // Roll forward one position (prepare hash for pos+1)
        if self.pos + self.w < self.seq.len() {
            let out = self.seq[self.pos];
            let in_ = self.seq[self.pos + self.w];
            if DNA_SEED[out as usize] == 0 { self.n_in_window -= 1; }
            if DNA_SEED[in_  as usize] == 0 { self.n_in_window += 1; }
            self.h = self.h.rotate_left(1)
                ^ DNA_SEED[out as usize].rotate_left(self.w as u32)
                ^ DNA_SEED[in_  as usize];
        }
        self.pos += 1;
        Some(h)
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// 2.  s-mer direct encoding  (A=0, C=1, G=2, T=3, base-4 positional)
// ─────────────────────────────────────────────────────────────────────────────

#[inline]
pub(crate) fn base4(b: u8) -> u64 {
    match b.to_ascii_uppercase() {
        b'A' => 0, b'C' => 1, b'G' => 2, b'T' => 3, _ => 0,
    }
}

#[inline]
pub(crate) fn encode_smer(smer: &[u8]) -> u64 {
    smer.iter().fold(0u64, |acc, &b| acc * 4 + base4(b))
}

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

/// Hash a block's values at the given 0-indexed level.
#[inline]
pub(crate) fn block_hash_for_level(values: &[u64], level_0idx: usize) -> u64 {
    block_hash_with_domain(values, level_domain(level_0idx))
}

// Backward-compatible aliases used by the demo/stats display code.

pub(crate) fn block_hash_with_domain(values: &[u64], domain: u64) -> u64 {
    block_hash_iter_domain(values.len(), values.iter().copied(), domain)
}

/// Hash `len` block values (yielded by `vals`) into the given domain. Taking an
/// iterator lets callers hash directly from an index list without materialising
/// a temporary `Vec` of the gathered values.
#[inline]
pub(crate) fn block_hash_iter_domain(len: usize, vals: impl Iterator<Item = u64>, domain: u64) -> u64 {
    // Seed mixes in both length and domain so that:
    //   • blocks of different lengths can't collide (even with same prefix)
    //   • L1 and L2 occupy disjoint hash spaces
    let mut h: u64 = (len as u64)
        .wrapping_mul(0xbf58476d1ce4e5b9)
        ^ domain;
    for v in vals {
        h = h.rotate_left(13).wrapping_add(v);
        // SplitMix64 avalanche step
        h ^= h >> 30;
        h = h.wrapping_mul(0xbf58476d1ce4e5b9);
        h ^= h >> 27;
        h = h.wrapping_mul(0x94d049bb133111eb);
        h ^= h >> 31;
    }
    h
}

/// Hash a block whose values are `src[i]` for each `i` in `indices`, at the given
/// 0-indexed level — without gathering them into a temporary `Vec` first.
#[inline]
pub(crate) fn block_hash_indices_for_level(indices: &[usize], src: &[u64], level_0idx: usize) -> u64 {
    block_hash_iter_domain(indices.len(), indices.iter().map(|&i| src[i]), level_domain(level_0idx))
}


// ─────────────────────────────────────────────────────────────────────────────
// 4.  Syncmer types
// ─────────────────────────────────────────────────────────────────────────────
