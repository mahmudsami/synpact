#![allow(unused_imports, dead_code)]
use crate::config::*;
use crate::hash::*;
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

/// Lightweight record for genome-scale processing.
#[derive(Clone)]
pub(crate) struct SyncmerLight {
    pub(crate) pos:   u32,
    pub(crate) value: u64,   // s-mer base-4 encoding — used for L1+ block hashing
}

pub(crate) fn select_syncmers_light(seq: &[u8], k: usize, s: usize, t: usize) -> Vec<SyncmerLight> {
    if seq.len() < k { return vec![]; }
    // Precompute all s-mer NT-hashes in one rolling pass.
    let smer_hashes: Vec<u64> = match DnaHashFwd::new(seq, s) {
        Some(iter) => iter.collect(),
        None => return vec![],
    };
    let mut out = Vec::new();
    for i in 0..=(seq.len() - k) {
        let kmer = &seq[i..i + k];
        if kmer.iter().any(|&b| matches!(b.to_ascii_uppercase(), b'N')) { continue; }
        let (_, min_j) = (0..=(k - s)).map(|j| (smer_hashes[i + j], j)).min().unwrap();
        if min_j == t {
            // Atom: full k-mer (high specificity, HiFi) or middle s-mer
            // (error-tolerant, ONT default).  k≤32 keeps the base-4 fit in u64.
            let atom: &[u8] = if kmer_atom() && k <= 32 { kmer } else { &kmer[t..t + s] };
            out.push(SyncmerLight { pos: i as u32, value: atom_value(atom) });
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

// ─────────────────────────────────────────────────────────────────────────────
// 5.  Locally consistent parsing — level-1 blocks
// ─────────────────────────────────────────────────────────────────────────────
