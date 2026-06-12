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
use std::collections::{HashSet, VecDeque};
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
    // s-mer NT-hashes are pulled on demand straight into the sliding-window min,
    // so the full L-length hash array is never materialised.
    let mut hashes = match DnaHashFwd::new(seq, s) {
        Some(it) => it,
        None => return vec![],
    };
    // Sliding-window minimum over the (k−s+1)-wide s-mer-hash window via a
    // monotone deque of (index, hash). Popping the back only on a strict `>`
    // keeps the leftmost index on ties — identical argmin to `.min()` on
    // `(hash, j)`, which breaks ties toward the smaller offset j.
    let w = k - s + 1;
    let n_pos = seq.len() - k + 1;
    let mut out = Vec::with_capacity(n_pos / w + 1);
    let mut dq: VecDeque<(usize, u64)> = VecDeque::new();
    let mut r = 0usize;          // next s-mer index to admit into the deque
    for i in 0..n_pos {
        let right = i + w - 1;
        while r <= right {
            let h = hashes.next().unwrap();
            while let Some(&(_, hb)) = dq.back() {
                if hb > h { dq.pop_back(); } else { break; }
            }
            dq.push_back((r, h));
            r += 1;
        }
        while let Some(&(f, _)) = dq.front() {
            if f < i { dq.pop_front(); } else { break; }
        }
        let min_j = dq.front().unwrap().0 - i;
        if min_j == t {
            let kmer = &seq[i..i + k];
            if kmer.iter().any(|&b| matches!(b.to_ascii_uppercase(), b'N')) { continue; }
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
