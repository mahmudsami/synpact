#![allow(unused_imports, dead_code)]
use crate::config::*;
use crate::hash::*;
use crate::syncmer::*;
use crate::lcp::*;
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

pub(crate) struct GIndex {
    /// levels[0] = L0 raw-syncmer entries (k-mer NT-hash, built sequentially).
    /// levels[1] = L1 block entries, levels[2] = L2, …
    /// Each level is stored struct-of-arrays (SoA): three parallel arrays sorted
    /// jointly by hash.  This is 13 bytes/entry (vs 16 for the AoS `(u64,u8,u32)`
    /// which pads to 16), an 18.75% memory saving, and — more importantly — the
    /// binary search touches only the contiguous `hashes` array (8 B/entry
    /// scanned instead of 16), roughly halving cache pressure during lookup.
    pub(crate) levels:    Vec<Level>,
    pub(crate) chr_names: Vec<String>,
}

/// One hierarchy level, struct-of-arrays.  `hashes` is sorted ascending; the
/// other two arrays are permuted in lockstep so index i refers to one entry.
#[derive(Default)]
pub(crate) struct Level {
    pub(crate) hashes: Vec<u64>,
    pub(crate) chrs:   Vec<u8>,
    pub(crate) poss:   Vec<u32>,
}

impl Level {
    #[inline] pub(crate) fn len(&self) -> usize { self.hashes.len() }
    #[inline] pub(crate) fn is_empty(&self) -> bool { self.hashes.is_empty() }

    /// Build a SoA level from an AoS entry vector (used during index build).
    pub(crate) fn from_aos(entries: Vec<(u64, u8, u32)>) -> Self {
        let n = entries.len();
        let mut hashes = Vec::with_capacity(n);
        let mut chrs   = Vec::with_capacity(n);
        let mut poss   = Vec::with_capacity(n);
        for (h, c, p) in entries {
            hashes.push(h); chrs.push(c); poss.push(p);
        }
        Level { hashes, chrs, poss }
    }
}

/// A view over the matching (chr, pos) entries for a looked-up hash.
pub(crate) struct Hits<'a> {
    pub(crate) chrs: &'a [u8],
    pub(crate) poss: &'a [u32],
}

impl<'a> Hits<'a> {
    #[inline] pub(crate) fn is_empty(&self) -> bool { self.chrs.is_empty() }
    /// Iterate (chr_id, genome_pos) pairs.
    #[inline]
    pub(crate) fn iter(&self) -> impl Iterator<Item = (u8, u32)> + '_ {
        self.chrs.iter().copied().zip(self.poss.iter().copied())
    }
}

impl GIndex {
    pub(crate) fn num_levels(&self) -> usize { self.levels.len() }

    /// Look up hash h. Returns (hits, too_repetitive):
    ///   - hits non-empty when found and occ ≤ max_occ
    ///   - too_repetitive=true means hash exists but occ > max_occ
    ///     (children are likely also repetitive — caller should NOT fall back)
    /// Out-of-range levels (e.g. a top-level block above the indexed depth)
    /// return empty hits, matching the original AoS `.get()` behaviour.
    /// The binary search touches only the contiguous `hashes` array.
    pub(crate) fn lookup_with_status(&self, level_0idx: usize, h: u64, max_occ: usize)
        -> (Hits<'_>, bool)
    {
        let empty = Hits { chrs: &[], poss: &[] };
        let Some(lvl) = self.levels.get(level_0idx) else { return (empty, false) };
        let lo = lvl.hashes.partition_point(|&x| x < h);
        let hi = lvl.hashes.partition_point(|&x| x <= h);
        if hi - lo > max_occ {
            (empty, true)
        } else {
            (Hits { chrs: &lvl.chrs[lo..hi], poss: &lvl.poss[lo..hi] }, false)
        }
    }

    /// Convenience wrapper — returns empty hits on not-found OR too-repetitive.
    pub(crate) fn lookup(&self, level_0idx: usize, h: u64, max_occ: usize) -> Hits<'_> {
        self.lookup_with_status(level_0idx, h, max_occ).0
    }
}

/// Read all chromosomes from a FASTA file (plain or .gz) into memory.
pub(crate) fn read_fasta_chrs(path: &str) -> Vec<(String, Vec<u8>)> {
    let file = File::open(path).unwrap_or_else(|e| panic!("Cannot open {path}: {e}"));
    let reader: Box<dyn BufRead> = if path.ends_with(".gz") {
        Box::new(BufReader::with_capacity(1 << 20, MultiGzDecoder::new(file)))
    } else {
        Box::new(BufReader::with_capacity(1 << 20, file))
    };
    let mut chrs: Vec<(String, Vec<u8>)> = Vec::new();
    let mut name = String::new();
    let mut seq: Vec<u8> = Vec::with_capacity(1 << 28);
    for line in reader.lines() {
        let line = line.expect("read error");
        let bytes = line.trim_end().as_bytes();
        if bytes.first() == Some(&b'>') {
            if !seq.is_empty() {
                chrs.push((name.clone(), std::mem::take(&mut seq)));
            }
            name = std::str::from_utf8(&bytes[1..])
                .unwrap_or("?").split_whitespace().next().unwrap_or("?")
                .to_string();
        } else {
            seq.extend(bytes.iter().map(|b| b.to_ascii_uppercase()));
        }
    }
    if !seq.is_empty() { chrs.push((name, seq)); }
    chrs
}

/// Extract entries for all N levels from one chromosome (N-split on 'N' runs).
/// Returns a Vec of N vecs: result[0] = L1 entries, result[1] = L2, …
pub(crate) fn chr_to_entries_n(chr_id: u8, seq: &[u8], k: usize, s: usize, t: usize, max_levels: usize,
                    idx_min_level: usize, mode: SeedMode) -> Vec<Vec<(u64, u8, u32)>>
{
    // slots 0..max_levels-1 = L1-L{max_levels} blocks. Levels below idx_min_level
    // (the GIndex level index, where L0 is the raw-syncmer level) are not stored —
    // per_level slot `li` holds L(li+1), so we keep slots with li+1 >= idx_min_level.
    let mut per_level: Vec<Vec<(u64, u8, u32)>> = vec![Vec::new(); max_levels];
    let mut seg_start: Option<usize> = None;

    let flush = |s0: usize, end: usize, per_level: &mut Vec<Vec<(u64, u8, u32)>>| {
        if end - s0 < k { return; }
        let levels = extract_all_levels(&seq[s0..end], k, s, t, max_levels, mode);
        for (li, entries) in levels.into_iter().enumerate() {
            if li < per_level.len() && (li + 1) >= idx_min_level {
                for (h, pos) in entries {
                    per_level[li].push((h, chr_id, s0 as u32 + pos));
                }
            }
        }
    };

    for (i, &b) in seq.iter().enumerate() {
        match (seg_start, matches!(b, b'N')) {
            (None, false) => seg_start = Some(i),
            (Some(s0), true) => { flush(s0, i, &mut per_level); seg_start = None; }
            _ => {}
        }
    }
    if let Some(s0) = seg_start { flush(s0, seq.len(), &mut per_level); }

    per_level
}

/// Build the L0 (raw-syncmer) index level sequentially in fixed-size chunks.
///
/// Processes one chromosome at a time, each chromosome in CHUNK_SIZE-bp slices.
/// Peak memory per chunk ≈ 2 × CHUNK_SIZE × 8 bytes (smer + kmer hash arrays) = 160 MB.
/// The growing `all_l0` accumulator is the dominant memory user (~700M × 13 bytes ≈ 9 GB
/// for the human genome before filtering); it is safe on a 32 GB machine.
///
/// Returns a sorted, occ-filtered slice ready for `GIndex.levels[0]`.
pub(crate) fn extract_l0_sequential(
    chrs: &[(String, Vec<u8>)],
    k: usize, s: usize, t: usize,
    mode: SeedMode,
    max_occ: usize,
) -> Vec<(u64, u8, u32)>
{
    const CHUNK_SIZE: usize = 10_000_000; // 10 Mbp → ~80 MB per hash array
    let overlap = k - 1; // left extension so syncmers near chunk boundaries aren't missed

    let mut all_l0: Vec<(u64, u8, u32)> = Vec::new();

    for (chr_id, (_name, seq)) in chrs.iter().enumerate() {
        if seq.len() < k { continue; }
        let chr_id = chr_id as u8;

        let mut owned_start = 0usize;
        while owned_start < seq.len() {
            // Extend chunk leftward so the full k-mer window is available for
            // syncmers whose left edge is just at `owned_start`.
            let chunk_start = owned_start.saturating_sub(overlap);
            let chunk_end   = (owned_start + CHUNK_SIZE).min(seq.len());

            let chunk = &seq[chunk_start..chunk_end];
            if chunk.len() < k { owned_start += CHUNK_SIZE; continue; }

            // Syncmer positions within this chunk (N-containing k-mers already skipped).
            let syncmers = select_seeds_light(chunk, k, s, t, mode);

            // k-mer NT-hashes for the chunk — the lookup keys stored in L0.
            let kmer_hashes: Vec<u64> = match DnaHashFwd::new(chunk, k) {
                Some(it) => it.collect(),
                None     => { owned_start += CHUNK_SIZE; continue; }
            };

            for sm in &syncmers {
                let abs_pos = chunk_start + sm.pos as usize;
                // Assign each syncmer to exactly one chunk (the one that "owns" its abs pos).
                if abs_pos >= owned_start && abs_pos < chunk_end {
                    let local = sm.pos as usize;
                    if local < kmer_hashes.len() {
                        all_l0.push((kmer_hashes[local], chr_id, abs_pos as u32));
                    }
                }
            }

            owned_start += CHUNK_SIZE;
        }
    }

    // Sort by hash for binary-search lookup, then filter repetitive k-mers.
    all_l0.sort_unstable_by_key(|e| e.0);
    let before = all_l0.len();
    let mut filtered: Vec<(u64, u8, u32)> = Vec::with_capacity(before);
    let mut i = 0;
    while i < all_l0.len() {
        let h  = all_l0[i].0;
        let mut j = i + 1;
        while j < all_l0.len() && all_l0[j].0 == h { j += 1 }
        if j - i <= max_occ { filtered.extend_from_slice(&all_l0[i..j]); }
        i = j;
    }
    eprintln!("      L0 raw syncmers: {} total → {} after occ≤{} filter",
              commas(before as u64), commas(filtered.len() as u64), max_occ);
    filtered
}

pub(crate) fn build_index(genome_path: &str, k: usize, s: usize, t: usize, max_levels: usize,
               idx_min_level: usize, mode: SeedMode) -> (GIndex, std::time::Duration)
{
    let t0 = Instant::now();

    print!("    reading sequences ... ");
    std::io::stdout().flush().ok();
    let chrs = read_fasta_chrs(genome_path);
    println!("{} chromosomes", chrs.len());

    let chr_names: Vec<String> = chrs.iter().map(|(n, _)| n.clone()).collect();

    print!("    extracting levels ≥ {idx_min_level} of {max_levels} (parallel) ... ");
    std::io::stdout().flush().ok();

    // Parallel per-chromosome extraction
    let all: Vec<Vec<Vec<(u64, u8, u32)>>> = chrs.par_iter().enumerate()
        .map(|(chr_id, (_, seq))| chr_to_entries_n(chr_id as u8, seq, k, s, t, max_levels, idx_min_level, mode))
        .collect();

    // Merge across chromosomes
    let mut merged: Vec<Vec<(u64, u8, u32)>> = vec![Vec::new(); max_levels];
    for chr_levels in all {
        for (li, entries) in chr_levels.into_iter().enumerate() {
            if li < merged.len() { merged[li].extend(entries); }
        }
    }
    // Drop trailing empty levels (can happen for short genomes / large k)
    while merged.last().map_or(false, |v| v.is_empty()) { merged.pop(); }

    print!("    sorting L1-L{max_levels} ... ");
    std::io::stdout().flush().ok();
    for level in &mut merged {
        level.sort_unstable_by_key(|e| e.0);
    }
    println!("done");

    // Build L0 (raw syncmers) sequentially in 10 Mbp chunks to stay within
    // memory budget — only when L0 is within the requested floor.  L1-L{max_levels}
    // are already sorted in `merged`.
    let l0 = if idx_min_level == 0 {
        print!("    building L0 raw-syncmer level (sequential, 10 Mbp chunks) ... ");
        std::io::stdout().flush().ok();
        let l0 = extract_l0_sequential(&chrs, k, s, t, mode, MAX_OCC_L1_DEFAULT);
        println!("done  ({} entries, occ≤{})", commas(l0.len() as u64), MAX_OCC_L1_DEFAULT);
        l0
    } else {
        Vec::new()  // L0 skipped; default mapping (--min-level 3) never reads it.
    };

    // Prepend L0 so GIndex.levels[0]=L0, levels[1]=L1, …
    merged.insert(0, l0);

    let level_counts: Vec<String> = merged.iter().enumerate()
        .map(|(li, v)| format!("L{}: {}", li, commas(v.len() as u64)))
        .collect();
    println!("    {}", level_counts.join("  "));

    // Convert AoS build buffers → SoA levels.
    let levels: Vec<Level> = merged.into_iter().map(Level::from_aos).collect();

    let elapsed = t0.elapsed();
    (GIndex { levels, chr_names }, elapsed)
}

// ── Index serialisation ───────────────────────────────────────────────────────
// Format v7  ("SYNCL2\x07\x00"):
//   8  bytes  magic
//   4  bytes  k  (u32 LE)
//   4  bytes  s  (u32 LE)
//   4  bytes  t  (u32 LE)
//   1  byte   seed_mode  (0 = Syncmer)
//   1  byte   atom flags (bit0 = k-mer atom, bit1 = canonical)
//   4  bytes  num_chrs  (u32 LE)
//   for each chr: 4-byte len + UTF-8 bytes
//   4  bytes  num_levels  (u32 LE)
//   for level 0 … num_levels-1:
//       8  bytes  num_entries  (u64 LE)
//       num_entries × 13 bytes  (8-byte hash + 1-byte chr_id + 4-byte pos)
//   Level 0 = L0 raw-syncmer fallback, level 1 = L1 blocks, …
// Older formats (v4–v6) still load (atom flags default off).

/// Write one level in the packed 13-bytes/entry on-disk format (unchanged from
/// the AoS era, so existing v6 index files remain compatible).
pub(crate) fn write_entries(w: &mut impl Write, lvl: &Level) {
    w.write_all(&(lvl.len() as u64).to_le_bytes()).unwrap();
    for i in 0..lvl.len() {
        w.write_all(&lvl.hashes[i].to_le_bytes()).unwrap();
        w.write_all(&[lvl.chrs[i]]).unwrap();
        w.write_all(&lvl.poss[i].to_le_bytes()).unwrap();
    }
}

pub(crate) fn save_index(idx: &GIndex, path: &str, k: usize, s: usize, t: usize, mode: SeedMode) {
    let file = File::create(path)
        .unwrap_or_else(|e| panic!("Cannot create {path}: {e}"));
    let mut w = BufWriter::with_capacity(1 << 23, file);
    // v7: like v6 but adds a 1-byte atom flag (1 = k-mer atom, 0 = s-mer atom)
    // immediately after the seed_mode byte.  v6 loads with atom defaulted to 0.
    w.write_all(b"SYNCL2\x07\x00").unwrap();
    w.write_all(&(k as u32).to_le_bytes()).unwrap();
    w.write_all(&(s as u32).to_le_bytes()).unwrap();
    w.write_all(&(t as u32).to_le_bytes()).unwrap();
    w.write_all(&[mode as u8]).unwrap();
    // v7 flag byte: bit0 = k-mer atom, bit1 = canonical atom.
    let flags = (kmer_atom() as u8) | ((canon_atom() as u8) << 1);
    w.write_all(&[flags]).unwrap();
    w.write_all(&(idx.chr_names.len() as u32).to_le_bytes()).unwrap();
    for name in &idx.chr_names {
        let b = name.as_bytes();
        w.write_all(&(b.len() as u32).to_le_bytes()).unwrap();
        w.write_all(b).unwrap();
    }
    w.write_all(&(idx.levels.len() as u32).to_le_bytes()).unwrap();
    for lvl in &idx.levels {
        write_entries(&mut w, lvl);
    }
    w.flush().unwrap();
}

#[inline]
pub(crate) fn read_u32le(r: &mut impl Read) -> u32 {
    let mut b = [0u8; 4]; r.read_exact(&mut b).unwrap(); u32::from_le_bytes(b)
}

pub(crate) fn load_index(path: &str) -> (GIndex, usize, usize, usize, SeedMode) {
    let file = File::open(path)
        .unwrap_or_else(|e| panic!("Cannot open {path}: {e}"));
    let mut r = BufReader::with_capacity(1 << 23, file);

    let mut magic = [0u8; 8];
    r.read_exact(&mut magic).unwrap();
    assert!(&magic[..6] == b"SYNCL2", "Not a valid SYNCL2 index file");

    let k = read_u32le(&mut r) as usize;
    let s = read_u32le(&mut r) as usize;
    let t = read_u32le(&mut r) as usize;

    let version = magic[6];
    let mode = if version >= 5 {
        let mut b = [0u8; 1];
        r.read_exact(&mut b).unwrap();   // seed-mode byte (always syncmer)
        SeedMode::Syncmer
    } else {
        SeedMode::Syncmer
    };
    // v7+: flag byte (bit0 = k-mer atom, bit1 = canonical).  Older = s-mer, fwd.
    if version >= 7 {
        let mut b = [0u8; 1];
        r.read_exact(&mut b).unwrap();
        KMER_ATOM.store(b[0] & 1 != 0, Ordering::Relaxed);
        CANON_ATOM.store(b[0] & 2 != 0, Ordering::Relaxed);
    } else {
        KMER_ATOM.store(false, Ordering::Relaxed);
        CANON_ATOM.store(false, Ordering::Relaxed);
    }

    let num_chrs = read_u32le(&mut r) as usize;
    let mut chr_names = Vec::with_capacity(num_chrs);
    for _ in 0..num_chrs {
        let len = read_u32le(&mut r) as usize;
        let mut buf = vec![0u8; len];
        r.read_exact(&mut buf).unwrap();
        chr_names.push(String::from_utf8(buf).unwrap());
    }

    // Read one level directly into SoA arrays from the packed 13-byte format.
    fn read_level(r: &mut impl Read) -> Level {
        let mut b8 = [0u8; 8];
        r.read_exact(&mut b8).unwrap();
        let num = u64::from_le_bytes(b8) as usize;
        let mut hashes = Vec::with_capacity(num);
        let mut chrs   = Vec::with_capacity(num);
        let mut poss   = Vec::with_capacity(num);
        const CHUNK: usize = 65536;
        let mut buf = vec![0u8; CHUNK * 13];
        let mut remaining = num;
        while remaining > 0 {
            let n = remaining.min(CHUNK);
            r.read_exact(&mut buf[..n * 13]).unwrap();
            for i in 0..n {
                let base = i * 13;
                hashes.push(u64::from_le_bytes(buf[base..base+8].try_into().unwrap()));
                chrs.push(buf[base+8]);
                poss.push(u32::from_le_bytes(buf[base+9..base+13].try_into().unwrap()));
            }
            remaining -= n;
        }
        Level { hashes, chrs, poss }
    }

    let levels: Vec<Level> = if version >= 6 {
        // v6: explicit num_levels + entries in level order (L0, L1, L2, …)
        let num_levels = read_u32le(&mut r) as usize;
        (0..num_levels).map(|_| read_level(&mut r)).collect()
    } else if version >= 4 {
        // v4/v5: no L0 level — prepend an empty L0 slot so existing level
        // numbering (L1 at GIndex.levels[1], etc.) remains consistent.
        let num_levels = read_u32le(&mut r) as usize;
        std::iter::once(Level::default())
            .chain((0..num_levels).map(|_| read_level(&mut r)))
            .collect()
    } else if version == 3 {
        // v3 legacy: L2, L1, L3 — prepend empty L0
        let l2 = read_level(&mut r);
        let l1 = read_level(&mut r);
        let l3 = read_level(&mut r);
        vec![Level::default(), l1, l2, l3]
    } else if version == 2 {
        // v2 legacy: L2, L1 — prepend empty L0
        let l2 = read_level(&mut r);
        let l1 = read_level(&mut r);
        vec![Level::default(), l1, l2]
    } else {
        // v1 legacy: L2 only — prepend empty L0 and empty L1
        let l2 = read_level(&mut r);
        vec![Level::default(), Level::default(), l2]
    };

    (GIndex { levels, chr_names }, k, s, t, mode)
}

// ─────────────────────────────────────────────────────────────────────────────
// 12.  FASTQ reader + reverse complement
// ─────────────────────────────────────────────────────────────────────────────

pub(crate) const MAX_OCC_DEFAULT:    usize = 500;

pub(crate) const MAX_OCC_L1_DEFAULT: usize = 200;

/// Maximum number of hierarchy levels to build (L1 … L_MAX_LEVELS).
/// MEASURED block spans (k=19,s=10): each level grows ~2.27× (≈avg children
/// per block), NOT 3×.  Actual mean spans:
///   L1≈40 bp → L2≈91 → L3≈207 → L4≈478 → L5≈1.1 kb → L6≈2.5 kb.
/// (The old "3×, up to 10 kb" estimate was a stale k=15 guess — overstated ~4×.)
/// Tested 7 levels + smaller k: same accuracy but 2× slower (denser index +
/// the L7 level is only ~1% usable, doesn't short-circuit) → kept at 6.
pub(crate) const MAX_LEVELS: usize = 6;

/// Vote weight for a match at hierarchy level L (0-indexed: L0=1, L1=1, L2=3, L3=9, …).
/// L0 (raw syncmer fallback) shares the same weight as L1 blocks.
/// Each level above L1 is 3× the one below it.
#[inline]
#[allow(dead_code)]
pub(crate) fn vote_weight(level: usize) -> u32 {
    if level == 0 { 1 } else { 3u32.pow((level - 1) as u32) }
}

/// Maximum occurrences allowed at each level.
/// Level is 0-indexed: 0 = L0 raw syncmers, 1 = L1 blocks, 2 = L2, …
///   L0/L1: base_l1   L2/L3: base   L4: base/5   L5: base/20   L6+: base/100
#[inline]
pub(crate) fn max_occ_for_level(level: usize, base: usize, base_l1: usize) -> usize {
    match level {
        0 | 1 => base_l1,
        2 | 3 => base,
        4 => (base / 5).max(10),
        5 => (base / 20).max(5),
        _ => (base / 100).max(2),
    }
}
