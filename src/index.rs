use crate::config::*;
use crate::syncmer::*;
use crate::lcp::*;
use std::io::{BufRead, BufReader, BufWriter, Read, Write};
use std::fs::File;
use std::time::Instant;
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

/// Stream a FASTA file (plain or .gz), invoking `f(name, seq)` once per
/// chromosome as it is parsed. Unlike loading the whole genome up front, this
/// lets the caller process and **drop each sequence before the next is read**,
/// so peak resident genome is bounded by what is in flight rather than the whole
/// file. Chromosomes are delivered in file order.
pub(crate) fn for_each_fasta_chr(path: &str, mut f: impl FnMut(String, Vec<u8>)) {
    let file = File::open(path).unwrap_or_else(|e| panic!("Cannot open {path}: {e}"));
    let reader: Box<dyn BufRead> = if path.ends_with(".gz") {
        Box::new(BufReader::with_capacity(1 << 20, MultiGzDecoder::new(file)))
    } else {
        Box::new(BufReader::with_capacity(1 << 20, file))
    };
    let mut name = String::new();
    let mut seq: Vec<u8> = Vec::new();
    for line in reader.lines() {
        let line = line.expect("read error");
        let bytes = line.trim_end().as_bytes();
        if bytes.first() == Some(&b'>') {
            if !seq.is_empty() {
                f(std::mem::take(&mut name), std::mem::take(&mut seq));
            }
            name = std::str::from_utf8(&bytes[1..])
                .unwrap_or("?").split_whitespace().next().unwrap_or("?")
                .to_string();
        } else {
            seq.extend(bytes.iter().map(|b| b.to_ascii_uppercase()));
        }
    }
    if !seq.is_empty() { f(name, seq); }
}

/// Extract entries for all N levels from one chromosome (N-split on 'N' runs).
/// Returns a Vec of N vecs: result[0] = L1 entries, result[1] = L2, …
pub(crate) fn chr_to_entries_n(chr_id: u8, seq: &[u8], k: usize, s: usize, t: usize, max_levels: usize,
                    mode: SeedMode) -> Vec<Vec<(u64, u8, u32)>>
{
    // per_level slot `li` holds L(li+1). Levels below --min-level are not stored
    // (kept as empty slots so a block's level still indexes its array).
    let min_level = min_lvl();
    let mut per_level: Vec<Vec<(u64, u8, u32)>> = vec![Vec::new(); max_levels];
    let mut seg_start: Option<usize> = None;

    let flush = |s0: usize, end: usize, per_level: &mut Vec<Vec<(u64, u8, u32)>>| {
        if end - s0 < k { return; }
        let levels = extract_all_levels(&seq[s0..end], k, s, t, max_levels, mode);
        for (li, entries) in levels.into_iter().enumerate() {
            if li < per_level.len() && (li + 1) >= min_level {
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

pub(crate) fn build_index(genome_path: &str, k: usize, s: usize, t: usize, max_levels: usize,
               mode: SeedMode) -> (GIndex, std::time::Duration)
{
    use std::sync::mpsc::sync_channel;
    use std::sync::Mutex;

    let t0 = Instant::now();
    let min_level = min_lvl();

    print!("    streaming sequences + extracting levels ≥ {min_level} of {max_levels} (parallel) ... ");
    std::io::stdout().flush().ok();

    // Per-level entry accumulators, written by all workers (entries land in
    // arbitrary completion order and are sorted afterwards). Each is a separate
    // mutex so workers appending to different levels never contend.
    let merged: Vec<Mutex<Vec<(u64, u8, u32)>>> =
        (0..max_levels).map(|_| Mutex::new(Vec::new())).collect();

    // A single reader thread parses the FASTA and feeds (chr_id, seq) into a
    // bounded channel; a rayon worker pool drains it, extracts entries, then
    // drops each sequence. The bound (2) caps *queued* sequences on top of the
    // ones being processed (≈ one per worker), so peak resident genome is a
    // handful of chromosomes rather than the whole file. Extraction dwarfs
    // gzip decode, so the reader stays comfortably ahead at this small depth.
    let (tx, rx) = sync_channel::<(u8, Vec<u8>)>(2);

    let chr_names: Vec<String> = std::thread::scope(|scope| {
        // Producer: assign chr_id in file order so it matches chr_names.
        let producer = scope.spawn(move || {
            let mut names: Vec<String> = Vec::new();
            let mut chr_id: usize = 0;
            for_each_fasta_chr(genome_path, |name, seq| {
                names.push(name);
                tx.send((chr_id as u8, seq)).expect("index build worker hung up");
                chr_id += 1;
            });
            names // tx dropped here → channel closes → consumers finish
        });

        // Consumers: one chromosome per task, sequence freed at task end.
        rx.into_iter().par_bridge().for_each(|(chr_id, seq)| {
            let per_level = chr_to_entries_n(chr_id, &seq, k, s, t, max_levels, mode);
            drop(seq);
            for (li, entries) in per_level.into_iter().enumerate() {
                if li < merged.len() && !entries.is_empty() {
                    merged[li].lock().unwrap().extend(entries);
                }
            }
        });

        producer.join().unwrap()
    });
    println!("{} chromosomes", chr_names.len());

    let mut merged: Vec<Vec<(u64, u8, u32)>> =
        merged.into_iter().map(|m| m.into_inner().unwrap()).collect();
    // Drop trailing empty levels (can happen for short genomes / large k)
    while merged.last().map_or(false, |v| v.is_empty()) { merged.pop(); }

    print!("    sorting L1-L{max_levels} ... ");
    std::io::stdout().flush().ok();
    // Full-tuple key: entries arrive in nondeterministic worker order, so sort
    // on (hash, chr, pos) — not hash alone — to keep the built index
    // reproducible regardless of thread scheduling.
    for level in &mut merged {
        level.sort_unstable();
    }
    println!("done");

    // levels[0] is an unused placeholder (the old L0 raw-syncmer slot) so block
    // levels still index directly: levels[1]=L1, levels[2]=L2, …
    merged.insert(0, Vec::new());

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
// Format v7 ("SYNCL2\x07\x00"):
//   8  bytes  magic
//   4  bytes  k  (u32 LE)   4 bytes  s   4 bytes  t
//   1  byte   seed_mode (0 = Syncmer)   1 byte  atom flag (always 1 = k-mer atom)
//   4  bytes  num_chrs (u32 LE);  per chr: 4-byte len + UTF-8 bytes
//   4  bytes  num_levels (u32 LE)
//   per level: 8-byte num_entries + entries × 13 bytes (8 hash | 1 chr | 4 pos)
//   Stored level 0 = unused placeholder, level 1 = L1 blocks, level 2 = L2, …
// v6/v7 load; older formats must be rebuilt.

/// Write one level in the packed 13-bytes/entry on-disk format.
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
    // Atom flag byte: bit0 = k-mer atom (always 1 now); kept for format stability.
    w.write_all(&[1u8]).unwrap();
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
    // v7+: 1-byte atom flag follows the seed-mode byte. The atom is always the
    // full k-mer now, so the flag is read and ignored.
    if version >= 7 {
        let mut b = [0u8; 1];
        r.read_exact(&mut b).unwrap();
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

    assert!(version >= 6, "index version {version} too old; rebuild with --build-index");
    // Explicit num_levels + entries in level order (level 0 = unused placeholder).
    let num_levels = read_u32le(&mut r) as usize;
    let levels: Vec<Level> = (0..num_levels).map(|_| read_level(&mut r)).collect();

    (GIndex { levels, chr_names }, k, s, t, mode)
}

// ─────────────────────────────────────────────────────────────────────────────
// 12.  FASTQ reader + reverse complement
// ─────────────────────────────────────────────────────────────────────────────

pub(crate) const MAX_OCC_DEFAULT: usize = 10_000;

/// Maximum number of hierarchy levels to build (L1 … L_MAX_LEVELS).
/// MEASURED block spans (k=19,s=10): each level grows ~2.27× (≈avg children
/// per block), NOT 3×.  Actual mean spans:
///   L1≈40 bp → L2≈91 → L3≈207 → L4≈478 → L5≈1.1 kb → L6≈2.5 kb.
/// (The old "3×, up to 10 kb" estimate was a stale k=15 guess — overstated ~4×.)
/// Tested 7 levels + smaller k: same accuracy but 2× slower (denser index +
/// the L7 level is only ~1% usable, doesn't short-circuit) → kept at 6.
pub(crate) const MAX_LEVELS: usize = 6;

/// Maximum occurrences allowed at each level (level index: 1 = L1, 2 = L2, …;
/// 0 is the unused placeholder). L1: base_l1  L2/L3: base  L4: /5  L5: /20  L6+: /100.
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
