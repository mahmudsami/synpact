use crate::config::*;
use crate::syncmer::*;
use crate::lcp::*;
use std::io::{BufRead, BufReader, BufWriter, Read, Seek, SeekFrom, Write};
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
    /// lookup touches only the contiguous `hashes` array (8 B/entry
    /// scanned instead of 16), roughly halving cache pressure during lookup.
    pub(crate) levels:    Vec<Level>,
    pub(crate) chr_names: Vec<String>,
}

/// One hierarchy level, struct-of-arrays.  `hashes` is sorted ascending; the
/// other two arrays are permuted in lockstep so index i refers to one entry.
pub(crate) struct Level {
    pub(crate) hashes: Vec<u64>,
    pub(crate) chrs:   Vec<u8>,
    pub(crate) poss:   Vec<u32>,
    /// Bucket table over the top bits of the hash (as in strobealign): the
    /// entries of bucket `b` are `buckets[b].start..buckets[b + 1].start`, so a
    /// lookup reads one table entry and searches a few hashes instead of
    /// binary-searching the whole level. Built at load time; not stored in the
    /// index file.
    pub(crate) buckets: Vec<Bucket>,
    pub(crate) shift:   u32,
}

/// One bucket-table entry: where the bucket starts, plus a 32-bit Bloom filter
/// (3 bits per hash) of the hashes it holds.
///
/// Most lookups miss: every block of the wrong strand, and nearly every L5/L6
/// block (a 1–2.5 kb block rarely survives 0.3% error unchanged); on CHM13 HiFi
/// reads ~87% of lookups find nothing. Without the filter a miss costs two
/// dependent cache misses — the table entry, then the bucket's hashes — on an
/// index far larger than the caches. The filter sits in the table entry that
/// must be read anyway, so all but ~8% of misses are settled by that one load
/// and never touch `hashes`. (2 bits per hash: ~9% get through, slightly slower.
/// A sparser table — ≈ 2–4 entries per bucket — lets fewer through and is a
/// little faster still, but doubles the table to ~116 MB on CHM13.)
#[derive(Clone, Copy, Default)]
#[repr(C, align(8))]
pub(crate) struct Bucket {
    pub(crate) start:  u32,
    pub(crate) filter: u32,
}

impl Level {
    #[inline] pub(crate) fn len(&self) -> usize { self.hashes.len() }

    /// Wrap sorted SoA arrays and build the bucket table (≈ 4–8 entries per bucket).
    pub(crate) fn new(hashes: Vec<u64>, chrs: Vec<u8>, poss: Vec<u32>) -> Self {
        assert!(hashes.len() <= u32::MAX as usize, "index level has more than 2^32 entries");
        let bits = (usize::BITS - hashes.len().leading_zeros()).saturating_sub(3).max(1);
        let shift = 64 - bits;
        let mut buckets = vec![Bucket::default(); (1usize << bits) + 1];
        for &h in &hashes {
            let b = (h >> shift) as usize;
            buckets[b + 1].start += 1;
            buckets[b].filter |= filter_bits(h, shift);
        }
        for b in 1..buckets.len() { buckets[b].start += buckets[b - 1].start; }
        Level { hashes, chrs, poss, buckets, shift }
    }

    /// The range of bucket `b` in the entry arrays, or None if the filter rules
    /// out `h` (which must lie in bucket `b`).
    #[inline]
    fn bucket_range(&self, b: usize, h: u64) -> Option<std::ops::Range<usize>> {
        let f = filter_bits(h, self.shift);
        if self.buckets[b].filter & f != f { return None; }
        Some(self.buckets[b].start as usize..self.buckets[b + 1].start as usize)
    }

    /// Build a SoA level from an AoS entry vector (used during index build).
    /// The split runs in parallel over fixed-size chunks (pure copy, so the
    /// result does not depend on the thread count).
    pub(crate) fn from_aos(entries: Vec<(u64, u8, u32)>) -> Self {
        const CHUNK: usize = 1 << 16;
        let n = entries.len();
        let mut hashes = vec![0u64; n];
        let mut chrs   = vec![0u8;  n];
        let mut poss   = vec![0u32; n];
        hashes.par_chunks_mut(CHUNK).zip(chrs.par_chunks_mut(CHUNK)).zip(poss.par_chunks_mut(CHUNK))
            .zip(entries.par_chunks(CHUNK))
            .for_each(|(((h, c), p), e)| {
                for (i, &(eh, ec, ep)) in e.iter().enumerate() {
                    h[i] = eh; c[i] = ec; p[i] = ep;
                }
            });
        drop(entries);
        Level::new(hashes, chrs, poss)
    }
}

/// The three Bloom-filter bits of `h` in its bucket's filter, from three 5-bit
/// fields of the hash bits just below the bucket number (`h >> shift`), which
/// are as uniform as the bucket number itself.
#[inline]
fn filter_bits(h: u64, shift: u32) -> u32 {
    let r = h << (64 - shift);  // drop the bucket bits (1 ≤ 64 − shift < 64)
    (1u32 << (r >> 59)) | (1u32 << ((r >> 54) & 31)) | (1u32 << ((r >> 49) & 31))
}

/// Hint the CPU to start loading the cache line holding `*p`. A prefetch never
/// faults and has no architectural effect, so any address is fine; unlike a
/// plain load it retires at once, so many can be in flight without stalling.
#[inline(always)]
fn prefetch<T>(p: *const T) {
    #[cfg(target_arch = "aarch64")]
    // SAFETY: PRFM is a hint: it does not access memory architecturally, cannot
    // fault, and touches no register or flag other than reading `p`.
    unsafe {
        std::arch::asm!("prfm pldl1keep, [{0}]", in(reg) p,
                        options(nostack, readonly, preserves_flags));
    }
    #[cfg(target_arch = "x86_64")]
    // SAFETY: as above — PREFETCHT0 is a hint that cannot fault.
    unsafe {
        std::arch::x86_64::_mm_prefetch(p as *const i8, std::arch::x86_64::_MM_HINT_T0);
    }
    #[cfg(not(any(target_arch = "aarch64", target_arch = "x86_64")))]
    let _ = p;
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
    /// Highest indexed level (levels[0] is the unused placeholder).
    pub(crate) fn max_level(&self) -> usize { self.levels.len().saturating_sub(1) }

    /// First prefetch stage for a later `lookup(level_0idx, h)`: start loading
    /// h's bucket-table entry. Needs no memory read, so it never stalls.
    #[inline]
    pub(crate) fn prefetch_bucket(&self, level_0idx: usize, h: u64) {
        if let Some(lvl) = self.levels.get(level_0idx) {
            prefetch(&lvl.buckets[(h >> lvl.shift) as usize]);
        }
    }

    /// Second prefetch stage, once the bucket entry has (mostly) arrived: if the
    /// bucket's filter admits h, start loading the bucket's first hash and the
    /// matching chr/pos entries (a hit range nearly always begins in the
    /// bucket's first line of each array).
    #[inline]
    pub(crate) fn prefetch_entries(&self, level_0idx: usize, h: u64) {
        if let Some(lvl) = self.levels.get(level_0idx) {
            if let Some(r) = lvl.bucket_range((h >> lvl.shift) as usize, h) {
                // Pointer arithmetic only (never dereferenced here); r.start ≤ len.
                prefetch(lvl.hashes.as_ptr().wrapping_add(r.start));
                prefetch(lvl.chrs.as_ptr().wrapping_add(r.start));
                prefetch(lvl.poss.as_ptr().wrapping_add(r.start));
            }
        }
    }

    /// Look up hash h. Returns empty hits when not found, when the level is
    /// out of range (e.g. a top-level block above the indexed depth), or when
    /// the hash is too repetitive (occ > max_occ). All entries with hash `h`
    /// lie in bucket `h >> shift`, so only that bucket is searched.
    pub(crate) fn lookup(&self, level_0idx: usize, h: u64, max_occ: usize) -> Hits<'_> {
        let empty = Hits { chrs: &[], poss: &[] };
        let Some(lvl) = self.levels.get(level_0idx) else { return empty };
        let Some(r) = lvl.bucket_range((h >> lvl.shift) as usize, h) else { return empty };
        let start = r.start;
        let bucket = &lvl.hashes[r];
        let lo = start + bucket.partition_point(|&x| x < h);
        let hi = start + bucket.partition_point(|&x| x <= h);
        if hi - lo > max_occ {
            empty
        } else {
            Hits { chrs: &lvl.chrs[lo..hi], poss: &lvl.poss[lo..hi] }
        }
    }
}

/// One line of FASTA as read by [`read_fasta_line`].
enum FastaLine {
    /// End of input.
    Eof,
    /// A sequence line, appended (trimmed, uppercased) to the caller's buffer.
    Seq,
    /// A '>' header line, with the record name (first word; "?" if none).
    Header(String),
}

/// Read one FASTA line from `r`. A sequence line is appended to the end of
/// `seq` with trailing whitespace ('\n', "\r\n", spaces, tabs — as
/// `str::trim_end`) removed and ASCII uppercased; a header line leaves `seq`
/// unchanged. Also returns the raw byte length consumed.
///
/// The line goes straight onto the end of `seq` (read_until scans for '\n'
/// with memchr and copies in bulk) and is trimmed/uppercased in place. This
/// avoids a String + UTF-8 validation per line and a byte-at-a-time extend: the
/// old `lines()` loop ran at ≈1 GB/s and, as the single producer, was the
/// critical path of the parallel index build.
fn read_fasta_line(r: &mut impl BufRead, seq: &mut Vec<u8>) -> (usize, FastaLine) {
    let start = seq.len();
    let n = r.read_until(b'\n', seq).expect("read error");
    if n == 0 { return (0, FastaLine::Eof); }
    let mut end = seq.len();
    while end > start && matches!(seq[end - 1], b'\t' | b'\n' | b'\x0b' | b'\x0c' | b'\r' | b' ') {
        end -= 1;
    }
    if end > start && seq[start] == b'>' {
        let name = std::str::from_utf8(&seq[start + 1..end])
            .unwrap_or("?").split_whitespace().next().unwrap_or("?")
            .to_string();
        seq.truncate(start);
        (n, FastaLine::Header(name))
    } else {
        seq[start..end].make_ascii_uppercase();
        seq.truncate(end);
        (n, FastaLine::Seq)
    }
}

/// Stream a FASTA file (plain or .gz), invoking `f(name, seq)` once per
/// chromosome as it is parsed. Unlike loading the whole genome up front, this
/// lets the caller process and **drop each sequence before the next is read**,
/// so peak resident genome is bounded by what is in flight rather than the whole
/// file. Chromosomes are delivered in file order; records with an empty
/// sequence are skipped, and sequence before the first header is delivered
/// under the name "".
pub(crate) fn for_each_fasta_chr(path: &str, mut f: impl FnMut(String, Vec<u8>)) {
    let file = File::open(path).unwrap_or_else(|e| panic!("Cannot open {path}: {e}"));
    let mut reader: Box<dyn BufRead> = if path.ends_with(".gz") {
        Box::new(BufReader::with_capacity(1 << 20, MultiGzDecoder::new(file)))
    } else {
        Box::new(BufReader::with_capacity(1 << 20, file))
    };
    let mut name = String::new();
    let mut seq: Vec<u8> = Vec::new();
    loop {
        match read_fasta_line(&mut reader, &mut seq).1 {
            FastaLine::Eof => break,
            FastaLine::Seq => {}
            FastaLine::Header(header) => {
                if !seq.is_empty() {
                    f(std::mem::replace(&mut name, header), std::mem::take(&mut seq));
                } else {
                    name = header;
                }
            }
        }
    }
    if !seq.is_empty() { f(name, seq); }
}

/// A FASTA record of an uncompressed file, located without parsing its
/// sequence: its name and the byte range of its sequence lines.
struct FastaRecord {
    name: String,
    body: std::ops::Range<u64>,
}

impl FastaRecord {
    /// Size of the sequence lines in bytes (newlines included): an upper bound
    /// on, and a close estimate of, the sequence length.
    fn body_len(&self) -> u64 { self.body.end - self.body.start }
}

/// Open `path` positioned at `range.start`, reading no further than `range.end`.
fn open_range(path: &str, range: std::ops::Range<u64>) -> BufReader<std::io::Take<File>> {
    let mut file = File::open(path).unwrap_or_else(|e| panic!("Cannot open {path}: {e}"));
    file.seek(SeekFrom::Start(range.start)).expect("seek error");
    BufReader::with_capacity(1 << 16, file.take(range.end - range.start))
}

/// Locate the records of an uncompressed FASTA file — exactly the (name, seq)
/// pairs, in the same order, that [`for_each_fasta_chr`] delivers — so each
/// can then be read independently with [`read_fasta_record`].
///
/// Why: with one streaming reader, chromosomes are processed in file order,
/// and a large one late in the file (chrX in CHM13) starts only after the
/// first wave of workers frees up and becomes the build's tail. Knowing every
/// record up front lets the build process the largest first, each worker
/// reading its own record, so nothing is read ahead and held in memory.
///
/// Header lines ('>' at the start of a line) are found by a parallel scan of
/// fixed-size blocks; '>' is rare, so a block is first tested with a
/// vectorisable OR-reduction and only searched when it contains one.
fn scan_fasta_records(path: &str) -> Vec<FastaRecord> {
    let len = std::fs::metadata(path)
        .unwrap_or_else(|e| panic!("Cannot open {path}: {e}")).len();
    const BLOCK: u64 = 1 << 22;
    let headers: Vec<u64> = (0..len.div_ceil(BLOCK)).into_par_iter().flat_map_iter(|bi| {
        let lo = bi * BLOCK;
        let hi = (lo + BLOCK).min(len);
        // Also read the byte before the block: a '>' is a header only at the
        // start of a line (file start or right after '\n').
        let from = lo.saturating_sub(1);
        let mut file = File::open(path).unwrap_or_else(|e| panic!("Cannot open {path}: {e}"));
        file.seek(SeekFrom::Start(from)).expect("seek error");
        let mut buf = vec![0u8; (hi - from) as usize];
        file.read_exact(&mut buf).expect("read error");
        let mut found = Vec::new();
        for (ci, chunk) in buf.chunks(64).enumerate() {
            if chunk.iter().fold(0u8, |acc, &b| acc | (b == b'>') as u8) == 0 { continue; }
            for (j, &b) in chunk.iter().enumerate() {
                let i = ci * 64 + j;
                let pos = from + i as u64;
                if b == b'>' && pos >= lo && (pos == 0 || buf[i - 1] == b'\n') {
                    found.push(pos);
                }
            }
        }
        found
    }).collect();

    // Whether the lines in `range` hold any sequence (a non-blank line).
    let has_seq = |range: std::ops::Range<u64>| -> bool {
        let mut r = open_range(path, range);
        let mut scratch = Vec::new();
        loop {
            match read_fasta_line(&mut r, &mut scratch).1 {
                FastaLine::Eof => return false,
                _ if !scratch.is_empty() => return true,
                _ => {}
            }
        }
    };

    // Sequence before the first header forms a record named "" (as when streaming).
    let first = headers.first().copied().unwrap_or(len);
    let lead = (first > 0 && has_seq(0..first))
        .then(|| FastaRecord { name: String::new(), body: 0..first });
    let named: Vec<Option<FastaRecord>> = (0..headers.len()).into_par_iter().map(|i| {
        let end = headers.get(i + 1).copied().unwrap_or(len);
        let mut r = open_range(path, headers[i]..end);
        let (n, line) = read_fasta_line(&mut r, &mut Vec::new());
        let FastaLine::Header(name) = line else { unreachable!("record does not start with a header") };
        let body = headers[i] + n as u64..end;
        has_seq(body.clone()).then(|| FastaRecord { name, body })
    }).collect();
    lead.into_iter().chain(named.into_iter().flatten()).collect()
}

/// Read one record located by [`scan_fasta_records`]: its sequence lines,
/// trimmed and uppercased exactly as [`for_each_fasta_chr`] does.
fn read_fasta_record(path: &str, rec: &FastaRecord) -> Vec<u8> {
    let mut r = open_range(path, rec.body.clone());
    let mut seq = Vec::with_capacity(rec.body_len() as usize);
    while !matches!(read_fasta_line(&mut r, &mut seq).1, FastaLine::Eof) {}
    seq
}

/// Position of the first byte at or after `from` that is 'N' (`want_n`) or is
/// not 'N' (`!want_n`); `seq.len()` if there is none. Whole 32-byte chunks
/// without a match are skipped with a vectorisable reduction, so splitting a
/// chromosome on N runs costs a fraction of a cycle per base instead of a
/// byte-at-a-time state machine.
fn next_n_or_not(seq: &[u8], from: usize, want_n: bool) -> usize {
    let is_match = |b: u8| (b == b'N') == want_n;
    let mut i = from;
    while let Some(chunk) = seq.get(i..i + 32) {
        if chunk.iter().fold(false, |acc, &b| acc | is_match(b)) { break; }
        i += 32;
    }
    seq[i..].iter().position(|&b| is_match(b)).map_or(seq.len(), |p| i + p)
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

    // Segments are the maximal runs of non-'N' bytes.
    let mut s0 = next_n_or_not(seq, 0, false);
    while s0 < seq.len() {
        let end = next_n_or_not(seq, s0, true);
        flush(s0, end, &mut per_level);
        s0 = next_n_or_not(seq, end, false);
    }

    per_level
}

/// Index build: cap on the sequence bases (per worker thread) whose entries are
/// being extracted at once from a plain FASTA. Extraction peaks at ≈5 bytes of
/// working memory per base, so 150 Mbp/thread ≈ 0.75 GB/thread — about what the
/// streaming build reached with its staggered starts (≈7 GB on CHM13 with 10
/// threads) while still keeping every thread busy.
const BUILD_BASES_PER_THREAD: u64 = 150_000_000;

pub(crate) fn build_index(genome_path: &str, k: usize, s: usize, t: usize, max_levels: usize,
               mode: SeedMode) -> (GIndex, std::time::Duration)
{
    use std::sync::mpsc::sync_channel;
    use std::sync::{Condvar, Mutex};

    let t0 = Instant::now();
    let min_level = min_lvl();

    print!("    streaming sequences + extracting levels ≥ {min_level} of {max_levels} (parallel) ... ");
    std::io::stdout().flush().ok();

    // Per-level entry accumulators, written by all workers (entries land in
    // arbitrary completion order and are sorted afterwards). Each is a separate
    // mutex so workers appending to different levels never contend.
    let merged: Vec<Mutex<Vec<(u64, u8, u32)>>> =
        (0..max_levels).map(|_| Mutex::new(Vec::new())).collect();

    // Extract one chromosome's entries and append them to `merged`.
    let process = |chr_id: u8, seq: Vec<u8>| {
        let per_level = chr_to_entries_n(chr_id, &seq, k, s, t, max_levels, mode);
        drop(seq);
        for (li, entries) in per_level.into_iter().enumerate() {
            if li < merged.len() && !entries.is_empty() {
                merged[li].lock().unwrap().extend(entries);
            }
        }
    };

    // The parallel record scan below seeks within the file, so anything that is
    // not a regular file (a pipe, `<(zcat genome.fa.gz)`, /dev/stdin) is streamed
    // like gzip. A path that cannot be opened also takes this branch and reports
    // the error when the reader opens it.
    let seekable = std::fs::metadata(genome_path).map(|m| m.is_file()).unwrap_or(false);
    let chr_names: Vec<String> = if genome_path.ends_with(".gz") || !seekable {
        // A gzip stream (or pipe) can only be read front to back: a single reader thread
        // parses the FASTA and feeds (chr_id, seq) into a bounded channel; a
        // rayon worker pool drains it, extracts entries, then drops each
        // sequence. The bound (2) caps *queued* sequences on top of the ones
        // being processed (≈ one per worker), so peak resident genome is a
        // handful of chromosomes rather than the whole file. Extraction dwarfs
        // gzip decode, so the reader stays comfortably ahead at this small depth.
        let (tx, rx) = sync_channel::<(u8, Vec<u8>)>(2);
        std::thread::scope(|scope| {
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
            rx.into_iter().par_bridge().for_each(|(chr_id, seq)| process(chr_id, seq));

            producer.join().unwrap()
        })
    } else {
        // Plain FASTA: locate every record first, then process them largest
        // first (LPT scheduling), each worker reading its own record when it
        // starts it, so a large chromosome late in the file no longer becomes
        // the tail. chr_id stays the record's rank in file order.
        //
        // Starting the largest records all at once would also line up their
        // extraction working sets (≈5 bytes/base at peak), which the streaming
        // reader used to stagger (+~1 GB peak RSS on CHM13, 10 threads). So the bases
        // in flight are capped at BUILD_BASES_PER_THREAD × threads: a worker
        // takes the largest pending record that fits under the cap (or any
        // record when nothing is in flight), and otherwise waits for one to
        // finish. Smaller records fill the gaps, so the cores stay busy.
        let records = scan_fasta_records(genome_path);
        let nthreads = rayon::current_num_threads();
        let budget = BUILD_BASES_PER_THREAD * nthreads as u64;
        // Ascending by size, so the largest pending record is popped off the end.
        let mut pending: Vec<usize> = (0..records.len()).collect();
        pending.sort_by_key(|&i| records[i].body_len());
        let sched = Mutex::new((pending, 0u64)); // (pending records, bases in flight)
        let freed = Condvar::new();
        // Returns a record's bases to the budget when its task ends, even by panic
        // (so the other workers are not left waiting forever).
        struct Release<'a> { sched: &'a Mutex<(Vec<usize>, u64)>, freed: &'a Condvar, len: u64 }
        impl Drop for Release<'_> {
            fn drop(&mut self) {
                if let Ok(mut g) = self.sched.lock() { g.1 -= self.len; }
                self.freed.notify_all();
            }
        }
        std::thread::scope(|scope| {
            for _ in 0..nthreads {
                scope.spawn(|| loop {
                    let (i, _release) = {
                        let mut g = sched.lock().unwrap();
                        loop {
                            let (pending, in_flight) = &mut *g;
                            if pending.is_empty() { return; }
                            let fits = |&i: &usize| *in_flight == 0 || *in_flight + records[i].body_len() <= budget;
                            if let Some(j) = pending.iter().rposition(fits) {
                                let i = pending.remove(j);
                                let len = records[i].body_len();
                                *in_flight += len;
                                break (i, Release { sched: &sched, freed: &freed, len });
                            }
                            g = freed.wait(g).unwrap();
                        }
                    };
                    process(i as u8, read_fasta_record(genome_path, &records[i]));
                });
            }
        });
        records.into_iter().map(|r| r.name).collect()
    };
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
    // Sorted in parallel (all levels at once, each with a parallel sort): the
    // key is a total order and equal tuples are indistinguishable, so the
    // result is identical to a sequential sort_unstable.
    merged.par_iter_mut().for_each(|level| level.par_sort_unstable());
    println!("done");

    // levels[0] is an unused placeholder (the old L0 raw-syncmer slot) so block
    // levels still index directly: levels[1]=L1, levels[2]=L2, …
    merged.insert(0, Vec::new());

    let level_counts: Vec<String> = merged.iter().enumerate()
        .map(|(li, v)| format!("L{}: {}", li, commas(v.len() as u64)))
        .collect();
    println!("    {}", level_counts.join("  "));

    // Convert AoS build buffers → SoA levels.
    let levels: Vec<Level> = merged.into_par_iter().map(Level::from_aos).collect();

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
    // Encode a chunk of records into one buffer and write it with a single call;
    // three tiny write_all calls per entry spent most of the save time in
    // BufWriter bookkeeping rather than in I/O.
    const CHUNK: usize = 1 << 16;
    let mut buf = vec![0u8; CHUNK * 13];
    for ((h, c), p) in lvl.hashes.chunks(CHUNK).zip(lvl.chrs.chunks(CHUNK)).zip(lvl.poss.chunks(CHUNK)) {
        let out = &mut buf[..h.len() * 13];
        for (i, rec) in out.chunks_exact_mut(13).enumerate() {
            rec[..8].copy_from_slice(&h[i].to_le_bytes());
            rec[8] = c[i];
            rec[9..].copy_from_slice(&p[i].to_le_bytes());
        }
        w.write_all(out).unwrap();
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
    // Only the header is read through this reader; the levels are read by
    // `read_level` with their own handles.
    let mut r = BufReader::with_capacity(1 << 16, file);

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
    // The entries are decoded in parallel: each fixed-size chunk opens its own
    // handle at its own offset and fills its slice of the pre-sized arrays, so
    // the read and the AoS→SoA decode both spread over the threads (on a
    // single thread this is the old chunked loop minus the per-entry pushes).
    // Returns the level and the file offset just past it.
    fn read_level(path: &str, offset: u64) -> (Level, u64) {
        let open_at = |at: u64| {
            let mut f = File::open(path).unwrap_or_else(|e| panic!("Cannot open {path}: {e}"));
            f.seek(SeekFrom::Start(at)).unwrap();
            f
        };
        let mut b8 = [0u8; 8];
        open_at(offset).read_exact(&mut b8).unwrap();
        let num = u64::from_le_bytes(b8) as usize;
        let data = offset + 8;
        let mut hashes = vec![0u64; num];
        let mut chrs   = vec![0u8;  num];
        let mut poss   = vec![0u32; num];
        const CHUNK: usize = 1 << 18;
        hashes.par_chunks_mut(CHUNK).zip(chrs.par_chunks_mut(CHUNK)).zip(poss.par_chunks_mut(CHUNK))
            .enumerate()
            .for_each(|(ci, ((h, c), p))| {
                let mut buf = vec![0u8; h.len() * 13];
                open_at(data + (ci * CHUNK * 13) as u64).read_exact(&mut buf).unwrap();
                for (i, rec) in buf.chunks_exact(13).enumerate() {
                    h[i] = u64::from_le_bytes(rec[..8].try_into().unwrap());
                    c[i] = rec[8];
                    p[i] = u32::from_le_bytes(rec[9..].try_into().unwrap());
                }
            });
        (Level::new(hashes, chrs, poss), data + num as u64 * 13)
    }

    assert!(version >= 6, "index version {version} too old; rebuild with --build-index");
    // Explicit num_levels + entries in level order (level 0 = unused placeholder).
    let num_levels = read_u32le(&mut r) as usize;
    let mut offset = r.stream_position().unwrap();
    drop(r);
    let levels: Vec<Level> = (0..num_levels).map(|_| {
        let (lvl, next) = read_level(path, offset);
        offset = next;
        lvl
    }).collect();

    (GIndex { levels, chr_names }, k, s, t, mode)
}

// ─────────────────────────────────────────────────────────────────────────────
// 12.  FASTQ reader + reverse complement
// ─────────────────────────────────────────────────────────────────────────────

pub(crate) const MAX_OCC_DEFAULT: usize = 10_000;

/// Default `--max-level`: the number of hierarchy levels built (L1 … L_MAX_LEVELS).
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
