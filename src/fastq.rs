use std::io::{BufRead, BufReader};
use std::fs::File;
use flate2::read::MultiGzDecoder;

/// Complement of each byte, uppercased; non-ACGT bytes map to their uppercase self.
const COMP: [u8; 256] = {
    let mut t = [0u8; 256];
    let mut i = 0;
    while i < 256 { t[i] = (i as u8).to_ascii_uppercase(); i += 1; }
    t[b'A' as usize] = b'T'; t[b'a' as usize] = b'T';
    t[b'C' as usize] = b'G'; t[b'c' as usize] = b'G';
    t[b'G' as usize] = b'C'; t[b'g' as usize] = b'C';
    t[b'T' as usize] = b'A'; t[b't' as usize] = b'A';
    t
};

pub(crate) fn revcomp(seq: &[u8]) -> Vec<u8> {
    seq.iter().rev().map(|&b| COMP[b as usize]).collect()
}

/// Yields raw sequence bytes for every FASTQ record in `path` (gzip or plain).
pub(crate) struct FastqReader {
    pub(crate) inner: Box<dyn BufRead + Send>,
    pub(crate) buf:   String,
}

impl FastqReader {
    pub(crate) fn open(path: &str) -> Self {
        let file = File::open(path)
            .unwrap_or_else(|e| panic!("Cannot open {path}: {e}"));
        let inner: Box<dyn BufRead + Send> = if path.ends_with(".gz") {
            Box::new(BufReader::with_capacity(1 << 20, MultiGzDecoder::new(file)))
        } else {
            Box::new(BufReader::with_capacity(1 << 20, file))
        };
        FastqReader { inner, buf: String::new() }
    }
}

impl Iterator for FastqReader {
    type Item = (String, Vec<u8>);   // (read_name, sequence)
    fn next(&mut self) -> Option<(String, Vec<u8>)> {
        // line 1: @header
        self.buf.clear();
        if self.inner.read_line(&mut self.buf).unwrap_or(0) == 0 { return None; }
        if !self.buf.starts_with('@') { return None; }
        // read name = first whitespace-delimited token after '@'
        let name = self.buf[1..].trim_end()
            .split_whitespace().next().unwrap_or("").to_string();
        // line 2: sequence
        self.buf.clear();
        self.inner.read_line(&mut self.buf).unwrap_or(0);
        let seq: Vec<u8> = self.buf.trim_end().as_bytes().iter().map(|b| b.to_ascii_uppercase()).collect();
        // line 3: + separator
        self.buf.clear();
        self.inner.read_line(&mut self.buf).unwrap_or(0);
        // line 4: quality
        self.buf.clear();
        self.inner.read_line(&mut self.buf).unwrap_or(0);
        Some((name, seq))
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// 13.  Read mapping
// ─────────────────────────────────────────────────────────────────────────────
