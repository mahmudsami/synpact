use std::io::{BufRead, BufReader};
use std::fs::File;
use flate2::read::MultiGzDecoder;

/// Complement of one byte, uppercased; non-ACGT bytes map to their uppercase self.
///
/// Branch-free (A↔T differ by 0x15, C↔G by 0x04) instead of a 256-entry table
/// lookup, so `revcomp` compiles to a 16-byte-wide NEON/SSE loop: a table gather
/// cannot be vectorised.
#[inline(always)]
fn complement(b: u8) -> u8 {
    let u = b.to_ascii_uppercase();
    let at = ((u == b'A') | (u == b'T')) as u8;
    let cg = ((u == b'C') | (u == b'G')) as u8;
    u ^ (at * 0x15) ^ (cg * 0x04)
}

pub(crate) fn revcomp(seq: &[u8]) -> Vec<u8> {
    seq.iter().rev().map(|&b| complement(b)).collect()
}

/// `char::is_whitespace` restricted to ASCII (`\t \n \v \f \r` and space): what
/// `str::trim_end` strips from an ASCII line.
#[inline]
fn is_ascii_ws(b: u8) -> bool { matches!(b, b'\t'..=b'\r' | b' ') }

/// Yields raw sequence bytes for every FASTQ record in `path` (gzip or plain).
///
/// Lines are read as bytes: the sequence line is read straight into the Vec that
/// is returned (one exact-size allocation, uppercased in place) and the quality
/// line is skipped without being copied or UTF-8 validated. Non-ASCII input
/// falls back to the `str` path so results match a `read_line`-based parser
/// byte for byte (an invalid-UTF-8 header ends the stream, an invalid-UTF-8
/// sequence line yields an empty sequence).
pub(crate) struct FastqReader {
    pub(crate) inner: Box<dyn BufRead + Send>,
    pub(crate) buf:   Vec<u8>,
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
        FastqReader { inner, buf: Vec::new() }
    }
}

impl Iterator for FastqReader {
    type Item = (String, Vec<u8>);   // (read_name, sequence)
    fn next(&mut self) -> Option<(String, Vec<u8>)> {
        // line 1: @header
        self.buf.clear();
        if self.inner.read_until(b'\n', &mut self.buf).unwrap_or(0) == 0 { return None; }
        let header = std::str::from_utf8(&self.buf).ok()?;
        if !header.starts_with('@') { return None; }
        // read name = first whitespace-delimited token after '@'
        let name = header[1..].trim_end()
            .split_whitespace().next().unwrap_or("").to_string();
        // line 2: sequence, trailing whitespace trimmed, uppercased
        let mut seq = Vec::new();
        let _ = self.inner.read_until(b'\n', &mut seq);
        if seq.is_ascii() {
            let end = seq.iter().rposition(|&b| !is_ascii_ws(b)).map_or(0, |p| p + 1);
            seq.truncate(end);
            seq.make_ascii_uppercase();
        } else {
            seq = match std::str::from_utf8(&seq) {
                Ok(s)  => s.trim_end().bytes().map(|b| b.to_ascii_uppercase()).collect(),
                Err(_) => Vec::new(),
            };
        }
        // line 3: + separator; line 4: quality
        let _ = self.inner.skip_until(b'\n');
        let _ = self.inner.skip_until(b'\n');
        Some((name, seq))
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// The branch-free complement equals the original table for every byte.
    #[test]
    fn complement_matches_table() {
        for i in 0..=255u8 {
            let want = match i.to_ascii_uppercase() {
                b'A' => b'T', b'C' => b'G', b'G' => b'C', b'T' => b'A', u => u,
            };
            assert_eq!(complement(i), want, "byte {i}");
        }
        assert_eq!(revcomp(b"ACGTNacgtnRy"), b"YRNACGTNACGT");
    }

    fn parse(text: &[u8]) -> Vec<(String, Vec<u8>)> {
        let mut r = FastqReader { inner: Box::new(std::io::Cursor::new(text.to_vec())), buf: Vec::new() };
        std::iter::from_fn(|| r.next()).collect()
    }

    #[test]
    fn parses_trims_and_uppercases() {
        let recs = parse(b"@r1 extra\nacgTN \r\n+\nIIIII\n@r2\tx\nAC\n+\nII");
        assert_eq!(recs, vec![("r1".to_string(), b"ACGTN".to_vec()),
                              ("r2".to_string(), b"AC".to_vec())]);
        // a line that is not a header ends the stream; invalid UTF-8 in the
        // sequence line yields an empty sequence (as read_line would)
        assert!(parse(b"\n@r1\nA\n+\nI\n").is_empty());
        assert_eq!(parse(b"@r1\nA\xffC\n+\nIII\n"), vec![("r1".to_string(), vec![])]);
    }
}
