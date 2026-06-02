# syncmer-hifi

A long-read mapper for **PacBio HiFi** reads built on a hierarchy of
**locally-consistent syncmer blocks**. Instead of seed-and-extend, it places a
read by chaining matches of multi-scale blocks that are reproducible between the
read and the reference even in the presence of sequencing errors.

On simulated T2T-CHM13 HiFi reads it reaches **99.43 % placement accuracy /
99.84 % precision**, matching `minimap2 -x map-hifi` precision (99.64 %) while
running CIGAR-free mapping at ~1,800 reads/s on 8 threads, with a ~17 % smaller
in-memory index.

---

## How it works

The method turns each sequence into a **multi-level hierarchy of blocks** and
maps by chaining block matches. Six steps:

### 1. Seeds — open syncmers
A k-mer (default `k = 19`) is selected as a *syncmer* iff the minimum-hash
s-mer (default `s = 10`) inside it lies at the central offset `t = (k−s)/2`.
This yields a sparse, reproducible ~`1/(k−s+1)` = 1/10 subset of positions.
Because selection is *local*, the read and the reference pick the **same**
positions, and a single base error only perturbs the few k-mers spanning it.

A fast forward-only rotation-XOR hash (no external dependency) drives selection.
The mapper processes the read and its reverse complement as two independent
passes, so exactly one pass lights up — that determines strand.

### 2. Locally-Consistent Parsing (LCP) → blocks
The stream of syncmer values is parsed into non-overlapping **L1 blocks** by
four content-only rules (local minimum, isolated local maximum, repetition run,
monotone run). Local minima are robust to mutation, so two copies of a region
parse into the *same* blocks. Long monotone runs are split deterministically
(Cole–Vishkin coin-tossing) into ≤3-unit blocks so they carry positional signal.

### 3. Recursive hierarchy (L1 … L6)
LCP is applied again to the sequence of L1 block hashes to form L2 blocks, then
L3 … up to **6 levels**. Each level’s blocks are ~2.3× longer than the level
below (measured spans: L1 ≈ 40 bp → L6 ≈ 2.5 kb). A long read is anchored by a
few high-level blocks if they are unique, and falls back to finer levels where
they are not.

### 4. Index
Each block is a 64-bit hash. The reference index stores, per level, three
parallel arrays `(hash, chr, pos)` sorted by hash (struct-of-arrays: 13 B/entry,
cache-friendly binary search). The **atom** fed into the hierarchy is the full
k-mer (not the s-mer) — at HiFi error rates this is rarely corrupted and gives
~97.5 % genome-unique L1 blocks, the key to minimap2-level specificity.

### 5. Anchoring + chaining
For each read, blocks are looked up; over-frequent blocks are skipped and finer
children tried. Each anchor’s weight is its **conservation-of-mass** score:
every syncmer carries mass 1, a block’s mass is the sum of its units’ masses, and
a unit shared by *m* blocks contributes 1/m to each. This automatically
down-weights repetitive blocks without any explicit repeat penalty. Anchors are
pruned to the densest diagonals, then a colinear chain-DP finds the best chain;
the gap to the second-best chain sets MAPQ.

### 6. Optional CIGAR (`--cigar`)
Base-level alignment reuses the chain as a scaffold: anchor spans are emitted
directly and only the short inter-anchor gaps are aligned, with an **affine-gap**
banded DP. This is ~5–8× faster than aligning the whole read and produces clean,
consolidated indel runs.

---

## Build

Requires a Rust toolchain (`cargo`).

```sh
cargo build --release
# binary: target/release/syncmer-hifi
```

## Usage

### 1. Build an index from a reference FASTA (once)
```sh
syncmer-hifi --build-index genome.fa.gz genome.idx --threads 8
```
Default preset is `k=19 s=10`. The index records its parameters, so mapping
needs no flags.

### 2. Map HiFi reads → PAF
```sh
syncmer-hifi --map reads.fq.gz genome.idx -o out.paf --threads 8
```

### 3. Map with base-level CIGAR (for variant calling)
The CIGAR path needs the reference sequence; pass it with `--ref`:
```sh
syncmer-hifi --map reads.fq.gz genome.idx --ref genome.fa.gz --cigar -o out.paf --threads 8
```

### Options
| Flag | Default | Meaning |
|------|---------|---------|
| `--k N` | 19 | k-mer length |
| `--s N` | 10 | syncmer s-mer length (density = 1/(k−s+1)) |
| `--threads N` | all cores | worker threads |
| `--cigar [BAND]` | off | emit `cg:Z:` CIGAR (needs reference); `BAND` = half-band in bp, auto if omitted |
| `--ref <fa>` | — | reference FASTA for `--cigar` against a `.idx` |
| `--max-occ N` | 500 | max genomic occurrences per L2+ block |
| `--rescue` | off | second relaxed-filter pass for reads that fail the first — recovers a few % more reads in repeat-rich regions, emitted at MAPQ 0 |

`--rescue` runs a second mapping pass (4× looser occurrence filter) **only** on
reads the default pass leaves unmapped. It never disturbs a confidently-mapped
read, and everything it recovers is reported at MAPQ 0 so downstream MAPQ
filtering can treat it as a flagged best-guess. On simulated HiFi it lifts the
mapping rate ~0.1 pp (≈40 reads / 50 k); use it when you want maximum recall.

You can also map directly against a FASTA (it is indexed in memory first):
```sh
syncmer-hifi --map reads.fq.gz genome.fa.gz -o out.paf
```

## Output

Standard [PAF](https://github.com/lh3/miniasm/blob/master/PAF.md). Each mapped
read yields one line: `qname qlen qstart qend strand target tlen tstart tend
matches alnlen mapq`, with a trailing `cg:Z:` tag when `--cigar` is used.
Unmapped reads get a `*` record.

## Evaluating accuracy

`eval.py` compares a PAF against a ground-truth TSV (`read_name  chr  start
end`, header optional):

```sh
python3 eval.py truth.tsv out.paf            # tolerance ±1000 bp
python3 eval.py truth.tsv out.paf 500        # custom tolerance
```
It reports mapping rate, accuracy (correct within tolerance), wrong-chromosome
count, and precision.

## Recommended settings

| Goal | Setting |
|------|---------|
| Default (best accuracy) | `k=19 s=10`, 6 levels |
| Slightly faster (~+10 %) | same, but build with 5 levels |
| Highest precision / fewest wrong-chr | `k=11 s=7` (≈2× slower) |

## Notes & limitations

- Designed for **HiFi** (≤ ~1 % error) long reads. The k-mer-atom specificity
  that makes it accurate relies on low error rates; it is not intended for ONT
  or short reads.
- The residual ~0.4 % unmapped and ~0.15 % mis-placed reads are dominated by
  true reference duplications (segmental duplications, acrocentric paralogs) and
  tandem-repeat arrays, where a single best locus is genuinely ambiguous; these
  are reported at low MAPQ.

## License

MIT
