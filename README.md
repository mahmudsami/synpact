# syncmer-hifi

A long-read mapper for **PacBio HiFi** reads built on a hierarchy of
**locally-consistent syncmer blocks**. Instead of seed-and-extend, it places a
read by voting matches of multi-scale blocks that are reproducible between the
read and the reference even in the presence of sequencing errors.

On simulated T2T-CHM13 HiFi reads (100 k × 24 kb, 0.1 % error) it reaches
**99.59 % placement accuracy / 99.97 % precision with zero wrong-chromosome
calls**, exceeding `minimap2 -x map-hifi` precision while running CIGAR-free
mapping at **~2,400 reads/s** on 8 threads. Because the default index stores only
the coarse blocks it actually uses (levels ≥ 3), the T2T-CHM13 index is just
**~560 MB** — roughly an order of magnitude smaller than the full hierarchy.

---

## How it works

The method turns each sequence into a **multi-level hierarchy of blocks** and
maps by voting on block matches. Six steps:

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
parse into the *same* blocks. Long monotone runs are split deterministically by
**Cole–Vishkin deterministic coin-tossing**: the run (repetition-free, since
strictly monotone) is reduced to a proper 3-colouring (iterated `coin' = 2·π + z`
down to ≤6 colours, then collapsed to {0,1,2}), and then parsed by the *same*
local-minimum and local-maximum rules used above. Over a 3-colour alphabet no two
slope points are adjacent, so every position is covered by a min/max triplet:
this provably yields blocks of size ≤3 with a bounded dependence window (locally
consistent), carrying positional signal without any heuristic length cap.

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

### 5. Anchoring + voting
For each read, blocks are looked up; over-frequent blocks are skipped and finer
children tried. Each anchor’s weight is its **conservation-of-mass** score:
every syncmer carries mass 1, a block’s mass is the sum of its units’ masses, and
a unit shared by *m* blocks contributes 1/m to each. This automatically
down-weights repetitive blocks without any explicit repeat penalty.

The read’s locus is then chosen by **diagonal voting** (default): anchors are
binned by diagonal (`r_pos − q_pos`) and the heaviest cluster wins, with the gap
to the next-best locus setting MAPQ. Voting is slightly more accurate than a
colinear chain-DP on contested (paralogous) reads and just as fast; pass
`--chaining` to use the gap-penalised chain-DP instead.

By default only blocks at **level ≥ 3** (≈ L3–L6, spans ≳ 150 bp) are indexed and
allowed to anchor a read (`--min-level`, default 3). The short L0–L2 blocks match
in many places across segmental duplications and paralogs, and at HiFi error rates
they add ambiguity rather than signal: dropping them raises accuracy, eliminates
wrong-chromosome calls, runs ~50 % faster, and shrinks the index ~10×. Build with
`--index-min-level 0` and map with `--min-level 0` for noisier (> 1 % error) data,
where the finer levels’ sensitivity is needed.

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
Default preset is `k=19 s=10`, indexing only the coarse blocks the default mapper
uses (levels ≥ 3) — the T2T-CHM13 index is ~560 MB and builds in well under a
minute. The index records its parameters, so mapping needs no flags. To also
support mapping noisy (> 1 % error) reads at `--min-level 0`, build the full
hierarchy with `--index-min-level 0`.

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
| `--min-level N` | 3 | lowest block level allowed to anchor a read; default 3 is tuned for HiFi, use 0 for >1 % error reads |
| `--index-min-level N` | 3 | *(build only)* lowest block level to index; match `--min-level`. Use 0 to support `--min-level 0` mapping |
| `--vote` | on | place reads by diagonal voting (default selector) |
| `--chaining` | off | place reads by colinear chain-DP instead of voting |
| `--batch-mb N` | 256 | per-batch read-sequence budget (MiB); keeps peak memory independent of read length. `0` = unbounded |
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
| Default (best HiFi accuracy + speed) | `k=19 s=10 --min-level 3` (the default) |
| Noisy reads (> 1 % error) | build `--index-min-level 0`, map `--min-level 0` |
| Highest precision / fewest wrong-chr | `k=11 s=7` (≈2× slower) |

## Source layout

The crate is split into one module per pipeline stage:

| Module | Responsibility |
|--------|----------------|
| `config.rs` | runtime flags (`--min-level`, `--vote`/`--chaining`, `--rescue`, …) and shared helpers |
| `hash.rs` | rolling DNA hash, atom encoding, per-level block hashing |
| `syncmer.rs` | open-syncmer selection, `SeedMode` |
| `lcp.rs` | locally-consistent parsing → `Block`/`HierBlock` hierarchy |
| `index.rs` | `GIndex` build / serialise / load (levels ≥ `--index-min-level`) |
| `fastq.rs` | FASTQ reader, reverse-complement |
| `align.rs` | banded affine-gap alignment and chain-guided CIGAR |
| `chain.rs` | anchor collection, diagonal voting, chain-DP |
| `map.rs` | per-read mapping, MAPQ, PAF output, the mapping driver |
| `main.rs` | CLI parsing and entry point |

## Notes & limitations

- Designed for **HiFi** (≤ ~1 % error) long reads. The k-mer-atom specificity
  that makes it accurate relies on low error rates; it is not intended for ONT
  or short reads. The default `--min-level 3` is tuned for this regime — above
  ~1 % error it trades too much sensitivity, so use `--min-level 0` there.
- The residual ~0.4 % unmapped and ~0.15 % mis-placed reads are dominated by
  true reference duplications (segmental duplications, acrocentric paralogs) and
  tandem-repeat arrays, where a single best locus is genuinely ambiguous; these
  are reported at low MAPQ.

## License

MIT
