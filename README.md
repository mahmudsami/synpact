# Synpact

A long-read mapper for **PacBio HiFi** reads built on a hierarchy of
**locally-consistent syncmer blocks**. Instead of seed-and-extend, it places a
read by voting matches of multi-scale blocks that are reproducible between the
read and the reference even in the presence of sequencing errors.

The output is standard PAF: one placement (chromosome, position, strand, MAPQ)
per read, no base-level alignment.

## How it works

The pipeline turns both the reference and each read into the same kind of
multi-level block hierarchy, then places a read by matching its blocks against
an index of the reference's blocks.

### 1. Seeds — open syncmers

A k-mer (default `k = 19`) is selected as a *syncmer* iff the minimum-hash
s-mer (default `s = 10`) inside it sits at the middle offset `t = (k − s) / 2`.
The s-mer hashes come from a forward-only rolling DNA hash, and the minimum over
each `(k − s + 1)`-wide window is tracked with a monotone deque, so selection is
linear in sequence length. This open-syncmer rule keeps roughly one seed every
`k − s + 1` bases and — crucially — is *locally consistent*: a sequencing error
only disturbs the seeds whose window it falls in, leaving the rest identical
between read and reference.

Each selected syncmer carries a **value**: the base-4 encoding of its full
k-mer (`k ≤ 32` keeps it in a `u64`). Using the whole k-mer this way is what
makes placement specific at HiFi error rates.

### 2. Locally-consistent parsing (LCP) into L1 blocks

The stream of syncmer values is parsed into non-overlapping **L1 blocks** by
a set of deterministic, locally-consistent rules (`locally_consistent_parsing`):

- a triplet around every **local minimum**,
- a triplet around every **local maximum** with no adjacent minimum,
- **repetition runs** (consecutive equal values) plus their neighbours,
- **monotone runs**, which are further split into ≤ 3-unit blocks by Cole–Vishkin
  deterministic coin-tossing (`dct_three_coloring`) so a long ramp doesn't become
  one giant block.

Each rule depends only on a bounded window of values, so the same content always
parses into the same blocks wherever it occurs — the read and the reference cut
at the same places. Adjacent blocks intentionally share boundary units.

### 3. Canonical block hashing and the level hierarchy

Each block is reduced to a single 64-bit **block hash** over its ordered s-mer
values (`block_hash_for_level`). The hash folds in block length and a per-level
domain constant so blocks of different sizes or different levels occupy disjoint
hash spaces, and applies a SplitMix64 avalanche step per value so reordering
changes the hash. Two blocks with identical content anywhere in the genome get
the same hash, making it a canonical block identifier.

The L1 block hashes are then themselves fed back through the *same* LCP rules to
form **L2 blocks**, and so on up to **L6** (`MAX_LEVELS = 6`). Higher levels
cover progressively longer spans and are more unique, so they place a read
faster and with less paralog ambiguity; lower levels are the fallback when a
coarse block isn't found.

### 4. Conservation of mass

Every syncmer carries **mass 1**. A unit shared by `m` blocks contributes `1/m`
to each, and a block's mass is the sum of its units' shares, so total mass is
conserved across every level. This mass becomes the anchor weight at mapping
time, giving larger / more-supported blocks proportionally more vote.

### 5. Index

For each stored level the reference's `(block hash, chromosome, position)` tuples
are sorted by hash into a struct-of-arrays layout and binary-searched at query
time. By default only levels ≥ `--min-level` (3) are stored — the coarse, unique
levels the mapper actually anchors from — which keeps the index compact. The
index is self-describing: it records `k`, `s`, `t`, and the seed mode, so mapping
needs no parameter flags.

### 6. Mapping a read

A read is turned into its own block hierarchy, on both the forward and
reverse-complement strands. Starting from the top level, each block is looked up
in the index; on a miss or a too-repetitive hit it recurses into its finer
children (`emit_anchors`). Every hit becomes a mass-weighted **anchor**
`(query pos, reference pos)`.

Placement is by **diagonal voting** (`vote_locus`): anchors are sorted by
diagonal (`r_pos − q_pos`) and a sliding window accumulates anchor weight; the
heaviest window is the locus. **MAPQ** is `(best − second) × 60 / best`, where
`second` is the heaviest window at a genuinely different locus — so a read with
one clear placement gets a high MAPQ and one with two equally-good loci gets a
low one.

## Build

Requires a Rust toolchain (`cargo`).

```sh
cargo build --release
# binary: target/release/synpact
```

## Usage

### 1. Build an index from a reference FASTA (once)
```sh
synpact --build-index genome.fa.gz genome.idx --threads 8
```
The reference FASTA may be plain or gzip-compressed. The index records its
parameters, so mapping against it needs no flags. The FASTA is streamed
chromosome-by-chromosome and processed in parallel, so peak memory stays bounded
by the chromosomes in flight rather than the whole genome.

### 2. Map HiFi reads → PAF
```sh
synpact --map reads.fq.gz genome.idx -o out.paf --threads 8
```

You can also map directly against a FASTA (it is indexed in memory first):
```sh
synpact --map reads.fq.gz genome.fa.gz -o out.paf
```

### Options
| Flag | Default | Meaning |
|------|---------|---------|
| `--k N` | 19 | k-mer (syncmer) length; must be ≤ 32 and > `--s` |
| `--s N` | 10 | syncmer s-mer length (density ≈ 1/(k−s+1)) |
| `--threads N` | all cores | worker threads |
| `--min-level N` | 3 | lowest block level both stored in the index and used to anchor reads (one unified floor) |
| `--max-occ N` | 10000 | max genomic occurrences for an L2/L3 block before it is skipped as too-repetitive; higher levels scale this down, the L1 cap scales it to 2/5 |
| `--batch-mb N` | 256 | per-batch read-sequence budget (MiB); keeps peak memory independent of read length. `0` = unbounded |

The middle s-mer offset `t = (k − s) / 2` is derived from `--k` and `--s`.

### Environment variables
| Variable | Effect |
|----------|--------|
| `PROFILE=1` | print a per-stage CPU-time breakdown (seeding, LCP, index lookups, voting, …) after mapping |
| `BATCH_MB=N` | fallback for `--batch-mb` when the flag is not given |

## Output

Standard [PAF](https://github.com/lh3/miniasm/blob/master/PAF.md). Each read
yields one line: `qname qlen qstart qend strand target tlen tstart tend matches
alnlen mapq`. There is no base-level alignment, so the placement spans the full
read length and no `cg:Z:` CIGAR tag is emitted. Unmapped reads get a `*`
record.

## Source layout

The crate is split into one module per pipeline stage:

| Module | Responsibility |
|--------|----------------|
| `main.rs` | CLI parsing and entry point |
| `config.rs` | runtime flags (`--min-level`, `--batch-mb`), per-stage profiling, shared helpers |
| `hash.rs` | rolling DNA hash, k-mer encoding, per-level block hashing |
| `syncmer.rs` | open-syncmer selection, `SeedMode` |
| `lcp.rs` | locally-consistent parsing → `Block` / `HierForest` hierarchy, mass accounting |
| `index.rs` | `GIndex` build / serialise / load, FASTA streaming, occurrence caps |
| `fastq.rs` | FASTQ reader, reverse-complement |
| `chain.rs` | anchor collection from the hierarchy, diagonal-voting locus selection |
| `map.rs` | per-read mapping, MAPQ, PAF output, the mapping driver |

## Notes & limitations

- Designed for **HiFi** (≤ ~1 % error) long reads. The k-mer specificity
  that makes placement accurate relies on low error rates; it is not intended
  for ONT or short reads.
- The default `--min-level 3` anchors only on the coarse, unique upper levels,
  which removes most paralog / segmental-duplication ambiguity at HiFi error
  rates. Lower `--min-level` values trade specificity for sensitivity.
- Placement is locus-only: the output is a position and MAPQ, not a base-level
  alignment. Reads in genuinely ambiguous regions (segmental duplications,
  acrocentric paralogs, tandem-repeat arrays) are reported at low MAPQ.
