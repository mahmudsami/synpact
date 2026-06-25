# Synpact

A long-read mapper for **PacBio HiFi** reads built on a hierarchy of
**locally-consistent syncmer blocks**. Instead of seed-and-extend, it places a
read by voting matches of multi-scale blocks that are reproducible between the
read and the reference even in the presence of sequencing errors.

The output is standard PAF: one placement (chromosome, position, strand, MAPQ)
per read, no base-level alignment.

## Build

Requires a Rust toolchain (`cargo`).

```sh
cargo build --release
# binary: target/release/synpact
```

## Usage

### 1. Build an index then map
```sh
synpact --build-index genome.fa.gz genome.idx --threads 8
synpact --map reads.fq.gz genome.idx -o out.paf --threads 8
```
The reference FASTA may be plain or gzip-compressed. The index records its
parameters, so mapping against it needs no flags. The FASTA is streamed
chromosome-by-chromosome and processed in parallel, so peak memory stays bounded
by the chromosomes in flight rather than the whole genome. Fewer threads at build time
lowers peak memory for the indexing step.

### 2. Index and map in one go
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


## Citation

Aydin M.S. and Sahlin, K. synpact: accurate, memory-light PacBio HiFi read mapping via a hierarchy of locally-consistent syncmer blocks, forthcoming.


## License
Synpact is available under the XXX.
