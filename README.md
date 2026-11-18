# Synpact

A long-read mapper for **PacBio HiFi** reads built on a hierarchy of
**locally-consistent syncmer blocks**. Instead of seed-and-extend, it places a
read by voting matches of multi-scale blocks that are reproducible between the
read and the reference even in the presence of sequencing errors.

The output is standard PAF: one placement (chromosome, position, strand, MAPQ)
per read, no base-level alignment.

## Installation

Requires a Rust toolchain (`cargo`).

```sh
git clone https://github.com/mahmudsami/synpact.git
cd synpact
cargo build --release
# binary is in: target/release/synpact
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

### 3. Presets
```sh
synpact --build-index genome.fa.gz genome.idx --preset accurate --threads 8
synpact --map reads.fq.gz genome.idx -o out.paf --preset accurate --threads 8
```
| Preset | Same as | Use |
|--------|---------|-----|
| `fast` | `--max-level 6 --min-level 3 --place-ties` | coarse levels only: fast, few confident errors; exact ties reported at MAPQ 0 |
| `accurate` | `--max-level 5 --min-level 2 --mapq0-margin 0.05 --span-switch 0.3,1.1` | one level finer; near-ties reported at MAPQ 0, and a read that a longer locus explains better than the vote-best is placed there |
| `optimized` | `--preset accurate --k 27 --s 24` | accurate with longer, denser syncmers (k 27, density 1/4): fewest wrong reads; about 2.4× the index (human 3.2 GB) and 3× the mapping time of `accurate` |

A preset only fills options that are not on the command line, so an explicit flag
overrides it (e.g. `--preset accurate --span-switch 0.3,1.25`). It also sets the
index's levels (and, for `optimized`, k and s), so build the index with the same preset
(or the same `--max-level` / `--min-level` / `--k` / `--s`) you map with; a prebuilt
index keeps its own k and s, and mapping warns if they differ from the ones asked for. Without a preset synpact runs as before:
L6→3, exact ties left unmapped.

### Options
| Flag | Default | Meaning |
|------|---------|---------|
| `--preset NAME` | none | `fast`, `accurate` or `optimized` (see Presets) |
| `--k N` | 19 | k-mer (syncmer) length; must be ≤ 32 and > `--s` |
| `--s N` | 10 | syncmer s-mer length (density ≈ 1/(k−s+1)) |
| `--threads N` | all cores | worker threads |
| `--min-level N` | 3 | lowest block level both stored in the index and used to anchor reads (one unified floor) |
| `--max-level N` | 6 | highest block level both built into the index and used to anchor reads (one unified ceiling, ≥ `--min-level`); mapping starts at min(N, the index's top level) and falls back down to `--min-level` |
| `--max-occ N` | 10000 | max genomic occurrences for an L2/L3 block before it is skipped as too-repetitive; higher levels scale this down, the L1 cap scales it to 2/5 |
| `--batch-mb N` | 128 | per-batch read-sequence budget (MiB; one batch is mapped while the next is read); keeps peak memory independent of read length. `0` = unbounded |
| `--place-ties` | off | report a read whose best locus ties with a different locus at MAPQ 0 (the first heaviest locus) instead of leaving it unmapped |
| `--mapq0-margin F` | 0 (off) | a read whose runner-up locus has ≥ (1 − F) × the best's votes is reported at MAPQ 0; implies `--place-ties`. `0 ≤ F < 1` |
| `--span-switch F,R` | off | among the loci with ≥ F × the best's votes, report the one whose anchors span ≥ R × more of the read than the vote-best's (MAPQ 0) |
| `--tiebreak R` | `votes` | experimental: which near-tied locus (within `--mapq0-margin`) is reported — `votes` (the heaviest), `floor` (most floor-level support), `chain` (best colinear chain of anchors) |

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
| `config.rs` | runtime flags (`--min-level`, `--max-level`, `--batch-mb`, the tie / margin / switch options) and presets, per-stage profiling, shared helpers |
| `hash.rs` | rolling DNA hash, k-mer encoding, per-level block hashing |
| `syncmer.rs` | open-syncmer selection, `SeedMode` |
| `lcp.rs` | locally-consistent parsing → `Block` / `HierForest` hierarchy, mass accounting |
| `index.rs` | `GIndex` build / serialise / load, FASTA streaming, occurrence caps |
| `fastq.rs` | FASTQ reader, reverse-complement |
| `chain.rs` | anchor collection from the hierarchy, diagonal-voting locus selection, near-tie candidates and their read span |
| `map.rs` | per-read mapping, MAPQ, PAF output, the mapping driver |

## Notes & limitations

- Designed for **HiFi** (≤ ~1 % error) long reads. The k-mer specificity
  that makes placement accurate relies on low error rates; it is not intended
  for ONT or short reads.
- The default `--min-level 3` anchors only on the coarse, unique upper levels,
  which removes most paralog / segmental-duplication ambiguity at HiFi error
  rates. Lower `--min-level` values trade specificity for sensitivity.
- `--max-level` sets the top of the range. An index built with `--max-level N`
  stores only levels up to N (smaller, faster to build) and caps mapping at N by
  itself; against a full index, `--max-level` at mapping time starts anchoring
  lower without rebuilding.
- Placement is locus-only: the output is a position and MAPQ, not a base-level
  alignment. Reads in genuinely ambiguous regions (segmental duplications,
  acrocentric paralogs, tandem-repeat arrays) are reported at low MAPQ; with
  `--place-ties` (every preset) an exactly tied read is reported at MAPQ 0 rather
  than left unmapped.


## Citation

Aydin M.S. and Sahlin, K. synpact: accurate, memory-light PacBio HiFi read mapping via a hierarchy of locally-consistent syncmer blocks, bioRxiv 2026.06.28.735066; doi: [https://doi.org/10.64898/2026.06.28.735066](https://doi.org/10.64898/2026.06.28.735066).


## License
Synpact is available under the XXX.
