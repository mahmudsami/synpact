#!/usr/bin/env python3
"""
Evaluate PAF mappings against a TSV truth file.
Truth TSV format: read_name\tchr\tstart\tend  (header line optional)

PAF format:
  readname qlen qst qen strand chrname rlen rstart rend matches alen mapq ...

Usage: python3 eval_tsv.py <truth.tsv> <mapped.paf> [tolerance=1000]
"""

import sys
import gzip


def load_truth_from_tsv(tsv_path: str) -> dict:
    """Return dict: read_name -> (chr, start, end)."""
    truth = {}
    opener = gzip.open if tsv_path.endswith('.gz') else open
    with opener(tsv_path, 'rt') as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith('#') or line.startswith('read_name'):
                continue
            parts = line.split('\t')
            if len(parts) < 4:
                continue
            rname, chrom, start, end = parts[0], parts[1], int(parts[2]), int(parts[3])
            truth[rname] = (chrom, start, end)
    return truth


def eval_paf(paf_file: str, truth: dict, tol: int = 1000):
    best = {}   # rname -> (target, rstart, mapq) or None

    opener = gzip.open if paf_file.endswith('.gz') else open
    with opener(paf_file, 'rt') as f:
        for line in f:
            if not line.strip():
                continue
            parts = line.split('\t')
            rname = parts[0]

            if rname not in truth:
                base = rname.rstrip('/12').rstrip('/')
                if base in truth:
                    rname = base
                else:
                    continue

            if len(parts) < 12:
                if rname not in best:
                    best[rname] = None
                continue

            target = parts[5]
            if target == '*':
                if rname not in best:
                    best[rname] = None
                continue

            rstart = int(parts[7])
            mapq   = int(parts[11])

            cur = best.get(rname)
            if cur is None or mapq > cur[2]:
                best[rname] = (target, rstart, mapq)

    # Stats
    mapped = unmapped = correct = wrong_chr = wrong_pos = low_mapq = 0

    for rname, entry in best.items():
        if rname not in truth:
            continue
        true_chr, true_start, true_end = truth[rname]
        if entry is None:
            unmapped += 1
            continue

        target, rstart, mapq = entry
        mapped += 1

        if mapq == 0:
            low_mapq += 1

        if target != true_chr:
            wrong_chr += 1
            continue

        if true_start - tol <= rstart <= true_end + tol:
            correct += 1
        else:
            wrong_pos += 1

    # Reads in truth not seen in PAF at all → unmapped
    for rname in truth:
        if rname not in best:
            unmapped += 1

    total = mapped + unmapped
    if total == 0:
        print("No reads found.")
        return

    accuracy      = 100.0 * correct  / total
    precision     = 100.0 * correct  / mapped if mapped else 0.0
    map_rate      = 100.0 * mapped   / total
    unmapped_pct  = 100.0 * unmapped / total
    wrong_chr_pct = 100.0 * wrong_chr / total

    print(f"  Total reads     : {total:>8,}")
    print(f"  Mapped          : {mapped:>8,}  ({map_rate:.2f}%)")
    print(f"  Unmapped        : {unmapped:>8,}  ({unmapped_pct:.2f}%)")
    print(f"  Correct (±{tol}bp): {correct:>8,}  ({accuracy:.2f}% accuracy)")
    print(f"  Wrong chr       : {wrong_chr:>8,}  ({wrong_chr_pct:.2f}%)")
    print(f"  Wrong position  : {wrong_pos:>8,}")
    print(f"  MAPQ=0 mapped   : {low_mapq:>8,}")
    print(f"  Precision       : {precision:.2f}%")


if __name__ == '__main__':
    if len(sys.argv) < 3:
        print(f"Usage: {sys.argv[0]} <truth.tsv> <mapped.paf> [tolerance=1000]")
        sys.exit(1)

    tsv_file = sys.argv[1]
    paf_file = sys.argv[2]
    tol = int(sys.argv[3]) if len(sys.argv) > 3 else 1000

    print(f"Loading truth from {tsv_file} …")
    truth = load_truth_from_tsv(tsv_file)
    print(f"  Loaded {len(truth):,} ground truth positions")

    print(f"\nEvaluating {paf_file} …")
    eval_paf(paf_file, truth, tol=tol)
