#!/usr/bin/env python3
import argparse
import gzip
from collections import Counter, defaultdict
from tqdm import tqdm


def open_maybe_gzip(path):
    return gzip.open(path, "rt") if path.endswith(".gz") else open(path, "r")


def bed_midpoints_by_chrom(path):
    out = defaultdict(list)
    with open_maybe_gzip(path) as f:
        for line in f:
            if not line or line.startswith("#") or line.startswith("track") or line.startswith("browser"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 3:
                continue
            try:
                start = int(parts[1])
                end = int(parts[2])
            except ValueError:
                continue
            mid = (start + end) // 2
            out[parts[0]].append(mid)
    return out


def nearest_signed_distances_counts_linear(smaller, larger, show_tqdm=False, desc=""):
    """
    smaller and larger sorted.
    Returns Counter of signed distances.
    """

    c = Counter()
    if not smaller or not larger:
        return c

    j = 0
    L = len(larger)

    iterator = smaller
    if show_tqdm:
        iterator = tqdm(smaller, desc=desc, leave=False)

    for x in iterator:
        while j < L and larger[j] < x:
            j += 1

        best_dist = None
        best_abs = None

        if j < L:
            d = larger[j] - x
            best_dist = d
            best_abs = abs(d)

        if j > 0:
            d2 = larger[j - 1] - x
            a2 = abs(d2)
            if best_dist is None or a2 < best_abs:
                best_dist = d2
                best_abs = a2
            elif a2 == best_abs:
                if best_dist < 0 and d2 > 0:
                    best_dist = d2

        if best_dist is not None:
            c[best_dist] += 1

    return c


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("-a", "--bed-a", required=True)
    ap.add_argument("-b", "--bed-b", required=True)
    ap.add_argument("-o", "--out", required=True)
    ap.add_argument("--no-auto-smaller", action="store_true")
    ap.add_argument("--require-same-chrom", action="store_true")
    ap.add_argument("--inner-tqdm-threshold", type=int, default=50000,
                    help="Show inner tqdm when queries per chrom exceed this (default 50k)")
    args = ap.parse_args()

    print("Loading BED A…")
    A = bed_midpoints_by_chrom(args.bed_a)

    print("Loading BED B…")
    B = bed_midpoints_by_chrom(args.bed_b)

    print("Sorting midpoints…")
    for chrom in A:
        A[chrom].sort()
    for chrom in B:
        B[chrom].sort()

    if args.require_same_chrom:
        chroms = sorted(set(A) & set(B))
    else:
        chroms = sorted(set(A) | set(B))

    print(f"Processing {len(chroms)} chromosomes")

    total = Counter()

    for chrom in tqdm(chroms, desc="Chromosomes"):
        a = A.get(chrom, [])
        b = B.get(chrom, [])
        if not a or not b:
            continue

        if args.no_auto_smaller:
            queries = a
            targets = b
        else:
            if len(a) <= len(b):
                queries = a
                targets = b
            else:
                queries = b
                targets = a

        show_inner = len(queries) >= args.inner_tqdm_threshold

        total += nearest_signed_distances_counts_linear(
            queries,
            targets,
            show_tqdm=show_inner,
            desc=f"{chrom} ({len(queries)} queries)"
        )

    print("Writing output…")

    with open(args.out, "w") as out:
        out.write("distance\tcount\n")
        for d in sorted(total):
            out.write(f"{d}\t{total[d]}\n")

    print("Done.")


if __name__ == "__main__":
    main()
