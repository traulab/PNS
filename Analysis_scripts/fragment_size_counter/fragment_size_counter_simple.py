#!/usr/bin/env python3
"""
Count paired-end fragment lengths from a BAM file.

This version is BAM-only:
- no BED file
- no chromatin states
- no human/autosome assumptions
- counts across all contigs in the BAM
- counts each paired fragment once using read1 and abs(TLEN)

Output TSV columns:
fragment_length    count
"""

import argparse
import os
from collections import Counter
import pysam


def count_fragment_lengths(
    bam_path: str,
    mapq: int = 0,
    include_duplicates: bool = False,
    require_proper_pair: bool = True,
    min_length: int = 1,
    max_length: int | None = None,
) -> Counter:
    counts = Counter()

    with pysam.AlignmentFile(bam_path, "rb") as bam:
        # until_eof=True makes this work without needing a BAM index.
        # It also naturally iterates over all contigs in the BAM.
        for read in bam.fetch(until_eof=True):
            if read.is_unmapped or read.mate_is_unmapped:
                continue
            if read.is_secondary or read.is_supplementary:
                continue
            if read.is_qcfail:
                continue
            if read.mapping_quality < mapq:
                continue
            if read.is_duplicate and not include_duplicates:
                continue
            if require_proper_pair and not read.is_proper_pair:
                continue

            # Count each physical fragment once.
            if not read.is_read1:
                continue

            frag_len = abs(read.template_length)
            if frag_len < min_length:
                continue
            if max_length is not None and frag_len > max_length:
                continue

            counts[frag_len] += 1

    return counts


def write_counts(counts: Counter, out_path: str) -> None:
    with open(out_path, "w") as out:
        out.write("fragment_length\tcount\n")
        for length in sorted(counts):
            out.write(f"{length}\t{counts[length]}\n")


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Count paired-end fragment lengths across all contigs in a BAM."
    )
    parser.add_argument(
        "-b",
        "--bam",
        required=True,
        help="Input paired-end BAM file.",
    )
    parser.add_argument(
        "-o",
        "--out",
        default=None,
        help=(
            "Output TSV file. "
            "Default: <BAM_basename>_fragment_length.tsv"
        ),
    )
    parser.add_argument(
        "--mapq",
        type=int,
        default=0,
        help="Minimum MAPQ. Default: 0.",
    )
    parser.add_argument(
        "--include-duplicates",
        action="store_true",
        help="Include duplicate-marked reads. Default: exclude duplicates.",
    )
    parser.add_argument(
        "--allow-improper-pairs",
        action="store_true",
        help="Do not require the proper-pair SAM flag. Default: require proper pairs.",
    )
    parser.add_argument(
        "--min-length",
        type=int,
        default=1,
        help="Minimum fragment length to count. Default: 1.",
    )
    parser.add_argument(
        "--max-length",
        type=int,
        default=None,
        help="Maximum fragment length to count. Default: no maximum.",
    )

    args = parser.parse_args()

    # Derive default output filename from BAM basename.
    if args.out is None:
        bam_basename = os.path.basename(args.bam)

        if bam_basename.lower().endswith(".bam"):
            bam_basename = bam_basename[:-4]
        else:
            bam_basename = os.path.splitext(bam_basename)[0]

        args.out = f"{bam_basename}_fragment_length.tsv"

    counts = count_fragment_lengths(
        bam_path=args.bam,
        mapq=args.mapq,
        include_duplicates=args.include_duplicates,
        require_proper_pair=not args.allow_improper_pairs,
        min_length=args.min_length,
        max_length=args.max_length,
    )

    write_counts(counts, args.out)

    total = sum(counts.values())
    print(f"Wrote: {args.out}")
    print(f"Fragments counted: {total:,}")
    print(f"Distinct fragment lengths: {len(counts):,}")


if __name__ == "__main__":
    main()

