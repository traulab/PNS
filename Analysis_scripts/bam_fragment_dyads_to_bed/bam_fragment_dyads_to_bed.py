#!/usr/bin/env python3

import argparse
import glob
import os
import re
import sys
from collections import Counter, defaultdict

import pysam


CIGAR_SOFT_HARD_OR_PAD = {4, 5, 6}  # S, H, P


def has_softclip_or_hardclip_or_padding(cigartuples):
    if not cigartuples:
        return False
    return any(op in CIGAR_SOFT_HARD_OR_PAD for op, _length in cigartuples)


def expand_bam_inputs(bam_inputs):
    bam_files = []

    for item in bam_inputs:
        matches = sorted(glob.glob(item))
        bam_files.extend(matches if matches else [item])

    bam_files = list(dict.fromkeys(bam_files))
    missing = [path for path in bam_files if not os.path.exists(path)]

    if missing:
        raise FileNotFoundError(
            "These BAM files or patterns were not found:\n  "
            + "\n  ".join(missing)
        )

    return bam_files


def matching_contig(chrom, references):
    if chrom in references:
        return chrom

    if chrom.startswith("chr"):
        without_chr = chrom[3:]
        if without_chr in references:
            return without_chr
    else:
        with_chr = "chr" + chrom
        if with_chr in references:
            return with_chr

    return None


def parse_region(region_text):
    """
    Parse a 1-based inclusive region such as:
        chr17:41,186,312-41,287,500

    Returns:
        chrom, start_0based, end_0based_exclusive
    """
    if region_text is None:
        return None

    cleaned = region_text.replace(",", "").replace(" ", "")
    match = re.fullmatch(r"([^:]+):(\d+)-(\d+)", cleaned)

    if not match:
        raise ValueError(
            "Invalid --region format. Use, for example, "
            "chr17:41,186,312-41,287,500"
        )

    chrom, start_text, end_text = match.groups()
    start_1based = int(start_text)
    end_1based = int(end_text)

    if start_1based < 1:
        raise ValueError("Region start must be at least 1.")

    if end_1based < start_1based:
        raise ValueError("Region end must be greater than or equal to region start.")

    return chrom, start_1based - 1, end_1based


def dyad_positions(frag_start, frag_len, mode):
    """
    Return one or two 0-based dyad positions.

    For an odd fragment length, all modes return the single central base.

    For an even fragment length:
      floor = left/lower middle base
      ceil  = right/higher middle base
      split = both middle bases
    """
    centre_2x = (2 * frag_start) + frag_len - 1
    floor_pos = centre_2x // 2
    ceil_pos = (centre_2x + 1) // 2

    if mode == "floor":
        return (floor_pos,)
    if mode == "ceil":
        return (ceil_pos,)
    if mode == "split":
        return (floor_pos,) if floor_pos == ceil_pos else (floor_pos, ceil_pos)

    raise ValueError(f"Unexpected dyad mode: {mode}")


def read_passes_filters(
    aln,
    mapq,
    include_duplicates,
    allow_improper_pairs,
    allow_softclipped,
):
    if aln.is_unmapped or aln.mate_is_unmapped:
        return False

    if aln.is_secondary or aln.is_supplementary:
        return False

    if aln.is_qcfail:
        return False

    if not include_duplicates and aln.is_duplicate:
        return False

    if not allow_improper_pairs and not aln.is_proper_pair:
        return False

    if not allow_softclipped and has_softclip_or_hardclip_or_padding(aln.cigartuples):
        return False

    if aln.mapping_quality < mapq:
        return False

    if aln.reference_id != aln.next_reference_id:
        return False

    # Count each paired-end fragment once.
    if aln.template_length <= 0:
        return False

    return True


def iter_fragments(
    bam,
    bam_chrom,
    fragment_length,
    fetch_start,
    fetch_end,
    mapq,
    include_duplicates,
    allow_improper_pairs,
    allow_softclipped,
):
    try:
        iterator = bam.fetch(bam_chrom, fetch_start, fetch_end)
    except ValueError:
        return

    for aln in iterator:
        if not read_passes_filters(
            aln=aln,
            mapq=mapq,
            include_duplicates=include_duplicates,
            allow_improper_pairs=allow_improper_pairs,
            allow_softclipped=allow_softclipped,
        ):
            continue

        frag_len = abs(aln.template_length)

        if frag_len != fragment_length:
            continue

        frag_start = aln.reference_start
        frag_end = frag_start + frag_len

        yield frag_start, frag_end


def main():
    parser = argparse.ArgumentParser(
        description=(
            "Extract dyad positions from exact-length paired-end BAM fragments and "
            "write a dense four-column BED file: chromosome, start, end, occurrence. "
            "Every genomic position in the requested region is written, including zero occurrences. "
            "Coordinate deduplication is performed separately within each BAM."
        )
    )

    parser.add_argument(
        "--bam",
        nargs="+",
        required=True,
        help="Input BAM file(s), explicit paths or wildcard patterns.",
    )

    parser.add_argument(
        "--length",
        type=int,
        required=True,
        help="Exact paired-end fragment length to retain.",
    )

    parser.add_argument(
        "--region",
        default=None,
        help=(
            "Optional 1-based inclusive region, for example "
            "chr17:41,186,312-41,287,500. Only dyads within the region are output."
        ),
    )

    parser.add_argument(
        "--out",
        required=True,
        help="Output four-column BED file.",
    )

    parser.add_argument(
        "--mapq",
        type=int,
        default=30,
        help="Minimum mapping quality. Default: 30.",
    )

    parser.add_argument(
        "--max-duplicates",
        type=int,
        default=0,
        help=(
            "Maximum additional fragments with identical start/end coordinates to "
            "retain within each BAM. Default 0 keeps one fragment per coordinate. "
            "Use -1 to disable coordinate deduplication. As in the source script, "
            "N retains up to N+1 fragments at the same coordinate."
        ),
    )

    parser.add_argument(
        "--include-duplicates",
        action="store_true",
        help="Include reads marked as duplicates in the BAM. Default: skip them.",
    )

    parser.add_argument(
        "--allow-improper-pairs",
        action="store_true",
        help="Allow pairs not marked as proper. Default: require proper pairs.",
    )

    parser.add_argument(
        "--allow-softclipped",
        action="store_true",
        help="Allow soft-clipped, hard-clipped, or padded alignments.",
    )

    parser.add_argument(
        "--dyad-mode",
        choices=("floor", "ceil", "split"),
        default="floor",
        help=(
            "How to represent the centre of even-length fragments. "
            "floor uses the left middle base, ceil the right middle base, and split "
            "writes both. Odd-length fragments always have one dyad. Default: floor."
        ),
    )

    args = parser.parse_args()

    if args.length < 1:
        sys.exit("ERROR: --length must be at least 1.")

    if args.mapq < 0:
        sys.exit("ERROR: --mapq must be at least 0.")

    try:
        region = parse_region(args.region)
        bam_files = expand_bam_inputs(args.bam)
    except (ValueError, FileNotFoundError) as exc:
        sys.exit(f"ERROR: {exc}")

    out_parent = os.path.dirname(args.out)
    if out_parent:
        os.makedirs(out_parent, exist_ok=True)

    total_occurrences = Counter()
    total_matching_fragments = 0
    total_retained_fragments = 0

    print("[INFO] BAM files:", file=sys.stderr)
    for bam_path in bam_files:
        print(f"  {bam_path}", file=sys.stderr)

    if region is not None:
        requested_chrom, region_start, region_end = region
        print(
            f"[INFO] region: {requested_chrom}:{region_start + 1:,}-{region_end:,} "
            "(1-based inclusive input)",
            file=sys.stderr,
        )

    for bam_path in bam_files:
        print(f"[INFO] Processing BAM: {bam_path}", file=sys.stderr)

        with pysam.AlignmentFile(bam_path, "rb") as bam:
            references = set(bam.references)

            if region is None:
                contigs = list(bam.references)
            else:
                requested_chrom, region_start, region_end = region
                bam_chrom = matching_contig(requested_chrom, references)

                if bam_chrom is None:
                    print(
                        f"[WARNING] Region chromosome {requested_chrom} was not found "
                        f"in {os.path.basename(bam_path)}; skipping this BAM.",
                        file=sys.stderr,
                    )
                    continue

                contigs = [bam_chrom]

            # This counter is intentionally reset for every BAM, preserving
            # independent observations at identical coordinates in different BAMs.
            coord_counts = defaultdict(int)

            for bam_chrom in contigs:
                chrom_len = bam.get_reference_length(bam_chrom)

                if region is None:
                    fetch_start = 0
                    fetch_end = chrom_len
                    region_start_for_filter = 0
                    region_end_for_filter = chrom_len
                else:
                    # Expand the fetch interval to capture fragments whose alignment
                    # begins before the requested region but whose dyad falls inside it.
                    fetch_start = max(0, region_start - args.length)
                    fetch_end = min(chrom_len, region_end + args.length)
                    region_start_for_filter = region_start
                    region_end_for_filter = region_end

                for frag_start, frag_end in iter_fragments(
                    bam=bam,
                    bam_chrom=bam_chrom,
                    fragment_length=args.length,
                    fetch_start=fetch_start,
                    fetch_end=fetch_end,
                    mapq=args.mapq,
                    include_duplicates=args.include_duplicates,
                    allow_improper_pairs=args.allow_improper_pairs,
                    allow_softclipped=args.allow_softclipped,
                ):
                    total_matching_fragments += 1

                    coord_key = (bam_chrom, frag_start, frag_end)

                    if args.max_duplicates >= 0:
                        if coord_counts[coord_key] > args.max_duplicates:
                            continue
                        coord_counts[coord_key] += 1

                    dyads = dyad_positions(
                        frag_start=frag_start,
                        frag_len=args.length,
                        mode=args.dyad_mode,
                    )

                    retained_any = False

                    for dyad in dyads:
                        if not (region_start_for_filter <= dyad < region_end_for_filter):
                            continue

                        total_occurrences[(bam_chrom, dyad)] += 1
                        retained_any = True

                    if retained_any:
                        total_retained_fragments += 1

    if not total_occurrences:
        print(
            "[WARNING] no dyad positions passed the requested filters; "
            "the output will contain zeros only.",
            file=sys.stderr,
        )

    contig_order = {}
    next_index = 0

    for bam_path in bam_files:
        with pysam.AlignmentFile(bam_path, "rb") as bam:
            for contig in bam.references:
                if contig not in contig_order:
                    contig_order[contig] = next_index
                    next_index += 1

    with open(args.out, "w") as out:
        if region is not None:
            requested_chrom, region_start, region_end = region

            # Use the contig naming style found in the first BAM containing the region.
            output_chrom = None
            for bam_path in bam_files:
                with pysam.AlignmentFile(bam_path, "rb") as bam:
                    output_chrom = matching_contig(requested_chrom, set(bam.references))
                    if output_chrom is not None:
                        break

            if output_chrom is None:
                sys.exit(f"ERROR: region chromosome {requested_chrom} was not found in any BAM.")

            for position in range(region_start, region_end):
                occurrence = total_occurrences.get((output_chrom, position), 0)
                out.write(f"{output_chrom}\t{position}\t{position + 1}\t{occurrence}\n")
        else:
            # Without --region, write every base of every BAM contig. This can produce
            # a very large file for a whole genome.
            written_contigs = set()

            for bam_path in bam_files:
                with pysam.AlignmentFile(bam_path, "rb") as bam:
                    for chrom, chrom_len in zip(bam.references, bam.lengths):
                        if chrom in written_contigs:
                            continue
                        written_contigs.add(chrom)

                        for position in range(chrom_len):
                            occurrence = total_occurrences.get((chrom, position), 0)
                            out.write(f"{chrom}\t{position}\t{position + 1}\t{occurrence}\n")

    print(f"[DONE] output: {args.out}", file=sys.stderr)
    print(
        f"[INFO] exact-length fragments seen before coordinate deduplication: "
        f"{total_matching_fragments:,}",
        file=sys.stderr,
    )
    print(
        f"[INFO] retained fragments contributing at least one dyad: "
        f"{total_retained_fragments:,}",
        file=sys.stderr,
    )
    print(
        f"[INFO] non-zero dyad positions: {len(total_occurrences):,}",
        file=sys.stderr,
    )
    print(
        f"[INFO] total BED occurrence count: {sum(total_occurrences.values()):,}",
        file=sys.stderr,
    )
    print(
        "[INFO] Deduplication was applied independently within each BAM. "
        "Therefore, the same fragment coordinates in different BAM files contribute "
        "separate occurrences.",
        file=sys.stderr,
    )


if __name__ == "__main__":
    main()
