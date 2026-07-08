#!/usr/bin/env python3

"""
Aggregate shorter-fragment start positions relative to longer-fragment starts.

Given a BAM file and two fragment lengths, this script:

1. Determines which length is longer and which is shorter.
2. Uses each longer fragment as an anchor.
3. Sets the start of each longer fragment to relative position 0.
4. Builds a relative window:

       -pad ... long_fragment_length + pad - 1

   For example, with long length 167 and pad 200:

       -200 ... 366

5. Adds +1 only at the start position of each longer fragment.
   Since longer fragments are the anchors, this is always relative position 0.

6. Finds shorter fragments whose starts fall within the relative window.
7. Adds +1 only at the start position of each shorter fragment,
   relative to the start of the longer fragment.

Optional coordinate restriction:

    --region chr1
    --region chr1:1000-2000

Regions are interpreted as 1-based inclusive coordinates.
Fragment inclusion is based on the fragment start position.

Output is a TSV with:

    relative_position    long_fragment_start_count    short_fragment_start_count
"""

import argparse
import bisect
import gzip
import re
import sys
from collections import defaultdict
from dataclasses import dataclass

try:
    import pysam
except ImportError:
    sys.exit(
        "ERROR: pysam is required.\n"
        "Install with:\n"
        "  pip install pysam\n"
        "or:\n"
        "  conda install -c bioconda pysam"
    )


@dataclass(frozen=True)
class Region:
    chrom: str
    start0: int | None
    end0: int | None
    original: str


def parse_args():
    parser = argparse.ArgumentParser(
        description="Aggregate shorter-fragment start positions relative to longer-fragment starts."
    )

    parser.add_argument(
        "-b", "--bam",
        required=True,
        help="Input BAM file."
    )

    parser.add_argument(
        "--length-a",
        type=int,
        required=True,
        help="First fragment length."
    )

    parser.add_argument(
        "--length-b",
        type=int,
        required=True,
        help="Second fragment length."
    )

    parser.add_argument(
        "-p", "--pad",
        type=int,
        default=200,
        help=(
            "Bases to extend before the long-fragment start and after the "
            "long-fragment end. Default: 200."
        )
    )

    parser.add_argument(
        "-o", "--out",
        required=True,
        help="Output TSV file. Use .gz for gzip-compressed output."
    )

    parser.add_argument(
        "--region",
        action="append",
        default=None,
        help=(
            "Restrict analysis to a chromosome or genomic interval. "
            "Examples: --region chr1 or --region chr1:1000-2000. "
            "Can be used multiple times. Coordinates are 1-based inclusive."
        )
    )

    parser.add_argument(
        "--chrom",
        action="append",
        default=None,
        help=(
            "Legacy shortcut for whole-chromosome restriction. "
            "Equivalent to --region chrN. Can be used multiple times."
        )
    )

    parser.add_argument(
        "--tolerance",
        type=int,
        default=0,
        help=(
            "Allow fragment lengths within +/- this many bp of the target lengths. "
            "Default: 0, meaning exact lengths only."
        )
    )

    parser.add_argument(
        "--min-mapq",
        type=int,
        default=0,
        help="Minimum MAPQ for the read used to represent the fragment. Default: 0."
    )

    parser.add_argument(
        "--no-proper-pair-required",
        action="store_true",
        help="Do not require reads to be marked as proper pairs."
    )

    parser.add_argument(
        "--include-duplicates",
        action="store_true",
        help="Include duplicate-marked reads. Default: duplicates are skipped."
    )

    parser.add_argument(
        "--include-qcfail",
        action="store_true",
        help="Include QC-fail reads. Default: QC-fail reads are skipped."
    )

    parser.add_argument(
        "--include-secondary",
        action="store_true",
        help="Include secondary alignments. Default: secondary alignments are skipped."
    )

    parser.add_argument(
        "--include-supplementary",
        action="store_true",
        help="Include supplementary alignments. Default: supplementary alignments are skipped."
    )

    return parser.parse_args()


def parse_region_string(region_string):
    """
    Parse region strings such as:

        chr1
        chr1:1000-2000
        chr1:1,000-2,000

    Coordinates are interpreted as 1-based inclusive and converted internally
    to 0-based half-open coordinates.
    """

    region_string = region_string.strip()

    if not region_string:
        raise ValueError("empty region string")

    if ":" not in region_string:
        return Region(
            chrom=region_string,
            start0=None,
            end0=None,
            original=region_string,
        )

    chrom, coord = region_string.split(":", 1)
    coord = coord.replace(",", "")

    match = re.fullmatch(r"(\d+)-(\d+)", coord)

    if not match:
        raise ValueError(
            f"invalid region '{region_string}'. Expected format like chr1 or chr1:1000-2000"
        )

    start1 = int(match.group(1))
    end1 = int(match.group(2))

    if start1 < 1:
        raise ValueError(f"invalid region '{region_string}': start must be >= 1")

    if end1 < start1:
        raise ValueError(f"invalid region '{region_string}': end must be >= start")

    return Region(
        chrom=chrom,
        start0=start1 - 1,
        end0=end1,
        original=region_string,
    )


def collect_regions(args):
    region_strings = []

    if args.region:
        region_strings.extend(args.region)

    if args.chrom:
        region_strings.extend(args.chrom)

    if not region_strings:
        return None

    regions = []

    for region_string in region_strings:
        try:
            regions.append(parse_region_string(region_string))
        except ValueError as e:
            sys.exit(f"ERROR: {e}")

    return regions


def region_contains_fragment_start(region, chrom, frag_start):
    if region is None:
        return True

    if chrom != region.chrom:
        return False

    if region.start0 is not None and frag_start < region.start0:
        return False

    if region.end0 is not None and frag_start >= region.end0:
        return False

    return True


def open_text(path, mode="wt"):
    if path.endswith(".gz"):
        return gzip.open(path, mode)
    return open(path, mode)


def length_matches(observed, target, tolerance):
    return abs(observed - target) <= tolerance


def read_passes_filters(read, args):
    if not read.is_paired:
        return False

    if read.is_unmapped or read.mate_is_unmapped:
        return False

    if read.reference_id != read.next_reference_id:
        return False

    if not args.no_proper_pair_required and not read.is_proper_pair:
        return False

    if read.mapping_quality < args.min_mapq:
        return False

    if read.is_duplicate and not args.include_duplicates:
        return False

    if read.is_qcfail and not args.include_qcfail:
        return False

    if read.is_secondary and not args.include_secondary:
        return False

    if read.is_supplementary and not args.include_supplementary:
        return False

    # Use only the positive-TLEN read so each paired-end fragment is counted once.
    if read.template_length <= 0:
        return False

    return True


def add_point_to_counts(counts, rel_pos, window_start, window_end_exclusive):
    """
    Add +1 at a single relative position if it falls inside the output window.
    """

    if rel_pos < window_start:
        return False

    if rel_pos >= window_end_exclusive:
        return False

    idx = rel_pos - window_start
    counts[idx] += 1
    return True


def get_fragment_key(read):
    """
    Used only when regions are supplied, to avoid double-counting fragments
    if the user gives overlapping regions.
    """

    return (
        read.query_name,
        read.reference_id,
        read.reference_start,
        read.next_reference_id,
        read.next_reference_start,
        read.template_length,
    )


def iter_reads_from_bam(bam, regions):
    """
    If no regions are supplied, stream the whole BAM.

    If regions are supplied, fetch each region. This requires a BAM index.
    """

    if regions is None:
        for read in bam.fetch(until_eof=True):
            yield read, None
        return

    for region in regions:
        try:
            if region.start0 is None:
                iterator = bam.fetch(region.chrom)
            else:
                iterator = bam.fetch(region.chrom, region.start0, region.end0)
        except ValueError as e:
            available = ", ".join(bam.references[:10])
            if len(bam.references) > 10:
                available += ", ..."
            sys.exit(
                f"ERROR: could not fetch region '{region.original}'.\n"
                f"Reason: {e}\n"
                f"Check chromosome naming. Example available contigs: {available}"
            )

        for read in iterator:
            yield read, region


def main():
    args = parse_args()

    if args.length_a <= 0 or args.length_b <= 0:
        sys.exit("ERROR: fragment lengths must be positive integers.")

    if args.length_a == args.length_b:
        sys.exit("ERROR: the two fragment lengths must be different.")

    if args.pad < 0:
        sys.exit("ERROR: --pad must be >= 0.")

    if args.tolerance < 0:
        sys.exit("ERROR: --tolerance must be >= 0.")

    regions = collect_regions(args)

    long_len = max(args.length_a, args.length_b)
    short_len = min(args.length_a, args.length_b)

    if abs(long_len - short_len) <= 2 * args.tolerance:
        sys.stderr.write(
            "WARNING: length ranges overlap because --tolerance is large. "
            "Fragments matching both lengths will be skipped as ambiguous.\n"
        )

    window_start = -args.pad
    window_end_exclusive = long_len + args.pad
    window_size = window_end_exclusive - window_start

    long_start_counts = [0] * window_size
    short_start_counts = [0] * window_size

    long_fragments_by_chrom = defaultdict(list)
    short_fragments_by_chrom = defaultdict(list)

    stats = {
        "reads_seen": 0,
        "reads_passing_filters": 0,
        "long_fragments": 0,
        "short_fragments": 0,
        "ambiguous_length_fragments": 0,
        "other_length_fragments": 0,
        "long_fragment_starts_counted": 0,
        "short_fragment_starts_near_long_fragments": 0,
        "short_fragment_starts_counted": 0,
        "duplicate_region_fragments_skipped": 0,
    }

    seen_fragments = set()

    with pysam.AlignmentFile(args.bam, "rb") as bam:
        if regions is not None and not bam.has_index():
            sys.exit(
                "ERROR: --region requires an indexed BAM.\n"
                "Create an index with:\n"
                f"  samtools index {args.bam}"
            )

        for read, region in iter_reads_from_bam(bam, regions):
            stats["reads_seen"] += 1

            if not read_passes_filters(read, args):
                continue

            chrom = bam.get_reference_name(read.reference_id)

            frag_len = abs(read.template_length)
            frag_start = read.reference_start
            frag_end = frag_start + frag_len

            if region is not None:
                if not region_contains_fragment_start(region, chrom, frag_start):
                    continue

                fragment_key = get_fragment_key(read)
                if fragment_key in seen_fragments:
                    stats["duplicate_region_fragments_skipped"] += 1
                    continue
                seen_fragments.add(fragment_key)

            stats["reads_passing_filters"] += 1

            is_long = length_matches(frag_len, long_len, args.tolerance)
            is_short = length_matches(frag_len, short_len, args.tolerance)

            if is_long and is_short:
                stats["ambiguous_length_fragments"] += 1
                continue

            if is_long:
                long_fragments_by_chrom[chrom].append((frag_start, frag_end, frag_len))
                stats["long_fragments"] += 1

                # Long fragments are the anchors, so their start is always relative position 0.
                counted = add_point_to_counts(
                    long_start_counts,
                    rel_pos=0,
                    window_start=window_start,
                    window_end_exclusive=window_end_exclusive,
                )

                if counted:
                    stats["long_fragment_starts_counted"] += 1

            elif is_short:
                short_fragments_by_chrom[chrom].append((frag_start, frag_end, frag_len))
                stats["short_fragments"] += 1

            else:
                stats["other_length_fragments"] += 1

    # Sort shorter fragments by genomic start for efficient lookup.
    short_index_by_chrom = {}

    for chrom, fragments in short_fragments_by_chrom.items():
        fragments.sort(key=lambda x: x[0])
        starts = [x[0] for x in fragments]
        short_index_by_chrom[chrom] = (starts, fragments)

    # For every longer fragment, find shorter-fragment starts in the relative window.
    for chrom, long_fragments in long_fragments_by_chrom.items():
        if chrom not in short_index_by_chrom:
            continue

        short_starts, short_fragments = short_index_by_chrom[chrom]

        for long_start, long_end, observed_long_len in long_fragments:
            abs_window_start = long_start + window_start
            abs_window_end_exclusive = long_start + window_end_exclusive

            # Since we are now counting only starts, query only shorter fragments
            # whose starts fall inside the window.
            left = bisect.bisect_left(short_starts, abs_window_start)
            right = bisect.bisect_left(short_starts, abs_window_end_exclusive)

            for short_start, short_end, observed_short_len in short_fragments[left:right]:
                rel_short_start = short_start - long_start

                stats["short_fragment_starts_near_long_fragments"] += 1

                counted = add_point_to_counts(
                    short_start_counts,
                    rel_pos=rel_short_start,
                    window_start=window_start,
                    window_end_exclusive=window_end_exclusive,
                )

                if counted:
                    stats["short_fragment_starts_counted"] += 1

    with open_text(args.out, "wt") as out:
        out.write(f"# bam: {args.bam}\n")
        out.write(f"# length_a: {args.length_a}\n")
        out.write(f"# length_b: {args.length_b}\n")
        out.write(f"# longer_fragment_length: {long_len}\n")
        out.write(f"# shorter_fragment_length: {short_len}\n")
        out.write(f"# pad: {args.pad}\n")
        out.write(f"# window_start_relative: {window_start}\n")
        out.write(f"# window_end_relative_exclusive: {window_end_exclusive}\n")
        out.write(f"# output_relative_positions: {window_start} to {window_end_exclusive - 1}\n")
        out.write(f"# count_mode: fragment_starts_only\n")
        out.write(f"# tolerance: {args.tolerance}\n")

        if regions is None:
            out.write("# regions: whole_bam\n")
        else:
            out.write("# regions: " + ",".join(r.original for r in regions) + "\n")
            out.write("# region_coordinate_system: 1-based inclusive\n")
            out.write("# region_filtering: fragment_start_must_fall_inside_region\n")

        for key, value in stats.items():
            out.write(f"# {key}: {value}\n")

        out.write(
            "relative_position\tlong_fragment_start_count\tshort_fragment_start_count\n"
        )

        for i in range(window_size):
            rel_pos = window_start + i
            out.write(
                f"{rel_pos}\t{long_start_counts[i]}\t{short_start_counts[i]}\n"
            )

    sys.stderr.write("Done.\n")
    sys.stderr.write(f"Longer fragment length: {long_len}\n")
    sys.stderr.write(f"Shorter fragment length: {short_len}\n")
    sys.stderr.write(f"Relative output range: {window_start} to {window_end_exclusive - 1}\n")

    if regions is None:
        sys.stderr.write("Regions: whole BAM\n")
    else:
        sys.stderr.write("Regions: " + ", ".join(r.original for r in regions) + "\n")

    sys.stderr.write(f"Long fragments counted: {stats['long_fragments']}\n")
    sys.stderr.write(f"Short fragments counted: {stats['short_fragments']}\n")
    sys.stderr.write(
        f"Short-fragment starts counted near long fragments: "
        f"{stats['short_fragment_starts_counted']}\n"
    )


if __name__ == "__main__":
    main()