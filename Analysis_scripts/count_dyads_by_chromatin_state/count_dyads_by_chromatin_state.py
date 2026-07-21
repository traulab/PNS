#!/usr/bin/env python3

"""
Count dyads within chromatin-state regions, genome-wide and per chromosome.

Input BED columns:

    1. chromosome
    2. start, 0-based
    3. end, exclusive
    4. chromatin-state name

The script produces three wide-format TSV files:

    PREFIX_dyad_counts.tsv
    PREFIX_state_lengths.tsv
    PREFIX_dyad_density_per_kb.tsv

Each output contains:

    column 1: Chromosome
    remaining columns: one column per chromatin state

Rows include:

    Total
    chr1
    chr2
    ...

If --output-prefix is omitted, the prefix is automatically generated from
the bigWig and BED basenames:

    BIGWIG_BASENAME__BED_BASENAME

Example:

    python3 count_dyads_by_chromatin_state.py \
        --bigwig BH01_chrAll_PNS_mode161_lower161_upper161_dyad.bw \
        --bed chromatin_states.bed
"""

import argparse
import gzip
import math
import re
import sys
from collections import OrderedDict
from pathlib import Path

try:
    import pyBigWig
except ImportError:
    sys.exit(
        "ERROR: pyBigWig is not installed.\n"
        "Install it with:\n"
        "    conda install -c bioconda pybigwig\n"
        "or:\n"
        "    pip install pyBigWig"
    )


def open_text_file(filename):
    """Open a plain-text or gzip-compressed text file."""
    if str(filename).endswith(".gz"):
        return gzip.open(filename, "rt")

    return open(filename, "r")


def remove_known_extensions(filename):
    """
    Remove common genomics file extensions.

    Examples:

        sample.bw       -> sample
        sample.bigWig   -> sample
        states.bed      -> states
        states.bed.gz   -> states
    """
    name = Path(filename).name

    known_extensions = (
        ".bedgraph.gz",
        ".bedGraph.gz",
        ".bed.gz",
        ".bigwig",
        ".bigWig",
        ".bedgraph",
        ".bedGraph",
        ".bw",
        ".bed",
        ".gz",
    )

    for extension in known_extensions:
        if name.endswith(extension):
            return name[:-len(extension)]

    return Path(name).stem


def chromosome_sort_key(chrom):
    """
    Sort chromosomes naturally.

    Examples:

        chr1, chr2, ..., chr22, chrX, chrY, chrM
    """
    chrom_without_prefix = re.sub(
        r"^chr",
        "",
        chrom,
        flags=re.IGNORECASE,
    )

    if chrom_without_prefix.isdigit():
        return 0, int(chrom_without_prefix)

    special_order = {
        "X": 23,
        "Y": 24,
        "M": 25,
        "MT": 25,
    }

    upper_name = chrom_without_prefix.upper()

    if upper_name in special_order:
        return 0, special_order[upper_name]

    return 1, chrom_without_prefix


def sum_bigwig_interval(bigwig, chrom, start, end):
    """
    Sum the bigWig signal across an interval.

    Each bigWig value is multiplied by the number of overlapping bases.

    For a 1-bp dyad entry:

        dyad value × 1 bp = dyad value

    This also correctly handles split dyads such as two adjacent positions
    each carrying a value of 0.5.
    """
    intervals = bigwig.intervals(chrom, start, end)

    if intervals is None:
        return 0.0

    total = 0.0

    for signal_start, signal_end, value in intervals:
        if not math.isfinite(value):
            continue

        overlap_start = max(start, signal_start)
        overlap_end = min(end, signal_end)

        if overlap_end > overlap_start:
            overlap_length = overlap_end - overlap_start
            total += value * overlap_length

    return total


def initialise_entry(summary, chrom, state):
    """Ensure that a chromosome/state summary entry exists."""
    if chrom not in summary:
        summary[chrom] = OrderedDict()

    if state not in summary[chrom]:
        summary[chrom][state] = {
            "length": 0,
            "dyad_count": 0.0,
        }


def add_to_summary(summary, chrom, state, length, dyad_count):
    """Add region length and dyad count to a summary entry."""
    initialise_entry(summary, chrom, state)

    summary[chrom][state]["length"] += length
    summary[chrom][state]["dyad_count"] += dyad_count


def format_count(value):
    """
    Format dyad counts.

    Whole-number counts are written as integers. Fractional values are
    retained when split dyads or weighted values are present.
    """
    if math.isclose(value, round(value), abs_tol=1e-9):
        return str(round(value))

    return f"{value:.10g}"


def format_density(value):
    """Format dyad-density values."""
    return f"{value:.10g}"


def write_wide_output(
    filename,
    chromosomes,
    states,
    summary,
    value_type,
):
    """
    Write a wide-format TSV file.

    value_type may be:

        dyad_count
        length
        density_per_kb
    """
    with open(filename, "w") as output_handle:
        output_handle.write(
            "Chromosome\t" + "\t".join(states) + "\n"
        )

        output_rows = ["Total"] + chromosomes

        for chrom in output_rows:
            values = []

            for state in states:
                entry = summary.get(chrom, {}).get(
                    state,
                    {
                        "length": 0,
                        "dyad_count": 0.0,
                    },
                )

                region_length = entry["length"]
                dyad_count = entry["dyad_count"]

                if value_type == "dyad_count":
                    value = format_count(dyad_count)

                elif value_type == "length":
                    value = str(region_length)

                elif value_type == "density_per_kb":
                    if region_length > 0:
                        density = (
                            dyad_count
                            / region_length
                            * 1000.0
                        )
                    else:
                        density = 0.0

                    value = format_density(density)

                else:
                    raise ValueError(
                        f"Unknown output value type: {value_type}"
                    )

                values.append(value)

            output_handle.write(
                chrom + "\t" + "\t".join(values) + "\n"
            )


def parse_arguments():
    parser = argparse.ArgumentParser(
        description=(
            "Calculate dyad counts, chromatin-state lengths and dyad "
            "densities genome-wide and per chromosome."
        )
    )

    parser.add_argument(
        "--bigwig",
        required=True,
        help="Input bigWig containing dyad positions and counts.",
    )

    parser.add_argument(
        "--bed",
        required=True,
        help=(
            "Input BED or BED.gz file. Columns 1-3 contain chromosome, "
            "start and end; column 4 contains the chromatin-state name."
        ),
    )

    parser.add_argument(
        "--output-prefix",
        default=None,
        help=(
            "Optional output prefix. If omitted, the prefix is generated "
            "automatically from the bigWig and BED basenames."
        ),
    )

    return parser.parse_args()


def main():
    args = parse_arguments()

    bigwig_path = Path(args.bigwig)
    bed_path = Path(args.bed)

    if not bigwig_path.is_file():
        sys.exit(
            f"ERROR: bigWig file not found: {bigwig_path}"
        )

    if not bed_path.is_file():
        sys.exit(
            f"ERROR: BED file not found: {bed_path}"
        )

    if args.output_prefix is None:
        bigwig_basename = remove_known_extensions(bigwig_path)
        bed_basename = remove_known_extensions(bed_path)

        output_prefix = (
            f"{bigwig_basename}__{bed_basename}"
        )
    else:
        output_prefix = args.output_prefix

    bw = pyBigWig.open(str(bigwig_path))

    if bw is None:
        sys.exit(
            f"ERROR: Could not open bigWig: {bigwig_path}"
        )

    bigwig_chroms = bw.chroms()

    # summary[chromosome][chromatin_state]
    summary = OrderedDict()

    # Preserve chromatin-state order from the BED file.
    states = []
    seen_states = set()

    chromosomes_seen = set()

    region_count = 0
    missing_chromosome_regions = 0
    clipped_regions = 0

    with open_text_file(bed_path) as bed_handle:
        for line_number, line in enumerate(
            bed_handle,
            start=1,
        ):
            line = line.strip()

            if not line or line.startswith("#"):
                continue

            fields = line.split()

            if len(fields) < 4:
                bw.close()
                sys.exit(
                    f"ERROR: BED line {line_number} contains fewer "
                    f"than four columns:\n{line}"
                )

            chrom = fields[0]
            state = fields[3]

            try:
                start = int(fields[1])
                end = int(fields[2])
            except ValueError:
                bw.close()
                sys.exit(
                    f"ERROR: Invalid BED coordinates on line "
                    f"{line_number}:\n{line}"
                )

            if start < 0:
                bw.close()
                sys.exit(
                    f"ERROR: BED start is negative on line "
                    f"{line_number}:\n{line}"
                )

            if end <= start:
                bw.close()
                sys.exit(
                    f"ERROR: BED end must be greater than start on "
                    f"line {line_number}:\n{line}"
                )

            if state not in seen_states:
                seen_states.add(state)
                states.append(state)

            chromosomes_seen.add(chrom)

            original_region_length = end - start
            dyad_count = 0.0
            region_count += 1

            if chrom not in bigwig_chroms:
                missing_chromosome_regions += 1

            else:
                chrom_length = bigwig_chroms[chrom]

                clipped_start = max(0, start)
                clipped_end = min(end, chrom_length)

                if (
                    clipped_start != start
                    or clipped_end != end
                ):
                    clipped_regions += 1

                if clipped_end > clipped_start:
                    dyad_count = sum_bigwig_interval(
                        bw,
                        chrom,
                        clipped_start,
                        clipped_end,
                    )

            # Retain the full BED interval length in the denominator.
            add_to_summary(
                summary,
                chrom,
                state,
                original_region_length,
                dyad_count,
            )

            # Genome-wide total.
            add_to_summary(
                summary,
                "Total",
                state,
                original_region_length,
                dyad_count,
            )

    bw.close()

    chromosomes = sorted(
        chromosomes_seen,
        key=chromosome_sort_key,
    )

    counts_output = (
        f"{output_prefix}_dyad_counts.tsv"
    )

    lengths_output = (
        f"{output_prefix}_state_lengths.tsv"
    )

    density_output = (
        f"{output_prefix}_dyad_density_per_kb.tsv"
    )

    write_wide_output(
        counts_output,
        chromosomes,
        states,
        summary,
        value_type="dyad_count",
    )

    write_wide_output(
        lengths_output,
        chromosomes,
        states,
        summary,
        value_type="length",
    )

    write_wide_output(
        density_output,
        chromosomes,
        states,
        summary,
        value_type="density_per_kb",
    )

    print(
        f"Processed BED regions: {region_count}",
        file=sys.stderr,
    )

    print(
        f"Chromatin states: {len(states)}",
        file=sys.stderr,
    )

    print(
        f"Chromosomes: {len(chromosomes)}",
        file=sys.stderr,
    )

    print(
        f"Dyad counts: {counts_output}",
        file=sys.stderr,
    )

    print(
        f"State lengths: {lengths_output}",
        file=sys.stderr,
    )

    print(
        f"Dyad density per kb: {density_output}",
        file=sys.stderr,
    )

    if missing_chromosome_regions:
        print(
            f"Warning: {missing_chromosome_regions} BED regions were "
            f"on chromosomes absent from the bigWig. Their lengths were "
            f"retained, but their dyad counts were set to zero.",
            file=sys.stderr,
        )

    if clipped_regions:
        print(
            f"Warning: {clipped_regions} BED regions extended beyond "
            f"the corresponding bigWig chromosome length and were "
            f"clipped when counting dyads.",
            file=sys.stderr,
        )


if __name__ == "__main__":
    main()