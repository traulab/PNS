#!/usr/bin/env python3
"""
Calculate distance autocorrelation (DAC) of dyad occurrence signal inside
strand-aware Alu-flanking regions, reading dyads directly from bigWig files.

Region definition from Alu BED/input:
  + strand: [Alu_start, Alu_start + extend)
  - strand: [Alu_end - extend, Alu_end)

Input dyad signal:
  --bigwig one_or_more_dyad_signal.bw

For every Alu-derived region, the script extracts bigWig signal with
bigWigToBedGraph, expands bedGraph intervals to per-base positions, and adds:
  dac[distance] += value_i * value_j
for every pair of non-zero positions separated by < dmax.

Notes:
  - DAC distances are orientation-independent, so reversing - strand intervals
    does not change the DAC. The strand controls which Alu-flanking interval is
    selected.
  - If you want only + or only - Alus, use --strands plus/minus.
  - Use --strand_col to specify which BED/input column contains + or -.
    Column numbering is 1-based, so standard BED strand is column 6.
  - Requires UCSC bigWigToBedGraph in PATH.
"""

import argparse
import glob
import os
import subprocess
import tempfile

import numpy as np
import pandas as pd
from tqdm import tqdm


def sanitize_filename(name):
    return str(name).replace("/", "_").replace("\\", "_").replace(" ", "_")


def standardize_chromosome_name(name):
    name = str(name).strip()

    if name.startswith("chr"):
        return name
    if name.isdigit():
        return f"chr{name}"
    if name in {"X", "Y"}:
        return f"chr{name}"
    if name in {"M", "MT"}:
        return "chrM"

    parts = name.split("_")
    if parts and parts[0].isdigit():
        return f"chr{parts[0]}"

    return name


def get_output_prefix(bigwig_patterns):
    files = bigwig_patterns.split()
    if len(files) == 1:
        base = os.path.basename(files[0].replace("*", "Combined"))
        return sanitize_filename(os.path.splitext(base)[0])
    return f"combined_{len(files)}_bigwig_files"


def expand_bigwig_patterns(bigwig_patterns):
    bigwig_files = []

    for pattern in bigwig_patterns.split():
        matches = sorted(glob.glob(pattern))
        if matches:
            bigwig_files.extend(matches)

    if not bigwig_files:
        raise FileNotFoundError(f"No bigWig files match: {bigwig_patterns}")

    return bigwig_files


def read_alu_bed(
    alu_bed,
    extend=2000,
    strands="both",
    min_region_length=1,
    strand_col=6,
):
    """
    Read Alu BED/input and convert each Alu to a strand-aware region.

    + strand: Start -> Start + extend
    - strand: End - extend -> End

    strand_col is 1-based, matching normal command-line column numbering.
    Standard BED strand is column 6.
    """
    if strand_col < 1:
        raise ValueError("--strand_col must be a 1-based column number, e.g. 6 for standard BED.")

    strand_index = strand_col - 1
    required_columns = max(3, strand_col)
    records = []

    with open(alu_bed, encoding="utf-8-sig") as f:
        for line_number, line in enumerate(f, start=1):
            if not line.strip() or line.startswith("#"):
                continue

            fields = line.rstrip("\n").split()
            if len(fields) < required_columns:
                raise ValueError(
                    f"Expected at least {required_columns} columns on line {line_number} "
                    f"because --strand_col {strand_col} was requested: {line!r}"
                )

            chrom = standardize_chromosome_name(fields[0])
            alu_start = int(fields[1])
            alu_end = int(fields[2])
            strand = fields[strand_index].strip()

            if strand not in {"+", "-"}:
                continue
            if strands == "plus" and strand != "+":
                continue
            if strands == "minus" and strand != "-":
                continue

            if strand == "+":
                region_start = alu_start
                region_end = alu_start + extend
            else:
                region_start = alu_end - extend
                region_end = alu_end

            if region_start < 0:
                region_start = 0

            region_length = region_end - region_start
            if region_length < min_region_length:
                continue

            records.append((chrom, region_start, region_end, strand))

    if not records:
        raise ValueError("No Alu-derived regions remained after filtering.")

    return records


def extract_bigwig_region_to_bedgraph(
    bigwig_file,
    chrom,
    start,
    end,
    output_bedgraph,
):
    cmd = [
        "bigWigToBedGraph",
        f"-chrom={chrom}",
        f"-start={start}",
        f"-end={end}",
        bigwig_file,
        output_bedgraph,
    ]

    try:
        subprocess.run(
            cmd,
            check=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
        )
        return True

    except FileNotFoundError:
        raise RuntimeError(
            "Could not find bigWigToBedGraph in PATH. "
            "Install UCSC tools or add them to PATH."
        )

    except subprocess.CalledProcessError as e:
        print(f"Warning: bigWigToBedGraph failed for {chrom}:{start}-{end}. Skipping.")
        if e.stderr:
            print(e.stderr)
        return False


def expand_bedgraph_to_positions(temp_bedgraph, value_limit=None):
    """
    Convert extracted bedGraph intervals into sorted per-base positions/values.
    Zero-valued intervals are ignored because they do not contribute to DAC.
    """
    if not os.path.exists(temp_bedgraph) or os.path.getsize(temp_bedgraph) == 0:
        return None, None

    bg = pd.read_csv(
        temp_bedgraph,
        sep="\t",
        header=None,
        names=["Chromosome", "Start", "End", "Value"],
        dtype={
            "Chromosome": str,
            "Start": int,
            "End": int,
            "Value": float,
        },
    )

    if bg.empty:
        return None, None

    pos_chunks = []
    val_chunks = []

    for start, end, value in zip(bg["Start"], bg["End"], bg["Value"]):
        if end <= start:
            continue

        if value_limit is not None:
            value = float(np.clip(value, -value_limit, value_limit))

        if value == 0:
            continue

        pos_chunks.append(np.arange(start, end, dtype=np.int64))
        val_chunks.append(np.full(end - start, value, dtype=float))

    if not pos_chunks:
        return None, None

    positions = np.concatenate(pos_chunks)
    values = np.concatenate(val_chunks)

    order = np.argsort(positions)
    return positions[order], values[order]


def update_opportunities(opportunities, region_length, dmax, cache):
    """
    Opportunity count for a dense region of length L:
      lag d has L - d possible pairs.

    This script keeps the original Alu-DAC convention:
      output distances are 1 to dmax - 1.
    """
    if region_length not in cache:
        arr = np.zeros(dmax, dtype=float)
        for lag in range(1, min(dmax, region_length)):
            arr[lag] = region_length - lag
        cache[region_length] = arr

    opportunities += cache[region_length]


def update_dac_from_sparse_region(dac, positions, values, dmax):
    """
    Add DAC contributions from one sparse dyad-signal region.
    positions should be sorted genomic coordinates inside one region.
    """
    if positions is None or values is None:
        return 0

    n = len(positions)
    if n < 2:
        return 0

    pair_count = 0

    for i in range(n - 1):
        pos1 = positions[i]
        val1 = values[i]

        # Positions are sorted, so stop as soon as distance reaches dmax.
        for j in range(i + 1, n):
            dist = int(positions[j] - pos1)
            if dist >= dmax:
                break
            if dist > 0:
                dac[dist] += val1 * values[j]
                pair_count += 1

    return pair_count


def save_dac_to_tsv(
    dac,
    raw_dac,
    opportunities,
    output_file,
    normalize_opportunities=False,
    total_signal_in_regions=0.0,
    cpm_scale=1_000_000,
):
    dac_values = dac[1:]
    total_dac = np.sum(dac_values)

    if total_dac != 0:
        dac_percent = (dac_values / total_dac) * 100
    else:
        dac_percent = np.zeros_like(dac_values, dtype=float)

    if total_signal_in_regions != 0:
        dac_per_million_signal_pairs = (
            raw_dac[1:] / (total_signal_in_regions ** 2)
        ) * cpm_scale
    else:
        dac_per_million_signal_pairs = np.zeros_like(raw_dac[1:], dtype=float)

    df = pd.DataFrame({
        "Distance": np.arange(1, len(dac), dtype=int),
        "DAC Value": dac_values,
        "DAC Value Percent": dac_percent,
        "Raw DAC Value": raw_dac[1:],
        "DAC per million signal-pairs": dac_per_million_signal_pairs,
    })

    if normalize_opportunities:
        df["Opportunities"] = opportunities[1:]

    df.to_csv(output_file, sep="\t", index=False)


def calculate_dac_for_bigwig_file(
    bigwig_file,
    alu_regions,
    output_file,
    dmax=1000,
    value_limit=None,
    normalize_opportunities=False,
    cpm_scale=1_000_000,
):
    print(f"Reading dyads from bigWig: {bigwig_file}")

    dac = np.zeros(dmax, dtype=float)
    opportunities = np.zeros(dmax, dtype=float)
    opportunity_cache = {}

    used_regions = 0
    skipped_bigwig_error = 0
    total_signal_positions_in_regions = 0
    total_signal_in_regions = 0.0
    total_pairs_used = 0

    with tempfile.TemporaryDirectory() as tmpdir:
        for region_index, (chrom, region_start, region_end, strand) in enumerate(
            tqdm(
                alu_regions,
                desc=f"DAC {os.path.basename(bigwig_file)}",
            )
        ):
            region_length = region_end - region_start

            temp_bedgraph = os.path.join(
                tmpdir,
                f"alu_{chrom}_{region_start}_{region_end}_{region_index}.bedgraph",
            )

            ok = extract_bigwig_region_to_bedgraph(
                bigwig_file=bigwig_file,
                chrom=chrom,
                start=region_start,
                end=region_end,
                output_bedgraph=temp_bedgraph,
            )

            if not ok:
                skipped_bigwig_error += 1
                continue

            positions, values = expand_bedgraph_to_positions(
                temp_bedgraph,
                value_limit=value_limit,
            )

            if os.path.exists(temp_bedgraph):
                os.remove(temp_bedgraph)

            used_regions += 1

            if normalize_opportunities:
                update_opportunities(
                    opportunities,
                    region_length,
                    dmax,
                    opportunity_cache,
                )

            if positions is None or values is None or len(positions) < 2:
                if positions is not None and values is not None:
                    total_signal_positions_in_regions += len(positions)
                    total_signal_in_regions += float(np.sum(values))
                continue

            total_signal_positions_in_regions += len(positions)
            total_signal_in_regions += float(np.sum(values))

            total_pairs_used += update_dac_from_sparse_region(
                dac,
                positions,
                values,
                dmax,
            )

    raw_dac = dac.copy()

    if normalize_opportunities:
        non_zero = opportunities != 0
        dac[non_zero] = dac[non_zero] / opportunities[non_zero]

    save_dac_to_tsv(
        dac=dac,
        raw_dac=raw_dac,
        opportunities=opportunities,
        output_file=output_file,
        normalize_opportunities=normalize_opportunities,
        total_signal_in_regions=total_signal_in_regions,
        cpm_scale=cpm_scale,
    )

    print(f"Wrote: {output_file}")
    print(f"  Alu regions used: {used_regions}")
    print(f"  Alu regions skipped after bigWig extraction errors: {skipped_bigwig_error}")
    print(f"  Non-zero bigWig signal positions inside selected regions, counted with overlap: {total_signal_positions_in_regions}")
    print(f"  Total bigWig signal inside selected regions, counted with overlap: {total_signal_in_regions:.10g}")
    print(f"  Signal pairs contributing to DAC: {total_pairs_used}")


def main():
    parser = argparse.ArgumentParser(
        description=(
            "Calculate DAC from dyad-occurrence bigWig files over strand-aware "
            "regions derived from an Alu BED/input file."
        )
    )

    parser.add_argument(
        "--alu_bed",
        required=True,
        help=(
            "Alu BED/input file. Columns 1-3 must be chrom start end. "
            "The strand column is set with --strand_col; default is 6 for standard BED."
        ),
    )

    parser.add_argument(
        "--bigwig",
        required=True,
        help=(
            "Space-separated dyad bigWig file paths or glob patterns, "
            "e.g. 'BH01_chr*.bw' or 'sample_dyad.bw'."
        ),
    )

    parser.add_argument(
        "--out_prefix",
        default=None,
        help="Output prefix. Default is based on the bigWig filename/pattern.",
    )

    parser.add_argument(
        "--extend",
        type=int,
        default=2000,
        help="Length of strand-aware region. +: start to start+extend; -: end-extend to end. Default: 2000.",
    )

    parser.add_argument(
        "--dmax",
        type=int,
        default=1000,
        help="Maximum DAC distance. Distances 1 to dmax-1 are output. Default: 1000.",
    )

    parser.add_argument(
        "--strands",
        choices=["plus", "minus", "both"],
        default="both",
        help="Which Alu strands to include. Default: both.",
    )

    parser.add_argument(
        "--strand_col",
        type=int,
        default=6,
        help=(
            "1-based column number containing strand information (+ or -) in --alu_bed. "
            "Default: 6, the standard BED strand column."
        ),
    )

    parser.add_argument(
        "--value_limit",
        type=float,
        default=None,
        help="Optional absolute cap for bigWig signal values before DAC calculation.",
    )

    parser.add_argument(
        "--min_region_length",
        type=int,
        default=1,
        help="Minimum Alu-derived region length to include after clipping start at zero. Default: 1.",
    )

    parser.add_argument(
        "--normalize_opportunities",
        action="store_true",
        help="Divide DAC at each distance by the number of possible base-pair opportunities in the selected regions.",
    )

    parser.add_argument(
        "--cpm_scale",
        type=float,
        default=1_000_000,
        help="Scale for the depth-normalized pairwise DAC column. Default: 1,000,000.",
    )

    args = parser.parse_args()

    print("Reading Alu BED and creating strand-aware DAC regions")
    alu_regions = read_alu_bed(
        args.alu_bed,
        extend=args.extend,
        strands=args.strands,
        min_region_length=args.min_region_length,
        strand_col=args.strand_col,
    )
    print(f"Alu-derived regions: {len(alu_regions)}")

    bigwig_files = expand_bigwig_patterns(args.bigwig)

    base_prefix = args.out_prefix if args.out_prefix else get_output_prefix(args.bigwig)
    suffix = "opportunity_normalized" if args.normalize_opportunities else "raw"

    for bigwig_file in bigwig_files:
        if len(bigwig_files) == 1:
            output_file = f"{base_prefix}_Alu_extend{args.extend}_{args.strands}_DAC_{suffix}.tsv"
        else:
            bigwig_base = sanitize_filename(os.path.splitext(os.path.basename(bigwig_file))[0])
            output_file = f"{base_prefix}_{bigwig_base}_Alu_extend{args.extend}_{args.strands}_DAC_{suffix}.tsv"

        calculate_dac_for_bigwig_file(
            bigwig_file=bigwig_file,
            alu_regions=alu_regions,
            output_file=output_file,
            dmax=args.dmax,
            value_limit=args.value_limit,
            normalize_opportunities=args.normalize_opportunities,
            cpm_scale=args.cpm_scale,
        )


if __name__ == "__main__":
    main()
