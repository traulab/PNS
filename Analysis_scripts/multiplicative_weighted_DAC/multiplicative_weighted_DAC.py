#!/usr/bin/env python3
"""
Streaming DAC from bigWig signal.

Modes:
1. BED regions:
   --chromatin_bed regions.bed

2. One chromosome only:
   --scope chromosome --chromosome chr1 --chrom_sizes genome.chrom.sizes

3. Whole genome:
   --scope genome --chrom_sizes genome.chrom.sizes

For chromosome/genome mode, regions are generated automatically and labelled "Genome"
unless --state_name is changed.
"""

import numpy as np
import pandas as pd
import argparse
import glob
import os
import subprocess
import tempfile
from tqdm import tqdm
import matplotlib
matplotlib.use("Agg")


def sanitize_filename(name):
    return str(name).replace("/", "_").replace("\\", "_").replace(" ", "_")


def standardize_chromosome_name(name):
    name = str(name)

    if name.startswith("chr"):
        return name
    if name.isdigit():
        return f"chr{name}"
    if name in {"X", "Y"}:
        return f"chr{name}"
    if name in {"M", "MT"}:
        return "chrM"

    parts = name.split("_")
    if parts[0].isdigit():
        return f"chr{parts[0]}"

    return name


def get_output_prefix(input_patterns):
    files = input_patterns.split()

    if len(files) == 1:
        base = os.path.basename(files[0].replace("*", "Combined"))
        return sanitize_filename(os.path.splitext(base)[0])

    return f"combined_{len(files)}"


def read_chrom_sizes(chrom_sizes_file):
    chrom_sizes = pd.read_csv(
        chrom_sizes_file,
        sep="\t",
        header=None,
        names=["Chromosome", "Size"],
        usecols=[0, 1],
        dtype={"Chromosome": str, "Size": int},
    )

    chrom_sizes["Chromosome"] = chrom_sizes["Chromosome"].apply(
        standardize_chromosome_name
    )

    return chrom_sizes


def make_windows_from_chrom_sizes(
    chrom_sizes_file,
    scope="genome",
    chromosome=None,
    window_size=100000,
    state_name="Genome",
):
    chrom_sizes = read_chrom_sizes(chrom_sizes_file)

    if scope == "chromosome":
        if chromosome is None:
            raise ValueError("--chromosome is required when --scope chromosome")

        chromosome = standardize_chromosome_name(chromosome)
        chrom_sizes = chrom_sizes[chrom_sizes["Chromosome"] == chromosome]

        if chrom_sizes.empty:
            raise ValueError(
                f"Chromosome {chromosome} was not found in {chrom_sizes_file}"
            )

    rows = []

    for _, row in chrom_sizes.iterrows():
        chrom = row["Chromosome"]
        size = int(row["Size"])

        for start in range(0, size, window_size):
            end = min(start + window_size, size)

            rows.append({
                "Chromosome": chrom,
                "Start": start,
                "End": end,
                "State": state_name,
                "Strand": "+",
            })

    return pd.DataFrame(rows)


def read_chromatin_states(
    chromatin_bed,
    convert_to_euchromatin=False,
    strand_column=None,
):
    print("Running read_chromatin_states")

    if strand_column is None:
        df = pd.read_csv(
            chromatin_bed,
            sep="\t",
            header=None,
            usecols=[0, 1, 2, 3],
            names=["Chromosome", "Start", "End", "State"],
            dtype={
                "Chromosome": str,
                "Start": int,
                "End": int,
                "State": str,
            },
        )

        df["Strand"] = "+"

    else:
        strand_col_index = strand_column - 1

        if strand_col_index < 0:
            raise ValueError(
                "--strand_column must be 1-based, e.g. 6 for BED column 6."
            )

        required_cols = [0, 1, 2, 3, strand_col_index]
        usecols = sorted(set(required_cols))

        raw = pd.read_csv(
            chromatin_bed,
            sep="\t",
            header=None,
            usecols=usecols,
            dtype=str,
        )

        col_map = {col: i for i, col in enumerate(usecols)}

        df = pd.DataFrame({
            "Chromosome": raw.iloc[:, col_map[0]].astype(str),
            "Start": raw.iloc[:, col_map[1]].astype(int),
            "End": raw.iloc[:, col_map[2]].astype(int),
            "State": raw.iloc[:, col_map[3]].astype(str),
            "Strand": raw.iloc[:, col_map[strand_col_index]].astype(str),
        })

        df["Strand"] = df["Strand"].fillna("+")
        df.loc[~df["Strand"].isin(["+", "-"]), "Strand"] = "+"

    df["Chromosome"] = df["Chromosome"].apply(standardize_chromosome_name)

    if convert_to_euchromatin:
        euchromatin_prefixes = {
            "1_", "2_", "4_", "5_", "6_", "7_", "9_", "10_", "11_"
        }

        df["State"] = df["State"].apply(
            lambda x: "Euchromatin"
            if any(str(x).startswith(prefix) for prefix in euchromatin_prefixes)
            else x
        )

    df = df.sort_values(
        ["State", "Chromosome", "Start", "End"]
    ).reset_index(drop=True)

    return df


def expand_bedgraph_to_positions(temp_bedgraph, value_limit=None):
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
            value = np.clip(value, -value_limit, value_limit)

        pos_chunks.append(np.arange(start, end, dtype=np.int64))
        val_chunks.append(np.full(end - start, value, dtype=float))

    if not pos_chunks:
        return None, None

    positions = np.concatenate(pos_chunks)
    values = np.concatenate(val_chunks)

    order = np.argsort(positions)

    return positions[order], values[order]


def reverse_interval_signal_if_needed(
    positions,
    values,
    interval_start,
    interval_end,
    strand,
):
    if positions is None or values is None:
        return positions, values

    if strand != "-":
        return positions, values

    reversed_positions = interval_start + interval_end - 1 - positions

    order = np.argsort(reversed_positions)

    return reversed_positions[order], values[order]


def update_opportunities(opportunities, region_length, dmax, opportunities_cache):
    if region_length not in opportunities_cache:
        opportunities_cache[region_length] = np.zeros(dmax + 1, dtype=float)

        max_lag = min(dmax, region_length - 1)

        for lag in range(1, max_lag + 1):
            opportunities_cache[region_length][lag] = region_length - lag

    opportunities += opportunities_cache[region_length]


def update_dac_from_region(dac, positions, values, dmax):
    if positions is None or values is None:
        return

    if len(positions) < 2:
        return

    grouped = (
        pd.DataFrame({"Position": positions, "Value": values})
        .groupby("Position", observed=False)["Value"]
        .sum()
        .reset_index()
        .sort_values("Position")
    )

    positions = grouped["Position"].to_numpy(dtype=np.int64)
    values = grouped["Value"].to_numpy(dtype=float)

    n = len(positions)

    for i in range(n):
        pos1 = positions[i]
        val1 = values[i]

        for j in range(i + 1, n):
            dist = positions[j] - pos1

            if dist > dmax:
                break

            dac[dist] += val1 * values[j]


def save_dac_to_tsv(dac, output_file):
    dac_values = dac[1:]
    total_dac = np.sum(dac_values)

    if total_dac != 0:
        dac_percent = (dac_values / total_dac) * 100
    else:
        dac_percent = np.zeros_like(dac_values, dtype=float)

    df = pd.DataFrame({
        "Distance": range(1, len(dac)),
        "DAC Value": dac_values,
        "DAC Value Percent": dac_percent,
    })

    df.to_csv(output_file, sep="\t", index=False)
    print(f"DAC values saved to {output_file}")


def calculate_dac_streaming_from_bigwig(
    bigwig_patterns,
    chromatin_df,
    output_prefix,
    dmax=1500,
    value_limit=None,
    min_region_length=2000,
    normalize_dac=True,
):
    print("Running calculate_dac_streaming_from_bigwig")

    bigwig_files = []

    for pattern in bigwig_patterns.split():
        bigwig_files.extend(glob.glob(pattern))

    if not bigwig_files:
        raise FileNotFoundError(f"No bigWig files match: {bigwig_patterns}")

    states = sorted(chromatin_df["State"].unique())

    opportunities_cache = {}
    completed_outputs = {}

    with tempfile.TemporaryDirectory() as tmpdir:
        for state in tqdm(states, desc="Processing states"):
            state_df = chromatin_df[chromatin_df["State"] == state]

            dac = np.zeros(dmax + 1, dtype=float)
            opportunities = np.zeros(dmax + 1, dtype=float)

            for bigwig_file in bigwig_files:
                print(f"Processing state {state}, bigWig {bigwig_file}")

                for region_index, row in tqdm(
                    state_df.iterrows(),
                    total=len(state_df),
                    desc=f"Extracting {state}",
                    leave=False,
                ):
                    chrom = str(row["Chromosome"])
                    start = int(row["Start"])
                    end = int(row["End"])
                    strand = str(row.get("Strand", "+"))

                    region_length = end - start

                    if region_length < min_region_length:
                        continue

                    temp_bedgraph = os.path.join(
                        tmpdir,
                        (
                            f"{sanitize_filename(state)}_"
                            f"{chrom}_{start}_{end}_{region_index}.bedgraph"
                        ),
                    )

                    cmd = [
                        "bigWigToBedGraph",
                        f"-chrom={chrom}",
                        f"-start={start}",
                        f"-end={end}",
                        bigwig_file,
                        temp_bedgraph,
                    ]

                    try:
                        subprocess.run(
                            cmd,
                            check=True,
                            stdout=subprocess.PIPE,
                            stderr=subprocess.PIPE,
                            text=True,
                        )

                    except FileNotFoundError:
                        raise RuntimeError(
                            "Could not find bigWigToBedGraph in PATH. "
                            "Install UCSC tools or add them to PATH."
                        )

                    except subprocess.CalledProcessError as e:
                        print(
                            f"Warning: bigWigToBedGraph failed for "
                            f"{chrom}:{start}-{end}. Skipping."
                        )
                        print(e.stderr)
                        continue

                    positions, values = expand_bedgraph_to_positions(
                        temp_bedgraph,
                        value_limit=value_limit,
                    )

                    positions, values = reverse_interval_signal_if_needed(
                        positions,
                        values,
                        interval_start=start,
                        interval_end=end,
                        strand=strand,
                    )

                    if normalize_dac:
                        update_opportunities(
                            opportunities,
                            region_length,
                            dmax,
                            opportunities_cache,
                        )

                    update_dac_from_region(
                        dac,
                        positions,
                        values,
                        dmax,
                    )

                    if os.path.exists(temp_bedgraph):
                        os.remove(temp_bedgraph)

            if normalize_dac:
                non_zero = opportunities != 0
                dac[non_zero] = dac[non_zero] / opportunities[non_zero]

            sanitized_state = sanitize_filename(state)
            suffix = "normalized" if normalize_dac else "raw"

            output_file = (
                f"{output_prefix}_{sanitized_state}_"
                f"bigwig_streaming_DAC_values_{suffix}.tsv"
            )

            save_dac_to_tsv(dac, output_file)

            completed_outputs[state] = output_file
            print(f"Completed and saved DAC for state: {state}")

    return completed_outputs


def main(args):
    if args.chromatin_bed is not None:
        chromatin_df = read_chromatin_states(
            args.chromatin_bed,
            convert_to_euchromatin=args.convert_to_euchromatin,
            strand_column=args.strand_column,
        )

    else:
        if args.chrom_sizes is None:
            raise ValueError(
                "If --chromatin_bed is not supplied, you must provide "
                "--chrom_sizes."
            )

        chromatin_df = make_windows_from_chrom_sizes(
            chrom_sizes_file=args.chrom_sizes,
            scope=args.scope,
            chromosome=args.chromosome,
            window_size=args.window_size,
            state_name=args.state_name,
        )

    output_prefix = get_output_prefix(args.bigwig)

    if args.scope == "chromosome" and args.chromosome is not None:
        output_prefix = f"{output_prefix}_{standardize_chromosome_name(args.chromosome)}"

    calculate_dac_streaming_from_bigwig(
        bigwig_patterns=args.bigwig,
        chromatin_df=chromatin_df,
        output_prefix=output_prefix,
        dmax=args.dmax,
        value_limit=args.value_limit,
        min_region_length=args.min_region_length,
        normalize_dac=not args.no_normalize_dac,
    )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description=(
            "Calculate streaming Distance Autocorrelation from bigWig signal. "
            "Can use BED intervals, one chromosome, or whole-genome windows."
        )
    )

    parser.add_argument(
        "--bigwig",
        required=True,
        help=(
            "Space-separated list of bigWig file patterns, "
            "e.g. 'CH01_chr*.bw' or 'sample.bw'."
        ),
    )

    parser.add_argument(
        "--chromatin_bed",
        default=None,
        help=(
            "Optional BED file containing regions/states. "
            "First four columns: chrom, start, end, state. "
            "If omitted, regions are generated from --chrom_sizes."
        ),
    )

    parser.add_argument(
        "--scope",
        choices=["genome", "chromosome"],
        default="genome",
        help=(
            "Used only when --chromatin_bed is omitted. "
            "'genome' uses all chromosomes in --chrom_sizes. "
            "'chromosome' uses only --chromosome."
        ),
    )

    parser.add_argument(
        "--chromosome",
        default=None,
        help=(
            "Chromosome to analyse when --scope chromosome, "
            "e.g. chr1, chrX, I, II, III."
        ),
    )

    parser.add_argument(
        "--chrom_sizes",
        default=None,
        help=(
            "Chromosome sizes file with two columns: chromosome and length. "
            "Required if --chromatin_bed is omitted."
        ),
    )

    parser.add_argument(
        "--window_size",
        type=int,
        default=100000,
        help=(
            "Window size used when generating genome/chromosome intervals. "
            "Default: 100000."
        ),
    )

    parser.add_argument(
        "--state_name",
        default="Genome",
        help=(
            "State/group name assigned to automatically generated intervals. "
            "Default: Genome."
        ),
    )

    parser.add_argument(
        "--dmax",
        type=int,
        default=1500,
        help="Maximum distance for DAC calculation. Inclusive. Default: 1500.",
    )

    parser.add_argument(
        "--value_limit",
        type=float,
        default=None,
        help="Optional absolute cap for bigWig signal values before DAC.",
    )

    parser.add_argument(
        "--min_region_length",
        type=int,
        default=1501,
        help=(
            "Minimum interval/window length to include. "
            "Default: 1501."
        ),
    )

    parser.add_argument(
        "--convert_to_euchromatin",
        action="store_true",
        help=(
            "Only used with --chromatin_bed. "
            "Collapse selected ChromHMM states into Euchromatin."
        ),
    )

    parser.add_argument(
        "--no_normalize_dac",
        action="store_true",
        help=(
            "Turn OFF opportunity-based DAC normalization. "
            "Default is ON."
        ),
    )

    parser.add_argument(
        "--strand_column",
        type=int,
        default=None,
        help=(
            "Only used with --chromatin_bed. "
            "Optional 1-based column number containing strand information. "
            "Use 6 for standard BED strand column."
        ),
    )

    args = parser.parse_args()

    main(args)