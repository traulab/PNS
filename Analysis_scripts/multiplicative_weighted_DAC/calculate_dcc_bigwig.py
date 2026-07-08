#!/usr/bin/env python3
"""
Streaming DCC from two bigWig dyad-occurrence signals.

DCC = Distance Cross-Correlation.

For two dyad occurrence tracks A and B, this script internally calculates
signed-lag cross-correlation:

    DCC[lag] = sum_x A[x] * B[x + lag]

where:

    lag = position_B - position_A

A positive lag means that signal in B is downstream/right of signal in A
in the coordinate system being analysed. If a BED file has a minus-strand
column and --strand_column is used, intervals on '-' are reversed first, so
positive lag means downstream in feature-oriented coordinates.

By default, the script collapses signed lags into absolute distances:

    Distance 0 = lag 0
    Distance d = lag -d + lag +d

This means the default output is directionless:

    Distance    DCC Value
    0
    1
    2
    ...
    dmax

To keep signed directional lags instead, use:

    --signed_lags

Typical default directionless use:

    python calculate_dcc_bigwig.py \
      --bigwig_a 'BH01_chr*_mode147_dyad.bw' \
      --bigwig_b 'BH01_chr*_mode167_dyad.bw' \
      --label_a 147bp \
      --label_b 167bp \
      --chrom_sizes /mnt/d/Snyder_bams/male/chrom.sizes \
      --scope genome \
      --dmax 50

This will report whether the strongest offset is 0 bp, 5 bp, 10 bp, etc.,
but not whether B is left or right of A.

Directional signed-lag use:

    python calculate_dcc_bigwig.py \
      --bigwig_a 'BH01_chr*_mode147_dyad.bw' \
      --bigwig_b 'BH01_chr*_mode167_dyad.bw' \
      --label_a 147bp \
      --label_b 167bp \
      --chrom_sizes /mnt/d/Snyder_bams/male/chrom.sizes \
      --scope genome \
      --dmax 50 \
      --signed_lags

If the strongest signed peak is at +5, then B dyads are typically +5 bp
relative to A dyads. If the strongest signed peak is at -5, then B dyads
are typically -5 bp relative to A dyads.
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


def compact_pattern_label(patterns, default_label):
    parts = str(patterns).split()
    if len(parts) == 1:
        base = os.path.basename(parts[0].replace("*", "Combined"))
        base = os.path.splitext(base)[0]
        return sanitize_filename(base)
    return default_label


def get_output_prefix(pattern_a, pattern_b, label_a=None, label_b=None):
    label_a = label_a or compact_pattern_label(pattern_a, "A")
    label_b = label_b or compact_pattern_label(pattern_b, "B")
    return f"{sanitize_filename(label_a)}_vs_{sanitize_filename(label_b)}"


def expand_bigwig_patterns(bigwig_patterns):
    bigwig_files = []

    for pattern in str(bigwig_patterns).split():
        matches = sorted(glob.glob(pattern))
        if matches:
            bigwig_files.extend(matches)
        elif os.path.exists(pattern):
            bigwig_files.append(pattern)

    # Preserve order but remove duplicates.
    seen = set()
    unique_files = []
    for path in bigwig_files:
        if path not in seen:
            unique_files.append(path)
            seen.add(path)

    if not unique_files:
        raise FileNotFoundError(f"No bigWig files match: {bigwig_patterns}")

    return unique_files


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

    return group_positions_values(positions, values)


def group_positions_values(positions, values):
    if positions is None or values is None or len(positions) == 0:
        return None, None

    grouped = (
        pd.DataFrame({"Position": positions, "Value": values})
        .groupby("Position", observed=False)["Value"]
        .sum()
        .reset_index()
        .sort_values("Position")
    )

    return (
        grouped["Position"].to_numpy(dtype=np.int64),
        grouped["Value"].to_numpy(dtype=float),
    )


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


def extract_positions_values_from_bigwigs(
    bigwig_files,
    chrom,
    start,
    end,
    tmpdir,
    temp_prefix,
    value_limit=None,
    verbose=False,
):
    all_positions = []
    all_values = []

    for file_index, bigwig_file in enumerate(bigwig_files):
        temp_bedgraph = os.path.join(
            tmpdir,
            f"{sanitize_filename(temp_prefix)}_{file_index}.bedgraph",
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
            # This commonly happens when chromosome-specific bigWigs are used and
            # the current bigWig does not contain the requested chromosome.
            if verbose:
                print(
                    f"Warning: bigWigToBedGraph failed for "
                    f"{bigwig_file} {chrom}:{start}-{end}. Skipping."
                )
                print(e.stderr)
            continue

        positions, values = expand_bedgraph_to_positions(
            temp_bedgraph,
            value_limit=value_limit,
        )

        if positions is not None and values is not None:
            all_positions.append(positions)
            all_values.append(values)

        if os.path.exists(temp_bedgraph):
            os.remove(temp_bedgraph)

    if not all_positions:
        return None, None

    positions = np.concatenate(all_positions)
    values = np.concatenate(all_values)

    return group_positions_values(positions, values)


def update_cross_opportunities(opportunities, region_length, dmax, opportunities_cache):
    if region_length not in opportunities_cache:
        lags = np.arange(-dmax, dmax + 1, dtype=np.int64)
        opp = region_length - np.abs(lags)
        opp[opp < 0] = 0
        opportunities_cache[region_length] = opp.astype(float)

    opportunities += opportunities_cache[region_length]


def collapse_signed_to_absolute(signed_values, dmax):
    """
    Collapse signed lags into absolute distances.

    signed_values is indexed as:
        -dmax ... 0 ... +dmax

    output is indexed as:
        0 ... dmax

    Distance 0 = lag 0
    Distance d = lag -d + lag +d
    """
    absolute_values = np.zeros(dmax + 1, dtype=float)

    absolute_values[0] = signed_values[dmax]

    for distance in range(1, dmax + 1):
        absolute_values[distance] = (
            signed_values[dmax - distance] +
            signed_values[dmax + distance]
        )

    return absolute_values


def update_dcc_from_region(dcc, positions_a, values_a, positions_b, values_b, dmax):
    if positions_a is None or values_a is None:
        return
    if positions_b is None or values_b is None:
        return
    if len(positions_a) == 0 or len(positions_b) == 0:
        return

    # DCC[lag] = sum_x A[x] * B[x + lag]
    # lag = position_B - position_A.
    for pos_a, val_a in zip(positions_a, values_a):
        left = np.searchsorted(positions_b, pos_a - dmax, side="left")
        right = np.searchsorted(positions_b, pos_a + dmax, side="right")

        if right <= left:
            continue

        lags = positions_b[left:right] - pos_a
        products = val_a * values_b[left:right]
        np.add.at(dcc, lags + dmax, products)


def save_dcc_to_tsv(dcc, output_file, dmax, signed_lags=False):
    dcc_values = dcc.astype(float)
    total_dcc = np.sum(dcc_values)

    if total_dcc != 0:
        dcc_percent = (dcc_values / total_dcc) * 100
    else:
        dcc_percent = np.zeros_like(dcc_values, dtype=float)

    if signed_lags:
        x_values = np.arange(-dmax, dmax + 1, dtype=int)
        x_column = "Lag"
    else:
        x_values = np.arange(0, dmax + 1, dtype=int)
        x_column = "Distance"

    df = pd.DataFrame({
        x_column: x_values,
        "DCC Value": dcc_values,
        "DCC Value Percent": dcc_percent,
    })

    df.to_csv(output_file, sep="\t", index=False)
    print(f"DCC values saved to {output_file}")


def save_shift_summary(
    dcc,
    output_file,
    dmax,
    summary_lag_window=50,
    signed_lags=False,
):
    values = dcc.astype(float)
    rows = []

    if signed_lags:
        lags = np.arange(-dmax, dmax + 1, dtype=int)

        if len(values) > 0:
            max_idx = int(np.argmax(values))
            rows.append({
                "Metric": "max_lag_all",
                "Lag_or_Distance": int(lags[max_idx]),
                "DCC Value": float(values[max_idx]),
            })

        in_window = np.abs(lags) <= summary_lag_window
        if np.any(in_window):
            window_lags = lags[in_window]
            window_values = values[in_window]
            max_idx = int(np.argmax(window_values))
            rows.append({
                "Metric": f"max_lag_within_{summary_lag_window}bp",
                "Lag_or_Distance": int(window_lags[max_idx]),
                "DCC Value": float(window_values[max_idx]),
            })

        for lag in [-10, -5, 0, 5, 10]:
            if -dmax <= lag <= dmax:
                idx = lag + dmax
                rows.append({
                    "Metric": f"lag_{lag}",
                    "Lag_or_Distance": lag,
                    "DCC Value": float(values[idx]),
                })

        if -dmax <= 5 <= dmax and -dmax <= 0 <= dmax:
            rows.append({
                "Metric": "lag_+5_minus_lag_0",
                "Lag_or_Distance": 5,
                "DCC Value": float(values[dmax + 5] - values[dmax]),
            })

        if -dmax <= -5 <= dmax and -dmax <= 0 <= dmax:
            rows.append({
                "Metric": "lag_-5_minus_lag_0",
                "Lag_or_Distance": -5,
                "DCC Value": float(values[dmax - 5] - values[dmax]),
            })

    else:
        distances = np.arange(0, dmax + 1, dtype=int)

        if len(values) > 0:
            max_idx = int(np.argmax(values))
            rows.append({
                "Metric": "max_distance_all",
                "Lag_or_Distance": int(distances[max_idx]),
                "DCC Value": float(values[max_idx]),
            })

        in_window = distances <= summary_lag_window
        if np.any(in_window):
            window_distances = distances[in_window]
            window_values = values[in_window]
            max_idx = int(np.argmax(window_values))
            rows.append({
                "Metric": f"max_distance_within_{summary_lag_window}bp",
                "Lag_or_Distance": int(window_distances[max_idx]),
                "DCC Value": float(window_values[max_idx]),
            })

        for distance in [0, 5, 10, 20]:
            if distance <= dmax:
                rows.append({
                    "Metric": f"distance_{distance}",
                    "Lag_or_Distance": distance,
                    "DCC Value": float(values[distance]),
                })

        if 5 <= dmax:
            rows.append({
                "Metric": "distance_5_minus_distance_0",
                "Lag_or_Distance": 5,
                "DCC Value": float(values[5] - values[0]),
            })

        if 10 <= dmax:
            rows.append({
                "Metric": "distance_10_minus_distance_0",
                "Lag_or_Distance": 10,
                "DCC Value": float(values[10] - values[0]),
            })

    pd.DataFrame(rows).to_csv(output_file, sep="\t", index=False)
    print(f"DCC shift summary saved to {output_file}")


def calculate_dcc_streaming_from_bigwig(
    bigwig_patterns_a,
    bigwig_patterns_b,
    chromatin_df,
    output_prefix,
    dmax=100,
    value_limit=None,
    min_region_length=None,
    normalize_dcc=True,
    normalize_by_signal_totals=False,
    summary_lag_window=50,
    signed_lags=False,
    verbose=False,
):
    print("Running calculate_dcc_streaming_from_bigwig")

    bigwig_files_a = expand_bigwig_patterns(bigwig_patterns_a)
    bigwig_files_b = expand_bigwig_patterns(bigwig_patterns_b)

    print(f"A bigWigs: {len(bigwig_files_a)}")
    print(f"B bigWigs: {len(bigwig_files_b)}")

    if min_region_length is None:
        min_region_length = dmax + 1

    states = sorted(chromatin_df["State"].unique())

    opportunities_cache = {}
    completed_outputs = {}

    with tempfile.TemporaryDirectory() as tmpdir:
        for state in tqdm(states, desc="Processing states"):
            state_df = chromatin_df[chromatin_df["State"] == state]

            # Internally always calculate signed DCC.
            # Output is later either kept signed or collapsed to absolute distances.
            dcc = np.zeros(2 * dmax + 1, dtype=float)
            opportunities = np.zeros(2 * dmax + 1, dtype=float)
            total_signal_a = 0.0
            total_signal_b = 0.0
            used_regions = 0

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

                temp_prefix_a = f"A_{sanitize_filename(state)}_{chrom}_{start}_{end}_{region_index}"
                temp_prefix_b = f"B_{sanitize_filename(state)}_{chrom}_{start}_{end}_{region_index}"

                positions_a, values_a = extract_positions_values_from_bigwigs(
                    bigwig_files=bigwig_files_a,
                    chrom=chrom,
                    start=start,
                    end=end,
                    tmpdir=tmpdir,
                    temp_prefix=temp_prefix_a,
                    value_limit=value_limit,
                    verbose=verbose,
                )

                positions_b, values_b = extract_positions_values_from_bigwigs(
                    bigwig_files=bigwig_files_b,
                    chrom=chrom,
                    start=start,
                    end=end,
                    tmpdir=tmpdir,
                    temp_prefix=temp_prefix_b,
                    value_limit=value_limit,
                    verbose=verbose,
                )

                positions_a, values_a = reverse_interval_signal_if_needed(
                    positions_a,
                    values_a,
                    interval_start=start,
                    interval_end=end,
                    strand=strand,
                )

                positions_b, values_b = reverse_interval_signal_if_needed(
                    positions_b,
                    values_b,
                    interval_start=start,
                    interval_end=end,
                    strand=strand,
                )

                if positions_a is None or positions_b is None:
                    continue

                if normalize_dcc:
                    update_cross_opportunities(
                        opportunities,
                        region_length,
                        dmax,
                        opportunities_cache,
                    )

                total_signal_a += float(np.sum(values_a))
                total_signal_b += float(np.sum(values_b))
                used_regions += 1

                update_dcc_from_region(
                    dcc,
                    positions_a,
                    values_a,
                    positions_b,
                    values_b,
                    dmax,
                )

            if signed_lags:
                dcc_out = dcc.copy()
                opportunities_out = opportunities.copy()

                if normalize_dcc:
                    non_zero = opportunities_out != 0
                    dcc_out[non_zero] = dcc_out[non_zero] / opportunities_out[non_zero]

            else:
                dcc_out = collapse_signed_to_absolute(dcc, dmax)
                opportunities_out = collapse_signed_to_absolute(opportunities, dmax)

                if normalize_dcc:
                    non_zero = opportunities_out != 0
                    dcc_out[non_zero] = dcc_out[non_zero] / opportunities_out[non_zero]

            if normalize_by_signal_totals:
                denom = total_signal_a * total_signal_b
                if denom != 0:
                    dcc_out = dcc_out / denom

            sanitized_state = sanitize_filename(state)

            suffix_parts = []
            suffix_parts.append("signed_lags" if signed_lags else "absolute_distances")
            suffix_parts.append("opportunity_normalized" if normalize_dcc else "raw")
            if normalize_by_signal_totals:
                suffix_parts.append("signal_total_normalized")
            suffix = "_".join(suffix_parts)

            output_file = (
                f"{output_prefix}_{sanitized_state}_"
                f"bigwig_streaming_DCC_values_{suffix}.tsv"
            )

            summary_file = (
                f"{output_prefix}_{sanitized_state}_"
                f"bigwig_streaming_DCC_shift_summary_{suffix}.tsv"
            )

            save_dcc_to_tsv(
                dcc_out,
                output_file,
                dmax,
                signed_lags=signed_lags,
            )

            save_shift_summary(
                dcc_out,
                summary_file,
                dmax,
                summary_lag_window=summary_lag_window,
                signed_lags=signed_lags,
            )

            completed_outputs[state] = {
                "dcc_values": output_file,
                "summary": summary_file,
                "used_regions": used_regions,
                "total_signal_a": total_signal_a,
                "total_signal_b": total_signal_b,
                "signed_lags": signed_lags,
            }

            print(
                f"Completed DCC for state: {state} "
                f"using {used_regions} regions; "
                f"sum(A)={total_signal_a:.6g}; sum(B)={total_signal_b:.6g}; "
                f"output_mode={'signed_lags' if signed_lags else 'absolute_distances'}"
            )

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

    output_prefix = args.output_prefix
    if output_prefix is None:
        output_prefix = get_output_prefix(
            args.bigwig_a,
            args.bigwig_b,
            args.label_a,
            args.label_b,
        )

    if args.scope == "chromosome" and args.chromosome is not None:
        output_prefix = f"{output_prefix}_{standardize_chromosome_name(args.chromosome)}"

    calculate_dcc_streaming_from_bigwig(
        bigwig_patterns_a=args.bigwig_a,
        bigwig_patterns_b=args.bigwig_b,
        chromatin_df=chromatin_df,
        output_prefix=output_prefix,
        dmax=args.dmax,
        value_limit=args.value_limit,
        min_region_length=args.min_region_length,
        normalize_dcc=not args.no_normalize_dcc,
        normalize_by_signal_totals=args.normalize_by_signal_totals,
        summary_lag_window=args.summary_lag_window,
        signed_lags=args.signed_lags,
        verbose=args.verbose,
    )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description=(
            "Calculate streaming Distance Cross-Correlation from two bigWig "
            "dyad-occurrence signals. Internally calculates signed "
            "DCC[lag] = sum_x A[x] * B[x + lag]. By default, output is "
            "collapsed to absolute distances. Use --signed_lags to keep "
            "directional signed lags."
        )
    )

    parser.add_argument(
        "--bigwig_a",
        required=True,
        help=(
            "Space-separated list of bigWig file patterns for signal A, "
            "e.g. 'CH01_chr*_147_dyad.bw' or 'sample_147.bw'."
        ),
    )

    parser.add_argument(
        "--bigwig_b",
        required=True,
        help=(
            "Space-separated list of bigWig file patterns for signal B, "
            "e.g. 'CH01_chr*_167_dyad.bw' or 'sample_167.bw'."
        ),
    )

    parser.add_argument(
        "--label_a",
        default=None,
        help="Optional short label for signal A, e.g. 147bp.",
    )

    parser.add_argument(
        "--label_b",
        default=None,
        help="Optional short label for signal B, e.g. 167bp.",
    )

    parser.add_argument(
        "--output_prefix",
        default=None,
        help="Optional output prefix. Default is inferred from labels/patterns.",
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
        default=100,
        help=(
            "Maximum signed lag/distance for DCC calculation. Inclusive. "
            "Default: 100. Use 50 for a focused 5-bp shift test, or 1500 "
            "for long-range phasing."
        ),
    )

    parser.add_argument(
        "--summary_lag_window",
        type=int,
        default=50,
        help=(
            "Lag/distance window used for max-lag or max-distance summary. "
            "Default: 50."
        ),
    )

    parser.add_argument(
        "--value_limit",
        type=float,
        default=None,
        help="Optional absolute cap for bigWig signal values before DCC.",
    )

    parser.add_argument(
        "--min_region_length",
        type=int,
        default=None,
        help=(
            "Minimum interval/window length to include. "
            "Default: dmax + 1."
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
        "--no_normalize_dcc",
        action="store_true",
        help=(
            "Turn OFF opportunity-based DCC normalization. "
            "Default is ON."
        ),
    )

    parser.add_argument(
        "--normalize_by_signal_totals",
        action="store_true",
        help=(
            "Also divide final DCC by total signal(A) * total signal(B) "
            "within the analysed regions. Default is OFF."
        ),
    )

    parser.add_argument(
        "--signed_lags",
        action="store_true",
        help=(
            "Keep signed DCC lags from -dmax to +dmax. "
            "By default, signed lags are collapsed to absolute distances, "
            "where distance d = lag -d + lag +d."
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

    parser.add_argument(
        "--verbose",
        action="store_true",
        help="Print warnings from failed per-chromosome bigWig extractions.",
    )

    args = parser.parse_args()

    main(args)