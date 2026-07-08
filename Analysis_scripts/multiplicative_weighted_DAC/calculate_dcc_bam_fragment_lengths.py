#!/usr/bin/env python3
"""
Streaming DCC directly from BAM files, using different fragment lengths for
signal A and signal B.

DCC = Distance Cross-Correlation.

For two fragment-derived occurrence tracks A and B, this script internally
calculates signed-lag cross-correlation:

    DCC[lag] = sum_x A[x] * B[x + lag]

where:

    lag = position_B - position_A

A positive lag means that B positions are downstream/right of A positions in the
genomic coordinate system being analysed. If a BED file has a minus-strand column and
--strand_column is used, intervals on '-' are reversed first, so positive lag
means downstream in feature-oriented coordinates.

By default, signed lags are collapsed into absolute distances:

    Distance 0 = lag 0
    Distance d = lag -d + lag +d

This means the default output is directionless:

    Distance    DCC Value
    0
    1
    2
    ...
    dmax

To keep directional lags instead, use:

    --signed_lags

Typical use, comparing 147-bp fragments in one BAM set against 167-bp fragments
in another BAM set:

    python calculate_dcc_bam_fragment_lengths.py \
      --bam_a '/mnt/c/Snyder_bams/CH01/BH01.bam' \
      --length_a 147 \
      --bam_b '/mnt/c/Snyder_bams/CH02/BH02.bam' \
      --length_b 167 \
      --label_a BH01_147 \
      --label_b BH02_167 \
      --chrom_sizes /mnt/d/Snyder_bams/male/chrom.sizes \
      --scope chromosome \
      --chromosome chr20 \
      --dmax 50 \
      --mapq 30

Same-length cross-sample comparison:

    python calculate_dcc_bam_fragment_lengths.py \
      --bam_a '/mnt/c/Snyder_bams/CH01/BH01.bam' \
      --length_a 147 \
      --bam_b '/mnt/c/Snyder_bams/CH02/BH02.bam' \
      --length_b 147 \
      --label_a BH01_147 \
      --label_b BH02_147 \
      --chrom_sizes /mnt/d/Snyder_bams/male/chrom.sizes \
      --scope chromosome \
      --chromosome chr20 \
      --dmax 50 \
      --mapq 30

Multiple BAMs can be supplied with quoted shell globs or a quoted space-separated
list, for example:

    --bam_a '/data/groupA/*.bam'
    --bam_b '/data/groupB/*.bam'

Position handling:
  - --position_a/--position_b dyad:
      odd fragment length, e.g. 147: one dyad base at start + length//2, weight 1
      even fragment length, e.g. 146: two middle bases, each weight 0.5
  - --position_a/--position_b left_end:
      genomic left fragment end = fragment_start, weight 1
  - --position_a/--position_b right_end:
      genomic right fragment end = fragment_end - 1, weight 1

Example comparing 147-bp left ends to 167-bp right ends:

    python calculate_dcc_bam_fragment_positions.py \
      --bam_a '/mnt/c/Snyder_bams/CH01/BH01.bam' \
      --length_a 147 \
      --position_a left_end \
      --bam_b '/mnt/c/Snyder_bams/CH01/BH01.bam' \
      --length_b 167 \
      --position_b right_end \
      --chrom_sizes /mnt/d/Snyder_bams/male/chrom.sizes \
      --scope chromosome \
      --chromosome chr20 \
      --dmax 50 \
      --mapq 30

BAM handling:
  - counts each paired-end fragment once using records with TLEN > 0
  - requires proper pairs by default
  - excludes unmapped, mate-unmapped, secondary, supplementary, QC-fail reads
  - excludes BAM-flagged duplicate reads by default
"""

import argparse
import glob
import os
from collections import defaultdict

import numpy as np
import pandas as pd
from tqdm import tqdm

try:
    import pysam
except ImportError as exc:
    raise SystemExit(
        "This script requires pysam. Install it with: pip install pysam"
    ) from exc


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


def expand_path_patterns(path_patterns, file_type_name="file"):
    files = []

    for pattern in str(path_patterns).split():
        matches = sorted(glob.glob(pattern))
        if matches:
            files.extend(matches)
        elif os.path.exists(pattern):
            files.append(pattern)

    seen = set()
    unique_files = []
    for path in files:
        if path not in seen:
            unique_files.append(path)
            seen.add(path)

    if not unique_files:
        raise FileNotFoundError(f"No {file_type_name}s match: {path_patterns}")

    return unique_files


def expand_bam_patterns(bam_patterns):
    return expand_path_patterns(bam_patterns, file_type_name="BAM")


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


def build_bam_chrom_lookup(bam):
    lookup = {}
    for ref_name in bam.references:
        lookup[ref_name] = ref_name
        lookup[standardize_chromosome_name(ref_name)] = ref_name

        if ref_name.startswith("chr"):
            no_chr = ref_name[3:]
            lookup[no_chr] = ref_name
            if no_chr == "M":
                lookup["MT"] = ref_name
        else:
            lookup[f"chr{ref_name}"] = ref_name
            if ref_name == "MT":
                lookup["chrM"] = ref_name

    return lookup


def resolve_bam_chrom(chrom, bam_chrom_lookup):
    if chrom in bam_chrom_lookup:
        return bam_chrom_lookup[chrom]

    std = standardize_chromosome_name(chrom)
    if std in bam_chrom_lookup:
        return bam_chrom_lookup[std]

    if std.startswith("chr") and std[3:] in bam_chrom_lookup:
        return bam_chrom_lookup[std[3:]]

    return None


def read_passes_filters(
    read,
    mapq=0,
    require_proper_pairs=True,
    include_duplicate_flag=False,
):
    if read.is_unmapped:
        return False
    if read.mate_is_unmapped:
        return False
    if read.is_secondary:
        return False
    if read.is_supplementary:
        return False
    if read.is_qcfail:
        return False
    if read.is_duplicate and not include_duplicate_flag:
        return False
    if require_proper_pairs and not read.is_proper_pair:
        return False
    if read.mapping_quality < mapq:
        return False
    if read.template_length <= 0:
        return False
    return True


def fragment_dyad_positions_and_weights(fragment_start, fragment_length):
    """
    Return dyad position(s) and weights for a fragment.

    For odd lengths, the exact centre base receives weight 1.
    For even lengths, the two middle bases each receive weight 0.5.
    """
    if fragment_length % 2 == 1:
        return [(fragment_start + fragment_length // 2, 1.0)]

    left_mid = fragment_start + fragment_length // 2 - 1
    right_mid = fragment_start + fragment_length // 2
    return [(left_mid, 0.5), (right_mid, 0.5)]


def fragment_positions_and_weights(fragment_start, fragment_length, position_type):
    """
    Return fragment-derived signal position(s) and weights.

    Coordinates are 0-based genomic base positions. The right end is therefore
    fragment_end - 1, not the half-open BED-style fragment_end coordinate.
    """
    fragment_end = fragment_start + fragment_length

    if position_type == "dyad":
        return fragment_dyad_positions_and_weights(fragment_start, fragment_length)

    if position_type == "left_end":
        return [(fragment_start, 1.0)]

    if position_type == "right_end":
        return [(fragment_end - 1, 1.0)]

    raise ValueError(
        f"Unknown position_type: {position_type}. "
        "Expected one of: dyad, left_end, right_end."
    )


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


def extract_fragment_signal_from_bams(
    bam_files,
    fragment_length,
    position_type,
    chrom,
    start,
    end,
    mapq=0,
    require_proper_pairs=True,
    include_duplicate_flag=False,
    max_duplicates=0,
    verbose=False,
):
    """
    Extract fragment-derived signal positions for one BAM set and one fragment length.

    max_duplicates means the maximum number of identical fragment coordinates
    retained per BAM. A value of 0 means no coordinate-count cap is applied.
    BAM-flagged duplicates are still excluded unless --include_duplicate_flag is
    used.
    """
    all_positions = []
    all_values = []

    fetch_start = max(0, start - fragment_length - 5)
    fetch_end = end + fragment_length + 5

    for bam_path in bam_files:
        try:
            bam = pysam.AlignmentFile(bam_path, "rb")
        except FileNotFoundError:
            raise FileNotFoundError(f"BAM not found: {bam_path}")

        with bam:
            bam_chrom_lookup = build_bam_chrom_lookup(bam)
            bam_chrom = resolve_bam_chrom(chrom, bam_chrom_lookup)

            if bam_chrom is None:
                if verbose:
                    print(f"Warning: {chrom} not found in {bam_path}. Skipping.")
                continue

            coordinate_counts = defaultdict(int)

            try:
                iterator = bam.fetch(bam_chrom, fetch_start, fetch_end)
            except ValueError as exc:
                if verbose:
                    print(
                        f"Warning: could not fetch {bam_chrom}:{fetch_start}-{fetch_end} "
                        f"from {bam_path}: {exc}"
                    )
                continue

            for read in iterator:
                if not read_passes_filters(
                    read,
                    mapq=mapq,
                    require_proper_pairs=require_proper_pairs,
                    include_duplicate_flag=include_duplicate_flag,
                ):
                    continue

                tlen = int(read.template_length)
                if tlen != fragment_length:
                    continue

                fragment_start = int(read.reference_start)
                fragment_end = fragment_start + tlen

                if fragment_end <= fragment_start:
                    continue

                # Count duplicate fragment coordinates within each BAM if requested.
                coord_key = (bam_chrom, fragment_start, fragment_end)
                if max_duplicates is not None and max_duplicates > 0:
                    coordinate_counts[coord_key] += 1
                    if coordinate_counts[coord_key] > max_duplicates:
                        continue

                for signal_pos, weight in fragment_positions_and_weights(
                    fragment_start,
                    fragment_length,
                    position_type,
                ):
                    if start <= signal_pos < end:
                        all_positions.append(signal_pos)
                        all_values.append(weight)

    if not all_positions:
        return None, None

    return group_positions_values(
        np.asarray(all_positions, dtype=np.int64),
        np.asarray(all_values, dtype=float),
    )


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

        for lag in [-20, -10, -5, 0, 5, 10, 20]:
            if -dmax <= lag <= dmax:
                idx = lag + dmax
                rows.append({
                    "Metric": f"lag_{lag}",
                    "Lag_or_Distance": lag,
                    "DCC Value": float(values[idx]),
                })

        for lag in [5, -5, 10, -10, 20, -20]:
            if -dmax <= lag <= dmax and -dmax <= 0 <= dmax:
                rows.append({
                    "Metric": f"lag_{lag}_minus_lag_0",
                    "Lag_or_Distance": lag,
                    "DCC Value": float(values[dmax + lag] - values[dmax]),
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

        for distance in [5, 10, 20]:
            if distance <= dmax:
                rows.append({
                    "Metric": f"distance_{distance}_minus_distance_0",
                    "Lag_or_Distance": distance,
                    "DCC Value": float(values[distance] - values[0]),
                })

    pd.DataFrame(rows).to_csv(output_file, sep="\t", index=False)
    print(f"DCC shift summary saved to {output_file}")


def calculate_dcc_streaming_from_bams(
    bam_patterns_a,
    length_a,
    position_a,
    bam_patterns_b,
    length_b,
    position_b,
    chromatin_df,
    output_prefix,
    dmax=100,
    mapq=0,
    require_proper_pairs=True,
    include_duplicate_flag=False,
    max_duplicates=0,
    min_region_length=None,
    normalize_dcc=True,
    normalize_by_signal_totals=False,
    summary_lag_window=50,
    signed_lags=False,
    verbose=False,
):
    print("Running calculate_dcc_streaming_from_bams")

    bam_files_a = expand_bam_patterns(bam_patterns_a)
    bam_files_b = expand_bam_patterns(bam_patterns_b)

    print(f"A BAMs: {len(bam_files_a)}; length A: {length_a}; position A: {position_a}")
    print(f"B BAMs: {len(bam_files_b)}; length B: {length_b}; position B: {position_b}")

    if min_region_length is None:
        min_region_length = dmax + 1

    states = sorted(chromatin_df["State"].unique())
    opportunities_cache = {}
    completed_outputs = {}

    for state in tqdm(states, desc="Processing states"):
        state_df = chromatin_df[chromatin_df["State"] == state]

        # Internally always calculate signed DCC.
        # Output is later either kept signed or collapsed to absolute distances.
        dcc = np.zeros(2 * dmax + 1, dtype=float)
        opportunities = np.zeros(2 * dmax + 1, dtype=float)
        total_signal_a = 0.0
        total_signal_b = 0.0
        used_regions = 0
        regions_with_a = 0
        regions_with_b = 0

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

            positions_a, values_a = extract_fragment_signal_from_bams(
                bam_files=bam_files_a,
                fragment_length=length_a,
                position_type=position_a,
                chrom=chrom,
                start=start,
                end=end,
                mapq=mapq,
                require_proper_pairs=require_proper_pairs,
                include_duplicate_flag=include_duplicate_flag,
                max_duplicates=max_duplicates,
                verbose=verbose,
            )

            positions_b, values_b = extract_fragment_signal_from_bams(
                bam_files=bam_files_b,
                fragment_length=length_b,
                position_type=position_b,
                chrom=chrom,
                start=start,
                end=end,
                mapq=mapq,
                require_proper_pairs=require_proper_pairs,
                include_duplicate_flag=include_duplicate_flag,
                max_duplicates=max_duplicates,
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

            if positions_a is not None:
                regions_with_a += 1
                total_signal_a += float(np.sum(values_a))
            if positions_b is not None:
                regions_with_b += 1
                total_signal_b += float(np.sum(values_b))

            if positions_a is None or positions_b is None:
                continue

            if normalize_dcc:
                update_cross_opportunities(
                    opportunities,
                    region_length,
                    dmax,
                    opportunities_cache,
                )

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
        suffix_parts.append(f"A{length_a}_{position_a}_B{length_b}_{position_b}")
        suffix_parts.append("signed_lags" if signed_lags else "absolute_distances")
        suffix_parts.append("opportunity_normalized" if normalize_dcc else "raw")
        if normalize_by_signal_totals:
            suffix_parts.append("signal_total_normalized")
        suffix = "_".join(suffix_parts)

        output_file = (
            f"{output_prefix}_{sanitized_state}_"
            f"bam_streaming_DCC_values_{suffix}.tsv"
        )

        summary_file = (
            f"{output_prefix}_{sanitized_state}_"
            f"bam_streaming_DCC_shift_summary_{suffix}.tsv"
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
            "used_regions_with_both_A_and_B": used_regions,
            "regions_with_A": regions_with_a,
            "regions_with_B": regions_with_b,
            "total_signal_a": total_signal_a,
            "total_signal_b": total_signal_b,
            "length_a": length_a,
            "position_a": position_a,
            "length_b": length_b,
            "position_b": position_b,
            "signed_lags": signed_lags,
        }

        print(
            f"Completed DCC for state: {state} "
            f"using {used_regions} regions with both A and B; "
            f"regions_with_A={regions_with_a}; regions_with_B={regions_with_b}; "
            f"sum(A)={total_signal_a:.6g}; sum(B)={total_signal_b:.6g}; "
            f"signals=A{length_a}:{position_a}/B{length_b}:{position_b}; "
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
        default_label_a = (
            f"{compact_pattern_label(args.bam_a, 'A')}_{args.length_a}bp_{args.position_a}"
        )
        default_label_b = (
            f"{compact_pattern_label(args.bam_b, 'B')}_{args.length_b}bp_{args.position_b}"
        )
        output_prefix = get_output_prefix(
            args.bam_a,
            args.bam_b,
            args.label_a or default_label_a,
            args.label_b or default_label_b,
        )

    if args.scope == "chromosome" and args.chromosome is not None:
        output_prefix = f"{output_prefix}_{standardize_chromosome_name(args.chromosome)}"

    calculate_dcc_streaming_from_bams(
        bam_patterns_a=args.bam_a,
        length_a=args.length_a,
        position_a=args.position_a,
        bam_patterns_b=args.bam_b,
        length_b=args.length_b,
        position_b=args.position_b,
        chromatin_df=chromatin_df,
        output_prefix=output_prefix,
        dmax=args.dmax,
        mapq=args.mapq,
        require_proper_pairs=not args.no_require_proper_pairs,
        include_duplicate_flag=args.include_duplicate_flag,
        max_duplicates=args.max_duplicates,
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
            "Calculate streaming Distance Cross-Correlation directly from two "
            "sets of BAM files, allowing signal A and signal B to use different "
            "fragment lengths and position types: dyad, left_end, or right_end. "
            "Internally calculates signed DCC[lag] = sum_x "
            "A[x] * B[x + lag]. By default, output is collapsed to absolute "
            "distances. Use --signed_lags to keep directional signed lags."
        )
    )

    parser.add_argument(
        "--bam_a",
        required=True,
        help=(
            "Space-separated list of BAM file patterns for signal A, "
            "e.g. '/data/A/*.bam' or '/data/A1.bam /data/A2.bam'. "
            "Quote shell globs/lists."
        ),
    )

    parser.add_argument(
        "--length_a",
        type=int,
        required=True,
        help="Exact paired-end fragment length to use for signal A, e.g. 147.",
    )

    parser.add_argument(
        "--position_a",
        choices=["dyad", "left_end", "right_end"],
        default="dyad",
        help=(
            "Fragment-derived position to use for signal A. "
            "dyad keeps the previous behaviour. left_end uses fragment_start. "
            "right_end uses fragment_end - 1. Default: dyad."
        ),
    )

    parser.add_argument(
        "--bam_b",
        required=True,
        help=(
            "Space-separated list of BAM file patterns for signal B, "
            "e.g. '/data/B/*.bam' or '/data/B1.bam /data/B2.bam'. "
            "Quote shell globs/lists."
        ),
    )

    parser.add_argument(
        "--length_b",
        type=int,
        required=True,
        help="Exact paired-end fragment length to use for signal B, e.g. 167.",
    )

    parser.add_argument(
        "--position_b",
        choices=["dyad", "left_end", "right_end"],
        default="dyad",
        help=(
            "Fragment-derived position to use for signal B. "
            "dyad keeps the previous behaviour. left_end uses fragment_start. "
            "right_end uses fragment_end - 1. Default: dyad."
        ),
    )

    parser.add_argument(
        "--label_a",
        default=None,
        help="Optional short label for signal A, e.g. BH01_147.",
    )

    parser.add_argument(
        "--label_b",
        default=None,
        help="Optional short label for signal B, e.g. BH02_167.",
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
            "e.g. chr1, 1, chr20, chrX."
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
            "Default: 100. Use 50 for a focused 5/10-bp shift test, or 1500 "
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
        "--mapq",
        type=int,
        default=0,
        help="Minimum read MAPQ for the counted TLEN-positive read. Default: 0.",
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
        "--max_duplicates",
        type=int,
        default=0,
        help=(
            "Maximum identical fragment coordinates retained per BAM. "
            "0 means no coordinate-count cap. BAM-flagged duplicates are still "
            "excluded unless --include_duplicate_flag is used. Default: 0."
        ),
    )

    parser.add_argument(
        "--include_duplicate_flag",
        action="store_true",
        help=(
            "Include reads marked as duplicate in the BAM flag. "
            "Default is to exclude BAM-flagged duplicates."
        ),
    )

    parser.add_argument(
        "--no_require_proper_pairs",
        action="store_true",
        help=(
            "Do not require the BAM proper-pair flag. Default requires proper pairs."
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
        help="Print warnings from missing chromosomes or fetch failures.",
    )

    args = parser.parse_args()
    main(args)
