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


def read_chromatin_states(chromatin_bed, convert_to_euchromatin=False):
    print("Running read_chromatin_states")

    df = pd.read_csv(
        chromatin_bed,
        sep="\t",
        header=None,
        usecols=[0, 1, 2, 3],
        names=["Chromosome", "Start", "End", "State"],
        dtype={"Chromosome": str, "Start": int, "End": int, "State": str},
    )

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

    df = df.sort_values(["State", "Chromosome", "Start", "End"]).reset_index(drop=True)
    return df


def expand_bedgraph_to_positions(temp_bedgraph, value_limit=None):
    """
    Read UCSC bigWigToBedGraph output and expand intervals to per-base values.
    """
    if not os.path.exists(temp_bedgraph) or os.path.getsize(temp_bedgraph) == 0:
        return None, None

    bg = pd.read_csv(
        temp_bedgraph,
        sep="\t",
        header=None,
        names=["Chromosome", "Start", "End", "Value"],
        dtype={"Chromosome": str, "Start": int, "End": int, "Value": float},
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


def update_opportunities(opportunities, region_length, dmax, opportunities_cache):
    if region_length not in opportunities_cache:
        opportunities_cache[region_length] = np.zeros(dmax, dtype=float)
        for lag in range(1, min(dmax, region_length)):
            opportunities_cache[region_length][lag] = region_length - lag + 1

    opportunities += opportunities_cache[region_length]


def update_dac_from_region(dac, positions, values, dmax):
    """
    Add one extracted interval to DAC.

    For each pair of bases separated by distance d:
        dac[d] += value1 * value2
    """
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

            if dist >= dmax:
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
    dmax=1000,
    value_limit=None,
    min_region_length=2000,
    normalize_dac=False,
):
    """
    Streaming DAC calculation.

    Workflow:
        for each chromatin state:
            for each interval in that state:
                run bigWigToBedGraph for that interval
                expand extracted signal
                immediately update DAC
                discard extracted signal

        after each chromatin state is fully processed:
            optionally normalize DAC
            save that state's DAC immediately
    """
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
        for state in tqdm(states, desc="Processing chromatin states"):
            state_df = chromatin_df[chromatin_df["State"] == state]

            dac = np.zeros(dmax, dtype=float)
            opportunities = np.zeros(dmax, dtype=float)

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
                    region_length = end - start

                    if region_length < min_region_length:
                        continue

                    temp_bedgraph = os.path.join(
                        tmpdir,
                        f"{sanitize_filename(state)}_{chrom}_{start}_{end}_{region_index}.bedgraph"
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

                if np.any(non_zero):
                    max_opportunities = np.max(opportunities)

                    scaling = np.zeros_like(opportunities, dtype=float)
                    scaling[non_zero] = max_opportunities / opportunities[non_zero]

                    dac *= scaling

            sanitized_state = sanitize_filename(state)
            suffix = "normalized" if normalize_dac else "raw"
            output_file = f"{output_prefix}_{sanitized_state}_bigwig_streaming_DAC_values_{suffix}.tsv"
            save_dac_to_tsv(dac, output_file)
            completed_outputs[state] = output_file
            print(f"Completed and saved DAC for state: {state}")

    return completed_outputs


def main(
    bigwig_patterns,
    chromatin_bed,
    dmax,
    value_limit,
    convert_to_euchromatin,
    min_region_length,
    normalize_dac,
):
    chromatin_df = read_chromatin_states(
        chromatin_bed,
        convert_to_euchromatin=convert_to_euchromatin,
    )

    output_prefix = get_output_prefix(bigwig_patterns)

    calculate_dac_streaming_from_bigwig(
        bigwig_patterns=bigwig_patterns,
        chromatin_df=chromatin_df,
        output_prefix=output_prefix,
        dmax=dmax,
        value_limit=value_limit,
        min_region_length=min_region_length,
        normalize_dac=normalize_dac,
    )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description=(
            "Calculate streaming Distance Autocorrelation from bigWig signal "
            "over chromatin-state BED intervals using UCSC bigWigToBedGraph."
        )
    )

    parser.add_argument(
        "--bigwig",
        required=True,
        help="Space-separated list of bigWig file patterns, e.g. 'CH01_chr*.bw'",
    )

    parser.add_argument(
        "--chromatin_bed",
        required=True,
        help="BED file containing chromatin states. First four columns: chrom, start, end, state.",
    )

    parser.add_argument(
        "--dmax",
        type=int,
        default=1000,
        help="Maximum distance for DAC calculation.",
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
        default=2000,
        help="Minimum chromatin interval length to include.",
    )

    parser.add_argument(
        "--convert_to_euchromatin",
        action="store_true",
        help="Collapse selected ChromHMM states into Euchromatin.",
    )

    parser.add_argument(
        "--normalize_dac",
        action="store_true",
        help="Apply opportunity-based DAC normalization. Default is off.",
    )

    args = parser.parse_args()

    main(
        bigwig_patterns=args.bigwig,
        chromatin_bed=args.chromatin_bed,
        dmax=args.dmax,
        value_limit=args.value_limit,
        convert_to_euchromatin=args.convert_to_euchromatin,
        min_region_length=args.min_region_length,
        normalize_dac=args.normalize_dac,
    )