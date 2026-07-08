#!/usr/bin/env python3

import argparse
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib.colors as mcolors
from matplotlib.ticker import MultipleLocator


def parse_dict(text, value_type=float):
    d = {}
    for item in text.split(","):
        k, v = item.split(":")
        d[int(k)] = value_type(v)
    return d


def normalize_probs(d):
    total = sum(d.values())
    return {k: v / total for k, v in d.items()}


def parse_float_list(text):
    return [float(x) for x in text.split(",")]


def generate_dyads(n_dyads, spacing_probs, seed=1, placement_sd=0):
    """
    Generate original dyad map.

    The intended spacing is first sampled from spacing_probs.
    Then placement error is added to that spacing before placing the next dyad.
    This placement error is used only once to build the original map.
    """

    rng = np.random.default_rng(seed)

    spacings = np.array(list(spacing_probs.keys()))
    probs = np.array(list(spacing_probs.values()))

    intended_spacings = rng.choice(spacings, size=n_dyads - 1, p=probs)

    if placement_sd > 0:
        placement_errors = rng.normal(0, placement_sd, size=n_dyads - 1)
        placement_errors = np.rint(placement_errors).astype(int)
    else:
        placement_errors = np.zeros(n_dyads - 1, dtype=int)

    placed_spacings = intended_spacings + placement_errors
    placed_spacings = np.maximum(placed_spacings, 1)

    dyads = np.zeros(n_dyads, dtype=int)
    dyads[1:] = np.cumsum(placed_spacings)

    return dyads, intended_spacings, placed_spacings, placement_errors


def apply_pairwise_shifts(dyads, shift_map, phase=0):
    shifted = dyads.copy()

    for i in range(phase, len(dyads) - 1, 2):
        old_spacing = dyads[i + 1] - dyads[i]

        if old_spacing in shift_map:
            new_spacing = shift_map[old_spacing]
            shifted[i + 1] = shifted[i] + new_spacing

    return shifted


def choose_non_adjacent_gaps(n_gaps, fraction, rng):
    max_selectable = (n_gaps + 1) // 2
    n_select = int(round(max_selectable * fraction))
    n_select = max(0, min(n_select, max_selectable))

    if n_select == 0:
        return np.array([], dtype=int)

    base = np.sort(
        rng.choice(
            n_gaps - n_select + 1,
            size=n_select,
            replace=False
        )
    )

    return base + np.arange(n_select)


def apply_random_embedded_shifts(dyads, shift_map, fraction, rng):
    shifted = dyads.copy()
    n_gaps = len(dyads) - 1
    selected_gaps = choose_non_adjacent_gaps(n_gaps, fraction, rng)

    for i in selected_gaps:
        old_spacing = dyads[i + 1] - dyads[i]

        if old_spacing in shift_map:
            new_spacing = shift_map[old_spacing]
            shifted[i + 1] = shifted[i] + new_spacing

    return shifted, selected_gaps


def call_positions(dyads, sd, rng):
    """
    Calling error is applied independently to each replicate.
    """

    if sd == 0:
        return dyads.copy()

    errors = rng.normal(0, sd, size=len(dyads))
    return np.rint(dyads + errors).astype(int)


def spacings_from_positions(pos):
    return np.diff(np.sort(pos))


def counts_table(spacings, calling_sd, compressed_fraction, label):
    distances, counts = np.unique(spacings, return_counts=True)
    total = counts.sum()

    return pd.DataFrame({
        "label": label,
        "compressed_fraction": compressed_fraction,
        "calling_sd": calling_sd,
        "distance": distances.astype(int),
        "count": counts.astype(int),
        "proportion": counts / total
    })


def apply_x_ticks(ax):
    ax.xaxis.set_major_locator(MultipleLocator(5))
    ax.xaxis.set_minor_locator(MultipleLocator(1))
    ax.tick_params(axis="x", which="major", length=6)
    ax.tick_params(axis="x", which="minor", length=3)


def print_debug_positions(original, compressed, n=12):
    print("\nDEBUG: first positions and spacings")
    print("----------------------------------")
    print(f"Original positions first {n}:")
    print(original[:n])
    print(f"Original spacings first {n - 1}:")
    print(np.diff(original[:n]))

    print(f"\nCompressed positions first {n}:")
    print(compressed[:n])
    print(f"Compressed spacings first {n - 1}:")
    print(np.diff(compressed[:n]))
    print("----------------------------------\n")


def get_replicate_types(frac, n_replicates):
    n_compressed = int(round(n_replicates * frac))
    n_original = n_replicates - n_compressed
    return ["original"] * n_original + ["compressed"] * n_compressed, n_original, n_compressed


def simulate_called_positions_for_fraction_sd(
    original,
    compressed,
    frac,
    calling_sd,
    n_replicates,
    seed,
    sd_i,
    compression_mode="embedded",
    shift_map=None
):
    rng = np.random.default_rng(
        seed
        + int(round(frac * 1000000))
        + 100000 * sd_i
    )

    called_position_arrays = []

    if compression_mode == "map_mix":
        rep_types, n_original, n_compressed = get_replicate_types(frac, n_replicates)

        for rep_type in rep_types:
            dyads = original if rep_type == "original" else compressed
            called = call_positions(dyads, calling_sd, rng)
            called_position_arrays.append(called)

        return called_position_arrays, n_original, n_compressed

    if shift_map is None:
        raise ValueError("shift_map is required when compression_mode='embedded'")

    selected_counts = []

    for _ in range(n_replicates):
        dyads, selected_gaps = apply_random_embedded_shifts(
            original,
            shift_map=shift_map,
            fraction=frac,
            rng=rng
        )

        selected_counts.append(len(selected_gaps))

        called = call_positions(dyads, calling_sd, rng)
        called_position_arrays.append(called)

    n_original = 0
    n_compressed = int(round(np.mean(selected_counts))) if selected_counts else 0

    return called_position_arrays, n_original, n_compressed


def calculate_dac_from_positions(positions, dmax, normalize_dac=True):
    positions = np.asarray(positions, dtype=int)
    positions = positions[positions >= 0]

    if len(positions) < 2:
        return np.zeros(dmax + 1, dtype=float)

    unique_pos, counts = np.unique(positions, return_counts=True)

    max_pos = int(unique_pos.max())
    region_length = max_pos + 1

    dac = np.zeros(dmax + 1, dtype=float)

    for i, p1 in enumerate(unique_pos):
        c1 = counts[i]

        max_target = p1 + dmax
        j_start = i + 1
        j_end = np.searchsorted(unique_pos, max_target, side="right")

        targets = unique_pos[j_start:j_end]
        target_counts = counts[j_start:j_end]

        distances = targets - p1
        dac[distances] += c1 * target_counts

    if normalize_dac:
        opportunities = np.zeros(dmax + 1, dtype=float)

        max_lag = min(dmax, region_length - 1)

        for lag in range(1, max_lag + 1):
            opportunities[lag] = region_length - lag

        non_zero = opportunities != 0

        if np.any(non_zero):
            max_opportunities = np.max(opportunities)
            scaling = np.zeros_like(opportunities)
            scaling[non_zero] = max_opportunities / opportunities[non_zero]
            dac *= scaling

    return dac


def dac_to_dataframe(dac, frac, calling_sd, label):
    dac_values = dac[1:]
    total = dac_values.sum()

    if total != 0:
        dac_percent = dac_values / total * 100
    else:
        dac_percent = np.zeros_like(dac_values)

    return pd.DataFrame({
        "label": label,
        "compressed_fraction": frac,
        "calling_sd": calling_sd,
        "Distance": np.arange(1, len(dac)),
        "DAC Value": dac_values,
        "DAC Value Percent": dac_percent
    })


def plot_error_distributions(error_sds, out_file):
    max_sd = max(error_sds)
    fig, ax = plt.subplots(figsize=(10, 5))

    norm = mcolors.Normalize(vmin=min(error_sds), vmax=max(error_sds))
    cmap = cm.coolwarm

    if max_sd == 0:
        ax.axvline(0, color=cmap(norm(0)), linewidth=1.2)
    else:
        x = np.linspace(-4 * max_sd, 4 * max_sd, 1000)

        for sd in error_sds:
            color = cmap(norm(sd))

            if sd == 0:
                ax.axvline(0, color=color, linewidth=1.2)
            else:
                y = (
                    1 / (sd * np.sqrt(2 * np.pi))
                    * np.exp(-0.5 * (x / sd) ** 2)
                )
                ax.plot(x, y, color=color, linewidth=1.2)

    sm = cm.ScalarMappable(norm=norm, cmap=cmap)
    sm.set_array([])
    cbar = fig.colorbar(sm, ax=ax)
    cbar.set_label("Calling error SD (bp)")

    ax.set_xlabel("Dyad calling error, bp")
    ax.set_ylabel("Density")
    ax.set_title("Dyad-calling error distributions")
    apply_x_ticks(ax)

    plt.tight_layout()
    plt.savefig(out_file, dpi=300)
    plt.close()


def simulate_counts_for_fraction(
    original,
    compressed,
    frac,
    n_replicates,
    error_sds,
    combine_method,
    seed,
    out_prefix,
    compression_mode,
    shift_map
):
    rep_types, n_original, n_compressed = get_replicate_types(frac, n_replicates)

    if compression_mode == "embedded":
        label = (
            f"{combine_method}_embedded_"
            f"frac_{frac:g}_"
            f"replicates_{n_replicates}"
        )
    else:
        label = (
            f"{combine_method}_mapmix_"
            f"frac_{frac:g}_"
            f"original_{n_original}_"
            f"compressed_{n_compressed}"
        )

    fraction_tables = []

    for sd_i, calling_sd in enumerate(error_sds):
        called_arrays, _, _ = simulate_called_positions_for_fraction_sd(
            original=original,
            compressed=compressed,
            frac=frac,
            calling_sd=calling_sd,
            n_replicates=n_replicates,
            seed=seed,
            sd_i=sd_i,
            compression_mode=compression_mode,
            shift_map=shift_map
        )

        if combine_method == "positions":
            final_spacings = spacings_from_positions(np.concatenate(called_arrays))
        else:
            final_spacings = np.concatenate([
                spacings_from_positions(x) for x in called_arrays
            ])

        count_df = counts_table(
            final_spacings,
            calling_sd=calling_sd,
            compressed_fraction=frac,
            label=label
        )

        fraction_tables.append(count_df)

    fraction_df = pd.concat(fraction_tables, ignore_index=True)

    fraction_out = f"{out_prefix}_{label}_distance_counts.tsv"
    fraction_df.to_csv(fraction_out, sep="\t", index=False)

    return fraction_df, label, fraction_out


def plot_fraction_errors(fraction_df, label, error_sds, x_min, x_max, out_prefix):
    fig, ax = plt.subplots(figsize=(12, 5))

    norm = mcolors.Normalize(vmin=min(error_sds), vmax=max(error_sds))
    cmap = cm.coolwarm

    for sd in error_sds:
        plot_df = fraction_df[
            (fraction_df["calling_sd"] == sd) &
            (fraction_df["distance"] >= x_min) &
            (fraction_df["distance"] <= x_max)
        ]

        color = cmap(norm(sd))

        ax.plot(
            plot_df["distance"],
            plot_df["proportion"],
            marker=".",
            linestyle="-",
            linewidth=0.8,
            markersize=3,
            color=color
        )

    sm = cm.ScalarMappable(norm=norm, cmap=cmap)
    sm.set_array([])
    cbar = fig.colorbar(sm, ax=ax)
    cbar.set_label("Calling error SD (bp)")

    ax.set_xlim(x_min, x_max)
    ax.set_xlabel("Called dyad spacing, bp")
    ax.set_ylabel("Proportion")
    ax.set_title(label)
    apply_x_ticks(ax)

    plt.tight_layout()

    plot_out = f"{out_prefix}_{label}_spacing_counts.png"
    plt.savefig(plot_out, dpi=300)
    plt.close()

    return plot_out


def plot_all_combinations(all_counts, compressed_fractions, error_sds, x_min, x_max, out_prefix):
    plot_outputs = []

    norm = mcolors.Normalize(
        vmin=min(compressed_fractions),
        vmax=max(compressed_fractions)
    )

    cmap = cm.coolwarm

    for sd in error_sds:
        fig, ax = plt.subplots(figsize=(12, 5))

        for frac in compressed_fractions:
            plot_df = all_counts[
                (all_counts["calling_sd"] == sd) &
                (all_counts["compressed_fraction"] == frac) &
                (all_counts["distance"] >= x_min) &
                (all_counts["distance"] <= x_max)
            ]

            color = cmap(norm(frac))

            ax.plot(
                plot_df["distance"],
                plot_df["proportion"],
                marker=".",
                linestyle="-",
                linewidth=0.8,
                markersize=3,
                color=color
            )

        sm = cm.ScalarMappable(norm=norm, cmap=cmap)
        sm.set_array([])
        cbar = fig.colorbar(sm, ax=ax)
        cbar.set_label("Compressed fraction")

        ax.set_xlim(x_min, x_max)
        ax.set_xlabel("Called dyad spacing, bp")
        ax.set_ylabel("Proportion")
        ax.set_title(f"All compressed fractions, calling SD = {sd:g} bp")
        apply_x_ticks(ax)

        plt.tight_layout()

        plot_out = f"{out_prefix}_all_fractions_calling_sd_{sd:g}_spacing_counts.png"
        plt.savefig(plot_out, dpi=300)
        plt.close()

        plot_outputs.append(plot_out)

    return plot_outputs


def plot_called_position_counts_region(
    original,
    compressed,
    compressed_fractions,
    error_sds,
    n_replicates,
    region_start,
    region_length,
    seed,
    out_prefix,
    compression_mode,
    shift_map
):
    region_end = region_start + region_length
    region_positions = np.arange(region_start, region_end + 1)

    outputs = []

    norm = mcolors.Normalize(
        vmin=min(compressed_fractions),
        vmax=max(compressed_fractions)
    )

    cmap = cm.coolwarm

    for sd_i, calling_sd in enumerate(error_sds):
        fig, ax = plt.subplots(figsize=(14, 4))

        for frac in compressed_fractions:
            called_arrays, _, _ = simulate_called_positions_for_fraction_sd(
                original=original,
                compressed=compressed,
                frac=frac,
                calling_sd=calling_sd,
                n_replicates=n_replicates,
                seed=seed,
                sd_i=sd_i,
                compression_mode=compression_mode,
                shift_map=shift_map
            )

            combined_positions = np.concatenate(called_arrays)

            region_called = combined_positions[
                (combined_positions >= region_start) &
                (combined_positions <= region_end)
            ]

            counts = np.zeros(region_length + 1, dtype=int)
            idx = region_called - region_start
            np.add.at(counts, idx, 1)

            label = (
                f"called_position_counts_"
                f"frac_{frac:g}_"
                f"calling_sd_{calling_sd:g}"
            )

            out_tsv = f"{out_prefix}_{label}.tsv"
            pd.DataFrame({
                "position": region_positions,
                "count": counts
            }).to_csv(out_tsv, sep="\t", index=False)

            color = cmap(norm(frac))

            ax.plot(
                region_positions,
                counts,
                linewidth=0.9,
                color=color
            )

            outputs.append(out_tsv)

        sm = cm.ScalarMappable(norm=norm, cmap=cmap)
        sm.set_array([])
        cbar = fig.colorbar(sm, ax=ax)
        cbar.set_label("Compressed fraction")

        ax.set_xlim(region_start, region_end)
        ax.set_xlabel("Position, bp")
        ax.set_ylabel("Called dyad count")
        ax.set_title(
            f"Combined called dyad counts, calling SD={calling_sd:g} bp"
        )

        ax.xaxis.set_major_locator(MultipleLocator(250))
        ax.xaxis.set_minor_locator(MultipleLocator(50))
        ax.tick_params(axis="x", which="major", length=6)
        ax.tick_params(axis="x", which="minor", length=3)

        plt.tight_layout()

        out_png = (
            f"{out_prefix}_called_position_counts_"
            f"all_fractions_calling_sd_{calling_sd:g}.png"
        )

        plt.savefig(out_png, dpi=300)
        plt.close()

        outputs.append(out_png)

    return outputs


def run_dac_analysis(
    original,
    compressed,
    compressed_fractions,
    error_sds,
    n_replicates,
    seed,
    dmax,
    normalize_dac,
    out_prefix,
    compression_mode,
    shift_map
):
    all_dac_tables = []

    norm = mcolors.Normalize(
        vmin=min(compressed_fractions),
        vmax=max(compressed_fractions)
    )

    cmap = cm.coolwarm

    for sd_i, calling_sd in enumerate(error_sds):
        fig, ax = plt.subplots(figsize=(12, 5))

        for frac in compressed_fractions:
            called_arrays, n_original, n_compressed = simulate_called_positions_for_fraction_sd(
                original=original,
                compressed=compressed,
                frac=frac,
                calling_sd=calling_sd,
                n_replicates=n_replicates,
                seed=seed,
                sd_i=sd_i,
                compression_mode=compression_mode,
                shift_map=shift_map
            )

            combined_positions = np.concatenate(called_arrays)

            label = (
                f"dac_frac_{frac:g}_"
                f"calling_sd_{calling_sd:g}_"
                f"original_{n_original}_"
                f"compressed_{n_compressed}"
            )

            dac = calculate_dac_from_positions(
                positions=combined_positions,
                dmax=dmax,
                normalize_dac=normalize_dac
            )

            dac_df = dac_to_dataframe(
                dac=dac,
                frac=frac,
                calling_sd=calling_sd,
                label=label
            )

            dac_out = f"{out_prefix}_{label}_DAC_values.tsv"
            dac_df.to_csv(dac_out, sep="\t", index=False)

            all_dac_tables.append(dac_df)

            color = cmap(norm(frac))

            ax.plot(
                dac_df["Distance"],
                dac_df["DAC Value Percent"],
                linewidth=0.4,
                color=color
            )

        sm = cm.ScalarMappable(norm=norm, cmap=cmap)
        sm.set_array([])

        cbar = fig.colorbar(sm, ax=ax)
        cbar.set_label("Compressed fraction")

        ax.set_xlabel("Distance, bp")
        ax.set_ylabel("DAC Value Percent")
        ax.set_title(f"DAC, calling SD = {calling_sd:g} bp")

        ax.xaxis.set_major_locator(MultipleLocator(50))
        ax.xaxis.set_minor_locator(MultipleLocator(10))

        plt.tight_layout()

        plot_out = f"{out_prefix}_DAC_all_fractions_calling_sd_{calling_sd:g}.png"
        plt.savefig(plot_out, dpi=300)
        plt.close()

        print(f"Saved combined DAC plot: {plot_out}")

    all_dac_df = pd.concat(all_dac_tables, ignore_index=True)

    all_dac_out = f"{out_prefix}_all_DAC_values.tsv"
    all_dac_df.to_csv(all_dac_out, sep="\t", index=False)

    print(f"Saved all DAC values: {all_dac_out}")


def main():
    parser = argparse.ArgumentParser(
        description="Simulate dyad spacing with separate placement and calling error."
    )

    parser.add_argument("--n-dyads", type=int, default=10000)

    parser.add_argument(
        "--spacing-probs",
        default="167:0.2,177:0.35,188:0.3,200:0.15"
    )

    parser.add_argument(
        "--pair-shifts",
        default="167:162,177:172,188:183,200:194,211:206"
    )

    parser.add_argument(
        "--placement-sd",
        type=float,
        default=0,
        help="SD of placement error applied to each sampled spacing while generating the original map."
    )

    parser.add_argument(
        "--compressed-fractions",
        default="0,0.1,0.2,0.3,0.5,0.7,1"
    )

    parser.add_argument("--n-replicates", type=int, default=10)

    parser.add_argument(
        "--shift-phase",
        choices=["0", "1", "random"],
        default="0"
    )

    parser.add_argument(
        "--compression-mode",
        choices=["embedded", "map_mix"],
        default="embedded"
    )

    parser.add_argument(
        "--debug-print",
        action="store_true"
    )

    parser.add_argument(
        "--combine-method",
        choices=["distances", "positions"],
        default="distances"
    )

    parser.add_argument(
        "--plot-mode",
        choices=["errors", "fractions", "both"],
        default="errors"
    )

    parser.add_argument(
        "--error-sds",
        default="0,1,2,3,5,10",
        help="Calling error SDs. This does not affect initial placement."
    )

    parser.add_argument("--x-min", type=int, default=0)
    parser.add_argument("--x-max", type=int, default=240)

    parser.add_argument("--run-dac", action="store_true")
    parser.add_argument("--dac-dmax", type=int, default=1500)

    parser.add_argument(
        "--no-normalize-dac",
        action="store_true"
    )

    parser.add_argument(
        "--plot-position-region",
        action="store_true"
    )

    parser.add_argument("--region-start", type=int, default=0)
    parser.add_argument("--region-length", type=int, default=3000)

    parser.add_argument("--seed", type=int, default=1)
    parser.add_argument("--out-prefix", default="dyad_spacing_simulation")

    args = parser.parse_args()

    spacing_probs = normalize_probs(parse_dict(args.spacing_probs, float))
    shift_map = parse_dict(args.pair_shifts, int)

    error_sds = parse_float_list(args.error_sds)
    compressed_fractions = parse_float_list(args.compressed_fractions)

    original, intended_spacings, placed_spacings, placement_errors = generate_dyads(
        args.n_dyads,
        spacing_probs,
        seed=args.seed,
        placement_sd=args.placement_sd
    )

    pd.DataFrame({
        "spacing_index": np.arange(len(placed_spacings)),
        "intended_spacing": intended_spacings,
        "placement_error": placement_errors,
        "placed_spacing": placed_spacings
    }).to_csv(
        f"{args.out_prefix}_placed_spacings.tsv",
        sep="\t",
        index=False
    )

    if args.shift_phase == "random":
        rng_phase = np.random.default_rng(args.seed + 99999)
        phase = int(rng_phase.integers(0, 2))
    else:
        phase = int(args.shift_phase)

    print(f"Using placement SD: {args.placement_sd}")
    print(f"Using compression mode: {args.compression_mode}")
    print(f"Using compression phase for map_mix/reference compressed map: {phase}")

    compressed = apply_pairwise_shifts(
        original,
        shift_map,
        phase=phase
    )

    if args.debug_print:
        print_debug_positions(original, compressed, n=12)

    all_count_tables = []

    for frac in compressed_fractions:
        fraction_df, label, fraction_out = simulate_counts_for_fraction(
            original=original,
            compressed=compressed,
            frac=frac,
            n_replicates=args.n_replicates,
            error_sds=error_sds,
            combine_method=args.combine_method,
            seed=args.seed,
            out_prefix=args.out_prefix,
            compression_mode=args.compression_mode,
            shift_map=shift_map
        )

        all_count_tables.append(fraction_df)

        if args.plot_mode in ["errors", "both"]:
            plot_fraction_errors(
                fraction_df,
                label,
                error_sds,
                args.x_min,
                args.x_max,
                args.out_prefix
            )

    all_counts = pd.concat(all_count_tables, ignore_index=True)

    all_counts.to_csv(
        f"{args.out_prefix}_all_fraction_distance_counts.tsv",
        sep="\t",
        index=False
    )

    if args.plot_mode in ["fractions", "both"]:
        plot_all_combinations(
            all_counts,
            compressed_fractions,
            error_sds,
            args.x_min,
            args.x_max,
            args.out_prefix
        )

    plot_error_distributions(
        error_sds,
        f"{args.out_prefix}_calling_error_distributions.png"
    )

    if args.plot_position_region:
        outputs = plot_called_position_counts_region(
            original=original,
            compressed=compressed,
            compressed_fractions=compressed_fractions,
            error_sds=error_sds,
            n_replicates=args.n_replicates,
            region_start=args.region_start,
            region_length=args.region_length,
            seed=args.seed,
            out_prefix=args.out_prefix,
            compression_mode=args.compression_mode,
            shift_map=shift_map
        )

        print(f"Saved {len(outputs)} called-position region outputs.")

    if args.run_dac:
        run_dac_analysis(
            original=original,
            compressed=compressed,
            compressed_fractions=compressed_fractions,
            error_sds=error_sds,
            n_replicates=args.n_replicates,
            seed=args.seed,
            dmax=args.dac_dmax,
            normalize_dac=not args.no_normalize_dac,
            out_prefix=args.out_prefix,
            compression_mode=args.compression_mode,
            shift_map=shift_map
        )


if __name__ == "__main__":
    main()