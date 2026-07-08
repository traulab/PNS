#!/usr/bin/env python3
"""
@author Andrew D Johnston

Compute adjacent-peak distance distributions (50–1000 bp) optionally stratified
by chromatin state, while supporting score-thresholding by percentile or by a
target number of peaks.

Key performance idea:
- Load peaks (and assign chromatin-state labels via interval trees),
  then iterate thresholds without re-reading BED files.

Input peak file (BED-like):
- Requires at least: chrom, start, end
- Uses either:
    * midpoint of (start,end) as peak position (default), OR
    * a user-specified column as the position
- Uses a user-specified column as the score for filtering

Chromatin state file (optional; BED-like):
- Requires at least: chrom, start, end, state_label
- State at a peak is defined by interval overlap at the peak position.
"""

import argparse
import os
import sys
from collections import defaultdict

import numpy as np
from intervaltree import IntervalTree
from scipy.signal import savgol_filter


# ----------------------------
# Helpers: weighted stats from counts (no expansion)
# ----------------------------
def weighted_mean_from_counts(distance_dict: dict[int, int]):
    """Weighted mean of distance keys using counts as weights."""
    if not distance_dict:
        return None

    total = 0
    s = 0
    for d, c in distance_dict.items():
        if c > 0:
            c = int(c)
            total += c
            s += int(d) * c

    return (s / total) if total > 0 else None


def weighted_median_from_counts(distance_dict: dict[int, int]):
    """Weighted median of distance keys using counts as weights."""
    if not distance_dict:
        return None

    items = [(d, c) for d, c in distance_dict.items() if c and c > 0]
    if not items:
        return None

    items.sort(key=lambda x: x[0])
    total = sum(int(c) for _, c in items)
    if total <= 0:
        return None

    target = (total + 1) // 2  # 1-indexed “middle”
    cum = 0
    for d, c in items:
        cum += int(c)
        if cum >= target:
            return d

    return items[-1][0]


# ----------------------------
# Smoothing + mode
# ----------------------------
def smooth_and_find_mode(distance_dict, window_length=21, polyorder=2):
    """
    Smooth counts by distance and return:
      - mode distance from SMOOTHED counts
      - weighted median and weighted mean from UNSMOOTHED counts
      - smoothed_count dict (rounded, non-negative ints)

    Notes:
    - Smoothing is applied to the *counts* sequence ordered by distance.
    - If smoothing parameters are invalid for the number of points, falls back to raw counts.
    """
    if not distance_dict:
        return None, None, None, defaultdict(int)

    # Sort distances so counts become an ordered 1D signal for smoothing
    distances = sorted(distance_dict.keys())
    counts = np.array([distance_dict[d] for d in distances], dtype=float)

    # Summary stats computed on raw counts (more interpretable as distribution properties)
    median_dist = weighted_median_from_counts(distance_dict)
    mean_dist = weighted_mean_from_counts(distance_dict)

    # Savitzky-Golay smoothing requires:
    # - odd window_length
    # - window_length < number of points
    # - polyorder < window_length
    if len(distances) > window_length and window_length % 2 == 1 and polyorder < window_length:
        try:
            smoothed_counts = savgol_filter(counts, window_length=window_length, polyorder=polyorder)
        except Exception:
            smoothed_counts = counts
    else:
        smoothed_counts = counts

    # Mode is the distance at which smoothed counts are maximal
    mode_value = distances[int(np.argmax(smoothed_counts))] if len(smoothed_counts) else None

    # Convert smoothed counts back into an int-valued dict (no negatives)
    smoothed_dict = defaultdict(int)
    for d, sc in zip(distances, smoothed_counts):
        smoothed_dict[d] = max(0, int(round(float(sc))))

    return mode_value, median_dist, mean_dist, smoothed_dict


def smooth_normalized_values(normalized_values, window_length=35, polyorder=3):
    """
    Smooth normalized (%-scaled) values for nicer output curves.
    Falls back to raw values if parameters are invalid.
    """
    if len(normalized_values) > window_length and window_length % 2 == 1 and polyorder < window_length:
        try:
            return savgol_filter(
                np.asarray(normalized_values, dtype=float),
                window_length=window_length,
                polyorder=polyorder,
            )
        except Exception:
            return normalized_values
    return normalized_values


# ----------------------------
# Percentile sweep helper
# ----------------------------
def frange_inclusive(start: float, stop: float, step: float) -> list[float]:
    """
    Floating-point range generator that includes the stop endpoint (within tolerance).
    Produces de-duplicated values after rounding noise.
    """
    if step <= 0:
        raise ValueError("--pct-step must be > 0")

    n = int(np.floor((stop - start) / step + 1e-12)) + 1
    vals = [start + i * step for i in range(n)]

    # Ensure stop is included even if step doesn't land exactly on it
    if vals and (vals[-1] < stop) and abs(vals[-1] - stop) > 1e-9:
        vals.append(stop)

    # Clamp to [start, stop]
    vals = [min(stop, max(start, v)) for v in vals]

    # De-duplicate nearly identical values
    out = []
    last = None
    for v in vals:
        if last is None or abs(v - last) > 1e-9:
            out.append(v)
        last = v
    return out


def pct_label(p: float) -> str:
    """Compact percentile label for filenames (e.g. 5, 5.5, 5.25)."""
    return f"{p:.2f}".rstrip("0").rstrip(".")


# ----------------------------
# Build chromatin interval trees
# ----------------------------
def build_chromatin_state_trees(chromatin_state_file: str):
    """
    Build an IntervalTree per chromosome from a chromatin-state BED-like file.

    Expected columns (>=4):
      chrom  start  end  state_label
    """
    trees = defaultdict(IntervalTree)
    with open(chromatin_state_file, "r") as csfile:
        for line in csfile:
            s = line.strip()
            if not s or s.startswith(("#", "track", "browser")):
                continue

            cols = s.split()
            if len(cols) < 4:
                continue

            chrom = cols[0]
            try:
                start = int(cols[1])
                end = int(cols[2])
            except ValueError:
                continue

            state = cols[3]
            trees[chrom].addi(start, end, state)

    return trees


# ----------------------------
# Load peaks (and annotate state)
# ----------------------------
def load_peaks(
    input_file: str,
    chromatin_state_trees,
    combine_euchromatin: bool = False,
    position_column: int | None = None,
    score_column: int = 3,
    store_line: bool = False,
):
    """
    One-pass loader that produces:
      - peaks_by_chrom: dict[chrom] -> list of (pos, score, state_label[, original_line])
      - scores: numpy array of scores (for percentile thresholding)
      - n_lines: total non-header/non-empty lines scanned
      - n_used: lines successfully parsed into peaks

    State label assignment:
      - If chromatin_state_trees is None, uses "All"
      - Else, overlaps at the peak position define states_here
      - If multiple states overlap, joins with '|'
      - If no overlaps, uses "NA"
      - Optional: collapse any overlap in selected prefixes to "Euchromatin"
    """
    euchromatin_prefixes = ("1_", "2_", "4_", "5_", "6_", "7_", "9_", "10_", "11_")

    peaks_by_chrom = defaultdict(list)
    scores = []
    n_lines = 0
    n_used = 0

    with open(input_file, "r") as f:
        for line in f:
            s = line.strip()
            if not s or s.startswith(("#", "track", "browser")):
                continue

            n_lines += 1
            cols = s.split()

            # Require enough columns to parse:
            # - chrom/start/end (cols[0..2]) always
            # - score_column
            # - optionally position_column
            min_len = 3
            if position_column is not None:
                min_len = max(min_len, position_column + 1)
            min_len = max(min_len, score_column + 1)
            if len(cols) < min_len:
                continue

            chrom = cols[0]

            # Parse score
            try:
                score = float(cols[score_column])
            except ValueError:
                continue

            # Parse position: either explicit column, or midpoint(start,end)
            try:
                if position_column is None:
                    pos = (int(cols[2]) + int(cols[1])) // 2
                else:
                    pos = int(cols[position_column])
            except (ValueError, IndexError):
                continue

            # Assign chromatin state label at this position
            if chromatin_state_trees is None:
                state_label = "All"
            else:
                overlaps = chromatin_state_trees[chrom][pos]
                states_here = {ov.data for ov in overlaps} if overlaps else {"NA"}

                # Optional: collapse “euchromatin-ish” states into a single label
                if combine_euchromatin and any(
                    str(st).startswith(pref) for st in states_here for pref in euchromatin_prefixes
                ):
                    state_label = "Euchromatin"
                else:
                    if len(states_here) == 1:
                        state_label = next(iter(states_here))
                    else:
                        state_label = "|".join(sorted(str(x) for x in states_here))

            n_used += 1
            scores.append(score)

            if store_line:
                peaks_by_chrom[chrom].append((pos, score, state_label, line))
            else:
                peaks_by_chrom[chrom].append((pos, score, state_label))

    # Sort peaks per chromosome by position so adjacency is meaningful
    for chrom in peaks_by_chrom:
        peaks_by_chrom[chrom].sort(key=lambda x: x[0])

    return peaks_by_chrom, np.asarray(scores, dtype=float), n_lines, n_used


# ----------------------------
# Compute distances for a given threshold
# ----------------------------
def compute_distances_for_threshold(peaks_by_chrom, threshold: float, min_distance: int, max_distance: int, max_order: int = 10):
    """
    For each chromosome, keep peaks with score >= threshold, then count distances
    between peak i and peak i+n for n = 1..max_order.

    Meaning of order:
      +1 = adjacent kept peaks
      +2 = next-nearest kept peaks, i.e. one nucleosome skipped
      +3 = third neighbour, etc.

    Distance counting rule:
    - Distances are only counted when both peaks have the SAME state label.
    - If no chromatin-state file was supplied, all peaks have state label "All".
    - Distances are only recorded in [min_distance, max_distance] bp.

    Outputs are dictionaries keyed by order:
      chrom_distances_by_order[order][chrom][state][distance] += 1
      chrom_all_distances_by_order[order][chrom][distance] += 1
      genome_distances_by_order[order][state][distance] += 1
      genome_all_distances_by_order[order][distance] += 1
    """
    if min_distance < 0 or max_distance < 0:
        raise ValueError("--min-distance and --max-distance must be >= 0")
    if max_distance < min_distance:
        raise ValueError("--max-distance must be >= --min-distance")
    if max_order < 1:
        raise ValueError("--max-order must be >= 1")

    chrom_distances_by_order = defaultdict(lambda: defaultdict(lambda: defaultdict(lambda: defaultdict(int))))
    chrom_all_distances_by_order = defaultdict(lambda: defaultdict(lambda: defaultdict(int)))
    genome_distances_by_order = defaultdict(lambda: defaultdict(lambda: defaultdict(int)))
    genome_all_distances_by_order = defaultdict(lambda: defaultdict(int))
    duplicate_positions = defaultdict(lambda: defaultdict(list))

    kept_count = 0

    for chrom, peaks in peaks_by_chrom.items():
        kept_peaks = []
        last_pos = None

        for rec in peaks:
            # Rec is either (pos, score, state) or (pos, score, state, line)
            pos, score, state = rec[0], rec[1], rec[2]

            if score < threshold:
                continue

            kept_count += 1

            # Track duplicates among kept peaks on the same chromosome
            if last_pos is not None and last_pos == pos:
                duplicate_positions[chrom][pos].append(state)
            last_pos = pos

            kept_peaks.append((pos, state))

        # Count +1, +2, +3 ... +max_order distances within this chromosome.
        # Because peaks were sorted by position during loading, pos2 >= pos1.
        n = len(kept_peaks)
        for i in range(n):
            pos1, state1 = kept_peaks[i]

            for order in range(1, max_order + 1):
                j = i + order
                if j >= n:
                    break

                pos2, state2 = kept_peaks[j]

                # Same logic as the original adjacent-only script: only compare matching states.
                if state1 != state2:
                    continue

                if pos2 == pos1:
                    continue

                d = pos2 - pos1
                if min_distance <= d <= max_distance:
                    chrom_distances_by_order[order][chrom][state1][d] += 1
                    chrom_all_distances_by_order[order][chrom][d] += 1
                    genome_distances_by_order[order][state1][d] += 1
                    genome_all_distances_by_order[order][d] += 1

    return (
        chrom_distances_by_order,
        chrom_all_distances_by_order,
        genome_distances_by_order,
        genome_all_distances_by_order,
        duplicate_positions,
        kept_count,
    )


def write_filtered_bed_from_loaded(peaks_by_chrom, threshold: float, out_path: str):
    """
    Write a filtered BED-like file containing only lines with score >= threshold.

    Requires load_peaks(..., store_line=True) so each record includes the original line.
    """
    with open(out_path, "w") as out:
        for chrom in sorted(peaks_by_chrom.keys()):
            for rec in peaks_by_chrom[chrom]:
                # Must have (pos, score, state, line)
                if len(rec) < 4:
                    raise ValueError(
                        "write_filtered_bed_from_loaded requires store_line=True when calling load_peaks()."
                    )
                score = rec[1]
                line = rec[3]
                if score >= threshold:
                    out.write(line if line.endswith("\n") else line + "\n")


# ----------------------------
# Threshold selection (percentile or target_peaks)
# ----------------------------
def choose_threshold(scores: np.ndarray, requested_percentile: float, target_peaks: int | None, n_used: int):
    """
    Returns (threshold, eff_percentile).

    If target_peaks is provided:
      - Choose an effective percentile so that approximately target_peaks are kept:
            eff_percentile = 100 * (1 - target_peaks / n_used)
      - Then threshold = percentile(scores, eff_percentile)

    Else:
      - threshold = percentile(scores, requested_percentile)
    """
    if scores.size == 0 or n_used <= 0:
        return float("-inf"), 0.0

    if target_peaks is not None and target_peaks > 0:
        k = min(int(target_peaks), int(n_used))
        frac_keep = float(k) / float(n_used)
        eff_percentile = max(0.0, min(100.0, 100.0 * (1.0 - frac_keep)))
        threshold = float(np.percentile(scores, eff_percentile))
        return threshold, eff_percentile

    p = max(0.0, min(100.0, float(requested_percentile)))
    threshold = float(np.percentile(scores, p))
    return threshold, p


# ----------------------------
# Write the standard output txt
# ----------------------------
def write_output_txt(
    output_file: str,
    score_column: int,
    requested_percentile: float,
    target_peaks: int | None,
    eff_percentile: float,
    threshold: float,
    n_used: int,
    n_lines: int,
    min_distance: int,
    max_distance: int,
    order: int,
    chrom_distances,
    chrom_all_distances,
    genome_distances,
    genome_all_distances,
    duplicate_positions,
):
    """
    Output format matches the previous layout you were using, including:
      - per-chrom per-state distance table
      - per-chrom all-state table
      - whole-genome per-state table
      - whole-genome all-state table
      - duplicate position report
    """
    with open(output_file, "w") as outfile:
        outfile.write(
            f"# Score filtering: score_column={score_column + 1}, "
            f"requested_percentile={requested_percentile:.2f}, "
            f"target_peaks={target_peaks}, effective_percentile={eff_percentile:.2f}, "
            f"threshold={threshold:.6g}, scored_lines={n_used}, total_lines_scanned={n_lines}\n"
        )
        outfile.write(f"# Distance filter: min_distance={min_distance}, max_distance={max_distance}, order=+{order}\n")

        # Chromosome Distances by Chromatin State
        outfile.write(f"Chromosome Distances by Chromatin State (+{order}):\n")
        for chrom in chrom_distances:
            outfile.write(f"\n{chrom}\n")
            for state in chrom_distances[chrom]:
                mode_dist, median_dist, mean_dist, smoothed_counts = smooth_and_find_mode(
                    chrom_distances[chrom][state]
                )
                if mode_dist is None:
                    continue

                total_count = sum(chrom_distances[chrom][state].values())
                if total_count <= 0:
                    continue

                dists_sorted = sorted(chrom_distances[chrom][state].keys())
                normalized_values = [(chrom_distances[chrom][state][d] / total_count) * 100.0 for d in dists_sorted]
                smoothed_normalized_values = smooth_normalized_values(normalized_values)

                outfile.write(
                    f"\nState: {state} (Mode Distance: {mode_dist}, "
                    f"Median Distance: {median_dist}, Mean Distance: {mean_dist})\n"
                )
                for i, d in enumerate(dists_sorted):
                    sc = smoothed_counts.get(d, 0)
                    c = chrom_distances[chrom][state][d]
                    outfile.write(
                        f"{chrom}\t{state}\t{d}\t{c}\t{sc}\t"
                        f"{normalized_values[i]:.3f}\t{float(smoothed_normalized_values[i]):.3f}\n"
                    )

        # Chromosome Distances (All States)
        outfile.write(f"\nChromosome Distances (All States, +{order}):\n")
        for chrom in chrom_all_distances:
            mode_dist, median_dist, mean_dist, smoothed_counts = smooth_and_find_mode(chrom_all_distances[chrom])
            if mode_dist is None:
                continue

            total_count = sum(chrom_all_distances[chrom].values())
            if total_count <= 0:
                continue

            dists_sorted = sorted(chrom_all_distances[chrom].keys())
            normalized_values = [(chrom_all_distances[chrom][d] / total_count) * 100.0 for d in dists_sorted]
            smoothed_normalized_values = smooth_normalized_values(normalized_values)

            outfile.write(
                f"\n{chrom} (Mode Distance: {mode_dist}, "
                f"Median Distance: {median_dist}, Mean Distance: {mean_dist})\n"
            )
            for i, d in enumerate(dists_sorted):
                sc = smoothed_counts.get(d, 0)
                c = chrom_all_distances[chrom][d]
                outfile.write(
                    f"{chrom}\tAll\t{d}\t{c}\t{sc}\t"
                    f"{normalized_values[i]:.3f}\t{float(smoothed_normalized_values[i]):.3f}\n"
                )

        # Whole Genome Distances by Chromatin State
        outfile.write(f"\nWhole Genome Distances by Chromatin State (+{order}):\n")
        for state in genome_distances:
            mode_dist, median_dist, mean_dist, smoothed_counts = smooth_and_find_mode(genome_distances[state])
            if mode_dist is None:
                continue

            total_count = sum(genome_distances[state].values())
            if total_count <= 0:
                continue

            dists_sorted = sorted(genome_distances[state].keys())
            normalized_values = [(genome_distances[state][d] / total_count) * 100.0 for d in dists_sorted]
            smoothed_normalized_values = smooth_normalized_values(normalized_values)

            outfile.write(
                f"\nState: {state} (Mode Distance: {mode_dist}, "
                f"Median Distance: {median_dist}, Mean Distance: {mean_dist})\n"
            )
            for i, d in enumerate(dists_sorted):
                sc = smoothed_counts.get(d, 0)
                c = genome_distances[state][d]
                outfile.write(
                    f"whole-genome\t{state}\t{d}\t{c}\t{sc}\t"
                    f"{normalized_values[i]:.3f}\t{float(smoothed_normalized_values[i]):.3f}\n"
                )

        # Whole Genome (All States)
        mode_dist, median_dist, mean_dist, smoothed_counts = smooth_and_find_mode(genome_all_distances)
        if mode_dist is not None:
            total_count = sum(genome_all_distances.values())
            if total_count > 0:
                dists_sorted = sorted(genome_all_distances.keys())
                normalized_values = [(genome_all_distances[d] / total_count) * 100.0 for d in dists_sorted]
                smoothed_normalized_values = smooth_normalized_values(normalized_values)

                outfile.write(
                    f"\nWhole Genome (All States) (Mode Distance: {mode_dist}, "
                    f"Median Distance: {median_dist}, Mean Distance: {mean_dist})\n"
                )
                for i, d in enumerate(dists_sorted):
                    sc = smoothed_counts.get(d, 0)
                    c = genome_all_distances[d]
                    outfile.write(
                        f"whole-genome\tAll\t{d}\t{c}\t{sc}\t"
                        f"{normalized_values[i]:.3f}\t{float(smoothed_normalized_values[i]):.3f}\n"
                    )

        # Duplicate Positions
        outfile.write("\nDuplicate Positions with Chromatin States:\n")
        for chrom in duplicate_positions:
            for position in duplicate_positions[chrom]:
                states = ", ".join(duplicate_positions[chrom][position])
                outfile.write(f"{chrom}\t{position}\t{states}\n")


# ----------------------------
# CLI / main
# ----------------------------
if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description=(
            "Loads peaks + chromatin-state labels, then can sweep percentiles "
            "without re-reading either BED file per percentile."
        )
    )
    parser.add_argument("input", help="Path to the input BED-like file.")
    parser.add_argument(
        "chromatin_state",
        nargs="?",
        default=None,
        help="Path to the chromatin state BED file (optional). If omitted, state is 'All'.",
    )
    parser.add_argument(
        "--combine-euchromatin",
        action="store_true",
        help=(
            "Collapse any overlap that includes 1_,2_,4_,5_,6_,7_,9_,10_,11_ to 'Euchromatin'. "
            "Ignored if no chromatin state file."
        ),
    )
    parser.add_argument(
        "--score-percentile",
        type=float,
        default=0.0,
        help="Single percentile (0-100). Ignored if --target-peaks is set or if --pct-range is used.",
    )
    parser.add_argument(
        "--target-peaks",
        type=int,
        default=None,
        help=(
            "If set, choose a score threshold so that approximately this many scored lines are kept "
            "(based on the chosen score column). Overrides --score-percentile and --pct-range."
        ),
    )
    parser.add_argument(
        "--write-filtered-bed",
        action="store_true",
        help="Write filtered BED(s) of kept peaks. In sweep mode, writes one per percentile.",
    )
    parser.add_argument(
        "--position-column",
        type=int,
        default=None,
        help="Column number to use directly as the position. If not set, uses midpoint of start/end.",
    )
    parser.add_argument(
        "--score-column",
        type=int,
        default=4,
        help="Column number to use as the score for filtering. Default 4.",
    )

    # NEW: distance bounds
    parser.add_argument(
        "--min-distance",
        type=int,
        default=1,
        help="Minimum adjacent-peak distance (bp) to include. Default 50.",
    )
    parser.add_argument(
        "--max-distance",
        type=int,
        default=10000,
        help="Maximum adjacent-peak distance (bp) to include. Default 10000.",
    )
    parser.add_argument(
        "--max-order",
        type=int,
        default=10,
        help="Maximum neighbour order to calculate. Default 10 gives separate +1 to +10 files.",
    )

    # Percentile sweep options
    parser.add_argument(
        "--pct-range",
        action="store_true",
        help="If set, run a sweep of percentiles from --pct-lower to --pct-upper (inclusive) with step --pct-step.",
    )
    parser.add_argument("--pct-lower", type=float, default=0.0, help="Lower percentile for --pct-range (default 0).")
    parser.add_argument("--pct-upper", type=float, default=99.0, help="Upper percentile for --pct-range (default 99).")
    parser.add_argument("--pct-step", type=float, default=1.0, help="Step size for --pct-range (default 1).")

    args = parser.parse_args()

    # Convert user-facing 1-based column numbers to internal 0-based indices
    if args.position_column is not None:
        if args.position_column <= 0:
            print("[ERROR] --position-column must be >= 1", file=sys.stderr)
            sys.exit(2)
        args.position_column -= 1

    if args.score_column is None or args.score_column <= 0:
        print("[ERROR] --score-column must be >= 1", file=sys.stderr)
        sys.exit(2)
    args.score_column -= 1

    # Basic validation for distance bounds
    if args.min_distance < 0 or args.max_distance < 0:
        print("[ERROR] --min-distance and --max-distance must be >= 0", file=sys.stderr)
        sys.exit(2)
    if args.max_distance < args.min_distance:
        print("[ERROR] --max-distance must be >= --min-distance", file=sys.stderr)
        sys.exit(2)
    if args.max_order < 1:
        print("[ERROR] --max-order must be >= 1", file=sys.stderr)
        sys.exit(2)

    # Build chromatin trees (optional)
    if args.chromatin_state is not None:
        print(f"[INFO] Loading chromatin states: {args.chromatin_state}", file=sys.stderr)
        chrom_trees = build_chromatin_state_trees(args.chromatin_state)
        state_base = os.path.splitext(os.path.basename(args.chromatin_state))[0]
    else:
        chrom_trees = None
        state_base = "whole_genome"

    # Load peaks
    store_line = bool(args.write_filtered_bed)
    print(f"[INFO] Loading peaks: {args.input}", file=sys.stderr)
    peaks_by_chrom, scores, n_lines, n_used = load_peaks(
        args.input,
        chrom_trees,
        combine_euchromatin=args.combine_euchromatin,
        position_column=args.position_column,
        score_column=args.score_column,
        store_line=store_line,
    )

    if scores.size == 0:
        print("[ERROR] No usable numeric scores found; nothing to do.", file=sys.stderr)
        sys.exit(1)

    input_base = os.path.splitext(os.path.basename(args.input))[0]
    in_base_path, in_ext = os.path.splitext(args.input)
    if not in_ext:
        in_ext = ".bed"

    # ----------------------------
    # Target peaks overrides everything
    # ----------------------------
    if args.target_peaks is not None and args.target_peaks > 0:
        thr, effp = choose_threshold(scores, args.score_percentile, args.target_peaks, n_used)

        print(
            f"[INFO] Target-peaks mode: target_peaks={args.target_peaks} "
            f"n_used={n_used} -> effective_percentile={effp:.2f} threshold={thr:.6g}",
            file=sys.stderr,
        )

        out_txt = f"{input_base}_{state_base}_scorepct{pct_label(effp)}_nuc_dis.txt"

        (
            chrom_distances_by_order,
            chrom_all_distances_by_order,
            genome_distances_by_order,
            genome_all_distances_by_order,
            duplicate_positions,
            kept_count,
        ) = compute_distances_for_threshold(peaks_by_chrom, thr, args.min_distance, args.max_distance, args.max_order)

        print(f"[INFO] threshold={thr:.6g} kept={kept_count}/{n_used}", file=sys.stderr)

        if args.write_filtered_bed:
            out_bed = f"{in_base_path}_scorepct{pct_label(effp)}{in_ext}"
            write_filtered_bed_from_loaded(peaks_by_chrom, thr, out_bed)
            print(f"[INFO] Filtered BED written: {out_bed}", file=sys.stderr)

        for order in range(1, args.max_order + 1):
            out_txt = f"{input_base}_{state_base}_scorepct{pct_label(effp)}_plus{order}_nuc_dis.txt"
            write_output_txt(
                out_txt,
                score_column=args.score_column,
                requested_percentile=float(args.score_percentile),
                target_peaks=args.target_peaks,
                eff_percentile=effp,
                threshold=thr,
                n_used=n_used,
                n_lines=n_lines,
                min_distance=args.min_distance,
                max_distance=args.max_distance,
                order=order,
                chrom_distances=chrom_distances_by_order[order],
                chrom_all_distances=chrom_all_distances_by_order[order],
                genome_distances=genome_distances_by_order[order],
                genome_all_distances=genome_all_distances_by_order[order],
                duplicate_positions=duplicate_positions,
            )
        sys.exit(0)

    # ----------------------------
    # Percentile sweep mode
    # ----------------------------
    if args.pct_range:
        lo = min(100.0, max(0.0, args.pct_lower))
        hi = min(100.0, max(0.0, args.pct_upper))
        if hi < lo:
            lo, hi = hi, lo
        step = float(args.pct_step)

        pcts = frange_inclusive(lo, hi, step)

        # Compute all thresholds in one percentile call (fast; no re-reading files)
        thresholds = np.percentile(scores, pcts)

        print(f"[INFO] Sweep: {len(pcts)} percentiles, thresholds computed.", file=sys.stderr)

        for p, thr in zip(pcts, thresholds):
            p = float(p)
            thr = float(thr)

            out_txt = f"{input_base}_{state_base}_scorepct{pct_label(p)}_nuc_dis.txt"

            (
                chrom_distances,
                chrom_all_distances,
                genome_distances,
                genome_all_distances,
                duplicate_positions,
                kept_count,
            ) = compute_distances_for_threshold(peaks_by_chrom, thr, args.min_distance, args.max_distance, args.max_order)

            print(f"[INFO] pct={p:.2f} threshold={thr:.6g} kept={kept_count}/{n_used}", file=sys.stderr)

            if args.write_filtered_bed:
                out_bed = f"{in_base_path}_scorepct{pct_label(p)}{in_ext}"
                write_filtered_bed_from_loaded(peaks_by_chrom, thr, out_bed)
                print(f"[INFO] Filtered BED written: {out_bed}", file=sys.stderr)

            for order in range(1, args.max_order + 1):
                out_txt = f"{input_base}_{state_base}_scorepct{pct_label(p)}_plus{order}_nuc_dis.txt"
                write_output_txt(
                    out_txt,
                    score_column=args.score_column,
                    requested_percentile=p,
                    target_peaks=None,
                    eff_percentile=p,
                    threshold=thr,
                    n_used=n_used,
                    n_lines=n_lines,
                    min_distance=args.min_distance,
                    max_distance=args.max_distance,
                    order=order,
                    chrom_distances=chrom_distances_by_order[order],
                    chrom_all_distances=chrom_all_distances_by_order[order],
                    genome_distances=genome_distances_by_order[order],
                    genome_all_distances=genome_all_distances_by_order[order],
                    duplicate_positions=duplicate_positions,
                )

    # ----------------------------
    # Single percentile mode
    # ----------------------------
    else:
        thr, effp = choose_threshold(scores, args.score_percentile, None, n_used)

        out_txt = f"{input_base}_{state_base}_scorepct{pct_label(effp)}_nuc_dis.txt"

        (
            chrom_distances_by_order,
            chrom_all_distances_by_order,
            genome_distances_by_order,
            genome_all_distances_by_order,
            duplicate_positions,
            kept_count,
        ) = compute_distances_for_threshold(peaks_by_chrom, thr, args.min_distance, args.max_distance, args.max_order)

        print(f"[INFO] pct={effp:.2f} threshold={thr:.6g} kept={kept_count}/{n_used}", file=sys.stderr)

        if args.write_filtered_bed:
            out_bed = f"{in_base_path}_scorepct{pct_label(effp)}{in_ext}"
            write_filtered_bed_from_loaded(peaks_by_chrom, thr, out_bed)
            print(f"[INFO] Filtered BED written: {out_bed}", file=sys.stderr)

        for order in range(1, args.max_order + 1):
            out_txt = f"{input_base}_{state_base}_scorepct{pct_label(effp)}_plus{order}_nuc_dis.txt"
            write_output_txt(
                out_txt,
                score_column=args.score_column,
                requested_percentile=effp,
                target_peaks=None,
                eff_percentile=effp,
                threshold=thr,
                n_used=n_used,
                n_lines=n_lines,
                min_distance=args.min_distance,
                max_distance=args.max_distance,
                order=order,
                chrom_distances=chrom_distances_by_order[order],
                chrom_all_distances=chrom_all_distances_by_order[order],
                genome_distances=genome_distances_by_order[order],
                genome_all_distances=genome_all_distances_by_order[order],
                duplicate_positions=duplicate_positions,
            )
