#!/usr/bin/env python3
"""
Aggregate aligned score vectors into a single mean profile.

Input:
  --scores : text file where each line is one region and contains space-separated
             numeric values (the aligned window).

Filtering:
  - Skip lines containing NaN/Inf, unless --nan_to_zero is used.
  - If --nan_to_zero is used, NaN values are converted to 0.
  - Skip lines containing >= N consecutive zeros (set --zero_thresh 0 to disable).
  - Skip lines containing any value > --max_score.
  - Skip lines with a different vector length than the first valid line.

Output:
  TSV with columns:
    relative_position   score
  relative_position is centered at 0 (midpoint of the vector).
"""

import argparse
import os
import numpy as np
from tqdm import tqdm


def parse_args():
    p = argparse.ArgumentParser(description="Aggregate aligned score vectors into a mean profile.")
    p.add_argument("--scores", required=True, help="Input score file (one line = one region).")
    p.add_argument("--output", default=None, help="Output TSV (default: <scores_basename>_agg.tsv).")
    p.add_argument(
        "--zero_thresh",
        type=int,
        default=5,
        help="Skip a line if it contains >= this many consecutive zeros (0 disables).",
    )
    p.add_argument(
        "--max_score",
        type=float,
        default=300.0,
        help="Skip a line if any value exceeds this (default: 300).",
    )
    p.add_argument(
        "--nan_to_zero",
        action="store_true",
        help="Convert NaN values to 0 instead of skipping lines containing NaN.",
    )

    args = p.parse_args()

    if args.output is None:
        base, _ = os.path.splitext(args.scores)
        args.output = f"{base}_agg.tsv"

    if args.zero_thresh < 0:
        p.error("--zero_thresh must be >= 0")

    return args


def has_consecutive_zeros(arr: np.ndarray, threshold: int) -> bool:
    if threshold <= 0 or arr.size < threshold:
        return False

    is_zero = (arr == 0.0).astype(np.int8)
    window = np.ones(threshold, dtype=np.int8)

    return np.any(np.convolve(is_zero, window, mode="valid") == threshold)


def load_and_aggregate(
    path: str,
    zero_thresh: int,
    max_score: float,
    nan_to_zero: bool,
) -> np.ndarray:
    running_sum = None
    kept = 0
    skipped = 0
    nan_converted_lines = 0

    with open(path, "r") as f:
        for line_num, line in enumerate(
            tqdm(f, desc="Aggregating", unit="line", dynamic_ncols=True),
            start=1,
        ):
            line = line.strip()

            if not line:
                skipped += 1
                continue

            values = np.fromstring(line, sep=" ", dtype=np.float64)

            if values.size == 0:
                skipped += 1
                continue

            if nan_to_zero:
                if np.any(np.isnan(values)):
                    nan_converted_lines += 1
                values = np.nan_to_num(values, nan=0.0, posinf=np.inf, neginf=-np.inf)

            if not np.all(np.isfinite(values)):
                skipped += 1
                continue

            if has_consecutive_zeros(values, zero_thresh):
                skipped += 1
                continue

            if max_score is not None and np.any(values > max_score):
                skipped += 1
                continue

            if running_sum is None:
                running_sum = np.zeros_like(values, dtype=np.float64)

            elif values.size != running_sum.size:
                skipped += 1
                continue

            running_sum += values
            kept += 1

    if kept == 0:
        raise RuntimeError("No valid lines found after filtering.")

    print(f"[INFO] Kept {kept} lines; skipped {skipped}.")

    if nan_to_zero:
        print(f"[INFO] Converted NaN values to 0 in {nan_converted_lines} lines.")

    return running_sum / kept


def write_aggregated(mean_vec: np.ndarray, out_tsv: str) -> None:
    center = mean_vec.size // 2

    with open(out_tsv, "w") as f:
        f.write("relative_position\tscore\n")

        for i, score in enumerate(mean_vec):
            f.write(f"{i - center}\t{score:.6f}\n")


def main():
    args = parse_args()

    mean_vec = load_and_aggregate(
        args.scores,
        args.zero_thresh,
        args.max_score,
        args.nan_to_zero,
    )

    write_aggregated(mean_vec, args.output)


if __name__ == "__main__":
    main()