#!/usr/bin/env python3
"""
Call nucleosome-protection and breakpoint peaks directly from an input bigWig.

Outputs:
  <out_prefix>_nucleosome_regions.bed
  <out_prefix>_breakpoint_peaks.bed

BED8 columns:
  chrom  start  end  name  score  strand  thickStart  thickEnd

Scoring:
  - nucleosome peaks: maximum signal within each called positive region
  - breakpoint peaks: minimum signal within each called negative region
    (implemented by calling peaks on the inverted signal, then flipping the peak
     score back to original-signal orientation)

Notes:
  - BED score is clamped to 0..1000
  - breakpoint peak raw scores will typically be negative in the original signal;
    for BED score, abs(score) is used before scaling/clamping.
"""

import argparse
import os
import sys
from typing import List, Tuple

import numpy as np
import pyBigWig
from scipy.signal import savgol_filter


def parse_region(region: str) -> Tuple[str, int, int]:
    """
    Parse region like:
      chr1
      chr1:1000-5000
    """
    if ":" not in region:
        return region, -1, -1

    chrom, rest = region.split(":", 1)
    start_s, end_s = rest.replace(",", "").split("-", 1)
    return chrom, int(start_s), int(end_s)


def smooth_signal(arr: np.ndarray, window: int, polyorder: int) -> np.ndarray:
    if len(arr) == 0:
        return arr.copy()

    if window < 3:
        return arr.copy()

    if window % 2 == 0:
        window += 1

    if len(arr) < window:
        return arr.copy()

    if polyorder >= window:
        polyorder = window - 1

    return savgol_filter(arr, window_length=window, polyorder=polyorder)


def get_bw_values(bw: pyBigWig.pyBigWig, chrom: str, start: int, end: int) -> np.ndarray:
    vals = bw.values(chrom, start, end, numpy=True)
    vals = np.array(vals, dtype=float)
    vals[np.isnan(vals)] = 0.0
    vals[np.isinf(vals)] = 0.0
    return vals


def find_peaks_and_regions(
    scores: np.ndarray,
    original_start: int,
    min_length: int = 50,
    max_neg_run: int = 5,
):
    """
    Identify positive regions and their max positions.
    Also identify the most negative point between adjacent positive regions.

    Returns:
      (
        (positive_peaks, positive_peak_scores),
        (negative_peaks, negative_peak_scores),
        (region_centres, positive_peak_regions),
      )
    """
    positive_regions = []
    current_region = None
    searching_for_positive = True
    neg_count = 0
    last_positive_end = None

    for i in range(len(scores)):
        score = scores[i]

        if searching_for_positive:
            if score > 0:
                if current_region is None:
                    current_region = [i, i]
                else:
                    current_region[1] = i
                neg_count = 0
            else:
                if current_region:
                    neg_count += 1
                    if neg_count == 1:
                        last_positive_end = i - 1
                    if neg_count >= max_neg_run:
                        if (
                            last_positive_end is not None
                            and last_positive_end - current_region[0] + 1 >= min_length
                        ):
                            positive_regions.append([current_region[0], last_positive_end])
                        current_region = None
                        searching_for_positive = False
                        neg_count = 0
        else:
            if score <= 0:
                if current_region is None:
                    current_region = [i, i]
                else:
                    current_region[1] = i
            else:
                current_region = [i, i]
                searching_for_positive = True
                neg_count = 0

    if current_region:
        if searching_for_positive and current_region[1] - current_region[0] + 1 >= min_length:
            positive_regions.append(current_region)

    negative_peaks = []
    negative_peak_scores = []
    for i in range(1, len(positive_regions)):
        prev_end = positive_regions[i - 1][1]
        next_start = positive_regions[i][0]
        inter_region_scores = scores[prev_end + 1 : next_start]
        if len(inter_region_scores) > 0:
            most_negative_index = int(np.argmin(inter_region_scores)) + prev_end + 1
            most_negative_score = float(scores[most_negative_index])
            negative_peaks.append(most_negative_index)
            negative_peak_scores.append(most_negative_score)

    positive_peaks = []
    positive_peak_scores = []
    for region in positive_regions:
        region_scores = scores[region[0] : region[1] + 1]
        peak_index = int(np.argmax(region_scores)) + region[0]
        positive_peaks.append(peak_index)
        positive_peak_scores.append(float(scores[peak_index]))

    positive_peak_regions = [
        (region[0] + original_start, region[1] + original_start)
        for region in positive_regions
    ]
    adjusted_positive_peaks = [
        ((region[0] + region[1]) // 2) + original_start
        for region in positive_regions
    ]

    positive_peaks = [p + original_start for p in positive_peaks]
    negative_peaks = [p + original_start for p in negative_peaks]

    return (
        (positive_peaks, positive_peak_scores),
        (negative_peaks, negative_peak_scores),
        (adjusted_positive_peaks, positive_peak_regions),
    )


def iter_peak_records(
    chrom: str,
    original_start: int,
    original_end: int,
    positive_peaks,
    region_centres,
    flip_scores: bool = False,
):
    """
    Output peak records using:
      - max score within region for nucleosome peaks
      - min score within region for breakpoint peaks

    The breakpoint peaks are called on the inverted signal, so when flip_scores=True
    the stored peak score is negated back into the original-signal orientation.
    """
    pos_peak_coords, pos_peak_scores = positive_peaks
    centres, regions = region_centres

    num_positive_peaks = len(centres)

    for i in range(num_positive_peaks):
        region_start = regions[i][0]
        region_end = regions[i][1]
        region_centre = centres[i]
        raw_peak = pos_peak_coords[i]

        if not (original_start <= region_centre < original_end):
            continue

        peak_score_callspace = float(pos_peak_scores[i])

        if flip_scores:
            # breakpoint peak: convert max(inverted signal) back to min(original signal)
            peak_score = -peak_score_callspace
        else:
            # nucleosome peak: max(original signal)
            peak_score = peak_score_callspace

        yield {
            "chrom": chrom,
            "region_start": int(region_start),
            "region_end": int(region_end),
            "region_centre": int(region_centre),
            "raw_peak": int(raw_peak),
            "peak_score": float(peak_score),
            "score": float(peak_score),
        }


def write_bed8(
    path: str,
    records: List[dict],
    label: str,
    score_scale: float = 1.0,
    append: bool = True,
):
    """
    Write BED8.

    BED score must be a non-negative integer. Since breakpoint peak scores are
    negative in the original signal, use abs(score) for the BED score field.
    """
    mode = "a" if append and os.path.exists(path) else "w"
    with open(path, mode) as out:
        for rec in records:
            score = int(round(abs(rec["score"]) * score_scale))

            if score < 0:
                score = 0
            if score > 1000:
                score = 1000

            name = f'{rec["chrom"]}:{rec["region_centre"]}_{label}'
            strand = "."
            thick_start = rec["region_centre"]
            thick_end = rec["region_centre"] + 1

            out.write(
                f'{rec["chrom"]}\t'
                f'{rec["region_start"]}\t'
                f'{rec["region_end"]}\t'
                f'{name}\t'
                f'{score}\t'
                f'{strand}\t'
                f'{thick_start}\t'
                f'{thick_end}\n'
            )


def process_interval(
    bw: pyBigWig.pyBigWig,
    chrom: str,
    start: int,
    end: int,
    smooth_window: int,
    smooth_polyorder: int,
    min_length: int,
    max_neg_run: int,
):
    signal = get_bw_values(bw, chrom, start, end)
    smoothed = smooth_signal(signal, smooth_window, smooth_polyorder)

    nuc_calls = find_peaks_and_regions(
        scores=smoothed,
        original_start=start,
        min_length=min_length,
        max_neg_run=max_neg_run,
    )

    brk_calls = find_peaks_and_regions(
        scores=(-1.0 * smoothed),
        original_start=start,
        min_length=min_length,
        max_neg_run=max_neg_run,
    )

    nuc_records = list(
        iter_peak_records(
            chrom=chrom,
            original_start=start,
            original_end=end,
            positive_peaks=nuc_calls[0],
            region_centres=nuc_calls[2],
            flip_scores=False,
        )
    )

    brk_records = list(
        iter_peak_records(
            chrom=chrom,
            original_start=start,
            original_end=end,
            positive_peaks=brk_calls[0],
            region_centres=brk_calls[2],
            flip_scores=True,
        )
    )

    return nuc_records, brk_records


def chunk_intervals(chrom: str, start: int, end: int, chunk_bp: int):
    cur = start
    while cur < end:
        chunk_end = min(cur + chunk_bp, end)
        yield chrom, cur, chunk_end
        cur = chunk_end


def main():
    parser = argparse.ArgumentParser(
        description="Call nucleosome and breakpoint peaks directly from a bigWig track."
    )
    parser.add_argument("-i", "--input_bw", required=True, help="Input bigWig file")
    parser.add_argument("-o", "--out_prefix", required=True, help="Output prefix")
    parser.add_argument(
        "-r",
        "--regions",
        nargs="*",
        default=None,
        help="Optional regions, e.g. chr1 chr2:100000-200000 chrX:0-500000",
    )
    parser.add_argument(
        "--chunk-bp",
        type=int,
        default=5000000,
        help="Process intervals in chunks of this size",
    )
    parser.add_argument(
        "--smooth-window",
        type=int,
        default=21,
        help="Savitzky-Golay window size (odd number preferred)",
    )
    parser.add_argument(
        "--smooth-polyorder",
        type=int,
        default=2,
        help="Savitzky-Golay polynomial order",
    )
    parser.add_argument(
        "--min-length",
        type=int,
        default=50,
        help="Minimum positive-region length to keep",
    )
    parser.add_argument(
        "--max-neg-run",
        type=int,
        default=5,
        help="Allowed negative run inside a positive region before closing it",
    )
    parser.add_argument(
        "--score-scale",
        type=float,
        default=1.0,
        help="Multiply abs(peak_score) by this before converting to BED score",
    )

    args = parser.parse_args()

    nuc_out = f"{args.out_prefix}_nucleosome_regions.bed"
    brk_out = f"{args.out_prefix}_breakpoint_peaks.bed"

    for path in (nuc_out, brk_out):
        if os.path.exists(path):
            os.remove(path)

    try:
        bw = pyBigWig.open(args.input_bw)
    except Exception as e:
        sys.stderr.write(f"ERROR: could not open bigWig: {e}\n")
        return 1

    chrom_sizes = bw.chroms()

    intervals = []
    if args.regions:
        for region in args.regions:
            chrom, start, end = parse_region(region)
            if chrom not in chrom_sizes:
                sys.stderr.write(f"WARNING: {chrom} not found in bigWig, skipping\n")
                continue

            if start == -1 and end == -1:
                start = 0
                end = chrom_sizes[chrom]

            start = max(0, start)
            end = min(end, chrom_sizes[chrom])

            if end <= start:
                continue

            intervals.extend(chunk_intervals(chrom, start, end, args.chunk_bp))
    else:
        for chrom, length in chrom_sizes.items():
            intervals.extend(chunk_intervals(chrom, 0, length, args.chunk_bp))

    first = True
    for chrom, start, end in intervals:
        nuc_records, brk_records = process_interval(
            bw=bw,
            chrom=chrom,
            start=start,
            end=end,
            smooth_window=args.smooth_window,
            smooth_polyorder=args.smooth_polyorder,
            min_length=args.min_length,
            max_neg_run=args.max_neg_run,
        )

        write_bed8(
            nuc_out,
            nuc_records,
            label="nuc",
            score_scale=args.score_scale,
            append=not first,
        )
        write_bed8(
            brk_out,
            brk_records,
            label="brk",
            score_scale=args.score_scale,
            append=not first,
        )

        first = False

    bw.close()
    return 0


if __name__ == "__main__":
    sys.exit(main())