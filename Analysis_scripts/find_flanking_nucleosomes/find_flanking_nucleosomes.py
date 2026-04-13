#!/usr/bin/env python3
"""
@author Andrew D Johnston

Given one or more TF binding-site BED files, this script finds the +/- Nth nucleosome
peak relative to each TF site and outputs per-offset files.

Key ideas
- Nucleosome peaks are loaded once (BED or BED.gz), grouped by chromosome and sorted by
  peak center for fast binary search (bisect).
- For each TF site, we compute its reference position ("TF center") and choose an orientation ("strand"):
    * TF center defaults to midpoint between BED start/end, but can be overridden with --tf_pos_col.
    * Strand can be forced for all sites with --force_strand (+, -, or random).
    * Otherwise strand uses a column from the TF BED if present (configurable; BED6 default col 6).
    * If missing/invalid, strand is inferred by looking at the nearest nucleosome center.
- "+N" and "-N" are defined relative to the TF strand:
    * On '+' strand: plus = downstream (increasing genomic coordinate),
                   minus = upstream (decreasing coordinate).
    * On '-' strand: plus = upstream (decreasing coordinate),
                   minus = downstream (increasing coordinate).

Outputs
For each TF BED and each offset k in [1..N], write four files:
    <tf_prefix>_plus<k>.tsv
    <tf_prefix>_minus<k>.tsv
    <tf_prefix>_plus<k>.bed
    <tf_prefix>_minus<k>.bed

TSV output columns:
    tf_index  tf_chrom  tf_start  tf_end  tf_center  strand  nuc_start  nuc_end  nuc_score

BED output columns:
    chrom  nuc_start  nuc_end  strand  nuc_score

Behavior
- One output row is written for every valid TF input row, in the same order as input.
- If the requested side runs off the chromosome's nucleosome list, the lookup is clamped
  to the first or last available nucleosome instead of returning NA.
- Only if a chromosome has no nucleosome entries at all do the nucleosome fields become NA in TSV.
- BED output is forced to remain BED-safe:
    * if no nucleosome exists for that chromosome, bed coords become tf_center and tf_center+1
    * score becomes 0 if unavailable
  This prevents downstream BED parsers from silently dropping rows.

Debugging / validation
- A skipped-line log is written per TF BED:
    <tf_prefix>_skipped_lines.tsv
- A summary file is written per TF BED:
    <tf_prefix>_summary.txt
- The script prints per-buffer row counts and asserts that every output buffer has exactly
  the same number of rows as the number of valid TF rows.
"""

import argparse
from bisect import bisect_left
import gzip
import glob
import os
import random
import sys
from typing import Optional

from tqdm import tqdm


# ----------------------------
# CLI
# ----------------------------
def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description="Find +N and -N nucleosome peaks for each TF binding site."
    )
    p.add_argument(
        "tf_beds",
        nargs="+",
        help="TF BED file(s) or glob pattern(s), e.g. *_hg19.bed or data/*.bed.gz",
    )
    p.add_argument(
        "-n",
        "--nuc_bed",
        required=True,
        help="Nucleosome peak BED file (.bed or .bed.gz).",
    )
    p.add_argument(
        "--num",
        type=int,
        default=100,
        help="Number of nucleosome peaks to find in each direction (default: 100).",
    )
    p.add_argument(
        "--tf_strand_col",
        type=int,
        default=6,
        help=(
            "1-based column in TF BED that contains strand (+/-). "
            "BED6 strand is column 6 (default: 6). "
            "If missing/invalid, strand is inferred unless --force_strand is used."
        ),
    )
    p.add_argument(
        "--force_strand",
        choices=["+", "-", "random"],
        default=None,
        help=(
            "Force all TF sites to strand '+' or '-', or assign each site a random strand. "
            "If not set, use TF BED strand or infer strand from nucleosomes."
        ),
    )
    p.add_argument(
        "--tf_pos_col",
        type=int,
        default=None,
        help=(
            "Optional 1-based column containing a single-base TF position "
            "(e.g. motif center). If set, TF center is taken from this column. "
            "Default: use midpoint of BED start/end."
        ),
    )
    p.add_argument(
        "--nuc_score_col",
        type=int,
        default=4,
        help=(
            "1-based column in nucleosome BED that contains peak score "
            "(default: 4). If missing, score is 'NA'."
        ),
    )
    p.add_argument(
        "--out_dir",
        default=".",
        help="Directory to write output files (default: current directory).",
    )
    p.add_argument(
        "--seed",
        type=int,
        default=None,
        help="Optional random seed for reproducible --force_strand random mode.",
    )
    return p.parse_args()


# ----------------------------
# I/O helpers
# ----------------------------
def smart_open(path: str):
    """Open a text file; use gzip if filename ends with .gz."""
    return gzip.open(path, "rt") if path.endswith(".gz") else open(path, "rt")


def expand_globs(patterns: list[str]) -> list[str]:
    """
    Expand glob patterns.

    - If a pattern matches files, include those.
    - If it doesn't match but is an existing path, include it as-is.
    - Otherwise warn.
    """
    out: list[str] = []
    for pat in patterns:
        matches = glob.glob(pat)
        if matches:
            out.extend(matches)
        else:
            if os.path.exists(pat):
                out.append(pat)
            else:
                print(f"[WARN] No files matched pattern: {pat}", file=sys.stderr)
    return out


# ----------------------------
# Nucleosome loading
# ----------------------------
def load_nucleosomes(nuc_bed_file: str, score_col_1based: int) -> dict:
    """
    Load nucleosome peaks into a dict keyed by chromosome.

    Returns:
      nuc_dict[chrom] = {
          'centers': [center0, center1, ...]  (sorted)
          'records': [(center, start, end, score), ...] (same order)
      }
    """
    score_idx = score_col_1based - 1
    tmp: dict[str, list[tuple[int, int, int, str]]] = {}

    with smart_open(nuc_bed_file) as f:
        for line in f:
            stripped = line.strip()
            if not stripped or stripped.startswith("#"):
                continue

            parts = stripped.split()
            if len(parts) < 3:
                continue

            try:
                chrom = parts[0]
                start = int(parts[1])
                end = int(parts[2])
            except ValueError:
                continue

            score = parts[score_idx] if len(parts) > score_idx else "NA"
            center = (start + end) // 2

            tmp.setdefault(chrom, []).append((center, start, end, score))

    out: dict[str, dict[str, list]] = {}
    for chrom, recs in tmp.items():
        recs.sort(key=lambda x: x[0])
        out[chrom] = {"centers": [r[0] for r in recs], "records": recs}

    return out


# ----------------------------
# Core lookup logic
# ----------------------------
def infer_strand(nuc_data: Optional[dict], tf_center: int) -> str:
    """
    Infer a strand for a TF site based on the nearest nucleosome center.

    Rule: return '+' if the nearest nucleosome center is left of tf_center, else '-'.

    If nuc_data is missing or empty, default to '+'.
    """
    if nuc_data is None or not nuc_data["centers"]:
        return "+"

    centers = nuc_data["centers"]
    idx = bisect_left(centers, tf_center)

    if idx == 0:
        nearest = centers[0]
    elif idx == len(centers):
        nearest = centers[-1]
    else:
        before = centers[idx - 1]
        after = centers[idx]
        nearest = before if abs(before - tf_center) <= abs(after - tf_center) else after

    return "+" if nearest < tf_center else "-"


def find_n_peak(
    nuc_data: Optional[dict],
    tf_center: int,
    strand: str,
    offset: int,
    direction: str,
):
    """
    Find the +N or -N nucleosome peak relative to tf_center and strand.

    direction: 'plus' or 'minus' (relative to TF strand)
    offset: 1-based index (1 = closest in that direction by center order)

    If the requested side runs out of peaks, clamp to the first/last nucleosome
    instead of returning None.

    Returns:
      (center, start, end, score) or None if chromosome has no nucleosomes
    """
    if nuc_data is None or not nuc_data["records"]:
        return None

    centers = nuc_data["centers"]
    recs = nuc_data["records"]

    idx = bisect_left(centers, tf_center)

    if strand == "+":
        target_idx = idx + (offset - 1) if direction == "plus" else idx - offset
    elif strand == "-":
        target_idx = idx - offset if direction == "plus" else idx + (offset - 1)
    else:
        return None

    if target_idx < 0:
        target_idx = 0
    elif target_idx >= len(recs):
        target_idx = len(recs) - 1

    return recs[target_idx]


def tf_strand_from_parts(parts: list[str], strand_col_1based: int) -> Optional[str]:
    """Return '+'/'-' if present and valid in the configured column, else None."""
    idx = strand_col_1based - 1
    if len(parts) > idx and parts[idx] in ("+", "-"):
        return parts[idx]
    return None


def tf_center_from_parts(parts: list[str], pos_col_1based: Optional[int]) -> Optional[int]:
    """
    Return TF reference position.

    - If pos_col_1based is provided, read TF position from that column (must be int).
    - Otherwise return None so caller can use BED midpoint.
    """
    if pos_col_1based is None:
        return None
    idx = pos_col_1based - 1
    if len(parts) <= idx:
        return None
    try:
        return int(parts[idx])
    except ValueError:
        return None


def assign_strand(
    parts: list[str],
    nuc_data: Optional[dict],
    tf_center: int,
    tf_strand_col_1based: int,
    force_strand: Optional[str],
) -> str:
    """
    Determine strand assignment priority:
      1. Forced strand if provided (+ or -)
      2. Random strand if force_strand == 'random'
      3. TF BED strand if valid
      4. Infer from nearest nucleosome
    """
    if force_strand is not None:
        if force_strand in ("+", "-"):
            return force_strand
        if force_strand == "random":
            return random.choice(["+", "-"])

    strand = tf_strand_from_parts(parts, tf_strand_col_1based)
    if strand is not None:
        return strand

    return infer_strand(nuc_data, tf_center)


def make_bed_safe_fields(
    chrom: str,
    tf_center: int,
    strand: str,
    peak,
):
    """
    Return BED-safe fields:
      chrom, nuc_start, nuc_end, strand, score

    If peak is None:
      - start = tf_center
      - end = tf_center + 1
      - score = 0

    This keeps every BED line valid so downstream tools do not silently drop rows.
    """
    if peak is not None:
        _, nuc_start, nuc_end, nuc_score = peak
        score_str = str(nuc_score)
        if score_str == "NA":
            score_str = "0"
        return chrom, str(nuc_start), str(nuc_end), strand, score_str

    return chrom, str(tf_center), str(tf_center + 1), strand, "0"


# ----------------------------
# Processing
# ----------------------------
def process_tf_bed(
    tf_bed: str,
    nuc_dict: dict,
    num: int,
    tf_strand_col_1based: int,
    force_strand: Optional[str],
    tf_pos_col_1based: Optional[int],
    out_dir: str,
):
    """
    Process one TF BED, writing per-offset TSV and BED files to out_dir.

    This version always writes one line per valid TF row, in the same order as input.
    If no chromosome-level nucleosome data exists, TSV nucleosome columns are written as NA.
    BED output is still forced to valid integer coordinates.
    Invalid input lines are logged to a skipped-lines file.
    """
    tsv_buffers = {f"plus{k}": [] for k in range(1, num + 1)}
    tsv_buffers.update({f"minus{k}": [] for k in range(1, num + 1)})

    bed_buffers = {f"plus{k}": [] for k in range(1, num + 1)}
    bed_buffers.update({f"minus{k}": [] for k in range(1, num + 1)})

    total_lines = 0
    written_rows = 0
    skipped_blank = 0
    skipped_comment = 0
    skipped_short = 0
    skipped_bad_coords = 0
    skipped_lines_log: list[str] = []

    base = os.path.basename(tf_bed)
    prefix = base.replace(".bed.gz", "").replace(".bed", "").replace(".gz", "")
    os.makedirs(out_dir, exist_ok=True)

    with smart_open(tf_bed) as f:
        for line_no, line in enumerate(
            tqdm(f, desc=os.path.basename(tf_bed), unit="line", dynamic_ncols=True),
            start=1,
        ):
            total_lines += 1
            raw = line.rstrip("\n")
            stripped = line.strip()

            if not stripped:
                skipped_blank += 1
                skipped_lines_log.append(f"{line_no}\tblank\t{raw}")
                continue

            if stripped.startswith("#"):
                skipped_comment += 1
                skipped_lines_log.append(f"{line_no}\tcomment\t{raw}")
                continue

            parts = stripped.split()
            if len(parts) < 3:
                skipped_short += 1
                skipped_lines_log.append(f"{line_no}\tlt_3_columns\t{raw}")
                continue

            chrom = parts[0]
            try:
                start = int(parts[1])
                end = int(parts[2])
            except ValueError:
                skipped_bad_coords += 1
                skipped_lines_log.append(f"{line_no}\tbad_start_end\t{raw}")
                continue

            written_rows += 1

            tf_center = tf_center_from_parts(parts, tf_pos_col_1based)
            if tf_center is None:
                tf_center = (start + end) // 2

            nuc_data = nuc_dict.get(chrom)

            strand = assign_strand(
                parts=parts,
                nuc_data=nuc_data,
                tf_center=tf_center,
                tf_strand_col_1based=tf_strand_col_1based,
                force_strand=force_strand,
            )

            for k in range(1, num + 1):
                for direction in ("plus", "minus"):
                    peak = find_n_peak(nuc_data, tf_center, strand, k, direction)

                    if peak is not None:
                        _, nuc_start, nuc_end, nuc_score = peak
                        nuc_start_str = str(nuc_start)
                        nuc_end_str = str(nuc_end)
                        nuc_score_str = str(nuc_score)
                    else:
                        nuc_start_str = "NA"
                        nuc_end_str = "NA"
                        nuc_score_str = "NA"

                    key = f"{direction}{k}"

                    tsv_buffers[key].append(
                        f"{line_no}\t{chrom}\t{start}\t{end}\t{tf_center}\t{strand}\t"
                        f"{nuc_start_str}\t{nuc_end_str}\t{nuc_score_str}"
                    )

                    bed_chrom, bed_start, bed_end, bed_strand, bed_score = make_bed_safe_fields(
                        chrom=chrom,
                        tf_center=tf_center,
                        strand=strand,
                        peak=peak,
                    )
                    bed_buffers[key].append(
                        f"{bed_chrom}\t{bed_start}\t{bed_end}\t{bed_strand}\t{bed_score}"
                    )

    # Assert all buffers have exactly the expected number of rows
    expected = written_rows
    buffer_counts: list[tuple[str, int, int]] = []

    for key, lines in tsv_buffers.items():
        count = len(lines)
        buffer_counts.append((f"{key}.tsv", count, expected))
        if count != expected:
            raise RuntimeError(
                f"[ERROR] Buffer {key}.tsv has {count} rows but expected {expected}."
            )

    for key, lines in bed_buffers.items():
        count = len(lines)
        buffer_counts.append((f"{key}.bed", count, expected))
        if count != expected:
            raise RuntimeError(
                f"[ERROR] Buffer {key}.bed has {count} rows but expected {expected}."
            )

    tsv_header = (
        "tf_index\ttf_chrom\ttf_start\ttf_end\ttf_center\tstrand\t"
        "nuc_start\tnuc_end\tnuc_score\n"
    )

    for key, lines in tsv_buffers.items():
        out_path = os.path.join(out_dir, f"{prefix}_{key}.tsv")
        with open(out_path, "w") as out_f:
            out_f.write(tsv_header)
            if lines:
                out_f.write("\n".join(lines) + "\n")

    for key, lines in bed_buffers.items():
        out_path = os.path.join(out_dir, f"{prefix}_{key}.bed")
        with open(out_path, "w") as out_f:
            if lines:
                out_f.write("\n".join(lines) + "\n")

    skipped_path = os.path.join(out_dir, f"{prefix}_skipped_lines.tsv")
    with open(skipped_path, "w") as skip_f:
        skip_f.write("line_number\treason\traw_line\n")
        if skipped_lines_log:
            skip_f.write("\n".join(skipped_lines_log) + "\n")

    summary_path = os.path.join(out_dir, f"{prefix}_summary.txt")
    with open(summary_path, "w") as summary_f:
        summary_f.write(f"TF BED: {tf_bed}\n")
        summary_f.write(f"Total physical lines read: {total_lines}\n")
        summary_f.write(f"Valid TF rows written per output file: {written_rows}\n")
        summary_f.write(f"Skipped blank lines: {skipped_blank}\n")
        summary_f.write(f"Skipped comment lines: {skipped_comment}\n")
        summary_f.write(f"Skipped lines with <3 columns: {skipped_short}\n")
        summary_f.write(f"Skipped lines with non-integer start/end: {skipped_bad_coords}\n")
        summary_f.write(f"Skipped-line log: {skipped_path}\n")
        summary_f.write("\nPer-buffer row counts:\n")
        for name, count, exp in buffer_counts:
            summary_f.write(f"{name}\t{count}\t(expected {exp})\n")

    print(f"[INFO] Finished: {tf_bed}")
    print(f"[INFO] Total physical lines read: {total_lines}")
    print(f"[INFO] Valid TF rows written per output file: {written_rows}")
    print(f"[INFO] Skipped blank lines: {skipped_blank}")
    print(f"[INFO] Skipped comment lines: {skipped_comment}")
    print(f"[INFO] Skipped lines with <3 columns: {skipped_short}")
    print(f"[INFO] Skipped lines with non-integer start/end: {skipped_bad_coords}")
    print(f"[INFO] Skipped-line log: {skipped_path}")
    print(f"[INFO] Summary: {summary_path}")
    print("[INFO] Per-buffer row counts:")
    for name, count, exp in buffer_counts:
        print(f"[INFO]   {name}: {count} (expected {exp})")


def main():
    args = parse_args()

    if args.seed is not None:
        random.seed(args.seed)

    if args.num <= 0:
        print("[ERROR] --num must be > 0", file=sys.stderr)
        sys.exit(1)

    tf_beds = expand_globs(args.tf_beds)
    if not tf_beds:
        print("[ERROR] No TF BED files found.", file=sys.stderr)
        sys.exit(1)

    print(f"[INFO] Loading nucleosomes: {args.nuc_bed}")
    nuc_dict = load_nucleosomes(args.nuc_bed, score_col_1based=args.nuc_score_col)

    for tf_bed in tf_beds:
        print(f"[INFO] Processing TF BED: {tf_bed}")
        process_tf_bed(
            tf_bed=tf_bed,
            nuc_dict=nuc_dict,
            num=args.num,
            tf_strand_col_1based=args.tf_strand_col,
            force_strand=args.force_strand,
            tf_pos_col_1based=args.tf_pos_col,
            out_dir=args.out_dir,
        )


if __name__ == "__main__":
    main()