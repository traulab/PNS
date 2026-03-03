#!/usr/bin/env python3
"""
@author Andrew D Johnston

Kircher-style WPS scoring + median-centering + Kircher-closely-matched peak calling

This version IMPORTS shared helper functions from:
  PNS_with_nucleosome_peak_calling.py  (located one directory down)

Shared helpers imported from PNS
  - resolve_fasta_contig
  - build_dinuc_index
  - uniform_randomize_fragments
  - dinuc_anchor_randomize_fragments
  - require_bam_indexes
  - split_into_regions
  - write_bedgraph
  - write_wig_gz_tracks
  - generate_paired_reads
  - generate_fragment_ranges

Randomization modes (same as PNS):
  - none (default)
  - uniform
  - dinuc_anchor (requires --fasta)
"""

import sys
import argparse
import os
import random
import math
from typing import List, Tuple, Optional
from collections import Counter

import numpy as np
import pysam
from tqdm import tqdm
from numpy.lib.stride_tricks import sliding_window_view
from scipy.signal import savgol_filter


# -------------------------------------------------------------------
# Import PNS script (located one directory down)
# -------------------------------------------------------------------
def _import_pns_helpers():
    """
    Import PNS_with_nucleosome_peak_calling.py by searching RELATIVE TO THIS WPS SCRIPT FILE,
    not the current working directory.

    Search order:
      1) same directory as WPS script
      2) immediate subdirectories of WPS script directory
      3) parent directory of WPS script
      4) immediate subdirectories of that parent directory
      5) recursive (depth-limited) search under WPS dir and its parent (fallback)
    """
    import importlib.util
    import glob

    here = os.path.dirname(os.path.abspath(__file__))
    parent = os.path.dirname(here)

    target_name = "PNS_with_nucleosome_peak_calling.py"

    candidates = []

    # 1) Same directory
    candidates.append(os.path.join(here, target_name))

    # 2) One directory down from WPS script dir (any subdir)
    candidates.extend(glob.glob(os.path.join(here, "*", target_name)))

    # 3) Parent directory
    candidates.append(os.path.join(parent, target_name))

    # 4) Siblings (subdirs of parent)
    candidates.extend(glob.glob(os.path.join(parent, "*", target_name)))

    # 5) Depth-limited recursive fallback (covers weird layouts)
    candidates.extend(glob.glob(os.path.join(here, "**", target_name), recursive=True))
    candidates.extend(glob.glob(os.path.join(parent, "**", target_name), recursive=True))

    # Deduplicate while preserving order
    seen = set()
    candidates = [p for p in candidates if not (p in seen or seen.add(p))]

    for path in candidates:
        if not os.path.isfile(path):
            continue

        spec = importlib.util.spec_from_file_location("PNS_with_nucleosome_peak_calling", path)
        if spec is None or spec.loader is None:
            continue
        mod = importlib.util.module_from_spec(spec)
        try:
            spec.loader.exec_module(mod)
            # Basic sanity check that we got the right module
            if not hasattr(mod, "generate_paired_reads"):
                continue
            print(f"[WPS] Imported PNS helpers from: {path}", file=sys.stderr)
            return mod
        except Exception:
            continue

    raise ImportError(
        "Could not locate/import PNS_with_nucleosome_peak_calling.py relative to the WPS script.\n"
        f"WPS script dir: {here}\n"
        f"Looked for: {target_name}\n"
        "Fix options:\n"
        "  - Ensure PNS_with_nucleosome_peak_calling.py exists somewhere under the WPS dir or its parent\n"
        "  - Or add its directory to PYTHONPATH\n"
    )


_pns = _import_pns_helpers()

# Shared helpers from PNS
resolve_fasta_contig = _pns.resolve_fasta_contig
build_dinuc_index = _pns.build_dinuc_index
uniform_randomize_fragments = _pns.uniform_randomize_fragments
dinuc_anchor_randomize_fragments = _pns.dinuc_anchor_randomize_fragments
require_bam_indexes = _pns.require_bam_indexes
split_into_regions = _pns.split_into_regions
write_bedgraph = _pns.write_bedgraph

# NEW: wig.gz writer (added in updated PNS)
write_wig_gz_tracks = getattr(_pns, "write_wig_gz_tracks", None)

# shared fragment extraction / filtering
generate_paired_reads = _pns.generate_paired_reads
generate_fragment_ranges = _pns.generate_fragment_ranges


# ----------------------------
# Kircher-equivalent WPS kernel (effective bx semantics)
# ----------------------------
def wps_kernel_kircher_exact(L_true: int, protection: int = 120) -> np.ndarray:
    """
    Kernel matching Kircher's effective computation with bx intervals.

    For protection=120 => half=60, effective window length behaves like 120 (2*half),
    and the interval semantics make the effective overlap/span logic equivalent to:

      total_len = L_true + 2*half - 2
      flank     = 2*half - 1
      mid       = L_true - 2*half

    If mid <= 0 => all -1.
    """
    half = protection // 2
    if L_true <= 0:
        return np.array([], dtype=np.int8)

    total_len = L_true + 2 * half - 2
    if total_len <= 0:
        return np.array([], dtype=np.int8)

    flank = 2 * half - 1
    mid = L_true - 2 * half

    if mid <= 0:
        return np.full(total_len, -1, dtype=np.int8)

    k = np.empty(total_len, dtype=np.int8)
    k[:flank] = -1
    k[flank : flank + mid] = +1
    k[flank + mid :] = -1
    return k


def precompute_distributions_kircher_exact(wps_frag_range, protection: int = 120):
    """Precompute kernels for each TRUE fragment length."""
    return {int(L): wps_kernel_kircher_exact(int(L), protection=protection) for L in wps_frag_range}


# ----------------------------
# Rolling median baseline (window=1000)
# ----------------------------
def rolling_median(x: np.ndarray, window: int = 1000) -> np.ndarray:
    """
    Rolling median for even window sizes.
    Places the median at the RIGHT-middle index (half = window//2).
    Returns NaN where the full window doesn't fit.
    """
    x = np.asarray(x, dtype=float)
    n = x.size
    if n < window:
        return np.full(n, np.nan, dtype=float)

    wins = sliding_window_view(x, window_shape=window)
    med = np.median(wins, axis=1)  # length n-window+1

    out = np.full(n, np.nan, dtype=float)
    half = window // 2
    out[half : n - half + 1] = med
    return out


# ----------------------------
# Scoring
# ----------------------------
def score_contig(
    bamfiles,
    contig: str,
    start: int,
    end: int,
    protection: int,
    wps_frag_range: set,
    max_duplicates: int,
    distributions: dict,
    subsample: Optional[float],
    baseline_window: int = 1000,
    sg_window: int = 21,
    sg_order: int = 2,
    # Randomization (same flags/behavior as PNS)
    randomize_mode: str = "none",  # none|uniform|dinuc_anchor
    fasta: Optional[pysam.FastaFile] = None,
    anchor_prob_start: float = 0.5,
    max_anchor_tries: int = 30,
    randomize_fallback: str = "uniform",  # uniform|keep|skip
):
    """
    Compute tracks over [start,end):
      - coverage: fragment overlap
      - dyad: fragment midpoint
      - wps: Kircher-equivalent kernel sum
      - wps_smoothed: Savitzky–Golay on raw wps
      - mWPS: wps - rolling median of raw wps
      - sm_mWPS: wps_smoothed - rolling median of raw wps

    Returns:
      scores_dict, fragments_filtered_post_randomize
    """
    ref_len = end - start
    coverage = np.zeros(ref_len, dtype=np.int32)
    dyad = np.zeros(ref_len, dtype=np.int32)
    wps = np.zeros(ref_len, dtype=np.float64)

    half = protection // 2

    # 1) Collect fragments for this window (shared PNS logic)
    fragments: List[Tuple[int, int]] = []
    for bamfile in bamfiles:
        for frag_start, frag_end in generate_fragment_ranges(
            bamfile, contig, start, end, max_duplicates, subsample
        ):
            if frag_end > frag_start:
                fragments.append((frag_start, frag_end))

    # 2) Randomize fragments (shared functions from PNS)
    if randomize_mode == "uniform" and fragments:
        fragments = uniform_randomize_fragments(fragments, start, end)

    elif randomize_mode == "dinuc_anchor" and fragments:
        if fasta is None:
            raise ValueError("randomize_mode=dinuc_anchor requires --fasta")
        fasta_contig = resolve_fasta_contig(fasta, contig)
        window_seq = fasta.fetch(fasta_contig, start, end).upper()
        if len(window_seq) != ref_len:
            window_seq = window_seq[:ref_len].ljust(ref_len, "N")
        dinuc_pos = build_dinuc_index(window_seq)

        fragments = dinuc_anchor_randomize_fragments(
            fragments=fragments,
            start=start,
            end=end,
            window_seq=window_seq,
            dinuc_pos=dinuc_pos,
            anchor_prob_start=anchor_prob_start,
            max_anchor_tries=max_anchor_tries,
            fallback=randomize_fallback,
        )

    # 3) Score fragments (possibly randomized)
    for frag_start, frag_end in fragments:
        L_true = frag_end - frag_start
        if L_true <= 0:
            continue

        # Coverage (true overlap)
        cov_s = max(frag_start, start)
        cov_e = min(frag_end, end)
        if cov_e > cov_s:
            coverage[cov_s - start : cov_e - start] += 1

        # Dyad (true midpoint)
        frag_center = frag_start + (L_true - 1) // 2
        if start <= frag_center < end:
            dyad[frag_center - start] += 1

        # WPS kernel sum
        if L_true not in wps_frag_range:
            continue

        kernel = distributions.get(L_true, None)
        if kernel is None or kernel.size == 0:
            continue

        # Effective placement matching Kircher bx overlap domain
        kernel_start_genome = frag_start - half + 1

        k_s = kernel_start_genome - start
        k_e = k_s + kernel.size

        arr_s = max(k_s, 0)
        arr_e = min(k_e, ref_len)
        if arr_e <= arr_s:
            continue

        ker_s = arr_s - k_s
        ker_e = ker_s + (arr_e - arr_s)

        wps[arr_s:arr_e] += kernel[ker_s:ker_e].astype(np.float64)

    # Smooth raw WPS
    if ref_len >= sg_window and sg_window % 2 == 1:
        wps_smoothed = savgol_filter(wps, sg_window, sg_order)
    else:
        wps_smoothed = wps.copy()

    # Baseline = rolling median of RAW WPS
    baseline = rolling_median(wps, window=baseline_window)
    if np.isnan(baseline).any():
        valid = np.where(~np.isnan(baseline))[0]
        if valid.size > 0:
            first, last = valid[0], valid[-1]
            baseline[:first] = baseline[first]
            baseline[last + 1 :] = baseline[last]
        else:
            baseline[:] = 0.0

    mWPS = wps - baseline
    sm_mWPS = wps_smoothed - baseline

    scores = {
        "coverage": [(contig, start, coverage)],
        "dyad": [(contig, start, dyad)],
        "wps": [(contig, start, wps)],
        "wps_smoothed": [(contig, start, wps_smoothed)],
        "mWPS": [(contig, start, mWPS)],
        "sm_mWPS": [(contig, start, sm_mWPS)],
    }

    return scores, fragments


# ----------------------------
# Kircher-matched peak calling
# ----------------------------
def kircher_median(values):
    """For even length, returns the UPPER middle element (Kircher behavior)."""
    helper = list(values)
    helper.sort()
    lVal = len(helper)
    point = 0.5
    if round(lVal * point) == int(lVal * point):
        return float(helper[int(lVal * point)])
    else:
        return float((helper[int(round(lVal * point))] + helper[int(lVal * point)]) * 0.5)


def kircher_continous_windows(region_pairs):
    """
    Faithful to Kircher continousWindows:
    when a gap occurs, it appends the current window and resets,
    but does NOT start a new window with the current point.
    """
    res = []
    cstart, cend = None, None
    cmax = None
    csum = 0.0
    for pos, val in region_pairs:
        if cmax is None:
            cend = pos
            cstart = pos
            cmax = val
            csum = val
        else:
            if cend + 1 == pos:
                cend = pos
                csum += val
                if cmax < val:
                    cmax = val
            else:
                res.append((csum, cstart, cend, cmax))
                cstart, cend = None, None
                cmax = None
                csum = 0.0
    if cmax is not None:
        res.append((csum, cstart, cend, cmax))
    return res


def round_half_up(x):
    return int(math.floor(float(x) + 0.5))


def evaluate_values_kircher(
    chrom_nochr,
    start1,
    end1,
    values,
    report=True,
    minlength=50,
    maxlength=150,
    vari_cutoff=5.0,
):
    if start1 is None or end1 is None or not values:
        return []

    L = len(values)
    calls = []

    if maxlength >= L >= minlength:
        cMed = kircher_median(values)
        region_pairs = [(pos, val) for pos, val in zip(range(start1, end1 + 1), values) if val >= cMed]
        windows = kircher_continous_windows(region_pairs)
        if not windows:
            return []
        windows.sort()
        _score_sum, cstart, cend, cval = windows[-1]
        if report and (cval > vari_cutoff):
            cmiddle = cstart + round_half_up((cend - cstart) * 0.5)
            chrom_out = f"chr{chrom_nochr}"
            bed_start = cstart - 1
            bed_end = cend
            name = f"{chrom_nochr}:{cstart}-{cend}"
            score_int = int(round(float(cval)))
            thick_start = cmiddle - 1
            thick_end = cmiddle
            calls.append((chrom_out, bed_start, bed_end, name, score_int, ".", thick_start, thick_end))
        return calls

    elif (3 * maxlength) >= L >= maxlength:
        cMed = kircher_median(values)
        region_pairs = [(pos, val) for pos, val in zip(range(start1, end1 + 1), values) if val >= cMed]
        windows = kircher_continous_windows(region_pairs)
        for _score_sum, cstart, cend, cval in windows:
            seg_len = (cend - cstart + 1)
            if maxlength >= seg_len >= minlength:
                if report and (cval > vari_cutoff):
                    cmiddle = cstart + round_half_up((cend - cstart) * 0.5)
                    chrom_out = f"chr{chrom_nochr}"
                    bed_start = cstart - 1
                    bed_end = cend
                    name = f"{chrom_nochr}:{cstart}-{cend}"
                    score_int = int(round(float(cval)))
                    thick_start = cmiddle - 1
                    thick_end = cmiddle
                    calls.append((chrom_out, bed_start, bed_end, name, score_int, ".", thick_start, thick_end))
        return calls

    return []


def call_peaks_kircher_matched(
    contig,
    adjusted_start0,
    track,
    merge_gap_bp=5,
    minlength=50,
    maxlength=150,
    vari_cutoff=5.0,
):
    chrom_nochr = contig.replace("chr", "")
    allowed = {str(i) for i in range(1, 23)} | {"X", "Y"}
    report = chrom_nochr in allowed

    calls = []
    cstart = None
    cend = None
    clist = []

    n = int(track.size)
    for i in range(n):
        ivalue = float(track[i])
        ipos = int(adjusted_start0 + i + 1)  # 1-based

        if ivalue > 0:
            if (cend is not None) and (ipos <= (cend + merge_gap_bp)):
                while (cend + 1) < ipos:
                    cend += 1
                    clist.append(0.0)
                clist.append(ivalue)
                cend = ipos
            else:
                if cstart is not None and clist and report:
                    calls.extend(
                        evaluate_values_kircher(
                            chrom_nochr,
                            cstart,
                            cend,
                            clist,
                            report=report,
                            minlength=minlength,
                            maxlength=maxlength,
                            vari_cutoff=vari_cutoff,
                        )
                    )
                clist = [ivalue]
                cstart = ipos
                cend = ipos

    if cstart is not None and clist and report:
        calls.extend(
            evaluate_values_kircher(
                chrom_nochr,
                cstart,
                cend,
                clist,
                report=report,
                minlength=minlength,
                maxlength=maxlength,
                vari_cutoff=vari_cutoff,
            )
        )

    return calls


# ----------------------------
# Output writers
# ----------------------------
def write_bed_rows(rows, path, mode):
    with open(path, mode) as f:
        for chrom, start, end, name, score, strand, thick_start, thick_end in rows:
            f.write(
                f"{chrom}\t{int(start)}\t{int(end)}\t{name}\t{int(score)}\t{strand}\t{int(thick_start)}\t{int(thick_end)}\n"
            )


def write_fragment_outputs(
    out_prefix: str,
    total_fragments_filtered_all: int,
    total_fragments_used_in_range: int,
    unique_bases_covered_by_used: int,
    length_counts: Counter,
):
    summary_path = f"{out_prefix}_fragment_summary.txt"
    lens_path = f"{out_prefix}_fragment_length_counts.tsv"

    with open(summary_path, "w") as f:
        f.write(f"total_fragments_filtered_all\t{total_fragments_filtered_all}\n")
        f.write(f"total_fragments_used_in_range\t{total_fragments_used_in_range}\n")
        f.write(f"unique_bases_covered_by_used_fragments\t{unique_bases_covered_by_used}\n")

    with open(lens_path, "w") as f:
        f.write("fragment_length\tcount\n")
        for L in sorted(length_counts.keys()):
            f.write(f"{int(L)}\t{int(length_counts[L])}\n")


# ----------------------------
# Main
# ----------------------------
def main():
    parser = argparse.ArgumentParser(
        description="Kircher-style WPS scoring + median-centering + Kircher-matched peak calling."
    )
    parser.add_argument("-b", "--bamfiles", nargs="+", required=True, help="BAM file(s) to process")
    parser.add_argument("-o", "--out_prefix", default=None, help="Output prefix")
    parser.add_argument(
        "-c",
        "--contigs",
        nargs="+",
        default=None,
        help='Limit to contig(s) and optional range, e.g. "12:51730340-52039340" or "12"',
    )
    parser.add_argument("--protection", type=int, default=120, help="Protection window (bp), default 120.")
    parser.add_argument("--frag-lower", type=int, default=127, help="Lower fragment size to include in WPS")
    parser.add_argument("--frag-upper", type=int, default=207, help="Upper fragment size to include in WPS")
    parser.add_argument("--max-duplicates", type=int, default=0, help="Maximum allowed duplicate fragments (same coords)")
    parser.add_argument("--subsample", type=float, default=None, help="Subsampling proportion (e.g. 0.5 keeps ~50%)")
    parser.add_argument("--chunk-bp", type=int, default=100000, help="Chunk size per contig")
    parser.add_argument("--overlap-bp", type=int, default=1000, help="Overlap padding for edge-safe scoring")
    parser.add_argument("--baseline-window", type=int, default=1000, help="Rolling median window for baseline subtraction")
    parser.add_argument("--sg-window", type=int, default=21, help="Savitzky-Golay window (odd)")
    parser.add_argument("--sg-order", type=int, default=2, help="Savitzky-Golay polynomial order")

    # Score-track output controls (NEW)
    parser.add_argument(
        "--score-format",
        choices=["bedgraph", "wiggz", "both", "none"],
        default="bedgraph",
        help="How to write per-base score tracks. 'wiggz' writes one <prefix>_<track>.wig.gz per track.",
    )
    parser.add_argument(
        "--score-tracks",
        nargs="*",
        default=["coverage", "sm_mWPS", "wps", "wps_smoothed", "mWPS", "dyad"],
        help=(
            "Which score tracks to output (space-separated). "
            "Valid: coverage dyad wps wps_smoothed mWPS sm_mWPS. "
            "Use '--score-tracks none' or '--score-format none' to disable."
        ),
    )

    # Peak calling controls
    parser.add_argument("--peak-minlen", type=int, default=50, help="Minimum length (bp) for candidate windows")
    parser.add_argument("--peak-maxlen", type=int, default=150, help="Maximum length (bp) for reported windows")
    parser.add_argument("--peak-maxregion", type=int, default=450, help="Reject merged regions longer than this (bp)")
    parser.add_argument("--peak-merge-gap", type=int, default=5, help="Merge positive runs if gap <= this many bp")
    parser.add_argument("--peak-varicutoff", type=float, default=5.0, help="Minimum max score to report a peak window")

    # Randomization controls (same as PNS)
    parser.add_argument("--seed", type=int, default=None, help="Random seed (for reproducibility)")
    parser.add_argument(
        "--randomize-mode",
        choices=["none", "uniform", "dinuc_anchor"],
        default="none",
        help="Fragment randomization mode within each processed window.",
    )
    parser.add_argument(
        "--fasta",
        default=None,
        help="Reference FASTA (required for --randomize-mode dinuc_anchor). Needs .fai index.",
    )
    parser.add_argument(
        "--anchor-prob-start",
        type=float,
        default=0.5,
        help="For dinuc_anchor: probability of anchoring on fragment START (otherwise anchor on END).",
    )
    parser.add_argument(
        "--max-anchor-tries",
        type=int,
        default=30,
        help="For dinuc_anchor: max attempts to find a valid placement per fragment.",
    )
    parser.add_argument(
        "--randomize-fallback",
        choices=["uniform", "keep", "skip"],
        default="uniform",
        help="For dinuc_anchor: what to do if no valid dinuc placement is found.",
    )

    args = parser.parse_args()

    # Normalize --score-tracks none + validate
    valid_tracks = {"coverage", "dyad", "wps", "wps_smoothed", "mWPS", "sm_mWPS"}
    if args.score_tracks and len(args.score_tracks) == 1 and args.score_tracks[0].lower() == "none":
        args.score_tracks = []
    else:
        bad = [t for t in args.score_tracks if t not in valid_tracks]
        if bad:
            parser.error(f"Unknown --score-tracks: {bad}. Valid: {sorted(valid_tracks)}")

    if args.score_format in ("wiggz", "both") and not callable(write_wig_gz_tracks):
        parser.error(
            "This WPS script was asked to write --score-format wiggz/both, but the imported PNS module "
            "does not provide write_wig_gz_tracks(). Please update/repoint PNS_with_nucleosome_peak_calling.py."
        )

    if args.seed is not None:
        random.seed(args.seed)
        np.random.seed(args.seed)

    require_bam_indexes(args.bamfiles, parser=parser)

    if args.randomize_mode == "dinuc_anchor" and not args.fasta:
        parser.error("--randomize-mode dinuc_anchor requires --fasta <ref.fa> (with .fai index).")

    bamfiles = []
    for p in args.bamfiles:
        try:
            bamfiles.append(pysam.AlignmentFile(p, "rb"))
        except Exception as e:
            parser.error(f"Unable to open BAM {p}: {e}")
            return 2

    fasta = None
    if args.fasta:
        try:
            fasta = pysam.FastaFile(args.fasta)
        except Exception as e:
            parser.error(f"Unable to open FASTA '{args.fasta}': {str(e)}")
            return 2

    if not args.out_prefix:
        bnames = [os.path.splitext(os.path.basename(b))[0] for b in args.bamfiles]
        args.out_prefix = "_".join(bnames)
        if args.contigs and len(args.contigs) == 1:
            args.out_prefix = f"{args.out_prefix}_{args.contigs[0].replace(':','_')}"

    args.out_prefix = (
        f"{args.out_prefix}_prot{args.protection}"
        f"_lower{args.frag_lower}_upper{args.frag_upper}"
        f"_maxdup{args.max_duplicates}"
    )
    if args.randomize_mode != "none":
        args.out_prefix += f"_rand{args.randomize_mode}"

    contigs = []
    if args.contigs:
        for contig_range in args.contigs:
            if ":" in contig_range:
                contig, positions = contig_range.split(":")
                start, end = map(int, positions.split("-"))
                contig_len = bamfiles[0].get_reference_length(contig)
                contigs.extend(
                    split_into_regions(
                        contig,
                        start,
                        end,
                        contig_len,
                        max_length=args.chunk_bp,
                        overlap=args.overlap_bp,
                    )
                )
            else:
                contig = contig_range
                start, end = 0, bamfiles[0].get_reference_length(contig)
                contig_len = bamfiles[0].get_reference_length(contig)
                contigs.extend(
                    split_into_regions(
                        contig,
                        start,
                        end,
                        contig_len,
                        max_length=args.chunk_bp,
                        overlap=args.overlap_bp,
                    )
                )
    else:
        for contig in bamfiles[0].references:
            start, end = 0, bamfiles[0].get_reference_length(contig)
            contig_len = bamfiles[0].get_reference_length(contig)
            contigs.extend(
                split_into_regions(
                    contig,
                    start,
                    end,
                    contig_len,
                    max_length=args.chunk_bp,
                    overlap=args.overlap_bp,
                )
            )

    wps_frag_range = set(range(args.frag_lower, args.frag_upper + 1))
    distributions = precompute_distributions_kircher_exact(wps_frag_range, protection=args.protection)

    combined_bedgraph = f"{args.out_prefix}_combined_scores.bedGraph"
    nuc_bed = f"{args.out_prefix}_nucleosome_regions.bed"
    brk_bed = f"{args.out_prefix}_breakpoint_peaks.bed"
    frag_summary = f"{args.out_prefix}_fragment_summary.txt"
    frag_lens = f"{args.out_prefix}_fragment_length_counts.tsv"

    wig_paths = [f"{args.out_prefix}_{t}.wig.gz" for t in args.score_tracks]

    # Remove old outputs (respecting format choices)
    to_remove = [nuc_bed, brk_bed, frag_summary, frag_lens]
    if args.score_format in ("bedgraph", "both"):
        to_remove.append(combined_bedgraph)
    if args.score_format in ("wiggz", "both"):
        to_remove.extend(wig_paths)

    for fn in to_remove:
        if os.path.exists(fn):
            os.remove(fn)

    total_fragments_filtered_all = 0
    total_fragments_used_in_range = 0
    unique_bases_covered_by_used = 0
    length_counts = Counter()

    first_region = True
    for contig, adjusted_start, adjusted_end, original_start, original_end in tqdm(contigs, desc="Scoring contigs"):
        scores, fragments_filtered = score_contig(
            bamfiles=bamfiles,
            contig=contig,
            start=adjusted_start,
            end=adjusted_end,
            protection=args.protection,
            wps_frag_range=wps_frag_range,
            max_duplicates=args.max_duplicates,
            distributions=distributions,
            subsample=args.subsample,
            baseline_window=args.baseline_window,
            sg_window=args.sg_window,
            sg_order=args.sg_order,
            randomize_mode=args.randomize_mode,
            fasta=fasta,
            anchor_prob_start=args.anchor_prob_start,
            max_anchor_tries=args.max_anchor_tries,
            randomize_fallback=args.randomize_fallback,
        )

        owned_fragments = [(fs, fe) for (fs, fe) in fragments_filtered if (original_start <= fs < original_end)]
        total_fragments_filtered_all += len(owned_fragments)

        if original_end > original_start:
            covered = np.zeros(original_end - original_start, dtype=bool)

            for fs, fe in owned_fragments:
                L = fe - fs
                if L not in wps_frag_range:
                    continue

                total_fragments_used_in_range += 1
                length_counts[L] += 1

                ov_s = max(fs, original_start)
                ov_e = min(fe, original_end)
                if ov_e > ov_s:
                    covered[ov_s - original_start : ov_e - original_start] = True

            unique_bases_covered_by_used += int(covered.sum())

        track = scores["sm_mWPS"][0][2]

        nuc_rows = call_peaks_kircher_matched(
            contig=contig,
            adjusted_start0=adjusted_start,
            track=track,
            merge_gap_bp=args.peak_merge_gap,
            minlength=args.peak_minlen,
            maxlength=args.peak_maxlen,
            vari_cutoff=args.peak_varicutoff,
        )
        brk_rows = call_peaks_kircher_matched(
            contig=contig,
            adjusted_start0=adjusted_start,
            track=(-1.0 * track),
            merge_gap_bp=args.peak_merge_gap,
            minlength=args.peak_minlen,
            maxlength=args.peak_maxlen,
            vari_cutoff=args.peak_varicutoff,
        )

        if args.peak_maxregion != 3 * args.peak_maxlen:
            nuc_rows = [r for r in nuc_rows if (r[2] - r[1]) <= args.peak_maxregion]
            brk_rows = [r for r in brk_rows if (r[2] - r[1]) <= args.peak_maxregion]

        def keep_core(rows):
            out = []
            for chrom, s, e, name, score, strand, ts, te in rows:
                if e <= original_start or s >= original_end:
                    continue
                out.append((chrom, s, e, name, score, strand, ts, te))
            return out

        nuc_rows = keep_core(nuc_rows)
        brk_rows = keep_core(brk_rows)

        write_bed_rows(nuc_rows, nuc_bed, mode=("w" if first_region else "a"))
        write_bed_rows(brk_rows, brk_bed, mode=("w" if first_region else "a"))

        # Per-base score tracks
        if args.score_format in ("bedgraph", "both"):
            write_bedgraph(scores, [(original_start, original_end)], args.out_prefix, first_region)

        if args.score_format in ("wiggz", "both") and args.score_tracks:
            # Use the shared PNS writer (one file per track, fixedStep, gzipped)
            write_wig_gz_tracks(
                scores=scores,
                contig=contig,
                adjusted_start=adjusted_start,
                original_start=original_start,
                original_end=original_end,
                out_prefix=args.out_prefix,
                tracks=args.score_tracks,
                first_region=first_region,
            )

        first_region = False

    write_fragment_outputs(
        out_prefix=args.out_prefix,
        total_fragments_filtered_all=total_fragments_filtered_all,
        total_fragments_used_in_range=total_fragments_used_in_range,
        unique_bases_covered_by_used=unique_bases_covered_by_used,
        length_counts=length_counts,
    )

    for b in bamfiles:
        try:
            b.close()
        except Exception:
            pass

    if fasta is not None:
        try:
            fasta.close()
        except Exception:
            pass

    return 0


if __name__ == "__main__":
    sys.exit(main())
