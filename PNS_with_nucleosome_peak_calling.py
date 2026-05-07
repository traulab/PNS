#!/usr/bin/env python3
"""
@author Andrew D Johnston
@author Fiach Antaw

Fragmentomics scoring + peak calling pipeline.

What this script does (high level):
1) Reads paired-end fragments from one or more BAMs in a contig or contig region.
2) Filters duplicates (same fragment coords) to max N and optionally subsamples.
3) For each fragment length in a specified range, adds a precomputed PNS score
   distribution across the fragment.
4) Also computes simple coverage and a dyad count (fragment centre) track.
5) Smooths the PNS track (Savitzky–Golay).
6) Calls:
   - positive peaks ("nucleosome regions") on smoothed PNS
   - negative peaks ("breakpoint peaks") by flipping the PNS sign and re-calling
7) Writes (configurable):
   - bedGraph (combined multi-track bedGraph; legacy behavior)
   - wig.gz (one file per track; fixedStep)
   - a BED file of nucleosome peaks
   - a BED file of breakpoint peaks

Randomization modes:
- none: no randomization (default)
- uniform: uniformly randomize fragment start positions within each processed window
- dinuc_anchor: for each fragment, choose to anchor on its start or end (probability
  controlled by --anchor-prob-start), then place the fragment so its anchored boundary
  dinucleotide matches a random occurrence of that dinucleotide in the reference sequence
  for that window (requires --fasta)

"""

import sys
from tqdm import tqdm
import argparse
import pysam
import numpy as np
import os
import gzip
from scipy.signal import savgol_filter
from collections import defaultdict, Counter
import random
from typing import Dict, List, Optional, Tuple

DINUCS = [a + b for a in "ACGT" for b in "ACGT"]


def resolve_fasta_contig(fasta: pysam.FastaFile, contig: str) -> str:
    """
    Resolve contig naming differences between BAM and FASTA.

    Tries:
      contig
      chr+contig
      contig without leading 'chr'
    """
    refs = set(fasta.references)
    if contig in refs:
        return contig
    if contig.startswith("chr"):
        alt = contig[3:]
        if alt in refs:
            return alt
    else:
        alt = f"chr{contig}"
        if alt in refs:
            return alt
    raise KeyError(
        f"Contig '{contig}' not found in FASTA. Tried '{contig}', 'chr{contig}', and "
        f"'{contig[3:]}' (if applicable)."
    )


def is_softclipped_or_padded(cigartuples):
    """Return True if CIGAR has S/H/P ops (soft clip, hard clip, pad)."""
    if not cigartuples:
        return False
    for op, _ln in cigartuples:
        if op in (4, 5, 6):  # 4=S, 5=H, 6=P
            return True
    return False


def generate_paired_reads(
    bamfile,
    contig=None,
    start=None,
    end=None,
    max_duplicates=0,
    subsample=None,
):
    """
    Behavior:
      - fetch(contig,start,end) with multiple_iterators=True
      - skip unmapped or mate-unmapped
      - skip duplicate reads and QC-fail reads
      - skip soft-clipped / hard-clipped / padded reads (CIGAR contains S/H/P)
      - pair mates by query_name (works on coordinate-sorted BAMs)
      - skip same-strand pairs
      - define fragment coords as (ref_name, min(start), max(end))
      - allow up to max_duplicates duplicates PER fragment coordinate (for non-deduplicated bams):
            if frag_counts[key] > max_duplicates: skip
        so max_duplicates=0 keeps exactly 1 instance, max_duplicates=1 keeps up to 2, etc.
      - optional subsampling: keep with probability=subsample
      - yields (fwd, rev) where fwd is the read on the forward strand
    """
    unpaired = {}
    frag_counts = defaultdict(int)

    try:
        it = bamfile.fetch(contig, start, end, multiple_iterators=True)
    except Exception:
        return

    for read in it:
        if read.is_unmapped or read.mate_is_unmapped:
            continue
        if read.is_duplicate or read.is_qcfail:
            continue
        if is_softclipped_or_padded(read.cigartuples):
            continue
        if read.reference_end is None or read.next_reference_start is None:
            continue

        qn = read.query_name
        if qn not in unpaired:
            unpaired[qn] = read
            continue

        mate = unpaired.pop(qn)

        # Mate filters too (defensive)
        if mate.is_unmapped or mate.mate_is_unmapped:
            continue
        if mate.is_duplicate or mate.is_qcfail:
            continue
        if is_softclipped_or_padded(mate.cigartuples):
            continue
        if mate.reference_end is None or mate.next_reference_start is None:
            continue

        # Must be opposite strands for proper PE
        if read.is_reverse == mate.is_reverse:
            continue

        # Optional subsampling
        if subsample is not None and random.random() > subsample:
            continue

        frag_contig = read.reference_name
        frag_start = min(read.reference_start, mate.reference_start)
        frag_end = max(read.reference_end, mate.reference_end)
        if frag_end <= frag_start:
            continue

        key = (frag_contig, frag_start, frag_end)
        if frag_counts[key] > max_duplicates:
            continue
        frag_counts[key] += 1

        # Yield (fwd, rev)
        if not read.is_reverse:
            yield read, mate
        else:
            yield mate, read


def generate_fragment_ranges(
    bamfile,
    contig,
    fetch_start,
    fetch_end,
    max_duplicates,
    subsample,
):
    """
    Convert paired reads into fragment genomic intervals (frag_start, frag_end).
    """
    for r_fwd, r_rev in generate_paired_reads(
        bamfile, contig, fetch_start, fetch_end, max_duplicates, subsample
    ):
        if r_fwd.is_reverse:
            r_fwd, r_rev = r_rev, r_fwd

        if r_fwd.reference_start > r_rev.reference_start:
            continue
        if r_rev.reference_end < r_fwd.reference_end:
            continue

        yield r_fwd.reference_start, r_rev.reference_end


def precompute_distributions(pns_frag_range, mode_DNA_length):
    """
    Precompute per-fragment-length score kernels used to add PNS signal.

    Returns:
      - pns_distributions: mean-normalized kernels
      - pospns_distributions: same kernels before mean normalization
    """
    pns_distributions = {}
    pospns_distributions = {}

    for fragment_length in pns_frag_range:

        if fragment_length < mode_DNA_length:
            total_length = mode_DNA_length + (mode_DNA_length - fragment_length)
        else:
            total_length = fragment_length

        midpoint = (mode_DNA_length - 1) // 2
        second_half_start = midpoint + 1

        scores = np.zeros(total_length)

        for i in range(total_length):
            if i <= midpoint:
                scores[i] = i / midpoint
            elif i <= mode_DNA_length - 1:
                scores[i] = 1 - (i - second_half_start) / midpoint

        end_scores = scores[::-1]
        combined_scores = scores + end_scores

        pospns_distributions[fragment_length] = combined_scores.copy()

        midpoint_val = np.mean(combined_scores)
        centered_scores = [x - midpoint_val for x in combined_scores]
        pns_distributions[fragment_length] = centered_scores

    return pns_distributions, pospns_distributions


def build_dinuc_index(seq: str) -> Dict[str, List[int]]:
    """
    Build index of dinucleotide start positions in a sequence (window-relative).

    Returns dict: dinuc -> list of i such that seq[i:i+2] == dinuc
    """
    idx = {d: [] for d in DINUCS}
    L = len(seq)
    for i in range(0, L - 1):
        d = seq[i : i + 2]
        if d in idx:
            idx[d].append(i)
    return idx


def uniform_randomize_fragments(
    fragments: List[Tuple[int, int]],
    start: int,
    end: int,
) -> List[Tuple[int, int]]:
    """
    Uniformly randomize fragment start positions within [start,end), preserving lengths.
    """
    ref_len = end - start
    randomized = []
    for frag_start, frag_end in fragments:
        L = frag_end - frag_start
        if L <= 0 or L > ref_len:
            continue
        max_start = end - L
        if max_start <= start:
            continue
        new_start = random.randint(start, max_start - 1)
        randomized.append((new_start, new_start + L))
    return randomized


def dinuc_anchor_randomize_fragments(
    fragments: List[Tuple[int, int]],
    start: int,
    end: int,
    window_seq: str,
    dinuc_pos: Dict[str, List[int]],
    anchor_prob_start: float = 0.5,
    max_anchor_tries: int = 30,
    fallback: str = "uniform",  # uniform|keep|skip
) -> List[Tuple[int, int]]:
    """
    Randomize fragments by anchoring on start or end dinucleotide.

    For each fragment:
      - compute start_dinuc and end_dinuc from window_seq using original coordinates
      - choose anchor side: start with prob=anchor_prob_start else end
      - pick a random occurrence of that dinuc in the window reference
      - place fragment so its anchored boundary matches that occurrence
      - require fragment fits fully inside [start,end)

    fallback:
      - uniform: if no placement found, place uniformly at random
      - keep: keep original coordinates
      - skip: drop fragment
    """
    ref_len = end - start
    randomized: List[Tuple[int, int]] = []

    for frag_start, frag_end in fragments:
        L = frag_end - frag_start
        if L <= 0 or L > ref_len:
            if fallback == "keep":
                randomized.append((frag_start, frag_end))
            continue

        s0 = frag_start - start
        e0 = frag_end - start

        if not (0 <= s0 <= ref_len - 2 and 2 <= e0 <= ref_len):
            if fallback == "keep":
                randomized.append((frag_start, frag_end))
            elif fallback == "uniform":
                randomized.extend(uniform_randomize_fragments([(frag_start, frag_end)], start, end))
            continue

        start_d = window_seq[s0 : s0 + 2]
        end_d = window_seq[e0 - 2 : e0]

        if start_d not in dinuc_pos or end_d not in dinuc_pos:
            if fallback == "keep":
                randomized.append((frag_start, frag_end))
            elif fallback == "uniform":
                randomized.extend(uniform_randomize_fragments([(frag_start, frag_end)], start, end))
            continue

        anchor_on_start = (random.random() < anchor_prob_start)

        placed = False
        if anchor_on_start:
            candidates = dinuc_pos.get(start_d, [])
            if candidates:
                for _ in range(max_anchor_tries):
                    i = random.choice(candidates)
                    if i > ref_len - L:
                        continue
                    new_start = start + i
                    new_end = new_start + L
                    if new_end <= end:
                        randomized.append((new_start, new_end))
                        placed = True
                        break
        else:
            candidates = dinuc_pos.get(end_d, [])
            if candidates:
                min_i = L - 2
                max_i = ref_len - 2
                for _ in range(max_anchor_tries):
                    i = random.choice(candidates)
                    if i < min_i or i > max_i:
                        continue
                    new_end = start + i + 2
                    new_start = new_end - L
                    if new_start >= start and new_end <= end:
                        randomized.append((new_start, new_end))
                        placed = True
                        break

        if not placed:
            if fallback == "keep":
                randomized.append((frag_start, frag_end))
            elif fallback == "uniform":
                randomized.extend(uniform_randomize_fragments([(frag_start, frag_end)], start, end))
            elif fallback == "skip":
                pass
            else:
                randomized.extend(uniform_randomize_fragments([(frag_start, frag_end)], start, end))

    return randomized


def score_contig(
    bamfiles,
    contig,
    start,
    end,
    mode_DNA_length,
    pns_frag_range,
    max_duplicates,
    pns_distributions,
    pospns_distributions,
    subsample,
    randomize_mode: str = "none",  # none|uniform|dinuc_anchor
    fasta: Optional[pysam.FastaFile] = None,
    anchor_prob_start: float = 0.5,
    max_anchor_tries: int = 30,
    randomize_fallback: str = "uniform",  # uniform|keep|skip
):
    """
    Score a genomic interval [start, end) on one contig.

    Randomization (within this window) is applied BEFORE scoring if requested.

    Returns:
      scores, pns_frag_range, fragments_filtered (post-filter/subsample/dup, post-randomize)
    """
    ref_len = end - start
    coverage = np.zeros(ref_len, dtype=int)
    dyad = np.zeros(ref_len, dtype=int)
    pns = np.zeros(ref_len, dtype=float)
    pospns = np.zeros(ref_len, dtype=float)
    fragment_ends = np.zeros(ref_len, dtype=int)

    # 1) Collect all fragments for this window (after filtering/subsampling/dup logic)
    fragments: List[Tuple[int, int]] = []
    for bamfile in bamfiles:
        for frag_start, frag_end in generate_fragment_ranges(
            bamfile, contig, start, end, max_duplicates, subsample
        ):
            fragments.append((frag_start, frag_end))

    # 2) Randomize fragments if requested
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

    # 3) Score fragments
    for frag_start, frag_end in fragments:
        frag_length = frag_end - frag_start
        pns_fragment_scores = pns_distributions.get(frag_length)
        pospns_fragment_scores = pospns_distributions.get(frag_length)

        # Add kernels if fragment length in allowed range
        if (
            frag_length in pns_frag_range
            and pns_fragment_scores is not None
            and pospns_fragment_scores is not None
        ):

            if frag_length < mode_DNA_length:
                total_length = mode_DNA_length + (mode_DNA_length - frag_length)
                fragment_center = frag_start + (frag_length // 2) - start
                start_pos = fragment_center - (total_length // 2)
                end_pos = start_pos + total_length
            else:
                start_pos = frag_start - start
                end_pos = frag_end - start

            pns_scores_to_add = np.array(pns_fragment_scores, dtype=float)
            pospns_scores_to_add = np.array(pospns_fragment_scores, dtype=float)

            if start_pos < 0:
                pns_scores_to_add = pns_scores_to_add[-start_pos:]
                pospns_scores_to_add = pospns_scores_to_add[-start_pos:]
                start_pos = 0
            if end_pos > ref_len:
                trim_len = ref_len - start_pos
                pns_scores_to_add = pns_scores_to_add[:trim_len]
                pospns_scores_to_add = pospns_scores_to_add[:trim_len]

            if 0 <= start_pos < ref_len:
                pns[start_pos : start_pos + len(pns_scores_to_add)] += pns_scores_to_add
                pospns[start_pos : start_pos + len(pospns_scores_to_add)] += pospns_scores_to_add

        # Coverage, dyad, and fragment ends
        if frag_start >= start and frag_end <= end:
            if frag_length in pns_frag_range:
                coverage[frag_start - start : frag_end - start] += 1

                # Dyad
                fragment_center = frag_start + (frag_length // 2) - start
                if 0 <= fragment_center < ref_len:
                    dyad[fragment_center] += 1

                # Fragment ends (both ends)
                left_end = frag_start - start
                right_end = frag_end - 1 - start

                if 0 <= left_end < ref_len:
                    fragment_ends[left_end] += 1
                if 0 <= right_end < ref_len:
                    fragment_ends[right_end] += 1

    # 4) Smooth PNS
    window_size = 21
    polyorder = 2
    if len(pns) >= window_size:
        pns_smoothed = savgol_filter(pns, window_size, polyorder)
    else:
        pns_smoothed = pns.copy()

    scores = {
        "coverage": [(contig, start, coverage)],
        "pns_smoothed": [(contig, start, pns_smoothed)],
        "pns": [(contig, start, pns)],
        "posPNS": [(contig, start, pospns)],
        "dyad": [(contig, start, dyad)],
        "fragment_ends": [(contig, start, fragment_ends)],
    }

    return scores, pns_frag_range, fragments


def find_peaks_and_regions(scores, original_start, min_length, max_neg_run):
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
                    if neg_count > max_neg_run:
                        if last_positive_end is not None and last_positive_end - current_region[0] + 1 >= min_length:
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
            most_negative_index = np.argmin(inter_region_scores) + prev_end + 1
            most_negative_score = scores[most_negative_index]
            negative_peaks.append(most_negative_index)
            negative_peak_scores.append(most_negative_score)

    positive_peaks = []
    positive_peak_scores = []
    for region in positive_regions:
        region_scores = scores[region[0] : region[1] + 1]
        peak_index = np.argmax(region_scores) + region[0]
        positive_peaks.append(peak_index)
        positive_peak_scores.append(scores[peak_index])

    # Convert to BED-style intervals
    positive_peak_regions = [
        (region[0] + original_start, region[1] + original_start + 1)
        for region in positive_regions
    ]

    adjusted_positive_peaks = [
        (start + end) // 2
        for start, end in positive_peak_regions
    ]

    positive_peaks = [p + original_start for p in positive_peaks]
    negative_peaks = [p + original_start for p in negative_peaks]

    return (
        (positive_peaks, positive_peak_scores),
        (negative_peaks, negative_peak_scores),
        (adjusted_positive_peaks, positive_peak_regions),
    )


# ----------------------------
# Score-track writers
# ----------------------------
def write_bedgraph(scores, contigs, out_prefix, first_region=False):
    """
    Legacy: combined bedGraph with multiple tracks as extra columns.
    """
    mode = "w" if first_region else "a"
    original_start, original_end = contigs[0]
    filename = f"{out_prefix}_combined_scores.bedGraph"

    with open(filename, mode) as f:
        first_score_type = list(scores.keys())[0]
        contig_data = scores[first_score_type]

        for contig, start, score_array in contig_data:
            if not contig.startswith("chr"):
                contig = f"chr{contig}"

            for i in range(len(score_array)):
                position = start + i
                if original_start <= position < original_end:
                    line = [contig, str(position), str(position + 1)]
                    for stype in scores.keys():
                        stype_score_array = scores[stype][0][2]
                        if i < len(stype_score_array):
                            v = stype_score_array[i]
                            if float(v).is_integer():
                                line.append(f"{int(v)}")
                            else:
                                line.append(f"{float(v):.6f}")
                        else:
                            line.append("0")
                    f.write("\t".join(line) + "\n")


def _wig_val_to_str(track: str, v) -> str:
    if track in ("coverage", "dyad", "fragment_ends"):
        return str(int(v))
    return f"{float(v):.6f}"


def write_wig_gz_tracks(
    scores: Dict[str, List[Tuple[str, int, np.ndarray]]],
    contig: str,
    adjusted_start: int,
    original_start: int,
    original_end: int,
    handles: Dict[str, gzip.GzipFile],
    tracks: List[str],
):
    """
    Write fixedStep WIG, gzipped, one file per track:
      <out_prefix>_<track>.wig.gz

    WIG fixedStep 'start' is 1-based. We write the core (non-overlap) region only.
    """
    if not tracks:
        return

    chrom = contig if contig.startswith("chr") else f"chr{contig}"
    core_len = max(0, original_end - original_start)
    if core_len == 0:
        return

    for track in tracks:
        if track not in scores:
            continue
        if track not in handles:
            continue

        f = handles[track]
        arr = scores[track][0][2]  # aligned to adjusted_start

        wig_start_1based = original_start + 1
        f.write(f"fixedStep chrom={chrom} start={wig_start_1based} step=1\n")

        for pos in range(original_start, original_end):
            i = pos - adjusted_start
            if 0 <= i < len(arr):
                f.write(_wig_val_to_str(track, arr[i]) + "\n")
            else:
                f.write("0\n")


# ----------------------------
# Peak writers
# ----------------------------
def write_bed8_rows(rows, path, mode):
    """
    BED8 writer (WPS-like), with strand provided by rows.
    Each row: (chrom, start, end, name, score, strand, thick_start, thick_end)
    """
    with open(path, mode) as f:
        for chrom, start, end, name, score, strand, thick_start, thick_end in rows:
            f.write(
                f"{chrom}\t{int(start)}\t{int(end)}\t{name}\t{int(score)}\t{strand}\t{int(thick_start)}\t{int(thick_end)}\n"
            )


def iter_peak_records(peaks, original_start, original_end, flip_scores=False):
    for (contig, _orig_start_key), peak_data in peaks.items():
        chrom = contig if contig.startswith("chr") else f"chr{contig}"

        num_positive_peaks = len(peak_data["region_centres"])
        num_negative_peaks = len(peak_data["negative_peaks"])

        for i in range(num_positive_peaks):
            region_start = peak_data["nucleosome_regions"][i][0]
            region_end = peak_data["nucleosome_regions"][i][1]

            region_centre = peak_data["region_centres"][i]
            raw_peak = peak_data["positive_peaks"][i]

            if not (original_start <= region_centre < original_end):
                continue

            upstream_index = None
            downstream_index = None

            for j in range(num_negative_peaks):
                if peak_data["negative_peaks"][j] < region_centre:
                    upstream_index = j
                else:
                    break

            for j in range(num_negative_peaks):
                if peak_data["negative_peaks"][j] > region_centre:
                    downstream_index = j
                    break

            if upstream_index is not None:
                upstream_negative_peak = peak_data["negative_peaks"][upstream_index]
                upstream_score = peak_data["negative_peak_scores"][upstream_index]
            else:
                upstream_negative_peak = region_centre
                upstream_score = peak_data["positive_peak_scores"][i]

            if downstream_index is not None:
                downstream_negative_peak = peak_data["negative_peaks"][downstream_index]
                downstream_score = peak_data["negative_peak_scores"][downstream_index]
            else:
                downstream_negative_peak = region_centre
                downstream_score = peak_data["positive_peak_scores"][i]

            peak_score = peak_data["positive_peak_scores"][i]
            if flip_scores:
                upstream_score *= -1
                downstream_score *= -1
                peak_score *= -1

            # prominence = float(peak_score) - float(np.mean([upstream_score, downstream_score]))
            peak_height = float(peak_score)

            yield {
                "chrom": chrom,
                "region_start": int(region_start),
                "region_end": int(region_end),
                "region_centre": int(region_centre),
                "raw_peak": int(raw_peak),
                "upstream_negative_peak": int(upstream_negative_peak),
                "downstream_negative_peak": int(downstream_negative_peak),
                "upstream_score": float(upstream_score),
                "downstream_score": float(downstream_score),
                "peak_score": float(peak_score),
                "prominence": float(peak_height),
                "max_coverage": peak_data["max_coverages"][i],
                "max_position": peak_data["max_positions"][i],
            }


def peaks_to_bed8_rows(peaks, original_start, original_end, label, flip_scores, score_scale):
    rows = []

    for rec in iter_peak_records(peaks, original_start, original_end, flip_scores):
        score_int = int(round(rec["prominence"] * score_scale))
        name = f'{rec["chrom"]}:{rec["region_centre"]}_{label}'
        strand = "."
        thick_start = rec["region_centre"]
        thick_end = rec["region_centre"] + 1

        rows.append((
            rec["chrom"],
            rec["region_start"],
            rec["region_end"],
            name,
            score_int,
            strand,
            thick_start,
            thick_end,
        ))

    return rows


def write_nucleosome_peaks_rich(peaks, contigs, out_prefix, first_region=False, flip_scores=False):
    mode = "w" if first_region else "a"
    original_start, original_end = contigs[0]

    nucleosome_filename = f"{out_prefix}.bed"
    with open(nucleosome_filename, mode) as f:
        for rec in iter_peak_records(peaks, original_start, original_end, flip_scores):
            f.write(
                f'{rec["chrom"]}\t{rec["region_start"]}\t{rec["region_end"]}\t'
                f'{rec["prominence"]:.2f}\t{rec["region_centre"]}\t'
                f'{rec["upstream_score"]:.2f}\t{rec["upstream_negative_peak"]}\t'
                f'{rec["downstream_score"]:.2f}\t{rec["downstream_negative_peak"]}\t'
                f'{rec["peak_score"]:.2f}\t{rec["raw_peak"]}\t'
                f'{rec["max_coverage"]}\t{rec["max_position"]}\n'
            )


def split_into_regions(contig, start, end, contig_len, max_length=100000, overlap=1000):
    regions = []
    current_start = start

    while current_start < end:
        original_start = current_start
        original_end = min(current_start + max_length, end)

        adjusted_start = max(0, original_start - overlap)
        adjusted_end = min(contig_len, original_end + overlap)

        regions.append((contig, adjusted_start, adjusted_end, original_start, original_end))
        current_start = original_end

    return regions


def call_and_write_peaks(
    scores,
    coverage_scores,
    adjusted_start,
    original_start,
    original_end,
    contig,
    out_prefix,
    first_region,
    peak_type_label,
    flip_scores,
    peak_format: str,
    peak_score_scale: float,
    min_region_length,
    max_neg_run,
):
    positive_peaks, negative_peaks, region_centres = find_peaks_and_regions(
        scores,
        adjusted_start,
        min_region_length,
        max_neg_run
    )

    max_coverages = []
    max_positions = []
    arr_len = coverage_scores.shape[0]

    for start_abs, end_abs in region_centres[1]:
        region_start_idx = max(0, start_abs - adjusted_start)
        region_end_idx = min(arr_len - 1, end_abs - adjusted_start)

        if region_end_idx < region_start_idx:
            max_coverages.append(0)
            max_positions.append(0)
            continue

        region_coverage = coverage_scores[region_start_idx : region_end_idx + 1]

        if region_coverage.size > 0:
            local_argmax = int(np.argmax(region_coverage))
            max_coverages.append(int(region_coverage[local_argmax]))
            max_positions.append(region_start_idx + local_argmax + adjusted_start)
        else:
            max_coverages.append(0)
            max_positions.append(0)

    peaks = {
        (contig, original_start): {
            "positive_peaks": positive_peaks[0],
            "positive_peak_scores": positive_peaks[1],
            "negative_peaks": negative_peaks[0],
            "negative_peak_scores": negative_peaks[1],
            "region_centres": region_centres[0],
            "nucleosome_regions": region_centres[1],
            "max_coverages": max_coverages,
            "max_positions": max_positions,
        }
    }

    if peak_format == "rich":
        write_nucleosome_peaks_rich(
            peaks,
            [(original_start, original_end)],
            out_prefix + peak_type_label,
            first_region,
            flip_scores,
        )
    elif peak_format == "bed8":
        mode = "w" if first_region else "a"
        out_path = f"{out_prefix}{peak_type_label}.bed"

        label = "nuc" if not flip_scores else "brk"
        rows = peaks_to_bed8_rows(
            peaks=peaks,
            original_start=original_start,
            original_end=original_end,
            label=label,
            flip_scores=flip_scores,
            score_scale=peak_score_scale,
        )
        write_bed8_rows(rows, out_path, mode=mode)
    else:
        raise ValueError(f"Unknown peak_format: {peak_format}")


def require_bam_indexes(bam_paths, parser=None):
    missing = []

    for bam in bam_paths:
        bam = os.path.abspath(bam)
        bam_dir = os.path.dirname(bam)

        idx1 = bam + ".bai"
        idx2 = os.path.join(bam_dir, os.path.splitext(os.path.basename(bam))[0] + ".bai")

        if not (os.path.exists(idx1) or os.path.exists(idx2)):
            missing.append((bam, idx1, idx2))

    if missing:
        msg_lines = [
            "ERROR: Missing BAM index (.bai) for the following BAM file(s).",
            "Each BAM input must have its .bai index in the SAME directory as the BAM:",
            "",
        ]
        for bam, idx1, idx2 in missing:
            msg_lines.append(f"  BAM: {bam}")
            msg_lines.append(f"    expected: {idx1}")
            msg_lines.append(f"         or : {idx2}")
            msg_lines.append("")
        msg_lines.append("Create indexes with:")
        msg_lines.append("  samtools index <file.bam>")
        msg = "\n".join(msg_lines)

        if parser is not None:
            parser.error(msg)
        else:
            raise FileNotFoundError(msg)


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


def main():
    parser = argparse.ArgumentParser(description="Score fragmentomics data.")
    parser.add_argument("-b", "--bamfiles", nargs="+", required=True, help="BAM file(s) to process")
    parser.add_argument("-o", "--out_prefix", help="prefix for output files (default: based on BAM names and contigs)")
    parser.add_argument("-c", "--contigs", nargs="+", help='limit to contig(s) and optional range, e.g. "2:100000-200000"')
    parser.add_argument("--mode-length", type=int, default=167, help="Mode fragment length (used in kernel geometry)")
    parser.add_argument("--frag-lower", type=int, default=137, help="Lower fragment length to include")
    parser.add_argument("--frag-upper", type=int, default=197, help="Upper fragment length to include")
    parser.add_argument("--max-duplicates", type=int, default=0, help="Max allowed duplicate fragments with same coords")
    parser.add_argument("--chunk-bp", type=int, default=100000, help="Chunk size for windowing")
    parser.add_argument("--overlap-bp", type=int, default=1000, help="Overlap padding on each side of chunk")
    parser.add_argument("--subsample", type=float, default=None, help="Subsampling proportion (e.g., 0.5 to subsample 50%% of the reads)")
    parser.add_argument("--seed", type=int, default=None, help="Random seed (for reproducibility)")

    # Randomization controls
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

    parser.add_argument(
        "--score-format",
        choices=["bedgraph", "wiggz", "both", "none"],
        default="wiggz",
        help="How to write per-base score tracks. 'wiggz' writes one <prefix>_<track>.wig.gz per track.",
    )
    parser.add_argument(
        "--score-tracks",
        nargs="*",
        default=["coverage", "pns_smoothed", "pns", "posPNS", "dyad", "fragment_ends"],
        help=(
            "Which score tracks to output (space-separated). "
            "Valid: coverage pns_smoothed pns posPNS dyad. "
            "Use '--score-tracks none' or '--score-format none' to disable."
        ),
    )

    parser.add_argument(
        "--peak-format",
        choices=["rich", "bed8"],
        default="bed8",
        help="Peak BED8 output format.",
    )
    parser.add_argument(
        "--peak-score-scale",
        type=float,
        default=1.0,
        help="Only for --peak-format bed8: score = round(prominence * scale) written as int.",
    )

    parser.add_argument(
        "--min-region-length",
        type=int,
        default=50,
        help="Minimum length (bp) of positive regions"
    )

    parser.add_argument(
        "--max-neg-run",
        type=int,
        default=5,
        help="Maximum consecutive non-positive bases allowed within a positive region"
    )

    args = parser.parse_args()

    valid_tracks = {"coverage", "pns_smoothed", "pns", "posPNS", "dyad", "fragment_ends"}
    if args.score_tracks and len(args.score_tracks) == 1 and args.score_tracks[0].lower() == "none":
        args.score_tracks = []
    else:
        bad = [t for t in args.score_tracks if t not in valid_tracks]
        if bad:
            parser.error(f"Unknown --score-tracks: {bad}. Valid: {sorted(valid_tracks)}")

    if args.seed is not None:
        random.seed(args.seed)
        np.random.seed(args.seed)

    require_bam_indexes(args.bamfiles, parser=parser)

    if args.randomize_mode == "dinuc_anchor" and not args.fasta:
        parser.error("--randomize-mode dinuc_anchor requires --fasta <ref.fa> (with .fai index).")

    if not args.out_prefix:
        bam_basenames = [os.path.splitext(os.path.basename(bam))[0] for bam in args.bamfiles]
        args.out_prefix = f"{'_'.join(bam_basenames)}"
        if args.contigs and len(args.contigs) == 1:
            safe_contig = args.contigs[0].replace(":", "_")
            args.out_prefix = f"{args.out_prefix}_{safe_contig}"

    args.out_prefix = f"{args.out_prefix}_mode{args.mode_length}_lower{args.frag_lower}_upper{args.frag_upper}"

    # Open BAMs
    bamfiles = []
    for bamfile_path in args.bamfiles:
        try:
            bamfile = pysam.AlignmentFile(bamfile_path, "rb")
            bamfiles.append(bamfile)
        except FileNotFoundError:
            parser.error(f"Unable to open bamfile {bamfile_path} (file not found)")
            return 2
        except Exception as e:
            parser.error(f"Unable to open bamfile {bamfile_path}: {str(e)}")
            return 2

    # Open FASTA if needed
    fasta = None
    if args.fasta:
        try:
            fasta = pysam.FastaFile(args.fasta)
        except Exception as e:
            parser.error(f"Unable to open FASTA '{args.fasta}': {str(e)}")
            return 2

    # Build list of regions to process
    contigs = []
    if args.contigs:
        for contig_range in args.contigs:
            if ":" in contig_range:
                contig, positions = contig_range.split(":")
                start, end = map(int, positions.split("-"))
                contig_len = bamfiles[0].get_reference_length(contig)
                contigs.extend(
                    split_into_regions(
                        contig, start, end, contig_len,
                        max_length=args.chunk_bp,
                        overlap=args.overlap_bp
                    )
                )
            else:
                contig = contig_range
                start, end = 0, bamfiles[0].get_reference_length(contig)
                contig_len = bamfiles[0].get_reference_length(contig)
                contigs.extend(
                    split_into_regions(
                        contig, start, end, contig_len,
                        max_length=args.chunk_bp,
                        overlap=args.overlap_bp
                    )
                )
    else:
        for contig in bamfiles[0].references:
            start, end = 0, bamfiles[0].get_reference_length(contig)
            contig_len = bamfiles[0].get_reference_length(contig)
            contigs.extend(
                split_into_regions(
                    contig, start, end, contig_len,
                    max_length=args.chunk_bp,
                    overlap=args.overlap_bp
                )
            )

    # Precompute kernels
    pns_frag_range = range(args.frag_lower, args.frag_upper + 1)
    pns_distributions, pospns_distributions = precompute_distributions(
        pns_frag_range, args.mode_length
    )

    # Remove old outputs
    combined_bedgraph = f"{args.out_prefix}_combined_scores.bedGraph"
    nucleosome_bed_rich = f"{args.out_prefix}_nucleosome_regions.bed"
    breakpoint_bed_rich = f"{args.out_prefix}_breakpoint_peaks.bed"
    nucleosome_bed8 = f"{args.out_prefix}_nucleosome_regions.bed"
    breakpoint_bed8 = f"{args.out_prefix}_breakpoint_peaks.bed"
    frag_summary = f"{args.out_prefix}_fragment_summary.txt"
    frag_lens = f"{args.out_prefix}_fragment_length_counts.tsv"

    # Track outputs
    wig_paths = [f"{args.out_prefix}_{t}.wig.gz" for t in args.score_tracks]

    to_remove = [frag_summary, frag_lens]
    if args.score_format in ("bedgraph", "both"):
        to_remove.append(combined_bedgraph)

    if args.peak_format == "rich":
        to_remove.extend([nucleosome_bed_rich, breakpoint_bed_rich])
    else:
        to_remove.extend([nucleosome_bed8, breakpoint_bed8])

    if args.score_format in ("wiggz", "both"):
        to_remove.extend(wig_paths)

    for fname in to_remove:
        if os.path.exists(fname):
            os.remove(fname)

    # Global fragment accounting
    total_fragments_filtered_all = 0
    total_fragments_used_in_range = 0
    unique_bases_covered_by_used = 0
    length_counts = Counter()

    # Open wig.gz handles ONCE for the whole run
    wig_handles = {}
    if args.score_format in ("wiggz", "both") and args.score_tracks:
        for track in args.score_tracks:
            path = f"{args.out_prefix}_{track}.wig.gz"
            wig_handles[track] = gzip.open(path, "wt")

    first_region = True
    for contig, adjusted_start, adjusted_end, original_start, original_end in tqdm(contigs, desc="Scoring contigs"):
        scores, pns_frag_range, fragments_filtered = score_contig(
            bamfiles=bamfiles,
            contig=contig,
            start=adjusted_start,
            end=adjusted_end,
            mode_DNA_length=args.mode_length,
            pns_frag_range=pns_frag_range,
            max_duplicates=args.max_duplicates,
            pns_distributions=pns_distributions,
            pospns_distributions=pospns_distributions,
            subsample=args.subsample,
            randomize_mode=args.randomize_mode,
            fasta=fasta,
            anchor_prob_start=args.anchor_prob_start,
            max_anchor_tries=args.max_anchor_tries,
            randomize_fallback=args.randomize_fallback,
        )

        # ---- Assign fragments to chunk by START in ORIGINAL window ----
        owned_fragments = []
        for frag_start, frag_end in fragments_filtered:
            if original_start <= frag_start < original_end:
                owned_fragments.append((frag_start, frag_end))

        total_fragments_filtered_all += len(owned_fragments)

        if original_end > original_start:
            covered = np.zeros(original_end - original_start, dtype=bool)

            for frag_start, frag_end in owned_fragments:
                L = frag_end - frag_start
                if L not in pns_frag_range:
                    continue

                total_fragments_used_in_range += 1
                length_counts[L] += 1

                ov_start = max(frag_start, original_start)
                ov_end = min(frag_end, original_end)
                if ov_end > ov_start:
                    covered[ov_start - original_start : ov_end - original_start] = True

            unique_bases_covered_by_used += int(covered.sum())

        pns_smoothed_scores = scores["pns_smoothed"][0][2]
        coverage_scores = scores["coverage"][0][2]

        # Nucleosome regions
        call_and_write_peaks(
            scores=pns_smoothed_scores,
            coverage_scores=coverage_scores,
            adjusted_start=adjusted_start,
            original_start=original_start,
            original_end=original_end,
            contig=contig,
            out_prefix=args.out_prefix,
            first_region=first_region,
            peak_type_label="_nucleosome_regions",
            flip_scores=False,
            peak_format=args.peak_format,
            peak_score_scale=args.peak_score_scale,
            min_region_length=args.min_region_length,
            max_neg_run=args.max_neg_run,
        )

        # Breakpoint peaks (flip sign)
        flipped_scores = -1 * pns_smoothed_scores
        call_and_write_peaks(
            scores=flipped_scores,
            coverage_scores=coverage_scores,
            adjusted_start=adjusted_start,
            original_start=original_start,
            original_end=original_end,
            contig=contig,
            out_prefix=args.out_prefix,
            first_region=first_region,
            peak_type_label="_breakpoint_peaks",
            flip_scores=True,
            peak_format=args.peak_format,
            peak_score_scale=args.peak_score_scale,
            min_region_length=args.min_region_length,
            max_neg_run=args.max_neg_run,
        )

        # Per-base scores, trimmed to non-overlap region
        if args.score_format in ("bedgraph", "both"):
            write_bedgraph(scores, [(original_start, original_end)], args.out_prefix, first_region)

        if args.score_format in ("wiggz", "both") and args.score_tracks:
            write_wig_gz_tracks(
                scores=scores,
                contig=contig,
                adjusted_start=adjusted_start,
                original_start=original_start,
                original_end=original_end,
                handles=wig_handles,
                tracks=args.score_tracks,
            )

        first_region = False

    write_fragment_outputs(
        out_prefix=args.out_prefix,
        total_fragments_filtered_all=total_fragments_filtered_all,
        total_fragments_used_in_range=total_fragments_used_in_range,
        unique_bases_covered_by_used=unique_bases_covered_by_used,
        length_counts=length_counts,
    )

    for f in wig_handles.values():
        try:
            f.close()
        except Exception:
            pass

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