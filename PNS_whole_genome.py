#!/usr/bin/env python3
"""
@author Andrew D Johnston
@author Fiach Antaw

Fragmentomics scoring + peak calling pipeline.

What this script does (high level):
1) Reads paired-end fragments from one or more BAMs in a contig or contig region.
2) Filters reads, optionally subsamples fragments, and limits the number of
   fragments retained with identical genomic coordinates.
3) For each fragment length in a specified range, adds a precomputed PNS score
   distribution across the fragment.
4) Also computes simple coverage, dyad count (fragment centre), and fragment-end tracks.
5) Optionally smooths the PNS track (Savitzky–Golay).
6) Optionally calls PNS-based peaks:
   - positive peaks ("nucleosome regions") on smoothed PNS
   - negative peaks ("breakpoint peaks") by flipping the PNS sign and re-calling
   Use --pns-mode off to skip PNS scoring, smoothing, and peak calling.
7) Optionally calculates pure observed dinucleotide profiles aligned to dyad 0.
8) Optionally classifies fragments into WW/SS Types 1-4 from a centred 147-bp core
   and writes separate outputs for each type plus the legacy all-fragment outputs.
9) Writes (configurable):
   - score tracks directly as BigWig files by default
   - optional wig.gz output instead of, or alongside, BigWig
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
import argparse
import pysam
import numpy as np
import os
import gzip

try:
    import pyBigWig
except ImportError:
    pyBigWig = None
from scipy.signal import savgol_filter
from collections import defaultdict, Counter
import random
from typing import Dict, List, Optional, Tuple

DINUCS = [a + b for a in "ACGT" for b in "ACGT"]
WW_DINUCS = {"AA", "AT", "TA", "TT"}
SS_DINUCS = {"CC", "CG", "GC", "GG"}
VALID_BASES = {"A", "C", "G", "T"}
WW_TYPE_GROUPS = ("type1", "type2", "type3", "type4")
ALL_OUTPUT_GROUPS = ("all",) + WW_TYPE_GROUPS
WW_MAJOR_SITE_COEFFICIENT = 36.0 / 32.0  # 1.125

# Base positions relative to the dyad (dyad = 0), copied from Supplementary
# Table S1. Each tuple is an inclusive range of bases belonging to one
# minor- or major-groove bending site. The SHL +/-0.5 sites are excluded,
# following the paper.
MINOR_GBS_RELATIVE = (
    (-69, -66),
    (-59, -56),
    (-48, -45),
    (-37, -34),
    (-27, -24),
    (-17, -14),
    (14, 17),
    (24, 27),
    (34, 37),
    (45, 48),
    (56, 59),
    (66, 69),
)

MAJOR_GBS_RELATIVE = (
    (-64, -61),
    (-53, -51),
    (-42, -40),
    (-32, -29),
    (-22, -19),
    (-12, -9),
    (9, 12),
    (19, 22),
    (29, 32),
    (40, 42),
    (51, 53),
    (61, 64),
)


def _dinuc_starts_from_base_ranges(base_ranges):
    """Convert inclusive base ranges into relative dinucleotide-start sites."""
    return tuple(
        pos
        for site_start, site_end in base_ranges
        for pos in range(site_start, site_end)
    )


MINOR_GBS_DINUC_STARTS = _dinuc_starts_from_base_ranges(MINOR_GBS_RELATIVE)
MAJOR_GBS_DINUC_STARTS = _dinuc_starts_from_base_ranges(MAJOR_GBS_RELATIVE)

if len(MINOR_GBS_DINUC_STARTS) != 36:
    raise RuntimeError("Expected 36 minor-groove dinucleotide positions.")
if len(MAJOR_GBS_DINUC_STARTS) != 32:
    raise RuntimeError("Expected 32 major-groove dinucleotide positions.")


def fragment_dyad(frag_start: int, frag_end: int) -> int:
    """
    Return the dyad base used for sequence alignment and WW-type assignment.

    Odd fragments use their true central base. For even fragments, this uses
    the right-hand central base, as requested:
        dyad = frag_start + floor(fragment_length / 2)
    """
    return frag_start + ((frag_end - frag_start) // 2)


def group_output_prefix(out_prefix: str, group: str) -> str:
    """Keep legacy unsuffixed output names for 'all'; suffix type-specific files."""
    if group == "all":
        return out_prefix
    return f"{out_prefix}_{group}"


def prepare_reference_context(
    fasta: pysam.FastaFile,
    contig: str,
    start: int,
    end: int,
):
    """Fetch one reference window and retain metadata for fast fragment slicing."""
    fasta_contig = resolve_fasta_contig(fasta, contig)
    fasta_length = fasta.get_reference_length(fasta_contig)
    seq = fasta.fetch(fasta_contig, start, end).upper()
    expected = end - start
    if len(seq) != expected:
        seq = seq[:expected].ljust(expected, "N")
    return {
        "contig": fasta_contig,
        "length": fasta_length,
        "start": start,
        "end": end,
        "seq": seq,
    }


def extract_reference_sequence(
    fasta: pysam.FastaFile,
    reference_context,
    seq_start: int,
    seq_end: int,
) -> Optional[str]:
    """Return an uppercase reference sequence or None when outside the contig."""
    if seq_start < 0 or seq_end <= seq_start:
        return None
    if seq_end > reference_context["length"]:
        return None

    window_start = reference_context["start"]
    window_end = reference_context["end"]
    if seq_start >= window_start and seq_end <= window_end:
        rel_start = seq_start - window_start
        rel_end = seq_end - window_start
        seq = reference_context["seq"][rel_start:rel_end]
    else:
        seq = fasta.fetch(reference_context["contig"], seq_start, seq_end).upper()

    if len(seq) != seq_end - seq_start:
        return None
    return seq


def sequence_is_acgt(seq: str) -> bool:
    return bool(seq) and all(base in VALID_BASES for base in seq)


def count_ww_ss_at_relative_sites(core_seq: str, relative_starts) -> Tuple[int, int]:
    """
    Count WW and SS dinucleotides at selected sites in a centred 147-bp core.

    core_seq[73] is the dyad base (relative position 0), so a relative
    dinucleotide start p maps to core_seq[p + 73 : p + 75].
    """
    ww = 0
    ss = 0
    for rel_start in relative_starts:
        i = rel_start + 73
        dinuc = core_seq[i : i + 2]
        if dinuc in WW_DINUCS:
            ww += 1
        elif dinuc in SS_DINUCS:
            ss += 1
    return ww, ss


def classify_fragment_ww_type(
    fasta: pysam.FastaFile,
    reference_context,
    frag_start: int,
    frag_end: int,
) -> Optional[str]:
    """
    Classify a fragment as WW/SS Type 1-4 using a centred 147-bp sequence.

    Every observed fragment length uses the same 147-bp reference template:
      core_start = dyad - 73
      core_end   = dyad + 74

    Returns None when the core crosses a reference boundary or contains a base
    other than A/C/G/T.
    """
    dyad = fragment_dyad(frag_start, frag_end)
    core_seq = extract_reference_sequence(
        fasta=fasta,
        reference_context=reference_context,
        seq_start=dyad - 73,
        seq_end=dyad + 74,
    )
    if core_seq is None or len(core_seq) != 147 or not sequence_is_acgt(core_seq):
        return None

    minor_ww, minor_ss = count_ww_ss_at_relative_sites(
        core_seq, MINOR_GBS_DINUC_STARTS
    )
    major_ww, major_ss = count_ww_ss_at_relative_sites(
        core_seq, MAJOR_GBS_DINUC_STARTS
    )

    ww_minor_enriched = minor_ww >= (major_ww * WW_MAJOR_SITE_COEFFICIENT)
    ss_minor_enriched = minor_ss > (major_ss * WW_MAJOR_SITE_COEFFICIENT)

    if ww_minor_enriched and not ss_minor_enriched:
        return "type1"
    if ww_minor_enriched and ss_minor_enriched:
        return "type2"
    if not ww_minor_enriched and not ss_minor_enriched:
        return "type3"
    return "type4"


def expected_dinuc_profile_positions(frag_lower: int, frag_upper: int) -> List[int]:
    """Return all possible dyad-relative dinucleotide-start positions in the range."""
    positions = set()
    for frag_length in range(frag_lower, frag_upper + 1):
        dyad_offset = frag_length // 2
        for i in range(frag_length - 1):
            positions.add(i - dyad_offset)
    return sorted(positions)


def new_dinuc_accumulator():
    return {
        "counts": defaultdict(Counter),
        "n_valid": Counter(),
        "fragments_used": 0,
        "fragments_skipped": 0,
    }


def add_fragment_to_dinuc_accumulator(
    accumulator,
    fragment_seq: Optional[str],
    frag_start: int,
    frag_end: int,
) -> bool:
    """Add a pure observed 16-dinucleotide profile aligned to dyad position 0."""
    if fragment_seq is None or len(fragment_seq) != frag_end - frag_start:
        accumulator["fragments_skipped"] += 1
        return False
    if not sequence_is_acgt(fragment_seq):
        accumulator["fragments_skipped"] += 1
        return False

    dyad = fragment_dyad(frag_start, frag_end)
    for i in range(len(fragment_seq) - 1):
        dinuc = fragment_seq[i : i + 2]
        rel_pos = (frag_start + i) - dyad
        accumulator["counts"][rel_pos][dinuc] += 1
        accumulator["n_valid"][rel_pos] += 1

    accumulator["fragments_used"] += 1
    return True


def write_dinuc_profile(
    out_path: str,
    accumulator,
    positions: List[int],
    fraction: bool = False,
):
    """Write pure observed dinucleotide fractions or percentages."""
    multiplier = 1.0 if fraction else 100.0
    suffix = "frac" if fraction else "pct"

    header = ["position", "n_valid"]
    header.extend([f"{dinuc}_{suffix}" for dinuc in DINUCS])
    header.extend([f"WW_{suffix}", f"SS_{suffix}"])

    with open(out_path, "w") as out:
        out.write("\t".join(header) + "\n")
        for rel_pos in positions:
            n_valid = int(accumulator["n_valid"].get(rel_pos, 0))
            counts = accumulator["counts"].get(rel_pos, Counter())
            row = [str(rel_pos), str(n_valid)]

            if n_valid == 0:
                row.extend(["NaN"] * (len(DINUCS) + 2))
            else:
                for dinuc in DINUCS:
                    row.append(f"{(counts[dinuc] / n_valid) * multiplier:.8g}")
                ww_count = sum(counts[d] for d in WW_DINUCS)
                ss_count = sum(counts[d] for d in SS_DINUCS)
                row.append(f"{(ww_count / n_valid) * multiplier:.8g}")
                row.append(f"{(ss_count / n_valid) * multiplier:.8g}")

            out.write("\t".join(row) + "\n")


def write_ww_type_summary(
    out_prefix: str,
    type_counts: Counter,
    total_in_range: int,
):
    """Write counts and percentages for Type 1-4 and unclassified fragments."""
    path = f"{out_prefix}_ww_type_summary.tsv"
    classified_total = sum(type_counts[t] for t in WW_TYPE_GROUPS)

    with open(path, "w") as out:
        out.write(
            "type\tfragment_count\tpercent_of_all_in_range"
            "\tpercent_of_classified\n"
        )
        for group in WW_TYPE_GROUPS:
            count = int(type_counts[group])
            pct_all = (100.0 * count / total_in_range) if total_in_range else float("nan")
            pct_classified = (
                100.0 * count / classified_total if classified_total else float("nan")
            )
            out.write(
                f"{group}\t{count}\t{pct_all:.8g}\t{pct_classified:.8g}\n"
            )

        unclassified = int(type_counts["unclassified"])
        pct_all = (
            100.0 * unclassified / total_in_range if total_in_range else float("nan")
        )
        out.write(f"unclassified\t{unclassified}\t{pct_all:.8g}\tNaN\n")
        out.write(f"all\t{int(total_in_range)}\t100\tNaN\n")


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
    max_duplicates=1,
    subsample=None,
    frag_counts=None,
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
      - coordinate-based deduplication limits the total number of fragments
        retained per coordinate:
            if max_duplicates > 0 and frag_counts[key] >= max_duplicates: skip
        max_duplicates=1 keeps exactly 1 fragment,
        max_duplicates=2 keeps up to 2 fragments, and
        max_duplicates=0 disables coordinate-based deduplication.
      - frag_counts may be shared across BAMs so coordinate-based deduplication
        is applied across all BAM inputs, or omitted to deduplicate each BAM
        independently.
      - optional subsampling is a separate operation and retains an accepted
        fragment with probability=subsample
      - yields (fwd, rev) where fwd is the read on the forward strand
    """
    unpaired = {}
    if frag_counts is None:
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

        frag_contig = read.reference_name
        frag_start = min(read.reference_start, mate.reference_start)
        frag_end = max(read.reference_end, mate.reference_end)
        if frag_end <= frag_start:
            continue

        key = (frag_contig, frag_start, frag_end)
        if max_duplicates > 0 and frag_counts[key] >= max_duplicates:
            continue
        frag_counts[key] += 1

        # Optional subsampling, applied after coordinate-based deduplication
        if subsample is not None and random.random() > subsample:
            continue

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
    frag_counts=None,
):
    """
    Convert paired reads into fragment genomic intervals (frag_start, frag_end).
    """
    for r_fwd, r_rev in generate_paired_reads(
        bamfile, contig, fetch_start, fetch_end, max_duplicates, subsample, frag_counts
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


def _new_track_arrays(ref_len: int, do_pns: bool):
    arrays = {
        "coverage": np.zeros(ref_len, dtype=int),
        "dyad": np.zeros(ref_len, dtype=float),
        "fragment_ends": np.zeros(ref_len, dtype=int),
        "fragment_left_ends": np.zeros(ref_len, dtype=int),
        "fragment_right_ends": np.zeros(ref_len, dtype=int),
    }
    if do_pns:
        arrays["pns"] = np.zeros(ref_len, dtype=float)
        arrays["posPNS"] = np.zeros(ref_len, dtype=float)
    return arrays


def _add_fragment_to_track_arrays(
    arrays,
    frag_start: int,
    frag_end: int,
    window_start: int,
    window_end: int,
    mode_DNA_length: int,
    pns_frag_range,
    pns_distributions,
    pospns_distributions,
    do_pns: bool,
):
    """Add one fragment to one output group's PNS/coverage/dyad/end arrays."""
    frag_length = frag_end - frag_start
    if frag_length not in pns_frag_range:
        return

    ref_len = window_end - window_start

    if do_pns:
        pns_fragment_scores = pns_distributions.get(frag_length)
        pospns_fragment_scores = pospns_distributions.get(frag_length)

        if pns_fragment_scores is not None and pospns_fragment_scores is not None:
            if frag_length < mode_DNA_length:
                total_length = mode_DNA_length + (mode_DNA_length - frag_length)
                centre = fragment_dyad(frag_start, frag_end) - window_start
                start_pos = centre - (total_length // 2)
                end_pos = start_pos + total_length
            else:
                start_pos = frag_start - window_start
                end_pos = frag_end - window_start

            pns_scores_to_add = np.asarray(pns_fragment_scores, dtype=float)
            pospns_scores_to_add = np.asarray(pospns_fragment_scores, dtype=float)

            if start_pos < 0:
                pns_scores_to_add = pns_scores_to_add[-start_pos:]
                pospns_scores_to_add = pospns_scores_to_add[-start_pos:]
                start_pos = 0
            if end_pos > ref_len:
                trim_len = max(0, ref_len - start_pos)
                pns_scores_to_add = pns_scores_to_add[:trim_len]
                pospns_scores_to_add = pospns_scores_to_add[:trim_len]

            if 0 <= start_pos < ref_len and len(pns_scores_to_add) > 0:
                arrays["pns"][start_pos : start_pos + len(pns_scores_to_add)] += (
                    pns_scores_to_add
                )
                arrays["posPNS"][start_pos : start_pos + len(pospns_scores_to_add)] += (
                    pospns_scores_to_add
                )

    # Preserve the original behavior: coverage, dyad, and end tracks require
    # the complete fragment to lie inside the adjusted scoring window.
    if frag_start < window_start or frag_end > window_end:
        return

    left = frag_start - window_start
    right = frag_end - window_start
    arrays["coverage"][left:right] += 1

    if frag_length % 2 == 1:
        centre = fragment_dyad(frag_start, frag_end) - window_start
        if 0 <= centre < ref_len:
            arrays["dyad"][centre] += 1.0
    else:
        right_centre = fragment_dyad(frag_start, frag_end) - window_start
        left_centre = right_centre - 1
        if 0 <= left_centre < ref_len:
            arrays["dyad"][left_centre] += 0.5
        if 0 <= right_centre < ref_len:
            arrays["dyad"][right_centre] += 0.5

    left_end = frag_start - window_start
    right_end = frag_end - 1 - window_start
    if 0 <= left_end < ref_len:
        arrays["fragment_ends"][left_end] += 1
        arrays["fragment_left_ends"][left_end] += 1
    if 0 <= right_end < ref_len:
        arrays["fragment_ends"][right_end] += 1
        arrays["fragment_right_ends"][right_end] += 1


def _arrays_to_scores(arrays, contig: str, start: int, do_pns: bool):
    scores = {
        "coverage": [(contig, start, arrays["coverage"])],
        "dyad": [(contig, start, arrays["dyad"])],
        "fragment_ends": [(contig, start, arrays["fragment_ends"])],
        "fragment_left_ends": [(contig, start, arrays["fragment_left_ends"])],
        "fragment_right_ends": [(contig, start, arrays["fragment_right_ends"])],
    }

    if do_pns:
        pns = arrays["pns"]
        if len(pns) >= 21:
            pns_smoothed = savgol_filter(pns, 21, 2)
        else:
            pns_smoothed = pns.copy()
        scores.update(
            {
                "pns_smoothed": [(contig, start, pns_smoothed)],
                "pns": [(contig, start, pns)],
                "posPNS": [(contig, start, arrays["posPNS"])],
            }
        )
    return scores


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
    pns_mode: str = "on",  # on|off
    randomize_mode: str = "none",  # none|uniform|dinuc_anchor
    fasta: Optional[pysam.FastaFile] = None,
    anchor_prob_start: float = 0.5,
    max_anchor_tries: int = 30,
    randomize_fallback: str = "uniform",  # uniform|keep|skip
    dedup_scope: str = "all_bams",  # all_bams|per_bam
    split_ww_types: bool = False,
    need_reference_sequence: bool = False,
):
    """
    Score one adjusted genomic window and optionally split fragments by WW type.

    Fragment type assignment uses a centred 147-bp reference sequence around
    dyad position 0, regardless of observed fragment length. Randomization is
    applied before sequence profiling and type assignment.

    Returns:
      scores_by_group:
        'all' always, plus type1-type4 when split_ww_types=True.
      pns_frag_range
      fragment_records:
        (frag_start, frag_end, ww_type_or_None) after all filtering and optional
        randomization.
      reference_context:
        cached FASTA window metadata when sequence access was requested.
    """
    do_pns = pns_mode == "on"
    ref_len = end - start
    groups = ALL_OUTPUT_GROUPS if split_ww_types else ("all",)
    arrays_by_group = {
        group: _new_track_arrays(ref_len, do_pns) for group in groups
    }

    fragments: List[Tuple[int, int]] = []
    shared_frag_counts = defaultdict(int) if dedup_scope == "all_bams" else None

    for bamfile in bamfiles:
        bam_frag_counts = shared_frag_counts if dedup_scope == "all_bams" else None
        for frag_start, frag_end in generate_fragment_ranges(
            bamfile, contig, start, end, max_duplicates, subsample, bam_frag_counts
        ):
            fragments.append((frag_start, frag_end))

    reference_context = None
    if need_reference_sequence:
        if fasta is None:
            raise ValueError("Reference-sequence processing requires --fasta")
        reference_context = prepare_reference_context(fasta, contig, start, end)

    if randomize_mode == "uniform" and fragments:
        fragments = uniform_randomize_fragments(fragments, start, end)

    elif randomize_mode == "dinuc_anchor" and fragments:
        if fasta is None or reference_context is None:
            raise ValueError("randomize_mode=dinuc_anchor requires --fasta")
        dinuc_pos = build_dinuc_index(reference_context["seq"])
        fragments = dinuc_anchor_randomize_fragments(
            fragments=fragments,
            start=start,
            end=end,
            window_seq=reference_context["seq"],
            dinuc_pos=dinuc_pos,
            anchor_prob_start=anchor_prob_start,
            max_anchor_tries=max_anchor_tries,
            fallback=randomize_fallback,
        )

    fragment_records = []
    for frag_start, frag_end in fragments:
        frag_length = frag_end - frag_start
        ww_type = None

        if split_ww_types and frag_length in pns_frag_range:
            ww_type = classify_fragment_ww_type(
                fasta=fasta,
                reference_context=reference_context,
                frag_start=frag_start,
                frag_end=frag_end,
            )

        fragment_records.append((frag_start, frag_end, ww_type))

        _add_fragment_to_track_arrays(
            arrays=arrays_by_group["all"],
            frag_start=frag_start,
            frag_end=frag_end,
            window_start=start,
            window_end=end,
            mode_DNA_length=mode_DNA_length,
            pns_frag_range=pns_frag_range,
            pns_distributions=pns_distributions,
            pospns_distributions=pospns_distributions,
            do_pns=do_pns,
        )

        if split_ww_types and ww_type in WW_TYPE_GROUPS:
            _add_fragment_to_track_arrays(
                arrays=arrays_by_group[ww_type],
                frag_start=frag_start,
                frag_end=frag_end,
                window_start=start,
                window_end=end,
                mode_DNA_length=mode_DNA_length,
                pns_frag_range=pns_frag_range,
                pns_distributions=pns_distributions,
                pospns_distributions=pospns_distributions,
                do_pns=do_pns,
            )

    scores_by_group = {
        group: _arrays_to_scores(arrays, contig, start, do_pns)
        for group, arrays in arrays_by_group.items()
    }

    return scores_by_group, pns_frag_range, fragment_records, reference_context

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
def write_bigwig_tracks(
    scores: Dict[str, List[Tuple[str, int, np.ndarray]]],
    contig: str,
    adjusted_start: int,
    original_start: int,
    original_end: int,
    handles,
    tracks: List[str],
):
    """Write completed chunk data directly to open BigWig handles.

    Dense tracks use fixed-step insertion. Sparse dyad and fragment-end tracks
    write only nonzero one-base intervals. Chunks must be supplied in BAM-header
    contig order and increasing coordinate order, as done by main().
    """
    sparse_tracks = {
        "dyad", "fragment_ends",
        "fragment_left_ends", "fragment_right_ends",
    }

    for track in tracks:
        if track not in scores or track not in handles:
            continue

        arr = scores[track][0][2]
        left = original_start - adjusted_start
        right = original_end - adjusted_start
        values = np.asarray(arr[left:right], dtype=np.float64)
        if values.size == 0:
            continue

        bw = handles[track]
        if track in sparse_tracks:
            nz = np.flatnonzero(values != 0)
            if nz.size == 0:
                continue
            starts = (nz + original_start).astype(np.int64)
            ends = starts + 1
            bw.addEntries(
                [contig] * len(starts),
                starts.tolist(),
                ends=ends.tolist(),
                values=values[nz].astype(float).tolist(),
            )
        else:
            # Fixed-step insertion avoids materialising chromosome and coordinate
            # lists for dense PNS/coverage arrays.
            bw.addEntries(
                contig,
                int(original_start),
                values=values.astype(float).tolist(),
                span=1,
                step=1,
            )


def _wig_val_to_str(track: str, v) -> str:
    if track in ("coverage", "fragment_ends", "fragment_left_ends", "fragment_right_ends"):
        return str(int(v))
    return f"{float(v):.1f}"


def write_wig_gz_tracks(
    scores: Dict[str, List[Tuple[str, int, np.ndarray]]],
    contig: str,
    adjusted_start: int,
    original_start: int,
    original_end: int,
    handles: Dict[str, gzip.GzipFile],
    tracks: List[str],
):
    if not tracks:
        return

    chrom = contig

    if not hasattr(write_wig_gz_tracks, "_last_varstep_chrom"):
        write_wig_gz_tracks._last_varstep_chrom = {}

    for track in tracks:
        if track not in scores or track not in handles:
            continue

        f = handles[track]
        arr = scores[track][0][2]

        if track in ("dyad", "fragment_ends", "fragment_left_ends", "fragment_right_ends"):

            state_key = (track, getattr(f, "name", id(f)))
            last_chrom = write_wig_gz_tracks._last_varstep_chrom.get(state_key)

            if last_chrom != chrom:
                f.write(f"variableStep chrom={chrom}\n")
                write_wig_gz_tracks._last_varstep_chrom[state_key] = chrom

            for pos in range(original_start, original_end):
                i = pos - adjusted_start
                if not (0 <= i < len(arr)):
                    continue

                value = arr[i]
                if value == 0:
                    continue

                f.write(f"{pos + 1}\t{_wig_val_to_str(track, value)}\n")

        else:
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
    dedup_scope: str,
    max_duplicates: int,
):
    summary_path = f"{out_prefix}_fragment_summary.txt"
    lens_path = f"{out_prefix}_fragment_length_counts.tsv"

    with open(summary_path, "w") as f:
        f.write(f"total_fragments_filtered_all\t{total_fragments_filtered_all}\n")
        f.write(f"total_fragments_used_in_range\t{total_fragments_used_in_range}\n")
        f.write(f"unique_bases_covered_by_used_fragments\t{unique_bases_covered_by_used}\n")
        f.write(f"dedup_scope\t{dedup_scope}\n")
        f.write(f"max_per_coordinate\t{max_duplicates}\n")

    with open(lens_path, "w") as f:
        f.write("fragment_length\tcount\n")
        for L in sorted(length_counts.keys()):
            f.write(f"{int(L)}\t{int(length_counts[L])}\n")



def resolve_bam_contig_name(requested: str, bam_references) -> str:
    refs = set(bam_references)
    if requested in refs:
        return requested
    if requested.startswith("chr") and requested[3:] in refs:
        return requested[3:]
    prefixed = f"chr{requested}"
    if prefixed in refs:
        return prefixed
    raise KeyError(f"Contig '{requested}' was not found in the BAM header.")


def expand_contig_tokens(raw_tokens, bam_references):
    """Expand comma lists, numeric ranges, autosomes, and all."""
    if not raw_tokens:
        return list(bam_references)

    tokens = []
    for raw in raw_tokens:
        tokens.extend(part.strip() for part in raw.split(",") if part.strip())

    requested = []
    for token in tokens:
        low = token.lower()
        if low == "all":
            requested.extend(bam_references)
            continue
        if low == "autosomes":
            requested.extend(str(i) for i in range(1, 23))
            continue
        if ":" not in token and "-" in token:
            left, right = token.split("-", 1)
            left_num = left[3:] if left.lower().startswith("chr") else left
            right_num = right[3:] if right.lower().startswith("chr") else right
            if left_num.isdigit() and right_num.isdigit():
                a, b = int(left_num), int(right_num)
                step = 1 if a <= b else -1
                requested.extend(str(i) for i in range(a, b + step, step))
                continue
        requested.append(token)

    resolved = []
    seen = set()
    for item in requested:
        if ":" in item:
            contig_part, region = item.split(":", 1)
            contig = resolve_bam_contig_name(contig_part, bam_references)
            value = f"{contig}:{region}"
        else:
            value = resolve_bam_contig_name(item, bam_references)
        if value not in seen:
            resolved.append(value)
            seen.add(value)
    return resolved

def main():
    parser = argparse.ArgumentParser(
        description=(
            "Score fragmentomics data, optionally output pure observed "
            "dinucleotide profiles and split fragments into WW/SS Types 1-4."
        )
    )
    parser.add_argument("-b", "--bamfiles", nargs="+", required=True, help="BAM file(s) to process")
    parser.add_argument("-o", "--out_prefix", help="prefix for output files (default: based on BAM names and contigs)")
    parser.add_argument(
        "-c", "--contigs", nargs="+",
        help=(
            "Contigs/regions to process. Supports comma lists and keywords, e.g. "
            "1,2,3; chr1-22,chrX,chrY; autosomes; all; or "
            "chr2:100000-200000. Default: all BAM contigs."
        ),
    )
    parser.add_argument("--mode-length", type=int, default=167, help="Mode fragment length (used in kernel geometry)")
    parser.add_argument("--frag-lower", type=int, default=137, help="Lower fragment length to include")
    parser.add_argument("--frag-upper", type=int, default=197, help="Upper fragment length to include")
    parser.add_argument(
        "--max-per-coordinate",
        "--max-duplicates",
        dest="max_duplicates",
        type=int,
        default=1,
        help=(
            "Maximum total number of fragments retained with identical coordinates. "
            "Use 1 for coordinate-based deduplication and 0 to disable it."
        ),
    )
    parser.add_argument(
        "--dedup-scope",
        choices=["all_bams", "per_bam"],
        default="all_bams",
        help=(
            "Coordinate-based deduplication scope. 'all_bams' (default) shares "
            "coordinate counts across every input BAM; 'per_bam' applies the "
            "limit independently within each BAM."
        ),
    )
    parser.add_argument("--chunk-bp", type=int, default=100000, help="Chunk size for windowing")
    parser.add_argument("--overlap-bp", type=int, default=1000, help="Overlap padding on each side of chunk")
    parser.add_argument("--subsample", type=float, default=None, help="Subsampling proportion (e.g., 0.5 to subsample 50%% of the reads)")
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
        help=(
            "Reference FASTA with .fai index. Required for dinuc_anchor, "
            "--dinuc-profile, and --split-ww-types."
        ),
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
        "--dinuc-profile",
        action="store_true",
        help=(
            "Output a pure observed 16-dinucleotide profile for the same in-range "
            "fragments used by the PNS/coverage/dyad tracks. Fragment sequences are "
            "aligned so dyad position is 0. Requires --fasta."
        ),
    )
    parser.add_argument(
        "--dinuc-fraction",
        action="store_true",
        help=(
            "With --dinuc-profile, output fractions from 0 to 1 instead of the "
            "default percentages from 0 to 100."
        ),
    )
    parser.add_argument(
        "--split-ww-types",
        action="store_true",
        help=(
            "Classify each in-range fragment as WW/SS type1-type4 using a centred "
            "147-bp reference sequence around dyad 0, then write separate PNS, dyad, "
            "coverage, fragment-end, peak, summary, and optional dinucleotide files "
            "for type1-type4 while retaining the legacy unsuffixed 'all' outputs. "
            "Requires --fasta."
        ),
    )

    parser.add_argument(
        "--pns-mode",
        choices=["on", "off"],
        default="on",
        help=(
            "Turn PNS scoring and PNS-based peak calling on or off. "
            "Use '--pns-mode off' to compute only non-PNS tracks such as coverage, dyad, "
            "fragment_ends, fragment_left_ends, and fragment_right_ends."
        ),
    )
    parser.add_argument(
        "--pns-format",
        choices=["bigwig", "wiggz", "both", "none"],
        default="bigwig",
        help=(
            "Output format for dense PNS tracks (pns_smoothed, pns, posPNS). "
            "Default: bigwig, written directly with pyBigWig."
        ),
    )
    parser.add_argument(
        "--other-format",
        choices=["bigwig", "wiggz", "both", "none"],
        default="bigwig",
        help=(
            "Output format for coverage, dyad, and fragment-end tracks. "
            "Default: bigwig, written directly with pyBigWig."
        ),
    )
    parser.add_argument(
        "--pns-tracks",
        nargs="*",
        default=["pns_smoothed", "pns", "posPNS"],
        help="PNS tracks to output. Valid: pns_smoothed pns posPNS. Use 'none' to disable.",
    )
    parser.add_argument(
        "--other-tracks",
        nargs="*",
        default=[
            "coverage", "dyad", "fragment_ends",
            "fragment_left_ends", "fragment_right_ends",
        ],
        help=(
            "Non-PNS tracks to output. Valid: coverage dyad fragment_ends "
            "fragment_left_ends fragment_right_ends. Use 'none' to disable."
        ),
    )
    parser.add_argument(
        "--peak-format",
        choices=["rich", "bed8"],
        default="bed8",
        help="Peak output format.",
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
        help="Minimum length (bp) of positive regions",
    )
    parser.add_argument(
        "--max-neg-run",
        type=int,
        default=5,
        help="Maximum consecutive non-positive bases allowed within a positive region",
    )

    args = parser.parse_args()

    if args.max_duplicates < 0:
        parser.error("--max-per-coordinate/--max-duplicates must be 0 or greater.")
    if args.frag_lower < 1 or args.frag_upper < args.frag_lower:
        parser.error("Require 1 <= --frag-lower <= --frag-upper.")
    if args.dinuc_profile and args.frag_lower < 2:
        parser.error("--dinuc-profile requires --frag-lower >= 2.")
    if args.subsample is not None and not (0.0 <= args.subsample <= 1.0):
        parser.error("--subsample must be between 0 and 1.")
    if not (0.0 <= args.anchor_prob_start <= 1.0):
        parser.error("--anchor-prob-start must be between 0 and 1.")

    pns_track_set = {"pns", "posPNS", "pns_smoothed"}
    other_track_set = {
        "coverage", "dyad", "fragment_ends",
        "fragment_left_ends", "fragment_right_ends",
    }

    if len(args.pns_tracks) == 1 and args.pns_tracks[0].lower() == "none":
        args.pns_tracks = []
    if len(args.other_tracks) == 1 and args.other_tracks[0].lower() == "none":
        args.other_tracks = []

    bad_pns = [t for t in args.pns_tracks if t not in pns_track_set]
    bad_other = [t for t in args.other_tracks if t not in other_track_set]
    if bad_pns:
        parser.error(f"Unknown --pns-tracks: {bad_pns}. Valid: {sorted(pns_track_set)}")
    if bad_other:
        parser.error(f"Unknown --other-tracks: {bad_other}. Valid: {sorted(other_track_set)}")

    if args.pns_mode == "off" and args.pns_tracks:
        print(
            "[INFO] --pns-mode off: PNS track output and peak calling are disabled.",
            file=sys.stderr,
        )
        args.pns_tracks = []
        args.pns_format = "none"

    if args.dinuc_fraction and not args.dinuc_profile:
        parser.error("--dinuc-fraction requires --dinuc-profile.")

    reference_required = (
        args.randomize_mode == "dinuc_anchor"
        or args.dinuc_profile
        or args.split_ww_types
    )
    if reference_required and not args.fasta:
        parser.error(
            "--fasta <ref.fa> is required for dinuc_anchor, --dinuc-profile, "
            "or --split-ww-types."
        )

    if args.seed is not None:
        random.seed(args.seed)
        np.random.seed(args.seed)

    require_bam_indexes(args.bamfiles, parser=parser)

    if not args.out_prefix:
        bam_basenames = [os.path.splitext(os.path.basename(bam))[0] for bam in args.bamfiles]
        args.out_prefix = "_".join(bam_basenames)
        if args.contigs and len(args.contigs) == 1:
            safe_contig = args.contigs[0].replace(":", "_")
            args.out_prefix = f"{args.out_prefix}_{safe_contig}"

    args.out_prefix = (
        f"{args.out_prefix}_mode{args.mode_length}_lower{args.frag_lower}"
        f"_upper{args.frag_upper}"
    )
    out_parent = os.path.dirname(os.path.abspath(args.out_prefix))
    if out_parent:
        os.makedirs(out_parent, exist_ok=True)

    bamfiles = []
    for bamfile_path in args.bamfiles:
        try:
            bamfiles.append(pysam.AlignmentFile(bamfile_path, "rb"))
        except FileNotFoundError:
            parser.error(f"Unable to open bamfile {bamfile_path} (file not found)")
        except Exception as exc:
            parser.error(f"Unable to open bamfile {bamfile_path}: {exc}")

    fasta = None
    if args.fasta:
        try:
            fasta = pysam.FastaFile(args.fasta)
        except Exception as exc:
            parser.error(f"Unable to open FASTA '{args.fasta}': {exc}")

    selected_contigs = expand_contig_tokens(args.contigs, bamfiles[0].references)

    # BigWig requires entries in chromosome/coordinate order. Sort all selected
    # contigs and regions according to the BAM header, regardless of the order
    # supplied on the command line.
    reference_order = {name: i for i, name in enumerate(bamfiles[0].references)}
    def _contig_spec_sort_key(spec):
        name = spec.split(":", 1)[0]
        region_start = 0
        if ":" in spec:
            region_start = int(spec.split(":", 1)[1].replace(",", "").split("-", 1)[0])
        return reference_order[name], region_start
    selected_contigs = sorted(selected_contigs, key=_contig_spec_sort_key)

    contigs = []
    selected_contig_names = []
    for contig_spec in selected_contigs:
        if ":" in contig_spec:
            contig, positions = contig_spec.split(":", 1)
            start, end = map(int, positions.replace(",", "").split("-", 1))
            contig_len = bamfiles[0].get_reference_length(contig)
            if not (0 <= start < end <= contig_len):
                parser.error(
                    f"Invalid region {contig_spec}; {contig} length is {contig_len:,}."
                )
        else:
            contig = contig_spec
            contig_len = bamfiles[0].get_reference_length(contig)
            start, end = 0, contig_len

        selected_contig_names.append(contig)
        contigs.extend(
            split_into_regions(
                contig, start, end, contig_len,
                max_length=args.chunk_bp, overlap=args.overlap_bp,
            )
        )

    n_distinct_contigs = len(set(selected_contig_names))
    print(
        f"[INFO] Selected {n_distinct_contigs} contig(s), split into "
        f"{len(contigs):,} processing chunk(s).",
        file=sys.stderr,
    )
    uses_bigwig = (
        args.pns_format in ("bigwig", "both")
        or args.other_format in ("bigwig", "both")
    )
    if uses_bigwig and pyBigWig is None:
        parser.error(
            "BigWig output requires pyBigWig. Install it with: "
            "conda install -c bioconda pybigwig  (or: pip install pyBigWig)"
        )
    if uses_bigwig:
        print(
            "[INFO] Score tracks will be written directly as BigWig files; "
            "no intermediate bedGraph or WIG conversion is required.",
            file=sys.stderr,
        )
    if args.pns_format in ("wiggz", "both") or args.other_format in ("wiggz", "both"):
        print(
            "[INFO] wig.gz output was explicitly requested and will also be written.",
            file=sys.stderr,
        )

    pns_frag_range = range(args.frag_lower, args.frag_upper + 1)
    if args.pns_mode == "on":
        pns_distributions, pospns_distributions = precompute_distributions(
            pns_frag_range, args.mode_length
        )
    else:
        pns_distributions, pospns_distributions = {}, {}

    output_groups = ALL_OUTPUT_GROUPS if args.split_ww_types else ("all",)

    # Remove stale outputs from previous runs.
    all_track_names = sorted(pns_track_set | other_track_set)
    to_remove = set()
    for group in ALL_OUTPUT_GROUPS:
        prefix = group_output_prefix(args.out_prefix, group)
        to_remove.update({
            f"{prefix}_fragment_summary.txt",
            f"{prefix}_fragment_length_counts.tsv",
            f"{prefix}_nucleosome_regions.bed",
            f"{prefix}_breakpoint_peaks.bed",
            f"{prefix}_dinuc_profile.tsv",
        })
        for track in all_track_names:
            to_remove.add(f"{prefix}_{track}.bw")
            to_remove.add(f"{prefix}_{track}.bedGraph")
            to_remove.add(f"{prefix}_{track}.wig.gz")
    to_remove.add(f"{args.out_prefix}_ww_type_summary.tsv")
    for fname in to_remove:
        if os.path.exists(fname):
            os.remove(fname)

    stats = {
        group: {
            "total_filtered": 0,
            "total_used": 0,
            "unique_bases": 0,
            "length_counts": Counter(),
        }
        for group in output_groups
    }

    type_counts = Counter()
    dinuc_positions = expected_dinuc_profile_positions(
        args.frag_lower, args.frag_upper
    )
    dinuc_accumulators = (
        {group: new_dinuc_accumulator() for group in output_groups}
        if args.dinuc_profile
        else {}
    )

    bigwig_handles = {group: {} for group in output_groups}
    wig_handles = {group: {} for group in output_groups}

    # BigWig headers use the exact BAM contig names and BAM-header order. Include
    # only selected contigs, but retain their original order and full lengths.
    selected_name_set = set(selected_contig_names)
    bigwig_header = [
        (name, int(length))
        for name, length in zip(bamfiles[0].references, bamfiles[0].lengths)
        if name in selected_name_set
    ]

    for group in output_groups:
        prefix = group_output_prefix(args.out_prefix, group)

        bigwig_tracks = []
        if args.pns_format in ("bigwig", "both"):
            bigwig_tracks.extend(args.pns_tracks)
        if args.other_format in ("bigwig", "both"):
            bigwig_tracks.extend(args.other_tracks)
        for track in dict.fromkeys(bigwig_tracks):
            handle = pyBigWig.open(f"{prefix}_{track}.bw", "w")
            handle.addHeader(bigwig_header)
            bigwig_handles[group][track] = handle

        wig_tracks = []
        if args.pns_format in ("wiggz", "both"):
            wig_tracks.extend(args.pns_tracks)
        if args.other_format in ("wiggz", "both"):
            wig_tracks.extend(args.other_tracks)
        for track in dict.fromkeys(wig_tracks):
            wig_handles[group][track] = gzip.open(
                f"{prefix}_{track}.wig.gz", "wt"
            )

    first_region = True
    need_reference_sequence = reference_required

    try:
        total_chunks = len(contigs)
        current_contig = None
        for chunk_index, (
            contig,
            adjusted_start,
            adjusted_end,
            original_start,
            original_end,
        ) in enumerate(contigs, start=1):
            if contig != current_contig:
                current_contig = contig
                print(
                    f"[INFO] Starting {contig} ({chunk_index:,}/{total_chunks:,} chunks).",
                    file=sys.stderr,
                )
            elif chunk_index % 100 == 0 or chunk_index == total_chunks:
                print(
                    f"[INFO] Processed {chunk_index:,}/{total_chunks:,} chunks; "
                    f"currently {contig}:{original_start:,}-{original_end:,}.",
                    file=sys.stderr,
                )
            (
                scores_by_group,
                pns_frag_range,
                fragment_records,
                reference_context,
            ) = score_contig(
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
                pns_mode=args.pns_mode,
                randomize_mode=args.randomize_mode,
                fasta=fasta,
                anchor_prob_start=args.anchor_prob_start,
                max_anchor_tries=args.max_anchor_tries,
                randomize_fallback=args.randomize_fallback,
                dedup_scope=args.dedup_scope,
                split_ww_types=args.split_ww_types,
                need_reference_sequence=need_reference_sequence,
            )

            owned_records = [
                record
                for record in fragment_records
                if original_start <= record[0] < original_end
            ]
            stats["all"]["total_filtered"] += len(owned_records)

            covered_by_group = {
                group: np.zeros(original_end - original_start, dtype=bool)
                for group in output_groups
            }

            for frag_start, frag_end, ww_type in owned_records:
                frag_length = frag_end - frag_start
                if frag_length not in pns_frag_range:
                    continue

                stats["all"]["total_used"] += 1
                stats["all"]["length_counts"][frag_length] += 1

                targets = ["all"]
                if args.split_ww_types:
                    if ww_type in WW_TYPE_GROUPS:
                        type_counts[ww_type] += 1
                        targets.append(ww_type)
                        stats[ww_type]["total_filtered"] += 1
                        stats[ww_type]["total_used"] += 1
                        stats[ww_type]["length_counts"][frag_length] += 1
                    else:
                        type_counts["unclassified"] += 1

                ov_start = max(frag_start, original_start)
                ov_end = min(frag_end, original_end)
                if ov_end > ov_start:
                    for group in targets:
                        covered_by_group[group][
                            ov_start - original_start : ov_end - original_start
                        ] = True

                if args.dinuc_profile:
                    fragment_seq = extract_reference_sequence(
                        fasta=fasta,
                        reference_context=reference_context,
                        seq_start=frag_start,
                        seq_end=frag_end,
                    )
                    for group in targets:
                        add_fragment_to_dinuc_accumulator(
                            accumulator=dinuc_accumulators[group],
                            fragment_seq=fragment_seq,
                            frag_start=frag_start,
                            frag_end=frag_end,
                        )

            for group in output_groups:
                stats[group]["unique_bases"] += int(covered_by_group[group].sum())

            for group in output_groups:
                scores = scores_by_group[group]
                prefix = group_output_prefix(args.out_prefix, group)
                coverage_scores = scores["coverage"][0][2]

                if args.pns_mode == "on":
                    pns_smoothed_scores = scores["pns_smoothed"][0][2]
                    call_and_write_peaks(
                        scores=pns_smoothed_scores,
                        coverage_scores=coverage_scores,
                        adjusted_start=adjusted_start,
                        original_start=original_start,
                        original_end=original_end,
                        contig=contig,
                        out_prefix=prefix,
                        first_region=first_region,
                        peak_type_label="_nucleosome_regions",
                        flip_scores=False,
                        peak_format=args.peak_format,
                        peak_score_scale=args.peak_score_scale,
                        min_region_length=args.min_region_length,
                        max_neg_run=args.max_neg_run,
                    )
                    call_and_write_peaks(
                        scores=-1 * pns_smoothed_scores,
                        coverage_scores=coverage_scores,
                        adjusted_start=adjusted_start,
                        original_start=original_start,
                        original_end=original_end,
                        contig=contig,
                        out_prefix=prefix,
                        first_region=first_region,
                        peak_type_label="_breakpoint_peaks",
                        flip_scores=True,
                        peak_format=args.peak_format,
                        peak_score_scale=args.peak_score_scale,
                        min_region_length=args.min_region_length,
                        max_neg_run=args.max_neg_run,
                    )

                bigwig_tracks = []
                if args.pns_format in ("bigwig", "both"):
                    bigwig_tracks.extend(args.pns_tracks)
                if args.other_format in ("bigwig", "both"):
                    bigwig_tracks.extend(args.other_tracks)
                if bigwig_tracks:
                    write_bigwig_tracks(
                        scores=scores,
                        contig=contig,
                        adjusted_start=adjusted_start,
                        original_start=original_start,
                        original_end=original_end,
                        handles=bigwig_handles[group],
                        tracks=list(dict.fromkeys(bigwig_tracks)),
                    )

                wig_tracks = []
                if args.pns_format in ("wiggz", "both"):
                    wig_tracks.extend(args.pns_tracks)
                if args.other_format in ("wiggz", "both"):
                    wig_tracks.extend(args.other_tracks)
                if wig_tracks:
                    write_wig_gz_tracks(
                        scores=scores,
                        contig=contig,
                        adjusted_start=adjusted_start,
                        original_start=original_start,
                        original_end=original_end,
                        handles=wig_handles[group],
                        tracks=wig_tracks,
                    )

            first_region = False

        for group in output_groups:
            prefix = group_output_prefix(args.out_prefix, group)
            write_fragment_outputs(
                out_prefix=prefix,
                total_fragments_filtered_all=stats[group]["total_filtered"],
                total_fragments_used_in_range=stats[group]["total_used"],
                unique_bases_covered_by_used=stats[group]["unique_bases"],
                length_counts=stats[group]["length_counts"],
                dedup_scope=args.dedup_scope,
                max_duplicates=args.max_duplicates,
            )

            if args.dinuc_profile:
                write_dinuc_profile(
                    out_path=f"{prefix}_dinuc_profile.tsv",
                    accumulator=dinuc_accumulators[group],
                    positions=dinuc_positions,
                    fraction=args.dinuc_fraction,
                )

        if args.split_ww_types:
            write_ww_type_summary(
                out_prefix=args.out_prefix,
                type_counts=type_counts,
                total_in_range=stats["all"]["total_used"],
            )

    finally:
        for group_handles in bigwig_handles.values():
            for handle in group_handles.values():
                try:
                    handle.close()
                except Exception:
                    pass
        for group_handles in wig_handles.values():
            for handle in group_handles.values():
                try:
                    handle.close()
                except Exception:
                    pass
        for bam in bamfiles:
            try:
                bam.close()
            except Exception:
                pass
        if fasta is not None:
            try:
                fasta.close()
            except Exception:
                pass

    if args.dinuc_profile:
        for group in output_groups:
            acc = dinuc_accumulators[group]
            print(
                f"[INFO] {group} dinucleotide fragments used: "
                f"{acc['fragments_used']:,}; skipped: {acc['fragments_skipped']:,}",
                file=sys.stderr,
            )

    return 0


if __name__ == "__main__":
    sys.exit(main())
