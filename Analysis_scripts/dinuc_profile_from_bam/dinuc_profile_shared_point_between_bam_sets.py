#!/usr/bin/env python3
"""
Generate observed dinucleotide profiles for fragments from two different BAM sets
that share a point: either a dyad/centre point or a fragment end point.

This is the between-BAM-set version of the low-memory shared-point profiler.
It treats signal A and signal B separately:

  A = fragments of --length_a from --bam_a
  B = fragments of --length_b from --bam_b

For each chromosome, it:
  1. scans all A BAMs and all B BAMs to build point->fragment-coordinate maps;
  2. intersects the A and B point keys;
  3. rescans A and B BAMs and profiles only fragments whose coordinates were
     matched through those shared points;
  4. clears chromosome-local match state before moving to the next chromosome.

Outputs are written separately for A and B, for example:

  <prefix>_<label_a>_len147_shared_dyad.tsv
  <prefix>_<label_b>_len167_shared_dyad.tsv
  <prefix>_summary.tsv

For --match-mode ends, outputs are shared_any/shared_left/shared_right for each
set. shared_any is the union of shared_left and shared_right.

Dinucleotide position is reported as the dinucleotide start relative to
frag_start + floor(fragment_length / 2), matching the previous profiler.
"""

import argparse
import glob
import os
import random
import sys
from collections import defaultdict

try:
    import pysam
except ImportError as exc:
    raise SystemExit("This script requires pysam. Install it with: pip install pysam") from exc


DINUCS = [
    "AA", "AC", "AG", "AT",
    "CA", "CC", "CG", "CT",
    "GA", "GC", "GG", "GT",
    "TA", "TC", "TG", "TT",
]

WW_DINUCS = {"AA", "AT", "TA", "TT"}
SS_DINUCS = {"CC", "CG", "GC", "GG"}
VALID_BASES = {"A", "C", "G", "T"}

CIGAR_SOFT_HARD_OR_PAD = {4, 5, 6}  # S, H, P
BASE_SUBSETS = ("any", "left", "right")
DYAD_SUBSET = "dyad"
SIGNALS = ("A", "B")


def sanitize_filename(name):
    return str(name).replace("/", "_").replace("\\", "_").replace(" ", "_")


def has_softclip_or_hardclip_or_padding(cigartuples):
    if not cigartuples:
        return False

    for op, _length in cigartuples:
        if op in CIGAR_SOFT_HARD_OR_PAD:
            return True

    return False


def expand_bam_inputs(bam_inputs):
    """
    Expand one or more BAM paths/globs.

    Works with either:
      --bam_a '/path/*.bam'
    or:
      --bam_a /path/a.bam /path/b.bam
    """
    bam_files = []

    for item in bam_inputs:
        # Also allow one quoted space-separated string, matching the DCC script style.
        for part in str(item).split():
            matches = sorted(glob.glob(part))
            if matches:
                bam_files.extend(matches)
            else:
                bam_files.append(part)

    bam_files = list(dict.fromkeys(bam_files))
    missing = [b for b in bam_files if not os.path.exists(b)]

    if missing:
        raise FileNotFoundError(
            "These BAM files/patterns were not found:\n  " + "\n  ".join(missing)
        )

    return bam_files


def expand_chroms(chrom_spec):
    chrom_spec = str(chrom_spec).replace(" ", "")

    if chrom_spec.lower() == "all":
        return [str(i) for i in range(1, 23)] + ["X", "Y"]

    chroms = []

    for token in chrom_spec.split(","):
        if not token:
            continue

        if token.startswith("chr"):
            token = token[3:]

        if "-" in token:
            start, end = token.split("-", 1)

            if not start.isdigit() or not end.isdigit():
                raise ValueError(f"Invalid chromosome range: {token}")

            chroms.extend([str(i) for i in range(int(start), int(end) + 1)])
        else:
            chroms.append(token)

    return chroms


def matching_contig(chrom, references):
    if chrom in references:
        return chrom

    if chrom.startswith("chr"):
        no_chr = chrom[3:]
        if no_chr in references:
            return no_chr
        if no_chr == "M" and "MT" in references:
            return "MT"
    else:
        with_chr = "chr" + chrom
        if with_chr in references:
            return with_chr
        if chrom == "MT" and "chrM" in references:
            return "chrM"
        if chrom == "M" and "MT" in references:
            return "MT"

    return None


def read_passes_filters(
    aln,
    mapq,
    include_duplicates,
    allow_improper_pairs,
    allow_softclipped,
):
    if aln.is_unmapped or aln.mate_is_unmapped:
        return False
    if aln.is_secondary or aln.is_supplementary:
        return False
    if aln.is_qcfail:
        return False
    if not include_duplicates and aln.is_duplicate:
        return False
    if not allow_improper_pairs and not aln.is_proper_pair:
        return False
    if not allow_softclipped and has_softclip_or_hardclip_or_padding(aln.cigartuples):
        return False
    if aln.mapping_quality < mapq:
        return False
    if aln.reference_id != aln.next_reference_id:
        return False
    if aln.template_length == 0:
        return False

    return True


def inferred_fragment_start_and_length(aln1, aln2):
    """Return the template start and length for a paired-end fragment."""
    if aln1.template_length > 0:
        frag_start = aln1.reference_start
        frag_len = aln1.template_length
    elif aln2.template_length > 0:
        frag_start = aln2.reference_start
        frag_len = aln2.template_length
    else:
        frag_start = min(aln1.reference_start, aln2.reference_start)
        frag_end = max(aln1.reference_end, aln2.reference_end)
        frag_len = frag_end - frag_start

    return int(frag_start), abs(int(frag_len))


def dyad_integer_values(frag_start, frag_len, dyad_match_mode):
    """
    Return dyad keys as integers for chromosome-local matching.

    Modes:
      exact: half-base-safe value 2*centre = 2*start + len - 1
      floor: integer floor of geometric centre
      ceil: integer ceil of geometric centre
      split: odd lengths get one centre base; even lengths get both middle bases
    """
    centre_2x = (2 * frag_start) + frag_len - 1

    if dyad_match_mode == "exact":
        return (centre_2x,)

    floor_key = centre_2x // 2
    ceil_key = (centre_2x + 1) // 2

    if dyad_match_mode == "floor":
        return (floor_key,)

    if dyad_match_mode == "ceil":
        return (ceil_key,)

    if dyad_match_mode == "split":
        if floor_key == ceil_key:
            return (floor_key,)
        return (floor_key, ceil_key)

    raise ValueError(f"unknown dyad match mode: {dyad_match_mode}")


def add_aligned_read_bases_to_fragment(aln, frag_start, frag_len, bases, base_quals, min_baseq):
    """
    Fill fragment bases from one aligned read.

    The resulting fragment sequence is in reference/genomic left-to-right
    orientation. Overlapping mates are resolved by keeping the base with the
    higher base quality.
    """
    seq = aln.query_sequence
    if seq is None:
        return

    quals = aln.query_qualities

    for query_pos, ref_pos in aln.get_aligned_pairs(matches_only=False):
        if query_pos is None or ref_pos is None:
            continue

        frag_i = ref_pos - frag_start
        if frag_i < 0 or frag_i >= frag_len:
            continue

        base = seq[query_pos].upper()
        if base not in VALID_BASES:
            continue

        qual = quals[query_pos] if quals is not None else 0
        if qual < min_baseq:
            continue

        if bases[frag_i] is None or qual > base_quals[frag_i]:
            bases[frag_i] = base
            base_quals[frag_i] = qual


def reconstruct_fragment_sequence(aln1, aln2, frag_start, frag_len, min_baseq):
    bases = [None] * frag_len
    base_quals = [-1] * frag_len

    add_aligned_read_bases_to_fragment(
        aln=aln1,
        frag_start=frag_start,
        frag_len=frag_len,
        bases=bases,
        base_quals=base_quals,
        min_baseq=min_baseq,
    )
    add_aligned_read_bases_to_fragment(
        aln=aln2,
        frag_start=frag_start,
        frag_len=frag_len,
        bases=bases,
        base_quals=base_quals,
        min_baseq=min_baseq,
    )

    return bases


def get_centered_window_offsets(fragment_length, window_length):
    """
    Return 0-based start/end offsets for a centred internal window.

    If fragment_length - window_length is odd, the extra base is left on the
    right side of the fragment. For example, length 168 with window 147 gives
    offsets 10:157.
    """
    if window_length < 2:
        raise ValueError("window length must be >= 2 for dinucleotide profiles")
    if window_length > fragment_length:
        raise ValueError("window length cannot be larger than fragment length")

    start = (fragment_length - window_length) // 2
    end = start + window_length

    return start, end


def make_profile(fragment_length, core_start_offset, core_end_offset, profile_start_offset, profile_end_offset):
    profile_length = profile_end_offset - profile_start_offset

    return {
        "fragment_length": fragment_length,
        "core_start_offset": core_start_offset,
        "core_end_offset": core_end_offset,
        "profile_start_offset": profile_start_offset,
        "profile_end_offset": profile_end_offset,
        "profile_length": profile_length,
        "dinuc_counts_by_position": [defaultdict(int) for _ in range(profile_length - 1)],
        "valid_dinuc_opportunities": [0] * (profile_length - 1),
        "candidate_fragments_seen": 0,
        "fragments_used": 0,
        "fragments_with_no_valid_dinucs": 0,
        "fragments_skipped_incomplete": 0,
        "fragments_skipped_reference_bounds": 0,
    }


def add_fragment_dinucs_to_profile(bases, profile, require_complete_fragment):
    """Add dinucleotide counts for one profiled interval."""
    if require_complete_fragment and any(base not in VALID_BASES for base in bases):
        profile["fragments_skipped_incomplete"] += 1
        return None

    added_positions = 0
    dinuc_counts_by_position = profile["dinuc_counts_by_position"]
    valid_dinuc_opportunities = profile["valid_dinuc_opportunities"]

    for i in range(len(bases) - 1):
        b1 = bases[i]
        b2 = bases[i + 1]

        if b1 not in VALID_BASES or b2 not in VALID_BASES:
            continue

        dinuc = b1 + b2
        dinuc_counts_by_position[i][dinuc] += 1
        valid_dinuc_opportunities[i] += 1
        added_positions += 1

    if added_positions == 0:
        profile["fragments_with_no_valid_dinucs"] += 1
    else:
        profile["fragments_used"] += 1

    return added_positions


def choose_shift_delta(shift_bp, shift_direction, rng):
    if shift_bp == 0:
        return 0, "none"
    if shift_direction == "plus":
        return shift_bp, "plus"
    if shift_direction == "minus":
        return -shift_bp, "minus"
    if shift_direction == "random":
        if rng.random() < 0.5:
            return -shift_bp, "minus"
        return shift_bp, "plus"

    raise ValueError(f"unknown shift direction: {shift_direction}")


def fetch_reference_bases(fasta, fasta_chrom, start_0based, end_0based):
    if start_0based < 0:
        return None

    chrom_len = fasta.get_reference_length(fasta_chrom)
    if end_0based > chrom_len:
        return None

    seq = fasta.fetch(fasta_chrom, start_0based, end_0based).upper()
    if len(seq) != end_0based - start_0based:
        return None

    return list(seq)


def iter_fragment_pairs_from_bam_chrom_length(
    bam,
    bam_chrom,
    fragment_length,
    mapq,
    include_duplicates,
    allow_improper_pairs,
    allow_softclipped,
    max_duplicates,
    coord_counts,
):
    """
    Iterate exact-length paired-end fragments from one BAM contig.

    This scans the BAM once and pairs reads by query name, avoiding a FASTA file
    and avoiding pysam.mate() calls for every read.
    """
    pending = {}

    try:
        iterator = bam.fetch(bam_chrom)
    except ValueError:
        return

    for aln in iterator:
        if not read_passes_filters(
            aln=aln,
            mapq=mapq,
            include_duplicates=include_duplicates,
            allow_improper_pairs=allow_improper_pairs,
            allow_softclipped=allow_softclipped,
        ):
            continue

        if abs(int(aln.template_length)) != fragment_length:
            continue

        qname = aln.query_name
        if qname not in pending:
            pending[qname] = aln
            continue

        mate = pending.pop(qname)
        if aln.reference_id != mate.reference_id:
            continue

        frag_start, frag_len = inferred_fragment_start_and_length(aln, mate)
        if frag_len != fragment_length:
            continue

        frag_end = frag_start + frag_len

        if max_duplicates >= 0:
            key = (bam_chrom, frag_start, frag_end)
            # max_duplicates=0 keeps one fragment per coordinate.
            # max_duplicates=N keeps N+1 observations, matching the earlier script.
            if coord_counts[key] > max_duplicates:
                continue
            coord_counts[key] += 1

        yield aln, mate, frag_start, frag_end, frag_len


def build_profile_specs_by_signal(length_a, length_b, window_length, extend_left_bp, extend_right_bp):
    specs = {}
    for signal, fragment_length in (("A", length_a), ("B", length_b)):
        core_profile_length = fragment_length if window_length is None else window_length

        try:
            core_start_offset, core_end_offset = get_centered_window_offsets(
                fragment_length=fragment_length,
                window_length=core_profile_length,
            )
        except ValueError as e:
            raise ValueError(f"for signal {signal} length {fragment_length}, --window-length invalid: {e}")

        profile_start_offset = core_start_offset - extend_left_bp
        profile_end_offset = core_end_offset + extend_right_bp
        profile_length = profile_end_offset - profile_start_offset

        if profile_length < 2:
            raise ValueError(
                f"for signal {signal} length {fragment_length}, profiled interval is < 2 bp after extension"
            )

        specs[signal] = {
            "fragment_length": fragment_length,
            "core_profile_length": core_profile_length,
            "core_start_offset": core_start_offset,
            "core_end_offset": core_end_offset,
            "profile_start_offset": profile_start_offset,
            "profile_end_offset": profile_end_offset,
            "profile_length": profile_length,
        }

    return specs


def make_profiles_by_signal(profile_specs, subset_names):
    profiles = {}

    for signal in SIGNALS:
        spec = profile_specs[signal]
        for subset in subset_names:
            profiles[(signal, subset)] = make_profile(
                fragment_length=spec["fragment_length"],
                core_start_offset=spec["core_start_offset"],
                core_end_offset=spec["core_end_offset"],
                profile_start_offset=spec["profile_start_offset"],
                profile_end_offset=spec["profile_end_offset"],
            )

    return profiles


def get_profile_bases_for_fragment(
    aln1,
    aln2,
    frag_start,
    frag_len,
    spec,
    min_baseq,
    use_reference_sequence,
    fasta,
    fasta_chrom,
    shift_delta,
):
    if use_reference_sequence:
        ref_start = frag_start + spec["profile_start_offset"] + shift_delta
        ref_end = frag_start + spec["profile_end_offset"] + shift_delta
        return fetch_reference_bases(
            fasta=fasta,
            fasta_chrom=fasta_chrom,
            start_0based=ref_start,
            end_0based=ref_end,
        )

    bases = reconstruct_fragment_sequence(
        aln1=aln1,
        aln2=aln2,
        frag_start=frag_start,
        frag_len=frag_len,
        min_baseq=min_baseq,
    )

    return bases[spec["profile_start_offset"]:spec["profile_end_offset"]]


def empty_key_maps(match_mode):
    maps = {
        "coords": set(),
    }
    if match_mode == "dyad":
        maps["dyad_to_coords"] = defaultdict(set)
    elif match_mode == "ends":
        maps["left_to_coords"] = defaultdict(set)
        maps["right_to_coords"] = defaultdict(set)
    else:
        raise ValueError(f"unknown match mode: {match_mode}")
    return maps


def scan_bam_set_for_chrom_keys(
    signal,
    bam_files,
    requested_chrom,
    fragment_length,
    match_mode,
    dyad_match_mode,
    mapq,
    include_duplicates,
    allow_improper_pairs,
    allow_softclipped,
    max_duplicates,
    total_pairs_by_signal,
):
    """
    Scan one signal's BAM set for one requested chromosome and return maps of
    shared-point keys to fragment coordinates.

    Coordinates are stored without the chromosome because this is called for one
    chromosome at a time. Matching is performed only within that chromosome.
    """
    maps = empty_key_maps(match_mode)
    bam_chroms_for_files = []
    missing_files = []

    for bam_path in bam_files:
        with pysam.AlignmentFile(bam_path, "rb") as bam:
            bam_refs = set(bam.references)
            bam_chrom = matching_contig(requested_chrom, bam_refs)
            if bam_chrom is None:
                missing_files.append(os.path.basename(bam_path))
                continue

            bam_chroms_for_files.append((bam_path, bam_chrom))
            coord_counts = defaultdict(int)

            for _aln1, _aln2, frag_start, frag_end, _frag_len in iter_fragment_pairs_from_bam_chrom_length(
                bam=bam,
                bam_chrom=bam_chrom,
                fragment_length=fragment_length,
                mapq=mapq,
                include_duplicates=include_duplicates,
                allow_improper_pairs=allow_improper_pairs,
                allow_softclipped=allow_softclipped,
                max_duplicates=max_duplicates,
                coord_counts=coord_counts,
            ):
                coord_id = (frag_start, frag_end)
                maps["coords"].add(coord_id)
                total_pairs_by_signal[signal] += 1

                if match_mode == "dyad":
                    for dyad_key in dyad_integer_values(
                        frag_start=frag_start,
                        frag_len=fragment_length,
                        dyad_match_mode=dyad_match_mode,
                    ):
                        maps["dyad_to_coords"][dyad_key].add(coord_id)
                else:
                    maps["left_to_coords"][frag_start].add(coord_id)
                    # Use the half-open end coordinate. Equality is equivalent to sharing
                    # the rightmost aligned base because all right ends are shifted by 1.
                    maps["right_to_coords"][frag_end].add(coord_id)

                if total_pairs_by_signal[signal] % 1000000 == 0:
                    print(
                        f"[MATCH] signal {signal} length {fragment_length} fragments seen: "
                        f"{total_pairs_by_signal[signal]:,}",
                        file=sys.stderr,
                    )

    return maps, bam_chroms_for_files, missing_files


def match_key_maps_between_sets(maps_a, maps_b, match_mode):
    if match_mode == "dyad":
        subset_names = (DYAD_SUBSET,)
        matched = {signal: {DYAD_SUBSET: set()} for signal in SIGNALS}
        shared_dyads = set(maps_a["dyad_to_coords"]).intersection(maps_b["dyad_to_coords"])

        for dyad_key in shared_dyads:
            matched["A"][DYAD_SUBSET].update(maps_a["dyad_to_coords"][dyad_key])
            matched["B"][DYAD_SUBSET].update(maps_b["dyad_to_coords"][dyad_key])

        return subset_names, matched

    if match_mode == "ends":
        subset_names = BASE_SUBSETS
        matched = {signal: {subset: set() for subset in subset_names} for signal in SIGNALS}

        shared_lefts = set(maps_a["left_to_coords"]).intersection(maps_b["left_to_coords"])
        shared_rights = set(maps_a["right_to_coords"]).intersection(maps_b["right_to_coords"])

        for left in shared_lefts:
            matched["A"]["left"].update(maps_a["left_to_coords"][left])
            matched["B"]["left"].update(maps_b["left_to_coords"][left])

        for right in shared_rights:
            matched["A"]["right"].update(maps_a["right_to_coords"][right])
            matched["B"]["right"].update(maps_b["right_to_coords"][right])

        for signal in SIGNALS:
            matched[signal]["any"] = matched[signal]["left"] | matched[signal]["right"]

        return subset_names, matched

    raise ValueError(f"unknown match mode: {match_mode}")


def profile_signal_chrom_matched_fragments(
    signal,
    bam_chroms_for_files,
    fragment_length,
    subset_names,
    matched,
    profile_specs,
    profiles,
    mapq,
    min_baseq,
    include_duplicates,
    allow_improper_pairs,
    allow_softclipped,
    max_duplicates,
    require_complete_fragment,
    max_fragments,
    use_reference_sequence,
    fasta,
    fasta_refs,
    shift_bp,
    shift_direction,
    rng,
    shift_direction_counts,
    progress,
):
    """Rescan one signal's BAMs for one chromosome and profile matched coordinates."""
    spec = profile_specs[signal]

    for bam_path, bam_chrom in bam_chroms_for_files:
        fasta_chrom = None
        if use_reference_sequence:
            fasta_chrom = matching_contig(bam_chrom, fasta_refs)
            if fasta_chrom is None:
                print(
                    f"[WARNING] {bam_chrom} is present in {os.path.basename(bam_path)} "
                    "but not found in FASTA; skipped",
                    file=sys.stderr,
                )
                continue

        with pysam.AlignmentFile(bam_path, "rb") as bam:
            coord_counts = defaultdict(int)

            for aln1, aln2, frag_start, frag_end, frag_len in iter_fragment_pairs_from_bam_chrom_length(
                bam=bam,
                bam_chrom=bam_chrom,
                fragment_length=fragment_length,
                mapq=mapq,
                include_duplicates=include_duplicates,
                allow_improper_pairs=allow_improper_pairs,
                allow_softclipped=allow_softclipped,
                max_duplicates=max_duplicates,
                coord_counts=coord_counts,
            ):
                coord_id = (frag_start, frag_end)
                matched_subsets = [
                    subset for subset in subset_names
                    if coord_id in matched[signal][subset]
                ]

                if not matched_subsets:
                    continue

                shift_delta, shift_label = choose_shift_delta(
                    shift_bp=shift_bp,
                    shift_direction=shift_direction,
                    rng=rng,
                )
                shift_direction_counts[shift_label] += 1

                profile_bases = get_profile_bases_for_fragment(
                    aln1=aln1,
                    aln2=aln2,
                    frag_start=frag_start,
                    frag_len=frag_len,
                    spec=spec,
                    min_baseq=min_baseq,
                    use_reference_sequence=use_reference_sequence,
                    fasta=fasta,
                    fasta_chrom=fasta_chrom,
                    shift_delta=shift_delta,
                )

                progress["profiled_matching_fragments"] += 1

                for subset in matched_subsets:
                    profile = profiles[(signal, subset)]
                    profile["candidate_fragments_seen"] += 1

                    if profile_bases is None:
                        profile["fragments_skipped_reference_bounds"] += 1
                        continue

                    add_fragment_dinucs_to_profile(
                        bases=profile_bases,
                        profile=profile,
                        require_complete_fragment=require_complete_fragment,
                    )

                if progress["profiled_matching_fragments"] % 100000 == 0:
                    print(
                        f"[PROFILE] reconstructed/fetched matching fragments: "
                        f"{progress['profiled_matching_fragments']:,}",
                        file=sys.stderr,
                    )

                if max_fragments is not None and progress["profiled_matching_fragments"] >= max_fragments:
                    return True

    return False


def process_bam_sets_lowmem(
    bam_files_a,
    length_a,
    label_a,
    bam_files_b,
    length_b,
    label_b,
    requested_chroms,
    match_mode,
    dyad_match_mode,
    profile_specs,
    mapq,
    min_baseq,
    include_duplicates,
    allow_improper_pairs,
    allow_softclipped,
    max_duplicates,
    require_complete_fragment,
    max_fragments,
    fasta_path,
    shift_bp,
    shift_direction,
    random_seed,
):
    subset_names = (DYAD_SUBSET,) if match_mode == "dyad" else BASE_SUBSETS
    profiles = make_profiles_by_signal(profile_specs, subset_names)

    signal_labels = {"A": label_a, "B": label_b}
    signal_lengths = {"A": length_a, "B": length_b}
    signal_bams = {"A": bam_files_a, "B": bam_files_b}

    total_pairs_by_signal = {signal: 0 for signal in SIGNALS}
    unique_coords_count_by_signal = {signal: 0 for signal in SIGNALS}
    matched_coords_count = {
        signal: {subset: 0 for subset in subset_names}
        for signal in SIGNALS
    }

    use_reference_sequence = fasta_path is not None
    fasta = None
    fasta_refs = set()

    if use_reference_sequence:
        try:
            fasta = pysam.FastaFile(fasta_path)
        except Exception as e:
            sys.exit(f"ERROR: could not open --fasta {fasta_path!r}: {e}")

        fasta_refs = set(fasta.references)
        print(f"[INFO] reference FASTA: {fasta_path}", file=sys.stderr)
        print(
            "[INFO] using reference FASTA sequence for the profiled interval; "
            "--min-baseq does not apply in this mode",
            file=sys.stderr,
        )

    rng = random.Random(random_seed)
    shift_direction_counts = defaultdict(int)
    progress = {"profiled_matching_fragments": 0}

    try:
        for requested_chrom in requested_chroms:
            print(f"[MATCH] Scanning chromosome {requested_chrom}", file=sys.stderr)

            maps_by_signal = {}
            bam_chroms_by_signal = {}

            for signal in SIGNALS:
                maps, bam_chroms_for_files, missing_files = scan_bam_set_for_chrom_keys(
                    signal=signal,
                    bam_files=signal_bams[signal],
                    requested_chrom=requested_chrom,
                    fragment_length=signal_lengths[signal],
                    match_mode=match_mode,
                    dyad_match_mode=dyad_match_mode,
                    mapq=mapq,
                    include_duplicates=include_duplicates,
                    allow_improper_pairs=allow_improper_pairs,
                    allow_softclipped=allow_softclipped,
                    max_duplicates=max_duplicates,
                    total_pairs_by_signal=total_pairs_by_signal,
                )
                maps_by_signal[signal] = maps
                bam_chroms_by_signal[signal] = bam_chroms_for_files
                unique_coords_count_by_signal[signal] += len(maps["coords"])

                if missing_files:
                    print(
                        f"[WARNING] signal {signal} ({signal_labels[signal]}), chromosome {requested_chrom} "
                        f"not found in {len(missing_files)} BAM(s); skipped for those BAMs",
                        file=sys.stderr,
                    )

                print(
                    f"[MATCH] signal {signal} ({signal_labels[signal]}), chromosome {requested_chrom}: "
                    f"{len(maps['coords']):,} unique coordinates",
                    file=sys.stderr,
                )

            subset_names, matched = match_key_maps_between_sets(
                maps_a=maps_by_signal["A"],
                maps_b=maps_by_signal["B"],
                match_mode=match_mode,
            )

            for signal in SIGNALS:
                for subset in subset_names:
                    matched_coords_count[signal][subset] += len(matched[signal][subset])

            match_parts = []
            for signal in SIGNALS:
                for subset in subset_names:
                    match_parts.append(
                        f"{signal}_{subset}={len(matched[signal][subset]):,}"
                    )
            print(
                f"[MATCH] chromosome {requested_chrom}: " + "; ".join(match_parts),
                file=sys.stderr,
            )

            if not any(
                len(matched[signal][subset]) > 0
                for signal in SIGNALS
                for subset in subset_names
            ):
                continue

            print(f"[PROFILE] Profiling chromosome {requested_chrom}", file=sys.stderr)

            for signal in SIGNALS:
                reached_limit = profile_signal_chrom_matched_fragments(
                    signal=signal,
                    bam_chroms_for_files=bam_chroms_by_signal[signal],
                    fragment_length=signal_lengths[signal],
                    subset_names=subset_names,
                    matched=matched,
                    profile_specs=profile_specs,
                    profiles=profiles,
                    mapq=mapq,
                    min_baseq=min_baseq,
                    include_duplicates=include_duplicates,
                    allow_improper_pairs=allow_improper_pairs,
                    allow_softclipped=allow_softclipped,
                    max_duplicates=max_duplicates,
                    require_complete_fragment=require_complete_fragment,
                    max_fragments=max_fragments,
                    use_reference_sequence=use_reference_sequence,
                    fasta=fasta,
                    fasta_refs=fasta_refs,
                    shift_bp=shift_bp,
                    shift_direction=shift_direction,
                    rng=rng,
                    shift_direction_counts=shift_direction_counts,
                    progress=progress,
                )

                if reached_limit:
                    print(
                        f"[INFO] stopping early because --max-fragments {max_fragments} was reached",
                        file=sys.stderr,
                    )
                    return (
                        profiles,
                        shift_direction_counts,
                        subset_names,
                        total_pairs_by_signal,
                        unique_coords_count_by_signal,
                        matched_coords_count,
                        True,
                    )

    finally:
        if fasta is not None:
            fasta.close()

    return (
        profiles,
        shift_direction_counts,
        subset_names,
        total_pairs_by_signal,
        unique_coords_count_by_signal,
        matched_coords_count,
        False,
    )


def write_profile_tsv(profile, out_path, fraction):
    fragment_length = profile["fragment_length"]
    profile_length = profile["profile_length"]
    n_dinuc_positions = profile_length - 1
    centre_offset = fragment_length // 2
    profile_start_offset = profile["profile_start_offset"]

    multiplier = 1.0 if fraction else 100.0
    value_suffix = "frac" if fraction else "pct"

    header = ["position", "n_valid"]
    header.extend([f"{dinuc}_{value_suffix}" for dinuc in DINUCS])
    header.extend([f"WW_{value_suffix}", f"SS_{value_suffix}"])

    dinuc_counts_by_position = profile["dinuc_counts_by_position"]
    valid_dinuc_opportunities = profile["valid_dinuc_opportunities"]

    with open(out_path, "w") as out:
        out.write("\t".join(header) + "\n")

        for dinuc_i in range(n_dinuc_positions):
            rel_pos = (profile_start_offset + dinuc_i) - centre_offset
            n_valid = valid_dinuc_opportunities[dinuc_i]
            counts = dinuc_counts_by_position[dinuc_i]
            row = [str(rel_pos), str(n_valid)]

            if n_valid == 0:
                row.extend(["NaN"] * len(DINUCS))
                row.extend(["NaN", "NaN"])
            else:
                for dinuc in DINUCS:
                    value = (counts[dinuc] / n_valid) * multiplier
                    row.append(f"{value:.8g}")

                ww_count = sum(counts[dinuc] for dinuc in WW_DINUCS)
                ss_count = sum(counts[dinuc] for dinuc in SS_DINUCS)
                row.append(f"{(ww_count / n_valid) * multiplier:.8g}")
                row.append(f"{(ss_count / n_valid) * multiplier:.8g}")

            out.write("\t".join(row) + "\n")


def write_summary_tsv(
    out_path,
    signal_labels,
    signal_lengths,
    subset_names,
    total_pairs_by_signal,
    unique_coords_count_by_signal,
    matched_coords_count,
    profiles,
    output_paths,
    stopped_early,
):
    header = [
        "signal",
        "label",
        "length",
        "subset",
        "fragments_seen_after_filters",
        "unique_fragment_coordinates_in_signal",
        "unique_matched_coordinates",
        "candidate_fragments_seen_in_profile_pass",
        "fragments_used",
        "fragments_skipped_incomplete",
        "fragments_skipped_reference_bounds",
        "fragments_with_no_valid_dinucs",
        "core_start_offset_0based",
        "core_end_offset_0based_exclusive",
        "profile_start_offset_0based",
        "profile_end_offset_0based_exclusive",
        "profile_length",
        "stopped_early_by_max_fragments",
        "output_tsv",
    ]

    with open(out_path, "w") as out:
        out.write("\t".join(header) + "\n")

        for signal in SIGNALS:
            for subset in subset_names:
                profile = profiles[(signal, subset)]
                row = [
                    signal,
                    signal_labels[signal],
                    str(signal_lengths[signal]),
                    subset,
                    str(total_pairs_by_signal[signal]),
                    str(unique_coords_count_by_signal[signal]),
                    str(matched_coords_count[signal][subset]),
                    str(profile["candidate_fragments_seen"]),
                    str(profile["fragments_used"]),
                    str(profile["fragments_skipped_incomplete"]),
                    str(profile["fragments_skipped_reference_bounds"]),
                    str(profile["fragments_with_no_valid_dinucs"]),
                    str(profile["core_start_offset"]),
                    str(profile["core_end_offset"]),
                    str(profile["profile_start_offset"]),
                    str(profile["profile_end_offset"]),
                    str(profile["profile_length"]),
                    str(bool(stopped_early)),
                    output_paths[(signal, subset)],
                ]
                out.write("\t".join(row) + "\n")


def main():
    parser = argparse.ArgumentParser(
        description=(
            "Generate observed dinucleotide profiles for fragments from two different "
            "BAM sets that share either dyads or fragment ends. Signal A uses --bam_a "
            "+ --length_a, and signal B uses --bam_b + --length_b."
        )
    )

    parser.add_argument(
        "--bam_a",
        nargs="+",
        required=True,
        help="BAM file(s)/glob(s) for signal A. Quote wildcards, e.g. '/data/A/*.bam'.",
    )
    parser.add_argument(
        "--length_a",
        type=int,
        required=True,
        help="Exact paired-end fragment length for signal A, e.g. 147.",
    )
    parser.add_argument(
        "--label_a",
        default="A",
        help="Short output label for signal A. Default: A.",
    )

    parser.add_argument(
        "--bam_b",
        nargs="+",
        required=True,
        help="BAM file(s)/glob(s) for signal B. Quote wildcards, e.g. '/data/B/*.bam'.",
    )
    parser.add_argument(
        "--length_b",
        type=int,
        required=True,
        help="Exact paired-end fragment length for signal B, e.g. 167.",
    )
    parser.add_argument(
        "--label_b",
        default="B",
        help="Short output label for signal B. Default: B.",
    )

    parser.add_argument(
        "--match-mode",
        choices=["ends", "dyad"],
        default="ends",
        help=(
            "Which shared point to require between A and B. 'ends' profiles fragments "
            "sharing left/start or right/end coordinates and writes shared_any, "
            "shared_left, and shared_right outputs. 'dyad' profiles fragments sharing "
            "the same dyad/centre and writes shared_dyad outputs. Default: ends."
        ),
    )
    parser.add_argument(
        "--dyad-match-mode",
        choices=["exact", "floor", "ceil", "split"],
        default="split",
        help=(
            "How to define shared dyads when --match-mode dyad is used. 'split' "
            "assigns odd lengths to one centre base and even lengths to both middle bases. "
            "'exact' uses half-base-safe geometric centres. Default: split."
        ),
    )
    parser.add_argument(
        "--window-length",
        type=int,
        default=None,
        help=(
            "Optional centred window length to profile within each fragment before extension/shift. "
            "For example, with --length_a 147 --length_b 167 --window-length 147, "
            "147 bp fragments are profiled as-is and 167 bp fragments use their centred 147 bp window."
        ),
    )
    parser.add_argument(
        "--extend-bp",
        type=int,
        default=0,
        help="Extend the profiled interval by this many bp on both sides. Requires --fasta when >0.",
    )
    parser.add_argument(
        "--extend-left-bp",
        "--extend-start-bp",
        dest="extend_left_bp",
        type=int,
        default=0,
        help="Extend only the genomic-left/start side. Requires --fasta when >0.",
    )
    parser.add_argument(
        "--extend-right-bp",
        "--extend-end-bp",
        dest="extend_right_bp",
        type=int,
        default=0,
        help="Extend only the genomic-right/end side. Requires --fasta when >0.",
    )
    parser.add_argument(
        "--shift-bp",
        type=int,
        default=0,
        help="Shift the profiled interval before fetching reference bases. Requires --fasta when >0.",
    )
    parser.add_argument(
        "--shift-direction",
        choices=["plus", "minus", "random"],
        default="plus",
        help="Direction for --shift-bp. Default: plus.",
    )
    parser.add_argument(
        "--random-seed",
        type=int,
        default=1,
        help="Random seed used when --shift-direction random. Default: 1.",
    )
    parser.add_argument(
        "--fasta",
        default=None,
        help=(
            "Reference genome FASTA. Required for shift/extension. When supplied, "
            "profiled bases are fetched from the reference genome instead of read bases."
        ),
    )
    parser.add_argument(
        "--chroms",
        default="all",
        help="Chromosomes to analyse, e.g. all, 1-22,X,Y, 20, or chr20. Default: all.",
    )
    parser.add_argument(
        "--out-prefix",
        required=True,
        help="Output prefix for profile TSVs and summary TSV.",
    )
    parser.add_argument(
        "--mapq",
        type=int,
        default=30,
        help="Minimum MAPQ. Default: 30.",
    )
    parser.add_argument(
        "--min-baseq",
        type=int,
        default=0,
        help="Minimum base quality when reconstructing fragment bases from BAM reads. Default: 0.",
    )
    parser.add_argument(
        "--max-duplicates",
        type=int,
        default=0,
        help=(
            "Maximum duplicate fragments with the same coordinates to keep per BAM. "
            "Default 0 keeps one fragment per coordinate. Use -1 to disable coordinate deduplication."
        ),
    )
    parser.add_argument(
        "--include-duplicates",
        "--include_duplicate_flag",
        dest="include_duplicates",
        action="store_true",
        help="Include duplicate-marked reads. Default: skip duplicate-marked reads.",
    )
    parser.add_argument(
        "--allow-improper-pairs",
        "--no_require_proper_pairs",
        dest="allow_improper_pairs",
        action="store_true",
        help="Allow pairs not marked proper. Default: require proper pairs.",
    )
    parser.add_argument(
        "--allow-softclipped",
        action="store_true",
        help="Allow soft/hard-clipped or padded reads. Default: skip them.",
    )
    parser.add_argument(
        "--require-complete-fragment",
        action="store_true",
        help="Only use fragments where every base in the profiled interval is A/C/G/T.",
    )
    parser.add_argument(
        "--fraction",
        action="store_true",
        help="Output fractions from 0 to 1 instead of percentages from 0 to 100.",
    )
    parser.add_argument(
        "--max-fragments",
        type=int,
        default=None,
        help="Optional maximum number of matched fragments to reconstruct/fetch in the profile pass.",
    )

    args = parser.parse_args()

    if args.length_a < 2 or args.length_b < 2:
        sys.exit("ERROR: --length_a and --length_b must both be >= 2 for dinucleotide profiles.")
    if args.window_length is not None and args.window_length < 2:
        sys.exit("ERROR: --window-length must be >= 2 for dinucleotide profiles.")
    if args.extend_bp < 0:
        sys.exit("ERROR: --extend-bp must be >= 0.")
    if args.extend_left_bp < 0:
        sys.exit("ERROR: --extend-left-bp/--extend-start-bp must be >= 0.")
    if args.extend_right_bp < 0:
        sys.exit("ERROR: --extend-right-bp/--extend-end-bp must be >= 0.")
    if args.shift_bp < 0:
        sys.exit("ERROR: --shift-bp must be >= 0. Use --shift-direction minus for lower-coordinate shifts.")
    if args.min_baseq < 0:
        sys.exit("ERROR: --min-baseq must be >= 0.")

    extend_left_bp = args.extend_bp + args.extend_left_bp
    extend_right_bp = args.extend_bp + args.extend_right_bp

    if (extend_left_bp > 0 or extend_right_bp > 0 or args.shift_bp > 0) and args.fasta is None:
        sys.exit(
            "ERROR: --fasta is required when using --shift-bp, --extend-bp, "
            "--extend-left-bp/--extend-start-bp, or --extend-right-bp/--extend-end-bp."
        )

    try:
        profile_specs = build_profile_specs_by_signal(
            length_a=args.length_a,
            length_b=args.length_b,
            window_length=args.window_length,
            extend_left_bp=extend_left_bp,
            extend_right_bp=extend_right_bp,
        )
    except ValueError as e:
        sys.exit(f"ERROR: {e}.")

    bam_files_a = expand_bam_inputs(args.bam_a)
    bam_files_b = expand_bam_inputs(args.bam_b)
    requested_chroms = expand_chroms(args.chroms)

    label_a = sanitize_filename(args.label_a)
    label_b = sanitize_filename(args.label_b)

    print("[INFO] signal A BAM files:", file=sys.stderr)
    for b in bam_files_a:
        print(f"  {b}", file=sys.stderr)
    print("[INFO] signal B BAM files:", file=sys.stderr)
    for b in bam_files_b:
        print(f"  {b}", file=sys.stderr)

    print(
        f"[INFO] A={label_a} length={args.length_a}; "
        f"B={label_b} length={args.length_b}",
        file=sys.stderr,
    )
    print(f"[INFO] out prefix: {args.out_prefix}", file=sys.stderr)
    print(f"[INFO] match mode: {args.match_mode}", file=sys.stderr)
    if args.match_mode == "dyad":
        print(f"[INFO] dyad match mode: {args.dyad_match_mode}", file=sys.stderr)
    print("[INFO] low-memory mode: processing one chromosome at a time", file=sys.stderr)
    print(
        f"[INFO] shift: {args.shift_bp} bp direction={args.shift_direction}; "
        f"extend left/start={extend_left_bp} bp; extend right/end={extend_right_bp} bp",
        file=sys.stderr,
    )

    for signal in SIGNALS:
        spec = profile_specs[signal]
        print(
            f"[INFO] signal {signal}: length {spec['fragment_length']}; "
            f"core centred window offsets {spec['core_start_offset']}:{spec['core_end_offset']} "
            f"length={spec['core_profile_length']}; profiled offsets after extension "
            f"{spec['profile_start_offset']}:{spec['profile_end_offset']} "
            f"length={spec['profile_length']}; positions are relative to "
            f"frag_start + floor({spec['fragment_length']} / 2)",
            file=sys.stderr,
        )

    (
        profiles,
        shift_direction_counts,
        subset_names,
        total_pairs_by_signal,
        unique_coords_count_by_signal,
        matched_coords_count,
        stopped_early,
    ) = process_bam_sets_lowmem(
        bam_files_a=bam_files_a,
        length_a=args.length_a,
        label_a=label_a,
        bam_files_b=bam_files_b,
        length_b=args.length_b,
        label_b=label_b,
        requested_chroms=requested_chroms,
        match_mode=args.match_mode,
        dyad_match_mode=args.dyad_match_mode,
        profile_specs=profile_specs,
        mapq=args.mapq,
        min_baseq=args.min_baseq,
        include_duplicates=args.include_duplicates,
        allow_improper_pairs=args.allow_improper_pairs,
        allow_softclipped=args.allow_softclipped,
        max_duplicates=args.max_duplicates,
        require_complete_fragment=args.require_complete_fragment,
        max_fragments=args.max_fragments,
        fasta_path=args.fasta,
        shift_bp=args.shift_bp,
        shift_direction=args.shift_direction,
        random_seed=args.random_seed,
    )

    signal_labels = {"A": label_a, "B": label_b}
    signal_lengths = {"A": args.length_a, "B": args.length_b}

    print("[INFO] matching/profile pass complete", file=sys.stderr)
    for signal in SIGNALS:
        print(
            f"[INFO] signal {signal} ({signal_labels[signal]}): "
            f"{total_pairs_by_signal[signal]:,} fragments seen after filters; "
            f"{unique_coords_count_by_signal[signal]:,} unique coordinates after coordinate deduplication",
            file=sys.stderr,
        )
        parts = [
            f"{subset}-shared={matched_coords_count[signal][subset]:,}"
            for subset in subset_names
        ]
        print(f"[INFO] signal {signal} ({signal_labels[signal]}): " + "; ".join(parts), file=sys.stderr)

    if all(
        matched_coords_count[signal][subset] == 0
        for signal in SIGNALS
        for subset in subset_names
    ):
        sys.exit("ERROR: no matching fragments found between signal A and signal B.")

    output_paths = {}
    for signal in SIGNALS:
        for subset in subset_names:
            out_path = (
                f"{args.out_prefix}_{signal_labels[signal]}_"
                f"len{signal_lengths[signal]}_shared_{subset}.tsv"
            )
            write_profile_tsv(
                profile=profiles[(signal, subset)],
                out_path=out_path,
                fraction=args.fraction,
            )
            output_paths[(signal, subset)] = out_path
            print(f"[DONE] output: {out_path}", file=sys.stderr)

    summary_path = f"{args.out_prefix}_summary.tsv"
    write_summary_tsv(
        out_path=summary_path,
        signal_labels=signal_labels,
        signal_lengths=signal_lengths,
        subset_names=subset_names,
        total_pairs_by_signal=total_pairs_by_signal,
        unique_coords_count_by_signal=unique_coords_count_by_signal,
        matched_coords_count=matched_coords_count,
        profiles=profiles,
        output_paths=output_paths,
        stopped_early=stopped_early,
    )
    print(f"[DONE] summary: {summary_path}", file=sys.stderr)

    print(
        "[INFO] shared profiles are between BAM sets: signal A fragments are matched only "
        "against signal B fragments, not against other fragments within the same signal.",
        file=sys.stderr,
    )
    if args.match_mode == "dyad":
        print(
            "[INFO] shared_dyad means A and B fragments have at least one dyad/centre key in common. "
            "With --dyad-match-mode split, even-length fragments contribute both integer middle positions.",
            file=sys.stderr,
        )
    else:
        print(
            "[INFO] shared_left means A and B fragments share their left/start coordinate; "
            "shared_right means A and B fragments share their right/end coordinate; "
            "shared_any is the union of those two matched subsets.",
            file=sys.stderr,
        )
    print(
        "[INFO] Dinucleotide position is the start of the dinucleotide relative to "
        "frag_start + floor(length / 2), after applying any centred window and extension. "
        "A shift changes which reference bases are fetched, but output position labels "
        "remain relative to the original fragment coordinate system.",
        file=sys.stderr,
    )
    print(
        "[INFO] shift directions used: "
        + ", ".join(f"{k}={v:,}" for k, v in sorted(shift_direction_counts.items())),
        file=sys.stderr,
    )


if __name__ == "__main__":
    main()
