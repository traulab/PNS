#!/usr/bin/env python3

import argparse
import glob
import os
import random
import sys
from collections import defaultdict

import pysam


DINUCS = [
    "AA", "AC", "AG", "AT",
    "CA", "CC", "CG", "CT",
    "GA", "GC", "GG", "GT",
    "TA", "TC", "TG", "TT",
]

WW_DINUCS = {"AA", "AT", "TA", "TT"}  # A/T-only dinucleotides
SS_DINUCS = {"CC", "CG", "GC", "GG"}  # G/C-only dinucleotides
VALID_BASES = {"A", "C", "G", "T"}

CIGAR_SOFT_HARD_OR_PAD = {4, 5, 6}  # S, H, P
BASE_SUBSETS = ("any", "left", "right")
DYAD_SUBSET = "dyad"


def has_softclip_or_hardclip_or_padding(cigartuples):
    if not cigartuples:
        return False

    for op, _length in cigartuples:
        if op in CIGAR_SOFT_HARD_OR_PAD:
            return True

    return False


def expand_bam_inputs(bam_inputs):
    bam_files = []

    for item in bam_inputs:
        matches = sorted(glob.glob(item))

        if matches:
            bam_files.extend(matches)
        else:
            bam_files.append(item)

    bam_files = list(dict.fromkeys(bam_files))

    missing = [b for b in bam_files if not os.path.exists(b)]

    if missing:
        raise FileNotFoundError(
            "These BAM files/patterns were not found:\n  " + "\n  ".join(missing)
        )

    return bam_files


def expand_chroms(chrom_spec):
    chrom_spec = chrom_spec.replace(" ", "")

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
    else:
        with_chr = "chr" + chrom
        if with_chr in references:
            return with_chr

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
    """
    Return the template start and length for a paired-end fragment.

    For standard paired-end BAMs, one mate has positive TLEN and the other has
    negative TLEN. The positive-TLEN read should be the left-most read, and its
    reference_start plus TLEN defines the full template span.
    """
    if aln1.template_length > 0:
        frag_start = aln1.reference_start
        frag_len = aln1.template_length
    elif aln2.template_length > 0:
        frag_start = aln2.reference_start
        frag_len = aln2.template_length
    else:
        # Fallback for unusual BAMs. This should rarely be needed.
        frag_start = min(aln1.reference_start, aln2.reference_start)
        frag_end = max(aln1.reference_end, aln2.reference_end)
        frag_len = frag_end - frag_start

    return frag_start, abs(frag_len)


def dyad_key_2x(bam_index, chrom, frag_start, frag_len):
    """
    Return a half-base-safe dyad key.

    The geometric fragment centre is:
        frag_start + (frag_len - 1) / 2

    To avoid floating point values for even-length fragments, store 2 * centre:
        2 * frag_start + frag_len - 1

    This makes odd and even fragment dyads comparable without rounding.
    """
    return (bam_index, chrom, (2 * frag_start) + frag_len - 1)


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
    """
    Add dinucleotide counts for one profiled interval.

    Dinucleotide position i means the dinucleotide starting at profiled-interval
    base i. Relative position is added later as:

        rel_position = profile_start_offset + i - (fragment_length // 2)

    Returns:
      - None if the interval is skipped because --require-complete-fragment was used
      - otherwise, the number of valid dinucleotide positions added
    """
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
    """Return the genomic-coordinate shift to apply to this fragment/window."""
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
    """
    Fetch reference bases for [start, end).

    Returns None when the requested interval falls outside the reference contig.
    """
    if start_0based < 0:
        return None

    chrom_len = fasta.get_reference_length(fasta_chrom)

    if end_0based > chrom_len:
        return None

    seq = fasta.fetch(fasta_chrom, start_0based, end_0based).upper()

    if len(seq) != end_0based - start_0based:
        return None

    return list(seq)


def iter_fragment_pairs_from_bam_chrom_lengths(
    bam,
    bam_chrom,
    frag_lengths_interest,
    mapq,
    include_duplicates,
    allow_improper_pairs,
    allow_softclipped,
    max_duplicates,
    coord_counts,
):
    """
    Iterate exact-length paired-end fragments from one BAM contig for either of
    the requested fragment lengths.

    This scans the BAM once and pairs reads by query name, avoiding a FASTA file
    and avoiding pysam.mate() calls for every read.
    """
    wanted_lengths = set(frag_lengths_interest)
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

        if abs(aln.template_length) not in wanted_lengths:
            continue

        qname = aln.query_name

        if qname not in pending:
            pending[qname] = aln
            continue

        mate = pending.pop(qname)

        if aln.reference_id != mate.reference_id:
            continue

        frag_start, frag_len = inferred_fragment_start_and_length(aln, mate)

        if frag_len not in wanted_lengths:
            continue

        frag_end = frag_start + frag_len

        if max_duplicates >= 0:
            key = (bam_chrom, frag_start, frag_end)

            # max_duplicates=0 keeps one fragment per coordinate.
            if coord_counts[key] > max_duplicates:
                continue

            coord_counts[key] += 1

        yield aln, mate, frag_start, frag_end, frag_len


def find_matching_fragment_coordinates(
    bam_files,
    requested_chroms,
    lengths,
    include_shared_dyads,
    mapq,
    include_duplicates,
    allow_improper_pairs,
    allow_softclipped,
    max_duplicates,
):
    """
    First pass: collect coordinates for both fragment lengths, then find the
    coordinate subsets where a fragment of one length shares its left end,
    right end, and optionally dyad with a fragment of the other length.

    BAM index is included in the matching key, so multiple BAM inputs are
    aggregated in the final profiles but fragments are only matched within the
    same BAM file, not artificially across different samples/files.
    """
    coords_by_len = {length: set() for length in lengths}
    left_end_keys_by_len = {length: set() for length in lengths}
    right_end_keys_by_len = {length: set() for length in lengths}
    dyad_keys_by_len = {length: set() for length in lengths}
    total_pairs_by_len = {length: 0 for length in lengths}

    for bam_index, bam_path in enumerate(bam_files):
        print(f"[PASS1] Opening BAM: {bam_path}", file=sys.stderr)

        with pysam.AlignmentFile(bam_path, "rb") as bam:
            bam_refs = set(bam.references)

            bam_chroms = []
            missing_chroms = []

            for chrom in requested_chroms:
                bam_chrom = matching_contig(chrom, bam_refs)

                if bam_chrom is None:
                    missing_chroms.append(chrom)
                else:
                    bam_chroms.append(bam_chrom)

            if missing_chroms:
                print(
                    f"[WARNING] In {os.path.basename(bam_path)}, chromosomes not found and skipped: "
                    f"{','.join(missing_chroms)}",
                    file=sys.stderr,
                )

            coord_counts = defaultdict(int)

            for bam_chrom in bam_chroms:
                print(f"[PASS1] Scanning {os.path.basename(bam_path)} {bam_chrom}", file=sys.stderr)

                for _aln1, _aln2, frag_start, frag_end, frag_len in iter_fragment_pairs_from_bam_chrom_lengths(
                    bam=bam,
                    bam_chrom=bam_chrom,
                    frag_lengths_interest=lengths,
                    mapq=mapq,
                    include_duplicates=include_duplicates,
                    allow_improper_pairs=allow_improper_pairs,
                    allow_softclipped=allow_softclipped,
                    max_duplicates=max_duplicates,
                    coord_counts=coord_counts,
                ):
                    coord_key = (bam_index, bam_chrom, frag_start, frag_end)
                    left_key = (bam_index, bam_chrom, frag_start)
                    right_key = (bam_index, bam_chrom, frag_end)
                    dyad_key = dyad_key_2x(bam_index, bam_chrom, frag_start, frag_len)

                    coords_by_len[frag_len].add(coord_key)
                    left_end_keys_by_len[frag_len].add(left_key)
                    right_end_keys_by_len[frag_len].add(right_key)
                    dyad_keys_by_len[frag_len].add(dyad_key)
                    total_pairs_by_len[frag_len] += 1

                    if total_pairs_by_len[frag_len] % 1000000 == 0:
                        print(
                            f"[PASS1] length {frag_len} fragments seen: {total_pairs_by_len[frag_len]:,}",
                            file=sys.stderr,
                        )

    length_a, length_b = lengths
    subset_names = list(BASE_SUBSETS)

    if include_shared_dyads:
        subset_names.append(DYAD_SUBSET)

    match_coords = {
        length: {subset: set() for subset in subset_names}
        for length in lengths
    }

    for length, other_length in ((length_a, length_b), (length_b, length_a)):
        for coord_key in coords_by_len[length]:
            bam_index, chrom, frag_start, frag_end = coord_key
            frag_len = frag_end - frag_start

            if (bam_index, chrom, frag_start) in left_end_keys_by_len[other_length]:
                match_coords[length]["left"].add(coord_key)

            if (bam_index, chrom, frag_end) in right_end_keys_by_len[other_length]:
                match_coords[length]["right"].add(coord_key)

            if include_shared_dyads:
                if dyad_key_2x(bam_index, chrom, frag_start, frag_len) in dyad_keys_by_len[other_length]:
                    match_coords[length]["dyad"].add(coord_key)

        # When dyads are not requested, this exactly matches the original script:
        # shared_any = shared_left union shared_right.
        any_coords = match_coords[length]["left"] | match_coords[length]["right"]

        if include_shared_dyads:
            any_coords = any_coords | match_coords[length]["dyad"]

        match_coords[length]["any"] = any_coords

    return coords_by_len, match_coords, total_pairs_by_len, tuple(subset_names)


def build_profile_specs(lengths, window_length, extend_left_bp, extend_right_bp):
    specs = {}

    for fragment_length in lengths:
        core_profile_length = fragment_length if window_length is None else window_length

        try:
            core_start_offset, core_end_offset = get_centered_window_offsets(
                fragment_length=fragment_length,
                window_length=core_profile_length,
            )
        except ValueError as e:
            raise ValueError(f"for length {fragment_length}, --window-length invalid: {e}")

        profile_start_offset = core_start_offset - extend_left_bp
        profile_end_offset = core_end_offset + extend_right_bp
        profile_length = profile_end_offset - profile_start_offset

        if profile_length < 2:
            raise ValueError(
                f"for length {fragment_length}, profiled interval is < 2 bp after extension"
            )

        specs[fragment_length] = {
            "fragment_length": fragment_length,
            "core_profile_length": core_profile_length,
            "core_start_offset": core_start_offset,
            "core_end_offset": core_end_offset,
            "profile_start_offset": profile_start_offset,
            "profile_end_offset": profile_end_offset,
            "profile_length": profile_length,
        }

    return specs


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


def profile_matched_fragments(
    bam_files,
    requested_chroms,
    lengths,
    subset_names,
    match_coords,
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
    profiles = {}

    for length in lengths:
        spec = profile_specs[length]

        for subset in subset_names:
            profiles[(length, subset)] = make_profile(
                fragment_length=length,
                core_start_offset=spec["core_start_offset"],
                core_end_offset=spec["core_end_offset"],
                profile_start_offset=spec["profile_start_offset"],
                profile_end_offset=spec["profile_end_offset"],
            )

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
    reconstructed_or_fetched_matching_fragments = 0

    for bam_index, bam_path in enumerate(bam_files):
        print(f"[PASS2] Opening BAM: {bam_path}", file=sys.stderr)

        with pysam.AlignmentFile(bam_path, "rb") as bam:
            bam_refs = set(bam.references)

            bam_chroms = []
            missing_chroms = []

            for chrom in requested_chroms:
                bam_chrom = matching_contig(chrom, bam_refs)

                if bam_chrom is None:
                    missing_chroms.append(chrom)
                else:
                    bam_chroms.append(bam_chrom)

            if missing_chroms:
                print(
                    f"[WARNING] In {os.path.basename(bam_path)}, chromosomes not found and skipped: "
                    f"{','.join(missing_chroms)}",
                    file=sys.stderr,
                )

            coord_counts = defaultdict(int)

            for bam_chrom in bam_chroms:
                fasta_chrom = None

                if use_reference_sequence:
                    fasta_chrom = matching_contig(bam_chrom, fasta_refs)

                    if fasta_chrom is None:
                        print(
                            f"[WARNING] {bam_chrom} is present in BAM but not found in FASTA; skipped",
                            file=sys.stderr,
                        )
                        continue

                print(f"[PASS2] Profiling {os.path.basename(bam_path)} {bam_chrom}", file=sys.stderr)

                for aln1, aln2, frag_start, frag_end, frag_len in iter_fragment_pairs_from_bam_chrom_lengths(
                    bam=bam,
                    bam_chrom=bam_chrom,
                    frag_lengths_interest=lengths,
                    mapq=mapq,
                    include_duplicates=include_duplicates,
                    allow_improper_pairs=allow_improper_pairs,
                    allow_softclipped=allow_softclipped,
                    max_duplicates=max_duplicates,
                    coord_counts=coord_counts,
                ):
                    coord_key = (bam_index, bam_chrom, frag_start, frag_end)

                    matched_subsets = [
                        subset for subset in subset_names
                        if coord_key in match_coords[frag_len][subset]
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
                        spec=profile_specs[frag_len],
                        min_baseq=min_baseq,
                        use_reference_sequence=use_reference_sequence,
                        fasta=fasta,
                        fasta_chrom=fasta_chrom,
                        shift_delta=shift_delta,
                    )

                    reconstructed_or_fetched_matching_fragments += 1

                    for subset in matched_subsets:
                        profile = profiles[(frag_len, subset)]
                        profile["candidate_fragments_seen"] += 1

                        if profile_bases is None:
                            profile["fragments_skipped_reference_bounds"] += 1
                            continue

                        add_fragment_dinucs_to_profile(
                            bases=profile_bases,
                            profile=profile,
                            require_complete_fragment=require_complete_fragment,
                        )

                    if reconstructed_or_fetched_matching_fragments % 100000 == 0:
                        print(
                            f"[PASS2] reconstructed/fetched matching fragments: {reconstructed_or_fetched_matching_fragments:,}",
                            file=sys.stderr,
                        )

                    if max_fragments is not None and reconstructed_or_fetched_matching_fragments >= max_fragments:
                        if fasta is not None:
                            fasta.close()
                        return profiles, shift_direction_counts

    if fasta is not None:
        fasta.close()

    return profiles, shift_direction_counts


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

                ww_value = (ww_count / n_valid) * multiplier
                ss_value = (ss_count / n_valid) * multiplier

                row.append(f"{ww_value:.8g}")
                row.append(f"{ss_value:.8g}")

            out.write("\t".join(row) + "\n")


def write_summary_tsv(out_path, lengths, subset_names, coords_by_len, match_coords, profiles, output_paths):
    header = [
        "length",
        "subset",
        "unique_fragment_coordinates_in_length",
        "unique_matched_coordinates",
        "candidate_fragments_seen_in_pass2",
        "fragments_used",
        "fragments_skipped_incomplete",
        "fragments_skipped_reference_bounds",
        "fragments_with_no_valid_dinucs",
        "core_start_offset_0based",
        "core_end_offset_0based_exclusive",
        "profile_start_offset_0based",
        "profile_end_offset_0based_exclusive",
        "profile_length",
        "output_tsv",
    ]

    with open(out_path, "w") as out:
        out.write("\t".join(header) + "\n")

        for length in lengths:
            for subset in subset_names:
                profile = profiles[(length, subset)]
                row = [
                    str(length),
                    subset,
                    str(len(coords_by_len[length])),
                    str(len(match_coords[length][subset])),
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
                    output_paths[(length, subset)],
                ]
                out.write("\t".join(row) + "\n")


def main():
    parser = argparse.ArgumentParser(
        description=(
            "Generate observed dinucleotide percentage profiles from two fragment lengths, "
            "using only fragments whose left end, right end, and optionally dyad is shared "
            "with a fragment of the other length. By default, bases are reconstructed from "
            "the BAM read pairs. Optionally, a reference FASTA can be used to profile a "
            "centred window, extend that window, and/or shift it in genomic coordinates "
            "before dinucleotide counting."
        )
    )

    parser.add_argument(
        "--bam",
        nargs="+",
        required=True,
        help="Input BAM file(s). Can use one BAM, multiple BAMs, or wildcard *.bam.",
    )

    parser.add_argument(
        "--lengths",
        nargs=2,
        type=int,
        required=True,
        metavar=("LEN_A", "LEN_B"),
        help="The two exact fragment lengths to compare, e.g. --lengths 147 167.",
    )

    parser.add_argument(
        "--include-shared-dyads",
        action="store_true",
        help=(
            "Also find fragments whose dyad/centre is shared with a fragment of the other length. "
            "When enabled, an extra _shared_dyad.tsv is written for each length, and "
            "_shared_any.tsv is the union of shared_left, shared_right, and shared_dyad. "
            "When disabled, behaviour matches the original script: shared_any is left union right."
        ),
    )

    parser.add_argument(
        "--window-length",
        type=int,
        default=None,
        help=(
            "Optional centred window length to profile within each fragment before extension/shift. "
            "For example, with --lengths 147 167 --window-length 147, 147 bp fragments are profiled "
            "as-is, and 167 bp fragments are profiled over their centred 147 bp window. "
            "Default: profile the full fragment length for each length."
        ),
    )

    parser.add_argument(
        "--extend-bp",
        type=int,
        default=0,
        help=(
            "Extend the selected full-fragment or centred-window interval by this many bp "
            "on BOTH sides before profiling. This is added to any side-specific extension. "
            "Requires --fasta when > 0. Default: 0."
        ),
    )

    parser.add_argument(
        "--extend-left-bp",
        "--extend-start-bp",
        dest="extend_left_bp",
        type=int,
        default=0,
        help=(
            "Extend only the genomic-left/start side by this many bp before profiling. "
            "This is added on top of --extend-bp. Requires --fasta when > 0. Default: 0."
        ),
    )

    parser.add_argument(
        "--extend-right-bp",
        "--extend-end-bp",
        dest="extend_right_bp",
        type=int,
        default=0,
        help=(
            "Extend only the genomic-right/end side by this many bp before profiling. "
            "This is added on top of --extend-bp. Requires --fasta when > 0. Default: 0."
        ),
    )

    parser.add_argument(
        "--shift-bp",
        type=int,
        default=0,
        help=(
            "Shift the selected/profiled interval by this many bp before fetching bases from "
            "the reference genome. Requires --fasta when > 0. Default: 0."
        ),
    )

    parser.add_argument(
        "--shift-direction",
        choices=["plus", "minus", "random"],
        default="plus",
        help=(
            "Direction for --shift-bp. plus means higher genomic coordinates, minus means "
            "lower genomic coordinates, and random chooses plus/minus independently for "
            "each matched fragment. Default: plus."
        ),
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
            "Reference genome FASTA. Required for --shift-bp > 0 or any extension option. "
            "When supplied, the profiled interval is fetched from the reference genome "
            "rather than reconstructed from BAM read bases. The FASTA must be indexed."
        ),
    )

    parser.add_argument(
        "--chroms",
        default="all",
        help="Chromosomes to analyse, e.g. all, 1-22,X,Y, 1,2,3, or X. Default: all.",
    )

    parser.add_argument(
        "--out-prefix",
        required=True,
        help=(
            "Output prefix. The script writes <prefix>_len<L>_shared_any.tsv, "
            "<prefix>_len<L>_shared_left.tsv, <prefix>_len<L>_shared_right.tsv, "
            "optionally <prefix>_len<L>_shared_dyad.tsv, and <prefix>_summary.tsv."
        ),
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
        help="Minimum base quality to use when reconstructing fragment bases. Default: 0.",
    )

    parser.add_argument(
        "--max-duplicates",
        type=int,
        default=0,
        help=(
            "Maximum duplicate fragments with same coordinates to keep. "
            "Default 0 keeps one fragment per coordinate. Use -1 to disable coordinate deduplication."
        ),
    )

    parser.add_argument(
        "--include-duplicates",
        action="store_true",
        help="Include duplicate-marked reads. Default: skip duplicate-marked reads.",
    )

    parser.add_argument(
        "--allow-improper-pairs",
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
        help=(
            "Only use fragments where every base in the profiled interval is A/C/G/T. "
            "Default: use each dinucleotide position whenever both bases are valid."
        ),
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
        help=(
            "Optional maximum number of matched fragments to reconstruct/fetch in pass 2. "
            "Useful for testing only; pass 1 still scans all requested chromosomes so matches are found correctly."
        ),
    )

    args = parser.parse_args()

    lengths = list(args.lengths)

    if len(set(lengths)) != 2:
        sys.exit("ERROR: --lengths must contain two different fragment lengths.")

    for length in lengths:
        if length < 2:
            sys.exit("ERROR: all fragment lengths must be >= 2 for dinucleotide profiles.")

    if args.window_length is not None and args.window_length < 2:
        sys.exit("ERROR: --window-length must be >= 2 for dinucleotide profiles.")

    if args.extend_bp < 0:
        sys.exit("ERROR: --extend-bp must be >= 0.")

    if args.extend_left_bp < 0:
        sys.exit("ERROR: --extend-left-bp/--extend-start-bp must be >= 0.")

    if args.extend_right_bp < 0:
        sys.exit("ERROR: --extend-right-bp/--extend-end-bp must be >= 0.")

    extend_left_bp = args.extend_bp + args.extend_left_bp
    extend_right_bp = args.extend_bp + args.extend_right_bp

    if args.shift_bp < 0:
        sys.exit("ERROR: --shift-bp must be >= 0. Use --shift-direction minus for left/lower-coordinate shifts.")

    if args.min_baseq < 0:
        sys.exit("ERROR: --min-baseq must be >= 0.")

    if (extend_left_bp > 0 or extend_right_bp > 0 or args.shift_bp > 0) and args.fasta is None:
        sys.exit(
            "ERROR: --fasta is required when using --shift-bp, --extend-bp, "
            "--extend-left-bp/--extend-start-bp, or --extend-right-bp/--extend-end-bp."
        )

    try:
        profile_specs = build_profile_specs(
            lengths=lengths,
            window_length=args.window_length,
            extend_left_bp=extend_left_bp,
            extend_right_bp=extend_right_bp,
        )
    except ValueError as e:
        sys.exit(f"ERROR: {e}.")

    bam_files = expand_bam_inputs(args.bam)

    print("[INFO] BAM files:", file=sys.stderr)
    for b in bam_files:
        print(f"  {b}", file=sys.stderr)

    print(f"[INFO] fragment lengths: {lengths[0]}, {lengths[1]}", file=sys.stderr)
    print(f"[INFO] out prefix: {args.out_prefix}", file=sys.stderr)
    print(f"[INFO] include shared dyads: {args.include_shared_dyads}", file=sys.stderr)
    print(
        f"[INFO] shift: {args.shift_bp} bp direction={args.shift_direction}; "
        f"extend left/start={extend_left_bp} bp; extend right/end={extend_right_bp} bp",
        file=sys.stderr,
    )

    for length in lengths:
        spec = profile_specs[length]
        print(
            f"[INFO] length {length}: core centred window offsets {spec['core_start_offset']}:{spec['core_end_offset']} "
            f"length={spec['core_profile_length']}; profiled offsets after extension "
            f"{spec['profile_start_offset']}:{spec['profile_end_offset']} length={spec['profile_length']}; "
            f"positions are relative to frag_start + floor({length} / 2)",
            file=sys.stderr,
        )

    requested_chroms = expand_chroms(args.chroms)

    coords_by_len, match_coords, total_pairs_by_len, subset_names = find_matching_fragment_coordinates(
        bam_files=bam_files,
        requested_chroms=requested_chroms,
        lengths=lengths,
        include_shared_dyads=args.include_shared_dyads,
        mapq=args.mapq,
        include_duplicates=args.include_duplicates,
        allow_improper_pairs=args.allow_improper_pairs,
        allow_softclipped=args.allow_softclipped,
        max_duplicates=args.max_duplicates,
    )

    print("[INFO] pass 1 complete", file=sys.stderr)
    for length in lengths:
        print(
            f"[INFO] length {length}: {total_pairs_by_len[length]:,} fragments seen; "
            f"{len(coords_by_len[length]):,} unique coordinates after coordinate deduplication",
            file=sys.stderr,
        )

        parts = []
        for subset in subset_names:
            parts.append(f"{subset}-shared={len(match_coords[length][subset]):,}")
        print(f"[INFO] length {length}: " + "; ".join(parts), file=sys.stderr)

    if all(len(match_coords[length]["any"]) == 0 for length in lengths):
        sys.exit("ERROR: no matching fragments found between the two lengths.")

    profiles, shift_direction_counts = profile_matched_fragments(
        bam_files=bam_files,
        requested_chroms=requested_chroms,
        lengths=lengths,
        subset_names=subset_names,
        match_coords=match_coords,
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

    output_paths = {}

    for length in lengths:
        for subset in subset_names:
            out_path = f"{args.out_prefix}_len{length}_shared_{subset}.tsv"
            write_profile_tsv(
                profile=profiles[(length, subset)],
                out_path=out_path,
                fraction=args.fraction,
            )
            output_paths[(length, subset)] = out_path
            print(f"[DONE] output: {out_path}", file=sys.stderr)

    summary_path = f"{args.out_prefix}_summary.tsv"
    write_summary_tsv(
        out_path=summary_path,
        lengths=lengths,
        subset_names=subset_names,
        coords_by_len=coords_by_len,
        match_coords=match_coords,
        profiles=profiles,
        output_paths=output_paths,
    )

    print(f"[DONE] summary: {summary_path}", file=sys.stderr)
    print(
        "[INFO] shared_left means fragments whose left/start coordinate is also used by a fragment "
        "of the other length.",
        file=sys.stderr,
    )
    print(
        "[INFO] shared_right means fragments whose right/end coordinate is also used by a fragment "
        "of the other length.",
        file=sys.stderr,
    )

    if args.include_shared_dyads:
        print(
            "[INFO] shared_dyad means fragments whose dyad/centre coordinate is also used by a fragment "
            "of the other length. Dyads are matched with 2*centre = 2*start + length - 1, "
            "so even-length half-base dyads are handled without rounding.",
            file=sys.stderr,
        )
        print(
            "[INFO] shared_any is the union of shared_left, shared_right, and shared_dyad.",
            file=sys.stderr,
        )
    else:
        print(
            "[INFO] shared_any is the union of shared_left and shared_right. Use --include-shared-dyads "
            "to also output shared_dyad and include dyads in shared_any.",
            file=sys.stderr,
        )

    print(
        "[INFO] Dinucleotide position is the start of the dinucleotide relative to "
        "frag_start + floor(length / 2), after applying any centred window and extension. "
        "A shift changes which reference bases are fetched, but the output position labels "
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
