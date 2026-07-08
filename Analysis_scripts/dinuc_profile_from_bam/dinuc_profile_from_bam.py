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


def sequence_is_acgt(seq):
    seq = seq.upper()
    return all(base in VALID_BASES for base in seq)


def sanitize_filename_token(token):
    """Make a string safe and compact for use inside an output filename."""
    token = str(token)
    safe = []

    for char in token:
        if char.isalnum() or char in ("-", "_", "."):
            safe.append(char)
        elif char in ("/", "\\", ":", ",", " "):
            safe.append("_")
        else:
            safe.append("_")

    cleaned = "".join(safe).strip("._-")

    while "__" in cleaned:
        cleaned = cleaned.replace("__", "_")

    return cleaned or "NA"


def bam_name_token(bam_files):
    """Return a compact name derived from one or more BAM filenames."""
    names = []

    for bam_path in bam_files:
        name = os.path.basename(bam_path)
        for suffix in (".bam", ".cram", ".sam"):
            if name.lower().endswith(suffix):
                name = name[: -len(suffix)]
                break
        names.append(sanitize_filename_token(name))

    if len(names) == 1:
        return names[0]

    if len(names) <= 3:
        return "plus".join(names)

    return f"{names[0]}_plus{len(names) - 1}BAMs"


def chroms_token(chrom_spec):
    """Return a filename token for the chromosome argument without expanding huge lists."""
    spec = str(chrom_spec).replace(" ", "")

    if spec.lower() == "all":
        return "chrAll"

    if spec.startswith("chr"):
        return sanitize_filename_token(spec)

    return "chr" + sanitize_filename_token(spec)


def auto_output_path(args, bam_files, profile_len, profile_basis):
    """Build an output filename from the BAM basename and key analysis parameters."""
    tokens = [
        bam_name_token(bam_files),
        f"len{args.length}",
        chroms_token(args.chroms),
    ]

    if profile_basis == "flank":
        tokens.append(f"flank{args.flank}")
    else:
        tokens.append("fragmentCentered")

    if args.extend_bp:
        tokens.append(f"extBoth{args.extend_bp}")
    if args.extend_left_bp:
        tokens.append(f"extL{args.extend_left_bp}")
    if args.extend_right_bp:
        tokens.append(f"extR{args.extend_right_bp}")
    if args.extend_random_side_bp:
        tokens.append(f"extRandSide{args.extend_random_side_bp}")

    if args.shift_bp:
        direction_token = {
            "plus": "plus",
            "minus": "minus",
            "random": "rand",
        }[args.shift_direction]
        tokens.append(f"shift{direction_token}{args.shift_bp}")

    tokens.append(f"profile{profile_len}bp")

    if args.allow_n_in_profile:
        tokens.append("allowN")

    tokens.append("frac" if args.fraction else "pct")

    if args.seed is not None and (args.extend_random_side_bp or args.shift_direction == "random"):
        tokens.append(f"seed{args.seed}")

    if args.max_fragments is not None:
        tokens.append(f"max{args.max_fragments}")

    filename = "_".join(tokens) + "_dinuc_freq.tsv"

    if args.out_dir:
        return os.path.join(args.out_dir, filename)

    return filename


def choose_random_side_extension(base_left_ext, base_right_ext, random_side_bp):
    """
    Return left/right extension for one fragment.

    Fixed extension is always applied. If random_side_bp > 0, that many bases
    are added to either the genomic-left/start side or genomic-right/end side.
    The final profiled interval is then centred on this new interval, not on the
    original fragment.
    """
    left_ext = base_left_ext
    right_ext = base_right_ext

    if random_side_bp > 0:
        if random.choice(("left", "right")) == "left":
            left_ext += random_side_bp
        else:
            right_ext += random_side_bp

    return left_ext, right_ext


def choose_shift(shift_bp, shift_direction):
    if shift_bp == 0:
        return 0

    if shift_direction == "plus":
        return shift_bp

    if shift_direction == "minus":
        return -shift_bp

    if shift_direction == "random":
        return random.choice((-shift_bp, shift_bp))

    raise ValueError(f"Unexpected shift direction: {shift_direction}")


def interval_for_fragment_profile(frag_start, frag_end, left_ext, right_ext, shift):
    """
    Profile the fragment itself, plus any extension.

    Example:
      length 145 + random one-sided extension 17 -> final interval length 162.
      Each interval is centred on its own final 162-bp span before dinucleotide
      positions are assigned.
    """
    profile_start = frag_start - left_ext + shift
    profile_end = frag_end + right_ext + shift
    return profile_start, profile_end


def interval_for_flank_profile(frag_start, frag_len, flank, left_ext, right_ext, shift):
    """
    Optional backwards-compatible flank mode.

    This is not the default. It profiles a centre-based window around the
    original fragment centre, then extends and/or shifts that window.
    """
    centre = frag_start + (frag_len // 2)
    profile_start = centre - flank - left_ext + shift

    # +2 rather than +1 because a flank profile of dinucleotide starts
    # -flank..+flank needs one extra base to complete the last dinucleotide.
    profile_end = centre + flank + 2 + right_ext + shift
    return profile_start, profile_end


def get_profile_seq(fasta, fa_chrom, profile_start, profile_end, chrom_len):
    if profile_start < 0 or profile_end > chrom_len:
        return None

    seq = fasta.fetch(fa_chrom, profile_start, profile_end).upper()
    expected_len = profile_end - profile_start

    if len(seq) != expected_len:
        return None

    return seq


def profile_position_range(profile_len):
    """
    Return dinucleotide-start positions for a centred sequence of profile_len bases.

    There are profile_len - 1 dinucleotide start positions.
    Position is relative to:
        profile_start + floor(profile_len / 2)

    Example:
      162-bp profile window -> 161 dinucleotide positions: -81..+79.
    """
    if profile_len < 2:
        raise ValueError("profile_len must be >= 2 for dinucleotide profiling")

    centre_offset = profile_len // 2
    return [i - centre_offset for i in range(profile_len - 1)]


def add_dinuc_profile_from_seq(
    seq,
    dinuc_counts_by_position,
    valid_dinuc_opportunities,
    require_acgt_profile=True,
):
    """
    Add all 16 dinucleotide counts from one profile sequence.

    The sequence is already centred as one complete interval. The output index is
    simply the dinucleotide start index within that interval.
    """
    if require_acgt_profile and not sequence_is_acgt(seq):
        return None

    added_positions = 0

    for i in range(len(seq) - 1):
        dinuc = seq[i : i + 2]

        if len(dinuc) != 2:
            continue

        if dinuc[0] not in VALID_BASES or dinuc[1] not in VALID_BASES:
            continue

        dinuc_counts_by_position[i][dinuc] += 1
        valid_dinuc_opportunities[i] += 1
        added_positions += 1

    return added_positions


def iter_fragments_from_bam_chrom(
    bam,
    bam_chrom,
    frag_len_interest,
    mapq,
    include_duplicates,
    allow_improper_pairs,
    allow_softclipped,
    max_duplicates,
    coord_counts,
):
    """
    Iterate exact-length paired-end fragments.

    Uses positive template_length so each paired fragment is counted once.
    """
    try:
        iterator = bam.fetch(bam_chrom)
    except ValueError:
        return

    for aln in iterator:
        if aln.is_unmapped or aln.mate_is_unmapped:
            continue

        if aln.is_secondary or aln.is_supplementary:
            continue

        if aln.is_qcfail:
            continue

        if not include_duplicates and aln.is_duplicate:
            continue

        if not allow_improper_pairs and not aln.is_proper_pair:
            continue

        if not allow_softclipped and has_softclip_or_hardclip_or_padding(aln.cigartuples):
            continue

        if aln.mapping_quality < mapq:
            continue

        if aln.reference_id != aln.next_reference_id:
            continue

        # Count each paired fragment once only.
        if aln.template_length <= 0:
            continue

        frag_len = abs(aln.template_length)

        if frag_len != frag_len_interest:
            continue

        frag_start = aln.reference_start
        frag_end = frag_start + frag_len

        if max_duplicates >= 0:
            key = (bam_chrom, frag_start, frag_end)

            # max_duplicates=0 keeps one fragment per coordinate.
            if coord_counts[key] > max_duplicates:
                continue

            coord_counts[key] += 1

        yield frag_start, frag_end


def main():
    parser = argparse.ArgumentParser(
        description=(
            "Generate pure observed dinucleotide frequency profiles from exact-length BAM fragments. "
            "By default, the profiled interval is the fragment itself, plus any requested extension, "
            "and each final interval is centred on itself. No expected profile, no shuffling, and no "
            "log2 observed/expected calculation are used."
        )
    )

    parser.add_argument(
        "--bam",
        nargs="+",
        required=True,
        help="Input BAM file(s). Can use one BAM, multiple BAMs, or wildcard *.bam.",
    )

    parser.add_argument(
        "--fasta",
        required=True,
        help="Reference FASTA. Must be indexed with samtools faidx.",
    )

    parser.add_argument(
        "--length",
        type=int,
        required=True,
        help="Exact fragment length to analyse, e.g. 145.",
    )

    parser.add_argument(
        "--flank",
        type=int,
        default=None,
        help=(
            "Optional old-style flank mode. If omitted, profiles the fragment itself. "
            "If supplied, profiles dinucleotide-start positions -flank..+flank around "
            "the original fragment centre before extension."
        ),
    )

    parser.add_argument(
        "--extend-bp",
        type=int,
        default=0,
        help="Extend both sides of the profile interval by this many bp. Default: 0.",
    )

    parser.add_argument(
        "--extend-left-bp",
        "--extend-start-bp",
        type=int,
        default=0,
        help="Extend only the genomic-left/start side by this many bp. Default: 0.",
    )

    parser.add_argument(
        "--extend-right-bp",
        "--extend-end-bp",
        type=int,
        default=0,
        help="Extend only the genomic-right/end side by this many bp. Default: 0.",
    )

    parser.add_argument(
        "--extend-random-side-bp",
        "--extend-random-end-bp",
        "--extend-random-ends-bp",
        type=int,
        default=0,
        help=(
            "For each fragment, randomly extend either the genomic-left/start side or "
            "the genomic-right/end side by this many bp. Default: 0. The final interval "
            "is centred on itself, so length 145 with this set to 17 is treated as a "
            "162-bp profiled interval."
        ),
    )

    parser.add_argument(
        "--shift-bp",
        type=int,
        default=0,
        help=(
            "Shift the final profile interval by this many bp after extension. Default: 0. "
            "Position labels remain centred on the final shifted interval."
        ),
    )

    parser.add_argument(
        "--shift-direction",
        choices=("plus", "minus", "random"),
        default="plus",
        help="Direction for --shift-bp. plus=higher coordinates, minus=lower coordinates, random=per-fragment random. Default: plus.",
    )

    parser.add_argument(
        "--chroms",
        default="all",
        help="Chromosomes to analyse, e.g. all, 1-22,X,Y, 1,2,3, or X. Default: all.",
    )

    parser.add_argument(
        "--out",
        default=None,
        help=(
            "Output TSV file with pure observed dinucleotide frequencies. "
            "If omitted, a filename is generated automatically from the BAM basename "
            "and key parameters such as length, chromosomes, extension, shift, seed, "
            "and percent/fraction mode."
        ),
    )

    parser.add_argument(
        "--out-dir",
        default=".",
        help=(
            "Directory for automatically generated output names when --out is omitted. "
            "Default: current directory."
        ),
    )

    parser.add_argument(
        "--mapq",
        type=int,
        default=30,
        help="Minimum MAPQ. Default: 30.",
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
        "--allow-n-in-profile",
        action="store_true",
        help=(
            "Allow Ns in profile windows. Default: require full A/C/G/T profile windows. "
            "When enabled, positions containing Ns are skipped position-by-position."
        ),
    )

    parser.add_argument(
        "--fraction",
        action="store_true",
        help="Output fractions from 0 to 1 instead of percentages from 0 to 100.",
    )

    parser.add_argument(
        "--seed",
        type=int,
        default=None,
        help="Random seed for reproducible random-side extension and random shifting.",
    )

    parser.add_argument(
        "--max-fragments",
        type=int,
        default=None,
        help="Optional maximum number of observed fragments to process.",
    )

    # Accepted only so older log2OE command lines fail less abruptly. They are ignored.
    parser.add_argument("--shuffles-per-fragment", type=int, default=None, help=argparse.SUPPRESS)
    parser.add_argument("--pseudocount-per-dinuc", type=float, default=None, help=argparse.SUPPRESS)

    args = parser.parse_args()

    if args.length < 2:
        sys.exit("ERROR: --length must be >= 2.")

    if args.flank is not None and args.flank < 0:
        sys.exit("ERROR: --flank must be >= 0.")

    for name, value in (
        ("--extend-bp", args.extend_bp),
        ("--extend-left-bp", args.extend_left_bp),
        ("--extend-right-bp", args.extend_right_bp),
        ("--extend-random-side-bp", args.extend_random_side_bp),
        ("--shift-bp", args.shift_bp),
    ):
        if value < 0:
            sys.exit(f"ERROR: {name} must be >= 0.")

    if args.seed is not None:
        random.seed(args.seed)

    if args.shuffles_per_fragment is not None or args.pseudocount_per_dinuc is not None:
        print(
            "[WARNING] --shuffles-per-fragment/--pseudocount-per-dinuc were supplied but are ignored; "
            "this script outputs pure observed frequencies only.",
            file=sys.stderr,
        )

    if args.flank is not None:
        print(
            "[WARNING] --flank was supplied, so the script is using old-style centre-flank mode. "
            "Omit --flank to profile the fragment itself plus extension.",
            file=sys.stderr,
        )

    bam_files = expand_bam_inputs(args.bam)

    print("[INFO] BAM files:", file=sys.stderr)
    for b in bam_files:
        print(f"  {b}", file=sys.stderr)

    fasta = pysam.FastaFile(args.fasta)
    fasta_refs = set(fasta.references)

    requested_chroms = expand_chroms(args.chroms)

    base_left_ext = args.extend_bp + args.extend_left_bp
    base_right_ext = args.extend_bp + args.extend_right_bp

    # All fragments have the same final profiled length because random-side extension
    # changes which side is extended, not the total extension amount.
    if args.flank is None:
        profile_len = args.length + base_left_ext + base_right_ext + args.extend_random_side_bp
        profile_basis = "fragment"
    else:
        # Old flank mode: dinucleotide starts -flank..+flank require 2*flank+2 bases,
        # then extension adds bases to the profile interval.
        profile_len = (2 * args.flank + 2) + base_left_ext + base_right_ext + args.extend_random_side_bp
        profile_basis = "flank"

    if profile_len < 2:
        sys.exit("ERROR: final profile interval length must be >= 2 for dinucleotide profiles.")

    if args.out is None:
        args.out = auto_output_path(
            args=args,
            bam_files=bam_files,
            profile_len=profile_len,
            profile_basis=profile_basis,
        )
        print(f"[INFO] auto-generated output name: {args.out}", file=sys.stderr)
    else:
        print(f"[INFO] output name supplied by --out: {args.out}", file=sys.stderr)

    out_parent = os.path.dirname(args.out)
    if out_parent:
        os.makedirs(out_parent, exist_ok=True)

    positions = profile_position_range(profile_len)
    n_positions = len(positions)

    dinuc_counts_by_position = [defaultdict(int) for _ in range(n_positions)]
    valid_dinuc_opportunities = [0] * n_positions

    total_fragments_seen = 0
    total_fragments_used = 0
    total_boundary_skipped = 0
    total_missing_fasta_contig = 0
    total_profile_n_skipped = 0
    total_fragments_with_no_valid_dinucs = 0

    require_acgt_profile = not args.allow_n_in_profile

    print(f"[INFO] profile basis: {profile_basis}", file=sys.stderr)
    print(f"[INFO] requested fragment length: {args.length}", file=sys.stderr)
    print(
        f"[INFO] extension: both={args.extend_bp}; left/start={args.extend_left_bp}; "
        f"right/end={args.extend_right_bp}; random-side={args.extend_random_side_bp}",
        file=sys.stderr,
    )
    print(f"[INFO] shift: {args.shift_bp} bp direction={args.shift_direction}", file=sys.stderr)
    print(
        f"[INFO] final profiled interval length: {profile_len} bp; "
        f"dinucleotide start positions: {positions[0]}..{positions[-1]} ({n_positions} positions)",
        file=sys.stderr,
    )

    for bam_path in bam_files:
        print(f"[INFO] Opening BAM: {bam_path}", file=sys.stderr)

        bam = pysam.AlignmentFile(bam_path, "rb")
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
            fa_chrom = matching_contig(bam_chrom, fasta_refs)

            if fa_chrom is None:
                print(
                    f"[WARNING] No matching FASTA contig for BAM contig {bam_chrom}; skipping.",
                    file=sys.stderr,
                )
                total_missing_fasta_contig += 1
                continue

            chrom_len = min(
                bam.get_reference_length(bam_chrom),
                fasta.get_reference_length(fa_chrom),
            )

            print(f"[INFO] Processing {os.path.basename(bam_path)} {bam_chrom}", file=sys.stderr)

            for frag_start, frag_end in iter_fragments_from_bam_chrom(
                bam=bam,
                bam_chrom=bam_chrom,
                frag_len_interest=args.length,
                mapq=args.mapq,
                include_duplicates=args.include_duplicates,
                allow_improper_pairs=args.allow_improper_pairs,
                allow_softclipped=args.allow_softclipped,
                max_duplicates=args.max_duplicates,
                coord_counts=coord_counts,
            ):
                total_fragments_seen += 1

                left_ext, right_ext = choose_random_side_extension(
                    base_left_ext=base_left_ext,
                    base_right_ext=base_right_ext,
                    random_side_bp=args.extend_random_side_bp,
                )
                shift = choose_shift(args.shift_bp, args.shift_direction)

                if args.flank is None:
                    profile_start, profile_end = interval_for_fragment_profile(
                        frag_start=frag_start,
                        frag_end=frag_end,
                        left_ext=left_ext,
                        right_ext=right_ext,
                        shift=shift,
                    )
                else:
                    profile_start, profile_end = interval_for_flank_profile(
                        frag_start=frag_start,
                        frag_len=args.length,
                        flank=args.flank,
                        left_ext=left_ext,
                        right_ext=right_ext,
                        shift=shift,
                    )

                obs_profile_seq = get_profile_seq(
                    fasta=fasta,
                    fa_chrom=fa_chrom,
                    profile_start=profile_start,
                    profile_end=profile_end,
                    chrom_len=chrom_len,
                )

                if obs_profile_seq is None:
                    total_boundary_skipped += 1
                    continue

                if len(obs_profile_seq) != profile_len:
                    # This should not happen unless the interval construction is inconsistent.
                    print(
                        f"[WARNING] unexpected profile length {len(obs_profile_seq)}; expected {profile_len}; skipping.",
                        file=sys.stderr,
                    )
                    total_boundary_skipped += 1
                    continue

                added_positions = add_dinuc_profile_from_seq(
                    seq=obs_profile_seq,
                    dinuc_counts_by_position=dinuc_counts_by_position,
                    valid_dinuc_opportunities=valid_dinuc_opportunities,
                    require_acgt_profile=require_acgt_profile,
                )

                if added_positions is None:
                    total_profile_n_skipped += 1
                    continue

                if added_positions == 0:
                    total_fragments_with_no_valid_dinucs += 1
                    continue

                total_fragments_used += 1

                if total_fragments_used % 100000 == 0:
                    print(
                        f"[INFO] observed fragments used: {total_fragments_used:,}",
                        file=sys.stderr,
                    )

                if args.max_fragments is not None and total_fragments_used >= args.max_fragments:
                    break

            if args.max_fragments is not None and total_fragments_used >= args.max_fragments:
                break

        bam.close()

        if args.max_fragments is not None and total_fragments_used >= args.max_fragments:
            break

    fasta.close()

    if total_fragments_used == 0:
        sys.exit("ERROR: no observed fragments were successfully processed.")

    multiplier = 1.0 if args.fraction else 100.0
    value_suffix = "frac" if args.fraction else "pct"

    header = ["position", "n_valid"]
    header.extend([f"{dinuc}_{value_suffix}" for dinuc in DINUCS])
    header.extend([f"WW_{value_suffix}", f"SS_{value_suffix}"])

    with open(args.out, "w") as out:
        out.write("\t".join(header) + "\n")

        for i, rel_pos in enumerate(positions):
            n_valid = valid_dinuc_opportunities[i]
            counts = dinuc_counts_by_position[i]

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

    print(f"[DONE] observed fragments used: {total_fragments_used:,}", file=sys.stderr)
    print(f"[DONE] output: {args.out}", file=sys.stderr)
    print(f"[INFO] matching fragments seen before profile filters: {total_fragments_seen:,}", file=sys.stderr)
    print(f"[INFO] boundary-skipped observed fragments: {total_boundary_skipped:,}", file=sys.stderr)
    print(f"[INFO] profiles skipped due to N: {total_profile_n_skipped:,}", file=sys.stderr)
    print(f"[INFO] fragments with no valid dinucleotide positions: {total_fragments_with_no_valid_dinucs:,}", file=sys.stderr)
    print(f"[INFO] missing FASTA contigs skipped: {total_missing_fasta_contig}", file=sys.stderr)
    print(
        "[INFO] Dinucleotide position is the start of the dinucleotide relative to "
        "profile_start + floor(profile_interval_length / 2).",
        file=sys.stderr,
    )
    print(
        f"[INFO] output relative-position range: {positions[0]}..{positions[-1]}",
        file=sys.stderr,
    )


if __name__ == "__main__":
    main()
