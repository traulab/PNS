#!/usr/bin/env bash
set -euo pipefail

usage() {
  cat <<'USAGE'
Usage:
  run_pns_dac_multi.sh \
    SCRIPT_DIR \
    BAM_SPEC \
    OUT_PREFIX \
    MODE \
    FRAG_LOWER \
    FRAG_UPPER \
    JOBS \
    CHROMS \
    [CHROM_SIZES] \
    [DMAX] \
    [VALUE_LIMIT] \
    [MIN_REGION_LENGTH] \
    [CHROM_PREFIX]

Arguments:

  SCRIPT_DIR
      Base BeadsOnASpring script directory.

  BAM_SPEC
      Single BAM, quoted glob, or comma-separated BAM list.

      All resolved BAM files are treated as one combined BAM set.

      Examples:

        /mnt/c/Snyder_bams/CH01/BH01.bam

        '/mnt/d/Snyder_bams/Gaffney2012/bam/*.bam'

        '/path/a.bam,/path/b.bam'

  OUT_PREFIX
      Output prefix for the combined BAM set.

      Examples:

        Gaffney2012_allBams_147
        CH01_167
        Combined_samples_147

  MODE
      PNS mode length.

  FRAG_LOWER
      Lower fragment-length limit.

  FRAG_UPPER
      Upper fragment-length limit.

  JOBS
      Number of chromosome or region jobs to run in parallel.

  CHROMS
      Chromosomes or genomic regions to process.

      Supported forms:

        all
        autosomes
        1,2,3,4,5,X,Y
        1-22,X,Y
        chr17
        chr17:41186312-41287499
        17:41186312-41287499
        chr1:100000-200000,chr17:41186312-41287499

  CHROM_SIZES
      Default:

        /mnt/d/Snyder_bams/male/chrom.sizes

  DMAX
      Default: 3000

  VALUE_LIMIT
      Default: 50

  MIN_REGION_LENGTH
      Default: 200

  CHROM_PREFIX
      Prefix added to chromosome names when passing them to PNS
      and DAC.

      Default:

        chr

      Use:

        none

      if BAM, bigWig, and chrom.sizes contigs are named:

        1,2,...,X,Y

      rather than:

        chr1,chr2,...,chrX,chrY


Examples
========

All Gaffney BAMs combined, whole genome:

  run_pns_dac_multi.sh \
    /mnt/c/git/beadsOnSpring/BeadsOnASpring \
    '/mnt/d/Snyder_bams/Gaffney2012/bam/*.bam' \
    Gaffney2012_allBams_147 \
    147 147 147 \
    10 \
    all


All Gaffney BAMs combined, autosomes only:

  run_pns_dac_multi.sh \
    /mnt/c/git/beadsOnSpring/BeadsOnASpring \
    '/mnt/d/Snyder_bams/Gaffney2012/bam/*.bam' \
    Gaffney2012_allBams_147_autosomes \
    147 147 147 \
    10 \
    autosomes


All Gaffney BAMs combined, chromosome 20 only:

  run_pns_dac_multi.sh \
    /mnt/c/git/beadsOnSpring/BeadsOnASpring \
    '/mnt/d/Snyder_bams/Gaffney2012/bam/*.bam' \
    Gaffney2012_allBams_147_chr20 \
    147 147 147 \
    10 \
    chr20


All Gaffney BAMs combined, one genomic region:

  run_pns_dac_multi.sh \
    /mnt/c/git/beadsOnSpring/BeadsOnASpring \
    '/mnt/d/Snyder_bams/Gaffney2012/bam/*.bam' \
    Gaffney2012_allBams_147_BRCA1 \
    147 147 147 \
    1 \
    chr17:41186312-41287499


Comma-separated BAM list:

  run_pns_dac_multi.sh \
    /mnt/c/git/beadsOnSpring/BeadsOnASpring \
    '/path/sample1.bam,/path/sample2.bam' \
    sample1_sample2_167 \
    167 166 168 \
    10 \
    1-22,X,Y


BAM contigs are not chr-prefixed:

  run_pns_dac_multi.sh \
    /mnt/c/git/beadsOnSpring/BeadsOnASpring \
    '/mnt/d/Snyder_bams/Gaffney2012/bam/*.bam' \
    Gaffney2012_allBams_147 \
    147 147 147 \
    10 \
    1-22,X,Y \
    /mnt/d/Snyder_bams/Gaffney2012/ref/hg19/chrom.sizes \
    3000 \
    50 \
    200 \
    none

USAGE
}


if [[ $# -lt 8 ]]; then
  usage
  exit 1
fi


SCRIPT_DIR="${1%/}"
BAM_SPEC="$2"
OUT_PREFIX="$3"
MODE="$4"
FRAG_LOWER="$5"
FRAG_UPPER="$6"
JOBS="$7"
CHROM_SPEC="$8"

CHROM_SIZES="${9:-/mnt/d/Snyder_bams/male/chrom.sizes}"
DMAX="${10:-3000}"
VALUE_LIMIT="${11:-50}"
MIN_REGION_LENGTH="${12:-200}"
CHROM_PREFIX="${13:-chr}"


PNS_SCRIPT="${SCRIPT_DIR}/PNS_with_nucleosome_peak_calling.py"

DAC_SCRIPT="${SCRIPT_DIR}/Analysis_scripts/multiplicative_weighted_DAC/multiplicative_weighted_DAC.py"


command -v python3 >/dev/null 2>&1 ||
  {
    echo "ERROR: python3 not found"
    exit 1
  }

command -v zcat >/dev/null 2>&1 ||
  {
    echo "ERROR: zcat not found"
    exit 1
  }

command -v wigToBigWig >/dev/null 2>&1 ||
  {
    echo "ERROR: wigToBigWig not found in PATH"
    exit 1
  }


[[ -f "$PNS_SCRIPT" ]] ||
  {
    echo "ERROR: PNS script not found:"
    echo "  $PNS_SCRIPT"
    exit 1
  }

[[ -f "$DAC_SCRIPT" ]] ||
  {
    echo "ERROR: DAC script not found:"
    echo "  $DAC_SCRIPT"
    exit 1
  }

[[ -f "$CHROM_SIZES" ]] ||
  {
    echo "ERROR: chrom.sizes not found:"
    echo "  $CHROM_SIZES"
    exit 1
  }


if ! [[ "$JOBS" =~ ^[0-9]+$ ]] || (( JOBS < 1 )); then
  echo "ERROR: JOBS must be a positive integer."
  echo "Got: $JOBS"
  exit 1
fi


if ! [[ "$MODE" =~ ^[0-9]+$ ]]; then
  echo "ERROR: MODE must be a non-negative integer."
  echo "Got: $MODE"
  exit 1
fi


if ! [[ "$FRAG_LOWER" =~ ^[0-9]+$ ]]; then
  echo "ERROR: FRAG_LOWER must be a non-negative integer."
  echo "Got: $FRAG_LOWER"
  exit 1
fi


if ! [[ "$FRAG_UPPER" =~ ^[0-9]+$ ]]; then
  echo "ERROR: FRAG_UPPER must be a non-negative integer."
  echo "Got: $FRAG_UPPER"
  exit 1
fi


if (( FRAG_UPPER < FRAG_LOWER )); then
  echo "ERROR: FRAG_UPPER is smaller than FRAG_LOWER."
  echo "FRAG_LOWER: $FRAG_LOWER"
  echo "FRAG_UPPER: $FRAG_UPPER"
  exit 1
fi


case "${CHROM_PREFIX,,}" in
  none)
    CHROM_PREFIX=""
    ;;
  *)
    ;;
esac


resolve_bams() {
  local spec="$1"

  python3 - "$spec" <<'PY'
import glob
import os
import sys

spec = sys.argv[1]

seen = set()
resolved = []

for token in spec.split(","):
    token = token.strip()

    if not token:
        continue

    if any(character in token for character in "*?["):
        matches = sorted(glob.glob(token))

        if not matches:
            print(
                f"ERROR: BAM glob matched no files: {token}",
                file=sys.stderr,
            )
            sys.exit(1)

        for path in matches:
            if not os.path.isfile(path):
                continue

            if path not in seen:
                seen.add(path)
                resolved.append(path)

    else:
        if not os.path.isfile(token):
            print(
                f"ERROR: BAM not found: {token}",
                file=sys.stderr,
            )
            sys.exit(1)

        if token not in seen:
            seen.add(token)
            resolved.append(token)

if not resolved:
    print(
        f"ERROR: no BAM files resolved from: {spec}",
        file=sys.stderr,
    )
    sys.exit(1)

for path in resolved:
    print(path)
PY
}


# Output:
#
#   target<TAB>label
#
# target does not contain the CHROM_PREFIX.
#
# Examples:
#
#   17<TAB>chr17
#
#   17:41186312-41287499<TAB>chr17_41186312_41287499
#
expand_targets() {
  local spec="$1"
  local token
  local start
  local end
  local chrom
  local region_start
  local region_end
  local label

  spec="${spec// /}"

  if [[ "${spec,,}" == "all" ]]; then
    for chrom in {1..22} X Y; do
      printf '%s\t%s\n' \
        "$chrom" \
        "chr${chrom}"
    done

    return
  fi


  if [[ "${spec,,}" == "autosomes" ]]; then
    for chrom in {1..22}; do
      printf '%s\t%s\n' \
        "$chrom" \
        "chr${chrom}"
    done

    return
  fi


  IFS=',' read -ra TOKENS <<< "$spec"

  for token in "${TOKENS[@]}"; do
    token="${token// /}"

    [[ -n "$token" ]] || continue

    # Remove an optional chr prefix from user input.
    token="${token#chr}"
    token="${token#CHR}"


    # Genomic region:
    #
    #   17:41186312-41287499
    #
    if [[ "$token" =~ ^([0-9]+|X|Y):([0-9]+)-([0-9]+)$ ]]; then
      chrom="${BASH_REMATCH[1]}"
      region_start="${BASH_REMATCH[2]}"
      region_end="${BASH_REMATCH[3]}"

      if (( region_start < 1 )); then
        echo "ERROR: region start must be at least 1:" >&2
        echo "  $token" >&2
        exit 1
      fi

      if (( region_end < region_start )); then
        echo "ERROR: region end is before region start:" >&2
        echo "  $token" >&2
        exit 1
      fi

      if [[ "$chrom" =~ ^[0-9]+$ ]] &&
         (( chrom < 1 || chrom > 22 )); then
        echo "ERROR: invalid chromosome:" >&2
        echo "  $chrom" >&2
        exit 1
      fi

      label="chr${chrom}_${region_start}_${region_end}"

      printf '%s\t%s\n' \
        "${chrom}:${region_start}-${region_end}" \
        "$label"


    # Chromosome range:
    #
    #   1-22
    #
    elif [[ "$token" =~ ^([0-9]+)-([0-9]+)$ ]]; then
      start="${BASH_REMATCH[1]}"
      end="${BASH_REMATCH[2]}"

      if (( start < 1 || end > 22 )); then
        echo "ERROR: chromosome ranges must be within 1-22:" >&2
        echo "  $token" >&2
        exit 1
      fi

      if (( end < start )); then
        echo "ERROR: invalid chromosome range:" >&2
        echo "  $token" >&2
        exit 1
      fi

      for chrom in $(seq "$start" "$end"); do
        printf '%s\t%s\n' \
          "$chrom" \
          "chr${chrom}"
      done


    # Single chromosome.
    elif [[ "$token" =~ ^([0-9]+|X|Y)$ ]]; then
      chrom="$token"

      if [[ "$chrom" =~ ^[0-9]+$ ]] &&
         (( chrom < 1 || chrom > 22 )); then
        echo "ERROR: invalid chromosome:" >&2
        echo "  $chrom" >&2
        exit 1
      fi

      printf '%s\t%s\n' \
        "$chrom" \
        "chr${chrom}"


    else
      echo "ERROR: invalid chromosome or region specification:" >&2
      echo "  $token" >&2
      exit 1
    fi
  done
}


add_prefix_to_target() {
  local target="$1"
  local chrom
  local coordinates

  chrom="${target%%:*}"

  if [[ "$target" == *:* ]]; then
    coordinates="${target#*:}"

    printf '%s:%s\n' \
      "${CHROM_PREFIX}${chrom}" \
      "$coordinates"
  else
    printf '%s\n' \
      "${CHROM_PREFIX}${chrom}"
  fi
}


mapfile -t BAMS < <(
  resolve_bams "$BAM_SPEC"
)

mapfile -t TARGET_LINES < <(
  expand_targets "$CHROM_SPEC"
)


if [[ ${#TARGET_LINES[@]} -eq 0 ]]; then
  echo "ERROR: no valid chromosomes or regions selected from:"
  echo "  $CHROM_SPEC"
  exit 1
fi


# Remove duplicate targets while preserving input order.
declare -A SEEN_TARGETS=()

TARGETS=()
LABELS=()

for line in "${TARGET_LINES[@]}"; do
  IFS=$'\t' read -r target label <<< "$line"

  if [[ -z "${SEEN_TARGETS[$target]:-}" ]]; then
    SEEN_TARGETS["$target"]=1

    TARGETS+=("$target")
    LABELS+=("$label")
  fi
done


echo "Resolved BAM files for one combined analysis:"
echo "Number of BAM files: ${#BAMS[@]}"

printf '  %s\n' "${BAMS[@]}"

echo
echo "Output prefix: $OUT_PREFIX"
echo "Parallel jobs: $JOBS"
echo "Mode: $MODE"
echo "Fragment lower: $FRAG_LOWER"
echo "Fragment upper: $FRAG_UPPER"

if [[ -n "$CHROM_PREFIX" ]]; then
  echo "Chromosome prefix passed to PNS and DAC: $CHROM_PREFIX"
else
  echo "Chromosome prefix passed to PNS and DAC: none"
fi

echo
echo "Targets:"

for target in "${TARGETS[@]}"; do
  pns_target="$(add_prefix_to_target "$target")"

  echo "  $pns_target"
done

echo


TRACKS=(
  dyad
  fragment_right_ends
  fragment_left_ends
)


BAM_LIST_FILE="$(mktemp)"

cleanup() {
  rm -f "$BAM_LIST_FILE"
}

trap cleanup EXIT


printf '%s\n' "${BAMS[@]}" > "$BAM_LIST_FILE"


export PNS_SCRIPT
export BAM_LIST_FILE
export OUT_PREFIX
export MODE
export FRAG_LOWER
export FRAG_UPPER
export CHROM_PREFIX


echo "============================================================"
echo "Running combined PNS analysis"
echo "Output prefix: $OUT_PREFIX"
echo "============================================================"
echo


# Pass each target and safe filename label to xargs.
for i in "${!TARGETS[@]}"; do
  printf '%s\t%s\n' \
    "${TARGETS[$i]}" \
    "${LABELS[$i]}"
done |
  xargs -P "$JOBS" -n 2 bash -c '
    set -euo pipefail

    target="$1"
    label="$2"

    mapfile -t bam_args < "$BAM_LIST_FILE"

    chrom="${target%%:*}"

    if [[ "$target" == *:* ]]; then
      coordinates="${target#*:}"
      pns_target="${CHROM_PREFIX}${chrom}:${coordinates}"
    else
      pns_target="${CHROM_PREFIX}${chrom}"
    fi

    out="${OUT_PREFIX}_${label}_PNS"

    echo "[START] ${OUT_PREFIX} ${pns_target}"
    echo "        BAMs passed to PNS: ${#bam_args[@]}"

    python3 "$PNS_SCRIPT" \
      -b "${bam_args[@]}" \
      --mode-length "$MODE" \
      --frag-lower "$FRAG_LOWER" \
      --frag-upper "$FRAG_UPPER" \
      -c "$pns_target" \
      --max-duplicates 0 \
      -o "$out" \
      --pns-mode off \
      --score-format wiggz \
      --score-tracks \
        fragment_left_ends \
        fragment_right_ends \
        dyad

    echo "[DONE]  ${OUT_PREFIX} ${pns_target}"
  ' _


echo
echo "All PNS jobs finished for combined BAM set:"
echo "  $OUT_PREFIX"
echo


# Use the region label when exactly one genomic region was requested.
#
# Otherwise use chrAll, including when exactly one whole chromosome
# was requested, to retain the existing naming convention.
if [[ ${#TARGETS[@]} -eq 1 &&
      "${TARGETS[0]}" == *:* ]]; then

  MERGED_LABEL="${LABELS[0]}"
else
  MERGED_LABEL="chrAll"
fi


for track in "${TRACKS[@]}"; do
  merged_wig="${OUT_PREFIX}_${MERGED_LABEL}_PNS_mode${MODE}_lower${FRAG_LOWER}_upper${FRAG_UPPER}_${track}.wig"

  merged_bw="${OUT_PREFIX}_${MERGED_LABEL}_PNS_mode${MODE}_lower${FRAG_LOWER}_upper${FRAG_UPPER}_${track}.bw"

  files=()

  for label in "${LABELS[@]}"; do
    f="${OUT_PREFIX}_${label}_PNS_mode${MODE}_lower${FRAG_LOWER}_upper${FRAG_UPPER}_${track}.wig.gz"

    if [[ ! -f "$f" ]]; then
      echo "ERROR: expected output file not found:"
      echo "  $f"
      exit 1
    fi

    files+=("$f")
  done


  echo "Merging ${track} WIG files:"
  echo "  ${merged_wig}"

  zcat "${files[@]}" > "$merged_wig"


  echo "Converting ${track} WIG to bigWig:"
  echo "  ${merged_bw}"

  wigToBigWig \
    "$merged_wig" \
    "$CHROM_SIZES" \
    "$merged_bw"


  echo "Running DAC for ${track}:"
  echo "  ${merged_bw}"


  DAC_SCOPE_ARGS=(
    --scope genome
  )


  # When exactly one chromosome or one genomic region was requested,
  # run DAC using chromosome scope.
  if [[ ${#TARGETS[@]} -eq 1 ]]; then
    target="${TARGETS[0]}"

    # Remove region coordinates:
    #
    #   17:41186312-41287499
    #
    # becomes:
    #
    #   17
    #
    dac_chrom="${target%%:*}"
    dac_chrom="${CHROM_PREFIX}${dac_chrom}"

    DAC_SCOPE_ARGS=(
      --scope chromosome
      --chromosome "$dac_chrom"
    )
  fi


  python3 "$DAC_SCRIPT" \
    --bigwig "$merged_bw" \
    --chrom_sizes "$CHROM_SIZES" \
    "${DAC_SCOPE_ARGS[@]}" \
    --dmax "$DMAX" \
    --value_limit "$VALUE_LIMIT" \
    --min_region_length "$MIN_REGION_LENGTH" \
    --no_normalize_dac


  echo "[DONE] DAC for ${OUT_PREFIX} ${track}"
  echo
done


echo "Finished combined PNS, bigWig conversion, and DAC analyses for:"
echo "  $OUT_PREFIX"