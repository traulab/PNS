#!/usr/bin/env bash
set -euo pipefail

usage() {
  cat <<'EOF'
Usage:
  run_pns_dac.sh SCRIPT_DIR BAM MODE FRAG_LOWER FRAG_UPPER JOBS CHROMS [CHROM_SIZES] [DMAX] [VALUE_LIMIT] [MIN_REGION_LENGTH]

CHROMS may be:

  all
  autosomes
  1,2,3,4,5,X,Y
  1-22,X,Y
  chr17:41186312-41287499
  17:41186312-41287499
  chr1:100000-200000,chr17:41186312-41287499

Whole-genome example:

  run_pns_dac.sh \
    /mnt/c/git/beadsOnSpring/BeadsOnASpring \
    /mnt/c/Snyder_bams/CH01/BH01.bam \
    167 166 168 \
    10 \
    all

Selected chromosomes:

  run_pns_dac.sh \
    /mnt/c/git/beadsOnSpring/BeadsOnASpring \
    /mnt/c/Snyder_bams/CH01/BH01.bam \
    167 166 168 \
    10 \
    1,2,3,4,5,X,Y

Chromosome range:

  run_pns_dac.sh \
    /mnt/c/git/beadsOnSpring/BeadsOnASpring \
    /mnt/c/Snyder_bams/CH01/BH01.bam \
    167 166 168 \
    10 \
    1-22,X,Y

Single genomic region:

  run_pns_dac.sh \
    /mnt/c/git/beadsOnSpring/BeadsOnASpring \
    /mnt/c/Snyder_bams/CH01/BH01.bam \
    167 166 168 \
    1 \
    chr17:41186312-41287499
EOF
}

if [[ $# -lt 7 ]]; then
  usage
  exit 1
fi

SCRIPT_DIR="${1%/}"
BAM="$2"
MODE="$3"
FRAG_LOWER="$4"
FRAG_UPPER="$5"
JOBS="$6"
CHROM_SPEC="$7"

CHROM_SIZES="${8:-/mnt/d/Snyder_bams/male/chrom.sizes}"
DMAX="${9:-3000}"
VALUE_LIMIT="${10:-50}"
MIN_REGION_LENGTH="${11:-200}"

PNS_SCRIPT="${SCRIPT_DIR}/PNS_with_nucleosome_peak_calling.py"
DAC_SCRIPT="${SCRIPT_DIR}/Analysis_scripts/multiplicative_weighted_DAC/multiplicative_weighted_DAC.py"

SAMPLE="$(basename "$BAM")"
SAMPLE="${SAMPLE%.bam}"

command -v python3 >/dev/null 2>&1 ||
  { echo "ERROR: python3 not found"; exit 1; }

command -v zcat >/dev/null 2>&1 ||
  { echo "ERROR: zcat not found"; exit 1; }

command -v wigToBigWig >/dev/null 2>&1 ||
  { echo "ERROR: wigToBigWig not found in PATH"; exit 1; }

[[ -f "$PNS_SCRIPT" ]] ||
  { echo "ERROR: PNS script not found: $PNS_SCRIPT"; exit 1; }

[[ -f "$DAC_SCRIPT" ]] ||
  { echo "ERROR: DAC script not found: $DAC_SCRIPT"; exit 1; }

[[ -f "$BAM" ]] ||
  { echo "ERROR: BAM not found: $BAM"; exit 1; }

[[ -f "$CHROM_SIZES" ]] ||
  { echo "ERROR: chrom.sizes not found: $CHROM_SIZES"; exit 1; }


# Output:
#
#   target<TAB>label
#
# Examples:
#
#   17<TAB>chr17
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
      printf '%s\t%s\n' "$chrom" "chr${chrom}"
    done
    return
  fi

  if [[ "${spec,,}" == "autosomes" ]]; then
    for chrom in {1..22}; do
      printf '%s\t%s\n' "$chrom" "chr${chrom}"
    done
    return
  fi  

  IFS=',' read -ra TOKENS <<< "$spec"

  for token in "${TOKENS[@]}"; do
    token="${token// /}"

    [[ -n "$token" ]] || continue

    # Remove an optional chr prefix.
    token="${token#chr}"

    # Genomic region, for example:
    #
    #   17:41186312-41287499
    #
    if [[ "$token" =~ ^([0-9]+|X|Y):([0-9]+)-([0-9]+)$ ]]; then
      chrom="${BASH_REMATCH[1]}"
      region_start="${BASH_REMATCH[2]}"
      region_end="${BASH_REMATCH[3]}"

      if (( region_start < 1 )); then
        echo "ERROR: region start must be at least 1: $token" >&2
        exit 1
      fi

      if (( region_end < region_start )); then
        echo "ERROR: region end is before region start: $token" >&2
        exit 1
      fi

      label="chr${chrom}_${region_start}_${region_end}"

      printf '%s\t%s\n' \
        "${chrom}:${region_start}-${region_end}" \
        "$label"

    # Chromosome range, for example:
    #
    #   1-22
    #
    elif [[ "$token" =~ ^([0-9]+)-([0-9]+)$ ]]; then
      start="${BASH_REMATCH[1]}"
      end="${BASH_REMATCH[2]}"

      if (( end < start )); then
        echo "ERROR: invalid chromosome range: $token" >&2
        exit 1
      fi

      for chrom in $(seq "$start" "$end"); do
        printf '%s\t%s\n' "$chrom" "chr${chrom}"
      done

    # Single chromosome.
    elif [[ "$token" =~ ^([0-9]+|X|Y)$ ]]; then
      printf '%s\t%s\n' "$token" "chr${token}"

    else
      echo "ERROR: invalid chromosome or region specification: $token" >&2
      exit 1
    fi
  done
}


mapfile -t TARGET_LINES < <(expand_targets "$CHROM_SPEC")

if [[ ${#TARGET_LINES[@]} -eq 0 ]]; then
  echo "ERROR: no valid chromosomes or regions selected from:"
  echo "  $CHROM_SPEC"
  exit 1
fi


# Remove duplicate entries while preserving their original order.
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


echo "Running PNS for sample: $SAMPLE"
echo "Parallel jobs: $JOBS"
echo "Mode: $MODE"
echo "Fragment lower: $FRAG_LOWER"
echo "Fragment upper: $FRAG_UPPER"
echo
echo "Targets:"

for i in "${!TARGETS[@]}"; do
  echo "  ${TARGETS[$i]}"
done

echo


export PNS_SCRIPT
export BAM
export MODE
export FRAG_LOWER
export FRAG_UPPER
export SAMPLE


# Pass each target and its safe output label to xargs.
for i in "${!TARGETS[@]}"; do
  printf '%s\t%s\n' "${TARGETS[$i]}" "${LABELS[$i]}"
done |
  xargs -P "$JOBS" -n 2 bash -c '
    set -euo pipefail

    target="$1"
    label="$2"

    out="${SAMPLE}_${label}_PNS"

    echo "[START] ${target}"

    python3 "$PNS_SCRIPT" \
      -b "$BAM" \
      --mode "$MODE" \
      --frag-lower "$FRAG_LOWER" \
      --frag-upper "$FRAG_UPPER" \
      -c "$target" \
      --max-duplicates 0 \
      -o "$out" \
      --pns-mode off \
      --score-format wiggz \
      --score-tracks \
        fragment_left_ends \
        fragment_right_ends \
        dyad

    echo "[DONE]  ${target}"
  ' _


echo
echo "All PNS jobs finished."
echo


TRACKS=(
  dyad
  fragment_right_ends
  fragment_left_ends
)


# Use chrAll for standard whole-genome/whole-chromosome requests.
# For a single region, include the region in the merged filename.
if [[ ${#TARGETS[@]} -eq 1 && "${TARGETS[0]}" == *:* ]]; then
  MERGED_LABEL="${LABELS[0]}"
else
  MERGED_LABEL="chrAll"
fi


for track in "${TRACKS[@]}"; do
  merged_wig="${SAMPLE}_${MERGED_LABEL}_PNS_mode${MODE}_lower${FRAG_LOWER}_upper${FRAG_UPPER}_${track}.wig"

  merged_bw="${SAMPLE}_${MERGED_LABEL}_PNS_mode${MODE}_lower${FRAG_LOWER}_upper${FRAG_UPPER}_${track}.bw"

  files=()

  for label in "${LABELS[@]}"; do
    f="${SAMPLE}_${label}_PNS_mode${MODE}_lower${FRAG_LOWER}_upper${FRAG_UPPER}_${track}.wig.gz"

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

  DAC_SCOPE_ARGS=(--scope genome)

  # When only one chromosome or one chromosome region was requested,
  # run DAC using chromosome scope.
  if [[ ${#TARGETS[@]} -eq 1 ]]; then
    target="${TARGETS[0]}"

    # Remove any region coordinates:
    # 17:41186312-41287499 -> 17
    dac_chrom="${target%%:*}"

    # Add the chr prefix expected by the DAC script.
    dac_chrom="chr${dac_chrom#chr}"

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

  echo "[DONE] DAC for ${track}"
  echo
done

echo "Finished all PNS, bigWig conversion, and DAC analyses."