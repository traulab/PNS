#!/usr/bin/env bash
set -euo pipefail

usage() {
  cat <<'EOF'
Usage:
  run_pns_dac.sh SCRIPT_DIR BAM MODE FRAG_LOWER FRAG_UPPER JOBS CHROMS [CHROM_SIZES] [DMAX] [VALUE_LIMIT] [MIN_REGION_LENGTH]

Example:
  run_pns_dac.sh \
    /mnt/c/git/beadsOnSpring/BeadsOnASpring \
    /mnt/c/Snyder_bams/CH01/BH01.bam \
    167 167 167 \
    10 \
    all

Example with selected chromosomes:
  run_pns_dac.sh \
    /mnt/c/git/beadsOnSpring/BeadsOnASpring \
    /mnt/c/Snyder_bams/CH01/BH01.bam \
    167 167 167 \
    10 \
    1,2,3,4,5,X,Y

Example with range:
  run_pns_dac.sh \
    /mnt/c/git/beadsOnSpring/BeadsOnASpring \
    /mnt/c/Snyder_bams/CH01/BH01.bam \
    167 167 167 \
    10 \
    1-22,X,Y
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

command -v python3 >/dev/null 2>&1 || { echo "ERROR: python3 not found"; exit 1; }
command -v zcat >/dev/null 2>&1 || { echo "ERROR: zcat not found"; exit 1; }
command -v wigToBigWig >/dev/null 2>&1 || { echo "ERROR: wigToBigWig not found in PATH"; exit 1; }

[[ -f "$PNS_SCRIPT" ]] || { echo "ERROR: PNS script not found: $PNS_SCRIPT"; exit 1; }
[[ -f "$DAC_SCRIPT" ]] || { echo "ERROR: DAC script not found: $DAC_SCRIPT"; exit 1; }
[[ -f "$BAM" ]] || { echo "ERROR: BAM not found: $BAM"; exit 1; }
[[ -f "$CHROM_SIZES" ]] || { echo "ERROR: chrom.sizes not found: $CHROM_SIZES"; exit 1; }

expand_chroms() {
  local spec="$1"

  spec="${spec// /}"

  if [[ "$spec" == "all" || "$spec" == "ALL" ]]; then
    seq 1 22
    echo X
    echo Y
    return
  fi

  IFS=',' read -ra TOKENS <<< "$spec"

  for token in "${TOKENS[@]}"; do
    token="${token// /}"
    token="${token#chr}"

    if [[ -z "$token" ]]; then
      continue
    fi

    if [[ "$token" =~ ^[0-9]+-[0-9]+$ ]]; then
      start="${token%-*}"
      end="${token#*-}"
      seq "$start" "$end"
    else
      echo "$token"
    fi
  done
}

mapfile -t REQUESTED_CHROMS < <(expand_chroms "$CHROM_SPEC")

declare -A WANT_CHROM=()

for c in "${REQUESTED_CHROMS[@]}"; do
  c="${c#chr}"
  WANT_CHROM["$c"]=1
done

CHROMS=()

for c in {1..22} X Y; do
  if [[ -n "${WANT_CHROM[$c]:-}" ]]; then
    CHROMS+=("$c")
  fi
done

if [[ ${#CHROMS[@]} -eq 0 ]]; then
  echo "ERROR: no valid chromosomes selected from: $CHROM_SPEC"
  echo "Use something like: all"
  echo "or: 1-22,X,Y"
  echo "or: 1,2,3,4,5"
  exit 1
fi

echo "Running PNS for sample: $SAMPLE"
echo "Chromosomes: ${CHROMS[*]}"
echo "Parallel jobs: $JOBS"
echo "Mode: $MODE"
echo "Fragment lower: $FRAG_LOWER"
echo "Fragment upper: $FRAG_UPPER"
echo

export PNS_SCRIPT BAM MODE FRAG_LOWER FRAG_UPPER SAMPLE

printf "%s\n" "${CHROMS[@]}" | xargs -P "$JOBS" -I {} bash -c '
  set -euo pipefail

  chrom="$1"
  out="${SAMPLE}_chr${chrom}_PNS"

  echo "[START] chr${chrom}"

  python3 "$PNS_SCRIPT" \
    -b "$BAM" \
    --mode "$MODE" \
    --frag-lower "$FRAG_LOWER" \
    --frag-upper "$FRAG_UPPER" \
    -c "$chrom" \
    --max-duplicates 0 \
    -o "$out" \
    --pns-mode off \
    --score-format wiggz \
    --score-tracks fragment_left_ends fragment_right_ends dyad

  echo "[DONE]  chr${chrom}"
' _ {}

echo
echo "All chromosome PNS jobs finished."
echo

TRACKS=(
  dyad
  fragment_right_ends
  fragment_left_ends
)

for track in "${TRACKS[@]}"; do
  merged_wig="${SAMPLE}_chrAll_PNS_mode${MODE}_lower${FRAG_LOWER}_upper${FRAG_UPPER}_${track}.wig"
  merged_bw="${SAMPLE}_chrAll_PNS_mode${MODE}_lower${FRAG_LOWER}_upper${FRAG_UPPER}_${track}.bw"

  files=()

  for chrom in "${CHROMS[@]}"; do
    f="${SAMPLE}_chr${chrom}_PNS_mode${MODE}_lower${FRAG_LOWER}_upper${FRAG_UPPER}_${track}.wig.gz"

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

  wigToBigWig "$merged_wig" "$CHROM_SIZES" "$merged_bw"

  echo "Running DAC for ${track}:"
  echo "  ${merged_bw}"

  python3 "$DAC_SCRIPT" \
    --bigwig "$merged_bw" \
    --chrom_sizes "$CHROM_SIZES" \
    --scope genome \
    --dmax "$DMAX" \
    --value_limit "$VALUE_LIMIT" \
    --min_region_length "$MIN_REGION_LENGTH" \
    --no_normalize_dac

  echo "[DONE] DAC for ${track}"
  echo
done

echo "Finished all PNS, bigWig conversion, and DAC analyses."