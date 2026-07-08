#!/usr/bin/env bash
set -euo pipefail

usage() {
  cat <<'USAGE'
Usage:
  run_pns_dac_multi.sh SCRIPT_DIR BAM_SPEC OUT_PREFIX MODE FRAG_LOWER FRAG_UPPER JOBS CHROMS [CHROM_SIZES] [DMAX] [VALUE_LIMIT] [MIN_REGION_LENGTH] [CHROM_PREFIX]

Arguments:
  SCRIPT_DIR          Base BeadsOnASpring script directory.

  BAM_SPEC            Single BAM, quoted glob, or comma-separated BAM list.
                      This is treated as ONE combined BAM set.
                      Examples:
                        /mnt/c/Snyder_bams/CH01/BH01.bam
                        '/mnt/d/Snyder_bams/Gaffney2012/bam/*.bam'
                        '/path/a.bam,/path/b.bam'

  OUT_PREFIX          Short output prefix for the combined BAM set.
                      Example:
                        Gaffney2012_147
                        CH01_167
                        Gaffney2012_allBams_147

  MODE                PNS mode length.
  FRAG_LOWER          Lower fragment length.
  FRAG_UPPER          Upper fragment length.
  JOBS                Number of chromosomes to process in parallel.
  CHROMS              all, 1-22,X,Y, chr1,chr2,chrX, etc.
  CHROM_SIZES         Default: /mnt/d/Snyder_bams/male/chrom.sizes
  DMAX                Default: 3000
  VALUE_LIMIT         Default: 50
  MIN_REGION_LENGTH   Default: 200
  CHROM_PREFIX        Prefix used when passing chromosomes to PNS.
                      Default: chr
                      Use none if BAM contigs are named 1,2,...,X,Y.

Example, all Gaffney BAMs combined:
  run_pns_dac_multi.sh \
    /mnt/c/git/beadsOnSpring/BeadsOnASpring \
    '/mnt/d/Snyder_bams/Gaffney2012/bam/*.bam' \
    Gaffney2012_allBams_147 \
    147 147 147 \
    10 \
    1-22,X,Y

Example, chr20 only:
  run_pns_dac_multi.sh \
    /mnt/c/git/beadsOnSpring/BeadsOnASpring \
    '/mnt/d/Snyder_bams/Gaffney2012/bam/*.bam' \
    Gaffney2012_allBams_147_chr20 \
    147 147 147 \
    10 \
    chr20

Example, BAM contigs are not chr-prefixed:
  run_pns_dac_multi.sh \
    /mnt/c/git/beadsOnSpring/BeadsOnASpring \
    '/mnt/d/Snyder_bams/Gaffney2012/bam/*.bam' \
    Gaffney2012_allBams_147 \
    147 147 147 \
    10 \
    1-22,X,Y \
    /mnt/d/Snyder_bams/male/chrom.sizes \
    3000 50 200 none
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

command -v python3 >/dev/null 2>&1 || { echo "ERROR: python3 not found"; exit 1; }
command -v zcat >/dev/null 2>&1 || { echo "ERROR: zcat not found"; exit 1; }
command -v wigToBigWig >/dev/null 2>&1 || { echo "ERROR: wigToBigWig not found in PATH"; exit 1; }

[[ -f "$PNS_SCRIPT" ]] || { echo "ERROR: PNS script not found: $PNS_SCRIPT"; exit 1; }
[[ -f "$DAC_SCRIPT" ]] || { echo "ERROR: DAC script not found: $DAC_SCRIPT"; exit 1; }
[[ -f "$CHROM_SIZES" ]] || { echo "ERROR: chrom.sizes not found: $CHROM_SIZES"; exit 1; }

if ! [[ "$JOBS" =~ ^[0-9]+$ ]] || [[ "$JOBS" -lt 1 ]]; then
  echo "ERROR: JOBS must be a positive integer. Got: $JOBS"
  exit 1
fi

resolve_bams() {
  local spec="$1"

  python3 - "$spec" <<'PY'
import sys, glob, os

spec = sys.argv[1]
seen = set()
out = []

for token in spec.split(","):
    token = token.strip()
    if not token:
        continue

    if any(ch in token for ch in "*?["):
        matches = sorted(glob.glob(token))
        if not matches:
            print(f"ERROR: BAM glob matched no files: {token}", file=sys.stderr)
            sys.exit(1)

        for f in matches:
            if os.path.isfile(f) and f not in seen:
                seen.add(f)
                out.append(f)
    else:
        if not os.path.isfile(token):
            print(f"ERROR: BAM not found: {token}", file=sys.stderr)
            sys.exit(1)

        if token not in seen:
            seen.add(token)
            out.append(token)

if not out:
    print(f"ERROR: no BAM files resolved from: {spec}", file=sys.stderr)
    sys.exit(1)

for f in out:
    print(f)
PY
}

expand_chroms() {
  local spec="$1"
  local token=""
  local start=""
  local end=""

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
    token="${token#CHR}"

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

mapfile -t BAMS < <(resolve_bams "$BAM_SPEC")
mapfile -t REQUESTED_CHROMS < <(expand_chroms "$CHROM_SPEC")

declare -A WANT_CHROM=()

for c in "${REQUESTED_CHROMS[@]}"; do
  c="${c#chr}"
  c="${c#CHR}"
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
  echo "or: chr1,chr2,chrX"
  exit 1
fi

if [[ "$CHROM_PREFIX" == "none" || "$CHROM_PREFIX" == "NONE" ]]; then
  PNS_CHROM_EXAMPLE="${CHROMS[0]}"
else
  PNS_CHROM_EXAMPLE="${CHROM_PREFIX}${CHROMS[0]}"
fi

echo "Resolved BAM files for ONE combined analysis:"
echo "Number of BAM files: ${#BAMS[@]}"
printf '  %s\n' "${BAMS[@]}"
echo
echo "Output prefix: $OUT_PREFIX"
echo "Chromosomes requested: ${CHROMS[*]}"
echo "Chromosome argument passed to PNS will look like: ${PNS_CHROM_EXAMPLE}"
echo "Parallel chromosome jobs: $JOBS"
echo "Mode: $MODE"
echo "Fragment lower: $FRAG_LOWER"
echo "Fragment upper: $FRAG_UPPER"
echo

TRACKS=(
  dyad
  fragment_right_ends
  fragment_left_ends
)

BAM_LIST_FILE="$(mktemp)"
trap 'rm -f "$BAM_LIST_FILE"' EXIT

printf '%s\n' "${BAMS[@]}" > "$BAM_LIST_FILE"

export PNS_SCRIPT
export BAM_LIST_FILE
export OUT_PREFIX
export MODE
export FRAG_LOWER
export FRAG_UPPER
export CHROM_PREFIX

echo "============================================================"
echo "Running combined PNS/DAC analysis"
echo "Output prefix: $OUT_PREFIX"
echo "============================================================"
echo

printf "%s\n" "${CHROMS[@]}" | xargs -P "$JOBS" -I {} bash -c '
  set -euo pipefail

  chrom="$1"
  mapfile -t bam_args < "$BAM_LIST_FILE"

  out="${OUT_PREFIX}_chr${chrom}_PNS"

  if [[ "$CHROM_PREFIX" == "none" || "$CHROM_PREFIX" == "NONE" ]]; then
    pns_chrom="$chrom"
  else
    pns_chrom="${CHROM_PREFIX}${chrom}"
  fi

  echo "[START] ${OUT_PREFIX} ${pns_chrom}"
  echo "        BAMs passed to PNS: ${#bam_args[@]}"

  python3 "$PNS_SCRIPT" \
    -b "${bam_args[@]}" \
    --mode-length "$MODE" \
    --frag-lower "$FRAG_LOWER" \
    --frag-upper "$FRAG_UPPER" \
    -c "$pns_chrom" \
    --max-duplicates 0 \
    -o "$out" \
    --pns-mode off \
    --score-format wiggz \
    --score-tracks fragment_left_ends fragment_right_ends dyad

  echo "[DONE]  ${OUT_PREFIX} ${pns_chrom}"
' _ {}

echo
echo "All chromosome PNS jobs finished for combined BAM set: $OUT_PREFIX"
echo

for track in "${TRACKS[@]}"; do
  merged_wig="${OUT_PREFIX}_chrAll_PNS_mode${MODE}_lower${FRAG_LOWER}_upper${FRAG_UPPER}_${track}.wig"
  merged_bw="${OUT_PREFIX}_chrAll_PNS_mode${MODE}_lower${FRAG_LOWER}_upper${FRAG_UPPER}_${track}.bw"

  files=()

  for chrom in "${CHROMS[@]}"; do
    f="${OUT_PREFIX}_chr${chrom}_PNS_mode${MODE}_lower${FRAG_LOWER}_upper${FRAG_UPPER}_${track}.wig.gz"

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

  echo "[DONE] DAC for ${OUT_PREFIX} ${track}"
  echo
done

echo "Finished combined PNS, bigWig conversion, and DAC analyses for: $OUT_PREFIX"