#!/usr/bin/env bash
set -euo pipefail

usage() {
  cat <<'USAGE'
Usage:
  run_pns_dac_ww_types.sh SCRIPT_DIR BAM_OR_QUOTED_GLOB CHROMS [PNS arguments...]

The wrapper reserves only the first three arguments:

  1. SCRIPT_DIR
  2. BAM file or quoted BAM glob
  3. CHROMS

Every remaining argument is forwarded unchanged to
pns_with_dinuc_ww_types.py.

The wrapper itself supplies:

  -b / --bamfiles
  -c / --contigs
  -o / --out_prefix

Do not supply those options yourself.

CHROMS may be:

  all          Use every contig listed in the first BAM header.
  autosomes    Detect whether BAM contigs are 1..22 or chr1..chr22.
  1,2,3,4,5,X,Y
  1-22,X,Y
  chr17:41186312-41287499
  17:41186312-41287499

Example:

  bash run_pns_dac_ww_types.sh \
    /mnt/c/git/beadsOnSpring/BeadsOnASpring \
    "/mnt/d/Snyder_bams/CH01/*.bam" \
    autosomes \
    --mode 145 \
    --frag-lower 145 \
    --frag-upper 145 \
    --fasta /mnt/d/Snyder_bams/Gaffney2012/ref/hg19/hg19.fa \
    --split-ww-types \
    --dinuc-profile \
    --max-duplicates 1 \
    --dedup-scope all_bams \
    --pns-mode off \
    --score-format wiggz \
    --score-tracks dyad

Wrapper settings may be changed through environment variables:

  PNS_JOBS=5
  DAC_DMAX=3000
  DAC_VALUE_LIMIT=50
  DAC_MIN_REGION_LENGTH=200

Resume/force settings:

  The wrapper skips completed PNS targets, merged profiles, bigWigs, and DAC
  runs. Set any of these to 1 to rerun a stage:

  FORCE_ALL=1
  FORCE_PNS=1
  FORCE_MERGE=1
  FORCE_DAC=1
USAGE
}

if [[ $# -lt 4 ]]; then
  usage
  exit 1
fi

SCRIPT_DIR="${1%/}"
BAM_SPEC="$2"
CHROM_SPEC="$3"
shift 3

# Bash preserves each remaining command-line word as a separate array element.
# For example, --dedup-scope and all_bams remain two adjacent elements.
PNS_ARGS=("$@")

JOBS="${PNS_JOBS:-5}"
DMAX="${DAC_DMAX:-3000}"
VALUE_LIMIT="${DAC_VALUE_LIMIT:-50}"
MIN_REGION_LENGTH="${DAC_MIN_REGION_LENGTH:-200}"

FORCE_ALL="${FORCE_ALL:-0}"
FORCE_PNS="${FORCE_PNS:-0}"
FORCE_MERGE="${FORCE_MERGE:-0}"
FORCE_DAC="${FORCE_DAC:-0}"

PNS_SCRIPT="${SCRIPT_DIR}/pns_with_dinuc_ww_types.py"
DAC_SCRIPT="${SCRIPT_DIR}/Analysis_scripts/multiplicative_weighted_DAC/multiplicative_weighted_DAC.py"

# Expand either a literal BAM path or a quoted shell glob.
BAM_FILES=()

if [[ -f "$BAM_SPEC" ]]; then
  BAM_FILES=("$BAM_SPEC")
else
  while IFS= read -r bam_path; do
    BAM_FILES+=("$bam_path")
  done < <(compgen -G "$BAM_SPEC" | sort)
fi

if [[ ${#BAM_FILES[@]} -eq 0 ]]; then
  echo "ERROR: no BAM files matched:" >&2
  echo "  $BAM_SPEC" >&2
  exit 1
fi

for bam_path in "${BAM_FILES[@]}"; do
  [[ -f "$bam_path" ]] || {
    echo "ERROR: BAM not found: $bam_path" >&2
    exit 1
  }
done

if [[ ${#BAM_FILES[@]} -eq 1 ]]; then
  SAMPLE="$(basename "${BAM_FILES[0]}")"
  SAMPLE="${SAMPLE%.bam}"
else
  FIRST_BAM_DIR="$(dirname "${BAM_FILES[0]}")"
  SAMPLE="$(basename "$FIRST_BAM_DIR")"
fi

command -v python3 >/dev/null 2>&1 ||
  { echo "ERROR: python3 not found"; exit 1; }

command -v samtools >/dev/null 2>&1 ||
  { echo "ERROR: samtools not found"; exit 1; }

command -v sha256sum >/dev/null 2>&1 ||
  { echo "ERROR: sha256sum not found"; exit 1; }

command -v zcat >/dev/null 2>&1 ||
  { echo "ERROR: zcat not found"; exit 1; }

command -v wigToBigWig >/dev/null 2>&1 ||
  { echo "ERROR: wigToBigWig not found in PATH"; exit 1; }

[[ -f "$PNS_SCRIPT" ]] ||
  { echo "ERROR: PNS script not found: $PNS_SCRIPT"; exit 1; }

[[ -f "$DAC_SCRIPT" ]] ||
  { echo "ERROR: DAC script not found: $DAC_SCRIPT"; exit 1; }

if ! [[ "$JOBS" =~ ^[1-9][0-9]*$ ]]; then
  echo "ERROR: JOBS must be a positive integer: $JOBS"
  exit 1
fi

for force_var_name in FORCE_ALL FORCE_PNS FORCE_MERGE FORCE_DAC; do
  force_value="${!force_var_name}"
  if [[ "$force_value" != "0" && "$force_value" != "1" ]]; then
    echo "ERROR: ${force_var_name} must be 0 or 1, not: ${force_value}" >&2
    exit 1
  fi
done

if [[ "$FORCE_ALL" == "1" ]]; then
  FORCE_PNS=1
  FORCE_MERGE=1
  FORCE_DAC=1
fi

# The wrapper owns these arguments because they vary for each parallel target.
for arg in "${PNS_ARGS[@]}"; do
  case "$arg" in
    -b|--bamfiles|-b=*|--bamfiles=*|-c|--contigs|-c=*|--contigs=*|-o|--out_prefix|-o=*|--out_prefix=*)
      echo "ERROR: do not include BAM, contig, or output-prefix options in PNS_ARGUMENTS:" >&2
      echo "  $arg" >&2
      exit 1
      ;;
  esac
done

has_pns_flag() {
  local wanted="$1"
  local arg
  for arg in "${PNS_ARGS[@]}"; do
    [[ "$arg" == "$wanted" ]] && return 0
  done
  return 1
}

get_pns_value() {
  local wanted="$1"
  local default_value="$2"
  local i arg

  for ((i = 0; i < ${#PNS_ARGS[@]}; i++)); do
    arg="${PNS_ARGS[$i]}"
    if [[ "$arg" == "$wanted" ]]; then
      if (( i + 1 >= ${#PNS_ARGS[@]} )); then
        echo "ERROR: $wanted has no value in PNS_ARGUMENTS." >&2
        exit 1
      fi
      printf '%s\n' "${PNS_ARGS[$((i + 1))]}"
      return
    elif [[ "$arg" == "$wanted="* ]]; then
      printf '%s\n' "${arg#*=}"
      return
    fi
  done

  printf '%s\n' "$default_value"
}

# These values are only read so the wrapper can locate the filenames produced
# by the Python script. They are still forwarded unchanged to the Python script.
MODE="$(get_pns_value --mode "$(get_pns_value --mode-length 167)")"
FRAG_LOWER="$(get_pns_value --frag-lower 137)"
FRAG_UPPER="$(get_pns_value --frag-upper 197)"
FASTA="$(get_pns_value --fasta '')"
SCORE_FORMAT="$(get_pns_value --score-format wiggz)"

if ! has_pns_flag --split-ww-types; then
  echo "ERROR: this wrapper requires --split-ww-types in PNS_ARGUMENTS." >&2
  exit 1
fi

if ! has_pns_flag --dinuc-profile; then
  echo "ERROR: this wrapper requires --dinuc-profile in PNS_ARGUMENTS." >&2
  exit 1
fi

if [[ -z "$FASTA" ]]; then
  echo "ERROR: --fasta is required in PNS_ARGUMENTS." >&2
  exit 1
fi

[[ -f "$FASTA" ]] || {
  echo "ERROR: FASTA not found: $FASTA" >&2
  exit 1
}

[[ -f "${FASTA}.fai" ]] || {
  echo "ERROR: FASTA index not found: ${FASTA}.fai" >&2
  exit 1
}

if [[ "$SCORE_FORMAT" != "wiggz" && "$SCORE_FORMAT" != "both" ]]; then
  echo "ERROR: --score-format must be wiggz or both so dyad WIG.GZ files exist." >&2
  exit 1
fi

# If --score-tracks was supplied explicitly, ensure dyad is among its values.
score_tracks_seen=false
dyad_requested=false
for ((i = 0; i < ${#PNS_ARGS[@]}; i++)); do
  if [[ "${PNS_ARGS[$i]}" == "--score-tracks" ]]; then
    score_tracks_seen=true
    for ((j = i + 1; j < ${#PNS_ARGS[@]}; j++)); do
      [[ "${PNS_ARGS[$j]}" == --* ]] && break
      [[ "${PNS_ARGS[$j]}" == "dyad" ]] && dyad_requested=true
    done
  fi
done

if $score_tracks_seen && ! $dyad_requested; then
  echo "ERROR: --score-tracks was supplied but does not include dyad." >&2
  exit 1
fi


# Read the BAM contig dictionary as:
#
#   contig<TAB>length
read_bam_contigs() {
  local bam="$1"

  samtools view -H "$bam" |
    awk -F '\t' '
      $1 == "@SQ" {
        sn = ""
        ln = ""

        for (i = 2; i <= NF; i++) {
          if ($i ~ /^SN:/) {
            sn = substr($i, 4)
          } else if ($i ~ /^LN:/) {
            ln = substr($i, 4)
          }
        }

        if (sn != "" && ln != "") {
          print sn "\t" ln
        }
      }
    '
}


mapfile -t BAM_CONTIG_LINES < <(read_bam_contigs "${BAM_FILES[0]}")

if [[ ${#BAM_CONTIG_LINES[@]} -eq 0 ]]; then
  echo "ERROR: no @SQ contigs found in BAM header:" >&2
  echo "  ${BAM_FILES[0]}" >&2
  exit 1
fi

declare -A BAM_CONTIG_LENGTH=()
BAM_CONTIGS=()

for line in "${BAM_CONTIG_LINES[@]}"; do
  IFS=$'\t' read -r contig length <<< "$line"
  BAM_CONTIGS+=("$contig")
  BAM_CONTIG_LENGTH["$contig"]="$length"
done


# All BAMs must have exactly the same contig names, order, and lengths.
FIRST_CONTIG_DICTIONARY="$(printf '%s\n' "${BAM_CONTIG_LINES[@]}")"

for bam_path in "${BAM_FILES[@]:1}"; do
  CURRENT_CONTIG_DICTIONARY="$(read_bam_contigs "$bam_path")"

  if [[ "$CURRENT_CONTIG_DICTIONARY" != "$FIRST_CONTIG_DICTIONARY" ]]; then
    echo "ERROR: BAM contig dictionaries differ." >&2
    echo "All BAMs must use identical contig names, order, and lengths." >&2
    echo >&2
    echo "First BAM:" >&2
    echo "  ${BAM_FILES[0]}" >&2
    echo "Different BAM:" >&2
    echo "  $bam_path" >&2
    exit 1
  fi
done


# Detect whether a complete autosomal set is named 1..22 or chr1..chr22.
has_plain_autosomes=true
has_chr_autosomes=true

for chrom in {1..22}; do
  [[ -n "${BAM_CONTIG_LENGTH[$chrom]:-}" ]] || has_plain_autosomes=false
  [[ -n "${BAM_CONTIG_LENGTH[chr${chrom}]:-}" ]] || has_chr_autosomes=false
done

if $has_chr_autosomes && ! $has_plain_autosomes; then
  BAM_AUTOSOME_PREFIX="chr"
elif $has_plain_autosomes && ! $has_chr_autosomes; then
  BAM_AUTOSOME_PREFIX=""
elif $has_chr_autosomes && $has_plain_autosomes; then
  echo "ERROR: BAM contains both 1..22 and chr1..chr22." >&2
  echo "The autosomal naming convention is ambiguous." >&2
  exit 1
else
  BAM_AUTOSOME_PREFIX=""
fi


# The Python script writes output chromosome names with a "chr" prefix,
# even when the BAM itself uses unprefixed contigs such as 1, 2, ..., 22.
#
# BAM access still uses the original BAM contig names. However, chrom.sizes
# must match the names written into the WIG/bigWig outputs. Therefore:
#
#   BAM contig 1     -> output contig chr1
#   BAM contig chr1  -> output contig chr1
#
# The same rule is applied to every BAM-header contig.
CHROM_SIZES="${SAMPLE}_output_chr.chrom.sizes"
CHROM_SIZES_TMP="${CHROM_SIZES}.tmp.$$"

printf '%s\n' "${BAM_CONTIG_LINES[@]}" |
  awk -F '\t' 'BEGIN { OFS="\t" }
    {
      chrom = $1
      if (chrom !~ /^chr/) {
        chrom = "chr" chrom
      }
      print chrom, $2
    }
  ' > "$CHROM_SIZES_TMP"

if [[ -f "$CHROM_SIZES" ]] && cmp -s "$CHROM_SIZES_TMP" "$CHROM_SIZES"; then
  rm -f "$CHROM_SIZES_TMP"
else
  mv -f "$CHROM_SIZES_TMP" "$CHROM_SIZES"
fi


safe_label() {
  local value="$1"
  value="${value//[^A-Za-z0-9_.-]/_}"
  printf '%s\n' "$value"
}


# Resolve chr20 versus 20 using the BAM dictionary.
resolve_chrom_name() {
  local requested="$1"
  local bare="${requested#chr}"
  local candidate

  if [[ -n "${BAM_CONTIG_LENGTH[$requested]:-}" ]]; then
    printf '%s\n' "$requested"
    return
  fi

  if [[ -n "${BAM_CONTIG_LENGTH[chr${bare}]:-}" ]]; then
    printf 'chr%s\n' "$bare"
    return
  fi

  if [[ -n "${BAM_CONTIG_LENGTH[$bare]:-}" ]]; then
    printf '%s\n' "$bare"
    return
  fi

  echo "ERROR: chromosome or contig not found in BAM header: $requested" >&2
  exit 1
}


# Output:
#
#   target<TAB>label
#
# target uses the exact BAM contig spelling.
expand_targets() {
  local spec="$1"
  local token
  local start
  local end
  local chrom
  local resolved
  local region_start
  local region_end
  local label

  spec="${spec// /}"

  if [[ "${spec,,}" == "all" ]]; then
    for chrom in "${BAM_CONTIGS[@]}"; do
      printf '%s\t%s\n' "$chrom" "$(safe_label "$chrom")"
    done
    return
  fi

  if [[ "${spec,,}" == "autosomes" ]]; then
    if ! $has_plain_autosomes && ! $has_chr_autosomes; then
      echo "ERROR: a complete autosome set was not found." >&2
      echo "Expected either 1..22 or chr1..chr22 in the BAM header." >&2
      exit 1
    fi

    for chrom in {1..22}; do
      resolved="${BAM_AUTOSOME_PREFIX}${chrom}"
      printf '%s\t%s\n' "$resolved" "$(safe_label "$resolved")"
    done
    return
  fi

  IFS=',' read -ra TOKENS <<< "$spec"

  for token in "${TOKENS[@]}"; do
    [[ -n "$token" ]] || continue

    # Region, e.g. 17:41186312-41287499 or chr17:...
    if [[ "$token" =~ ^([^:]+):([0-9]+)-([0-9]+)$ ]]; then
      chrom="${BASH_REMATCH[1]}"
      region_start="${BASH_REMATCH[2]}"
      region_end="${BASH_REMATCH[3]}"
      resolved="$(resolve_chrom_name "$chrom")"

      if (( region_start < 1 )); then
        echo "ERROR: region start must be at least 1: $token" >&2
        exit 1
      fi

      if (( region_end < region_start )); then
        echo "ERROR: region end is before region start: $token" >&2
        exit 1
      fi

      if (( region_end > BAM_CONTIG_LENGTH[$resolved] )); then
        echo "ERROR: region end exceeds contig length:" >&2
        echo "  region: $token" >&2
        echo "  contig length: ${BAM_CONTIG_LENGTH[$resolved]}" >&2
        exit 1
      fi

      label="$(safe_label "${resolved}_${region_start}_${region_end}")"
      printf '%s\t%s\n' "${resolved}:${region_start}-${region_end}" "$label"

    # Numeric chromosome range, e.g. 1-22.
    elif [[ "$token" =~ ^([0-9]+)-([0-9]+)$ ]]; then
      start="${BASH_REMATCH[1]}"
      end="${BASH_REMATCH[2]}"

      if (( end < start )); then
        echo "ERROR: invalid chromosome range: $token" >&2
        exit 1
      fi

      for chrom in $(seq "$start" "$end"); do
        resolved="$(resolve_chrom_name "$chrom")"
        printf '%s\t%s\n' "$resolved" "$(safe_label "$resolved")"
      done

    # Single chromosome or arbitrary contig.
    else
      resolved="$(resolve_chrom_name "$token")"
      printf '%s\t%s\n' "$resolved" "$(safe_label "$resolved")"
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


echo "Running WW-type dyad analysis for sample: $SAMPLE"
echo "BAM files (${#BAM_FILES[@]}):"
for bam_path in "${BAM_FILES[@]}"; do
  echo "  $bam_path"
done
echo "Reference FASTA: $FASTA"
echo "Detected autosome prefix: ${BAM_AUTOSOME_PREFIX:-<none>}"
echo "Output-track chrom.sizes: $CHROM_SIZES"
echo "Parallel jobs: $JOBS"
echo "Resume mode: enabled"
echo "Force stages: PNS=$FORCE_PNS MERGE=$FORCE_MERGE DAC=$FORCE_DAC"
echo "PNS arguments:"
printf '  %q' "${PNS_ARGS[@]}"
echo
echo "Filename mode/lower/upper: $MODE/$FRAG_LOWER/$FRAG_UPPER"
echo
echo "Targets:"

for i in "${!TARGETS[@]}"; do
  echo "  ${TARGETS[$i]}"
done

echo


export PNS_SCRIPT
BAM_LIST="$(printf '%s\n' "${BAM_FILES[@]}")"
export BAM_LIST
export SAMPLE
export MODE FRAG_LOWER FRAG_UPPER
export FORCE_PNS

PNS_ARGS_FILE="$(mktemp)"
printf '%s\0' "${PNS_ARGS[@]}" > "$PNS_ARGS_FILE"
export PNS_ARGS_FILE

cleanup() {
  rm -f "$PNS_ARGS_FILE"
}
trap cleanup EXIT


# Run one PNS job per selected chromosome or region. Each job writes:
#   all dyads (legacy unsuffixed output)
#   type1 dyads
#   type2 dyads
#   type3 dyads
#   type4 dyads
# as well as the corresponding dinucleotide profiles and fragment summaries.
for i in "${!TARGETS[@]}"; do
  printf '%s\t%s\n' "${TARGETS[$i]}" "${LABELS[$i]}"
done |
  xargs -P "$JOBS" -n 2 bash -c '
    set -euo pipefail

    target="$1"
    label="$2"
    out="${SAMPLE}_${label}_PNS"

    mapfile -t bam_args <<< "$BAM_LIST"
    mapfile -d "" -t pns_args < "$PNS_ARGS_FILE"

    output_base="${out}_mode${MODE}_lower${FRAG_LOWER}_upper${FRAG_UPPER}"
    done_marker="${output_base}_ww_types.done"

    required_outputs=(
      "${output_base}_dyad.wig.gz"
      "${output_base}_dinuc_profile.tsv"
      "${output_base}_type1_dyad.wig.gz"
      "${output_base}_type1_dinuc_profile.tsv"
      "${output_base}_type2_dyad.wig.gz"
      "${output_base}_type2_dinuc_profile.tsv"
      "${output_base}_type3_dyad.wig.gz"
      "${output_base}_type3_dinuc_profile.tsv"
      "${output_base}_type4_dyad.wig.gz"
      "${output_base}_type4_dinuc_profile.tsv"
    )

    pns_signature="$({
      printf "%s\0" "$target" "${bam_args[@]}" "${pns_args[@]}"
      for bam in "${bam_args[@]}"; do
        stat -c "%n\0%s\0%Y\0" "$bam"
      done
    } | sha256sum | awk "{print \$1}")"

    outputs_complete=true
    for required in "${required_outputs[@]}"; do
      if [[ ! -s "$required" ]]; then
        outputs_complete=false
        break
      fi
    done

    marker_matches=false
    if [[ -s "$done_marker" ]] && [[ "$(head -n 1 "$done_marker")" == "$pns_signature" ]]; then
      marker_matches=true
    fi

    if [[ "$FORCE_PNS" != "1" ]] && $outputs_complete; then
      if [[ -s "$done_marker" ]]; then
        if $marker_matches; then
          echo "[SKIP]  ${target} (matching PNS run already complete)"
          exit 0
        fi
        echo "[RERUN] ${target} (PNS arguments or BAM inputs changed)"
      else
        # Supports resuming analyses completed by an older wrapper that did not
        # yet create done markers. The current signature is adopted once.
        printf "%s\n" "$pns_signature" > "${done_marker}.tmp.$$"
        mv -f "${done_marker}.tmp.$$" "$done_marker"
        echo "[SKIP]  ${target} (existing PNS outputs adopted as complete)"
        exit 0
      fi
    fi

    echo "[START] ${target}"
    rm -f "$done_marker"

    python3 "$PNS_SCRIPT" \
      -b "${bam_args[@]}" \
      -c "$target" \
      -o "$out" \
      "${pns_args[@]}"

    for required in "${required_outputs[@]}"; do
      if [[ ! -s "$required" ]]; then
        echo "ERROR: PNS completed but expected output is missing or empty:" >&2
        echo "  $required" >&2
        exit 1
      fi
    done

    printf "%s\n" "$pns_signature" > "${done_marker}.tmp.$$"
    mv -f "${done_marker}.tmp.$$" "$done_marker"

    echo "[DONE]  ${target}"
  ' _


echo
echo "All PNS/type-splitting jobs finished."
echo


# The all group uses the legacy unsuffixed filename. Type-specific groups use
# _type1 through _type4 immediately before _dyad.wig.gz.
DYAD_GROUPS=(
  all
  type1
  type2
  type3
  type4
)


# Use chrAll for standard whole-genome/whole-chromosome requests.
# For a single region, include the region in the merged filename.
if [[ ${#TARGETS[@]} -eq 1 && "${TARGETS[0]}" == *:* ]]; then
  MERGED_LABEL="${LABELS[0]}"
else
  MERGED_LABEL="chrAll"
fi


# DAC scope is the same for every WW group.
DAC_SCOPE_ARGS=(--scope genome)

# When only one chromosome or one chromosome region was requested,
# run DAC using chromosome scope.
if [[ ${#TARGETS[@]} -eq 1 ]]; then
  target="${TARGETS[0]}"

  # Remove region coordinates, then convert the BAM contig spelling to
  # the chromosome spelling used in the Python-generated WIG/bigWig output.
  dac_chrom="${target%%:*}"
  if [[ "$dac_chrom" != chr* ]]; then
    dac_chrom="chr${dac_chrom}"
  fi

  DAC_SCOPE_ARGS=(
    --scope chromosome
    --chromosome "$dac_chrom"
  )
fi



# Merge chromosome/region dinucleotide profiles.
#
# Each source TSV contains frequencies plus n_valid for every relative position.
# For each frequency column and position, the merged value is:
#
#   sum(source_frequency * source_n_valid) / sum(source_n_valid)
#
# The merged n_valid is the sum across source files. A temporary output is used
# so this also works safely when a single-region source and destination filename
# are identical.
merge_dinuc_profiles() {
  local output="$1"
  shift

  local tmp="${output}.tmp.$$"
  rm -f "$tmp"

  if ! python3 - "$tmp" "$@" <<'PY'
import csv
import math
import os
import sys
from collections import OrderedDict

if len(sys.argv) < 3:
    raise SystemExit(
        "ERROR: merge_dinuc_profiles requires an output path and at least one input TSV."
    )

output_path = sys.argv[1]
input_paths = sys.argv[2:]

header = None
position_order = []
position_seen = set()

# position -> {"n_valid": int, "weighted": [float, ...]}
merged = OrderedDict()

for input_path in input_paths:
    if not os.path.isfile(input_path):
        raise SystemExit(f"ERROR: dinucleotide profile not found: {input_path}")

    with open(input_path, "r", newline="") as handle:
        reader = csv.reader(handle, delimiter="\t")

        try:
            current_header = next(reader)
        except StopIteration:
            raise SystemExit(f"ERROR: empty dinucleotide profile: {input_path}")

        if len(current_header) < 3:
            raise SystemExit(
                f"ERROR: expected position, n_valid, and frequency columns in: {input_path}"
            )

        if current_header[0:2] != ["position", "n_valid"]:
            raise SystemExit(
                f"ERROR: first two columns must be position and n_valid in: {input_path}"
            )

        if header is None:
            header = current_header
            n_profile_columns = len(header) - 2
        elif current_header != header:
            raise SystemExit(
                "ERROR: dinucleotide-profile headers differ:\n"
                f"  expected: {' | '.join(header)}\n"
                f"  observed: {' | '.join(current_header)}\n"
                f"  file: {input_path}"
            )

        file_positions = []

        for line_number, row in enumerate(reader, start=2):
            if not row:
                continue

            if len(row) != len(header):
                raise SystemExit(
                    f"ERROR: wrong number of columns in {input_path}, line {line_number}"
                )

            try:
                position = int(row[0])
                n_valid = int(row[1])
            except ValueError as exc:
                raise SystemExit(
                    f"ERROR: invalid position or n_valid in {input_path}, "
                    f"line {line_number}: {exc}"
                )

            if n_valid < 0:
                raise SystemExit(
                    f"ERROR: negative n_valid in {input_path}, line {line_number}"
                )

            file_positions.append(position)

            if position not in merged:
                merged[position] = {
                    "n_valid": 0,
                    "weighted": [0.0] * n_profile_columns,
                }

            entry = merged[position]
            entry["n_valid"] += n_valid

            if n_valid == 0:
                continue

            for index, value_text in enumerate(row[2:]):
                try:
                    value = float(value_text)
                except ValueError as exc:
                    raise SystemExit(
                        f"ERROR: invalid frequency in {input_path}, "
                        f"line {line_number}, column {header[index + 2]}: {exc}"
                    )

                if not math.isfinite(value):
                    raise SystemExit(
                        f"ERROR: non-finite frequency with n_valid > 0 in "
                        f"{input_path}, line {line_number}, "
                        f"column {header[index + 2]}"
                    )

                entry["weighted"][index] += value * n_valid

        if not position_order:
            position_order = file_positions
            position_seen = set(file_positions)
        elif file_positions != position_order:
            raise SystemExit(
                "ERROR: dinucleotide-profile position rows differ between files.\n"
                f"  file: {input_path}"
            )

with open(output_path, "w", newline="") as out:
    writer = csv.writer(out, delimiter="\t", lineterminator="\n")
    writer.writerow(header)

    for position in position_order:
        entry = merged[position]
        total_n = entry["n_valid"]
        row = [position, total_n]

        if total_n == 0:
            row.extend(["NaN"] * (len(header) - 2))
        else:
            row.extend(
                f"{weighted_sum / total_n:.8g}"
                for weighted_sum in entry["weighted"]
            )

        writer.writerow(row)
PY
  then
    rm -f "$tmp"
    return 1
  fi

  mv -f "$tmp" "$output"
}

# Return success when OUTPUT exists, is non-empty, and is at least as new as
# every INPUT. This allows completed merged stages to be skipped while ensuring
# they are rebuilt after any chromosome-level input changes.
is_up_to_date() {
  local output="$1"
  shift
  local input

  [[ -s "$output" ]] || return 1

  for input in "$@"; do
    [[ -s "$input" ]] || return 1
    [[ "$output" -nt "$input" || "$output" -ef "$input" ]] || return 1
  done

  return 0
}

for group in "${DYAD_GROUPS[@]}"; do
  if [[ "$group" == "all" ]]; then
    group_suffix=""
    group_label="all"
  else
    group_suffix="_${group}"
    group_label="$group"
  fi

  merged_profile="${SAMPLE}_${MERGED_LABEL}_PNS_mode${MODE}_lower${FRAG_LOWER}_upper${FRAG_UPPER}${group_suffix}_dinuc_profile.tsv"
  profile_files=()

  for label in "${LABELS[@]}"; do
    profile_file="${SAMPLE}_${label}_PNS_mode${MODE}_lower${FRAG_LOWER}_upper${FRAG_UPPER}${group_suffix}_dinuc_profile.tsv"

    if [[ ! -f "$profile_file" ]]; then
      echo "ERROR: expected ${group_label} dinucleotide profile not found:"
      echo "  $profile_file"
      exit 1
    fi

    profile_files+=("$profile_file")
  done

  if [[ "$FORCE_MERGE" != "1" ]] && is_up_to_date "$merged_profile" "${profile_files[@]}"; then
    echo "[SKIP] merged dinucleotide profile for ${group_label}:"
    echo "  ${merged_profile}"
  else
    echo "Merging ${group_label} dinucleotide profiles with n_valid weighting:"
    echo "  ${merged_profile}"

    merge_dinuc_profiles "$merged_profile" "${profile_files[@]}"

    echo "[DONE] merged dinucleotide profile for ${group_label}"
  fi
  echo

  merged_wig="${SAMPLE}_${MERGED_LABEL}_PNS_mode${MODE}_lower${FRAG_LOWER}_upper${FRAG_UPPER}${group_suffix}_dyad.wig"
  merged_bw="${SAMPLE}_${MERGED_LABEL}_PNS_mode${MODE}_lower${FRAG_LOWER}_upper${FRAG_UPPER}${group_suffix}_dyad.bw"

  files=()

  for label in "${LABELS[@]}"; do
    f="${SAMPLE}_${label}_PNS_mode${MODE}_lower${FRAG_LOWER}_upper${FRAG_UPPER}${group_suffix}_dyad.wig.gz"

    if [[ ! -f "$f" ]]; then
      echo "ERROR: expected ${group_label} dyad output file not found:"
      echo "  $f"
      exit 1
    fi

    files+=("$f")
  done

  if [[ "$FORCE_MERGE" != "1" ]] && is_up_to_date "$merged_bw" "${files[@]}" "$CHROM_SIZES"; then
    echo "[SKIP] merged bigWig for ${group_label}:"
    echo "  ${merged_bw}"
  else
    echo "Merging ${group_label} dyad WIG files:"
    echo "  ${merged_wig}"

    wig_tmp="${merged_wig}.tmp.$$"
    bw_tmp="${merged_bw}.tmp.$$"
    rm -f "$wig_tmp" "$bw_tmp"

    zcat "${files[@]}" > "$wig_tmp"
    mv -f "$wig_tmp" "$merged_wig"

    echo "Converting ${group_label} dyad WIG to bigWig:"
    echo "  ${merged_bw}"

    wigToBigWig \
      "$merged_wig" \
      "$CHROM_SIZES" \
      "$bw_tmp"

    mv -f "$bw_tmp" "$merged_bw"
    echo "[DONE] merged bigWig for ${group_label}"
  fi

  dac_marker="${merged_bw%.bw}_DAC.done"
  dac_signature="$({
    printf "%s\0" \
      "$merged_bw" \
      "$CHROM_SIZES" \
      "${DAC_SCOPE_ARGS[@]}" \
      "$DMAX" \
      "$VALUE_LIMIT" \
      "$MIN_REGION_LENGTH" \
      "--no_normalize_dac"
    stat -c "%n\0%s\0%Y\0" "$merged_bw" "$CHROM_SIZES" "$DAC_SCRIPT"
  } | sha256sum | awk '{print $1}')"

  if [[ "$FORCE_DAC" != "1" && -s "$dac_marker" && "$(head -n 1 "$dac_marker")" == "$dac_signature" ]]; then
    echo "[SKIP] DAC for ${group_label} dyads (matching completed run)"
  else
    echo "Running DAC for ${group_label} dyads:"
    echo "  ${merged_bw}"
    rm -f "$dac_marker"

    python3 "$DAC_SCRIPT" \
      --bigwig "$merged_bw" \
      --chrom_sizes "$CHROM_SIZES" \
      "${DAC_SCOPE_ARGS[@]}" \
      --dmax "$DMAX" \
      --value_limit "$VALUE_LIMIT" \
      --min_region_length "$MIN_REGION_LENGTH" \
      --no_normalize_dac

    printf "%s\n" "$dac_signature" > "${dac_marker}.tmp.$$"
    mv -f "${dac_marker}.tmp.$$" "$dac_marker"
    echo "[DONE] DAC for ${group_label} dyads"
  fi
  echo
done

echo "Finished all WW-type dinucleotide-profile merges, dyad bigWig conversions, and DAC analyses (completed stages were resumed/skipped)."