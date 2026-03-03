# WPS with Nucleosome + Breakpoint Peak Calling (Kircher-style)

This script generates **Kircher-style WPS (Windowed Protection Score)** tracks from paired-end BAM files, applies **rolling-median baseline correction**, optional **Savitzky–Golay smoothing**, and calls:

- **Positive peaks** → **nucleosome regions**
- **Negative peaks (troughs)** → **breakpoint peaks** (called by flipping the sign and re-using the same peak caller)

It is designed to match **Martin Kircher’s 2015 WPS behavior** as closely as possible while operating on standard BAM inputs (0-based, half-open coordinates via `pysam`).

**Important implementation detail:** this WPS script **imports shared helper functions** from `PNS_with_nucleosome_peak_calling.py` (expected somewhere under the WPS script directory or its parent directory).

---

## Table of contents

- [Overview](#overview)
- [Requirements](#requirements)
- [How this implementation matches Kircher WPS](#how-this-implementation-matches-kircher-wps)
  - [Kircher-style overlap semantics](#kircher-style-overlap-semantics)
  - [Fragment size filtering](#fragment-size-filtering)
- [Tracks produced](#tracks-produced)
- [Baseline correction and smoothing](#baseline-correction-and-smoothing)
- [Peak calling (Kircher-matched)](#peak-calling-kircher-matched)
- [Fragment randomization modes](#fragment-randomization-modes)
- [Script usage](#script-usage)
  - [Basic examples](#basic-examples)
  - [Arguments](#arguments)
  - [How contigs / regions are processed](#how-contigs--regions-are-processed)
- [Outputs](#outputs)
  - [Score tracks](#score-tracks)
  - [Peak calls](#peak-calls)
  - [Fragment summaries](#fragment-summaries)
- [Converting to bigWig / bigBed](#converting-to-bigwig--bigbed)
- [Conda / Mamba environment setup](#conda--mamba-environment-setup)
- [Notes / gotchas](#notes--gotchas)

---

## Overview

High-level pipeline:

1. Reads **paired-end fragments** from one or more BAM files across contigs or user-specified regions.
2. Fragment extraction / filtering is performed by shared functions imported from **PNS** (`generate_paired_reads`, `generate_fragment_ranges`, duplicate handling, etc.).
3. Optionally randomizes fragment positions (see [Fragment randomization modes](#fragment-randomization-modes)).
4. For each fragment:
   - Updates **coverage** track (+1 for each covered base).
   - Updates **dyad** track (+1 at fragment midpoint).
   - If the fragment length is within `[--frag-lower, --frag-upper]`, adds a **Kircher-matched WPS kernel** to the WPS track.
5. Computes a **rolling-median baseline** (default 1000 bp) on the **raw WPS** track and subtracts it.
6. Optionally smooths **raw WPS** via Savitzky–Golay (default window/order 21/2).
7. Builds:
   - `mWPS = wps - rolling_median(wps)`
   - `sm_mWPS = wps_smoothed - rolling_median(wps)`
8. Calls peaks on **`sm_mWPS`**:
   - Positive runs → nucleosome calls
   - Negative runs (by sign flip) → breakpoint calls
9. Writes score tracks (bedGraph and/or `.wig.gz`) plus peak calls and fragment summaries.

---

## Requirements

### Python
- Python 3 (modern 3.x)

### Python libraries
- `numpy`
- `scipy`
- `pysam`
- `tqdm`

### External tools (optional but recommended)
- `samtools` (indexing / inspection)
- UCSC tools (for visualization-friendly formats):
  - `bedGraphToBigWig`
  - `bedToBigBed`

---

## How this implementation matches Kircher WPS

### Kircher-style overlap semantics

Kircher’s original implementation uses bx-python interval logic and 1-based coordinate handling that leads to subtle off-by-one effects in how “spanning the protection window” is evaluated.

This script reproduces that behavior by adding a **precomputed per-fragment kernel** whose effective length and placement match Kircher’s bx semantics:

- For protection window `--protection 120`, half-window is `60`.
- The WPS kernel is placed at `kernel_start = frag_start - half + 1` (genome coordinates).
- Each allowed fragment length has a specific kernel shape that yields Kircher-equivalent “+1 in the middle / -1 in the flanks” behavior.

### Fragment size filtering

Only fragments in the configured size range contribute to WPS:

- `--frag-lower` (default `120`)
- `--frag-upper` (default `180`)

---

## Tracks produced

The script computes the following tracks over each processed region:

- `coverage` — per-base fragment overlap count
- `dyad` — per-base fragment midpoint count
- `wps` — raw Kircher-style WPS kernel sum
- `wps_smoothed` — Savitzky–Golay smoothing of `wps` (if possible given window length)
- `mWPS` — `wps - rolling_median(wps)`
- `sm_mWPS` — `wps_smoothed - rolling_median(wps)`

**Peak calling uses `sm_mWPS`.**

---

## Baseline correction and smoothing

### Rolling-median baseline
A rolling median baseline is computed on the **raw** WPS track:

- `--baseline-window` (default `1000` bp)

The baseline is computed via a sliding median and extended at region edges by holding the first/last valid baseline value.

### Savitzky–Golay smoothing
Smoothing is applied to raw WPS (not to the baseline):

- `--sg-window` (default `21`, must be odd and ≤ region length)
- `--sg-order` (default `2`)

Resulting smoothed track:
- `wps_smoothed`

The baseline-subtracted smoothed track:
- `sm_mWPS = wps_smoothed - baseline`

---

## Peak calling (Kircher-matched)

Peak calling is designed to mirror Kircher’s `evaluateValues()` logic as closely as possible.

Process (on a given track):

1. Build runs where `value > 0`.
2. Merge runs if separated by gaps ≤ `--peak-merge-gap` (default `5` bp), filling gaps with zeros.
3. Evaluate merged runs:
   - If run length is in `[--peak-minlen, --peak-maxlen]` (default 50–150 bp), compute median-based windows and report best.
   - If run length is in `[--peak-maxlen, 3*--peak-maxlen]` (default 150–450 bp), split into median-threshold windows and report each segment that meets length constraints.
   - Longer regions are rejected (default behavior is >450 bp, but you can override with `--peak-maxregion`).
4. A window is reported only if its **maximum value** exceeds:
   - `--peak-varicutoff` (default `5.0`)

### What gets called
- **Nucleosome regions**: peak calling on `sm_mWPS`
- **Breakpoint peaks**: peak calling on `-1 * sm_mWPS`

---

## Fragment randomization modes

Randomization modes match those in the imported PNS module.

- `--randomize-mode none` (default)  
  Use observed fragment positions.

- `--randomize-mode uniform`  
  Uniformly randomize fragment placements within each processed window.

- `--randomize-mode dinuc_anchor` *(requires `--fasta`)*  
  Randomize fragment placements while anchoring fragment start or end on matching dinucleotide contexts.

Additional controls for dinucleotide-anchored randomization:
- `--fasta <ref.fa>` *(must have `.fai`)*
- `--anchor-prob-start` (default `0.5`) — probability of anchoring on fragment start vs end
- `--max-anchor-tries` (default `30`) — attempts per fragment
- `--randomize-fallback {uniform,keep,skip}` (default `uniform`) — behavior if no valid placement is found

---

## Script usage

### Basic examples

#### Whole genome (all contigs in BAM header)
```bash
python3 wps_with_nucleosome_peak_calling.py   -b sample.bam
```

#### Restrict to a contig
```bash
python3 wps_with_nucleosome_peak_calling.py   -b sample.bam   -c chr12
```

#### Restrict to an interval
Intervals are interpreted as **0-based, end-exclusive**: `contig:start-end`
```bash
python3 wps_with_nucleosome_peak_calling.py   -b sample.bam   -c chr12:52621135-52641135
```

#### Write both bedGraph and `.wig.gz`, and only selected tracks
```bash
python3 wps_with_nucleosome_peak_calling.py   -b sample.bam   -c chr12:52621135-52641135   --score-format both   --score-tracks sm_mWPS wps coverage
```

#### Disable all score tracks (peaks + summaries only)
```bash
python3 wps_with_nucleosome_peak_calling.py   -b sample.bam   -c chr12:52621135-52641135   --score-format none
```

#### Randomize fragments uniformly
```bash
python3 wps_with_nucleosome_peak_calling.py   -b sample.bam   -c chr12:52621135-52641135   --randomize-mode uniform   --seed 1
```

#### Dinucleotide-anchored randomization
```bash
python3 wps_with_nucleosome_peak_calling.py   -b sample.bam   -c chr12:52621135-52641135   --randomize-mode dinuc_anchor   --fasta /path/to/ref.fa   --seed 1
```

---

## Arguments

### Core inputs
- `-b/--bamfiles` **(required)**: one or more paired-end BAMs
- `-o/--out_prefix`: output prefix (default derived from BAM basename(s) and region)
- `-c/--contigs`: contig(s) and optional interval(s), e.g. `chr12` or `chr12:51730340-52039340`

### WPS scoring
- `--protection` (default `120`): protection window size (bp)
- `--frag-lower` (default `127`): minimum fragment length contributing to WPS
- `--frag-upper` (default `207`): maximum fragment length contributing to WPS
- `--max-duplicates` (default `0`): maximum allowed duplicate fragments with identical `(start,end)` within the window
- `--subsample` (default `None`): subsample proportion (e.g. `0.5` keeps ~50%)

### Chunking / overlap
- `--chunk-bp` (default `100000`): chunk size per contig
- `--overlap-bp` (default `1000`): overlap padding for edge-safe scoring

### Baseline / smoothing
- `--baseline-window` (default `1000`): rolling median window (bp)
- `--sg-window` (default `21`): Savitzky–Golay window (odd)
- `--sg-order` (default `2`): Savitzky–Golay polynomial order

### Score-track output controls
- `--score-format {bedgraph,wiggz,both,none}` (default `wiggz`)
  - `bedgraph`: write `<prefix>_combined_scores.bedGraph`
  - `wiggz`: write one `<prefix>_<track>.wig.gz` per selected track
  - `both`: write both bedGraph and `.wig.gz`
  - `none`: do not write any per-base score tracks
- `--score-tracks` (default: `coverage sm_mWPS wps wps_smoothed mWPS dyad`)
  - Valid: `coverage dyad wps wps_smoothed mWPS sm_mWPS`
  - Use `--score-tracks none` to disable track writing (equivalent to empty list; still allows peaks/summaries).

### Peak calling
- `--peak-minlen` (default `50`): minimum window length to report
- `--peak-maxlen` (default `150`): maximum window length to report
- `--peak-maxregion` (default `450`): reject merged positive regions longer than this (bp)
- `--peak-merge-gap` (default `5`): merge positive runs across gaps ≤ this many bp
- `--peak-varicutoff` (default `5.0`): minimum max score within a reported window

### Randomization
- `--seed` (default `None`): random seed (reproducibility)
- `--randomize-mode {none,uniform,dinuc_anchor}` (default `none`)
- `--fasta <ref.fa>`: required for `dinuc_anchor` (must have `.fai`)
- `--anchor-prob-start` (default `0.5`)
- `--max-anchor-tries` (default `30`)
- `--randomize-fallback {uniform,keep,skip}` (default `uniform`)

---

## How contigs / regions are processed

To scale to genome-wide BAMs efficiently, contigs are processed in windows:

- window length: `--chunk-bp` (default 100,000 bp)
- overlap padding: `--overlap-bp` (default 1,000 bp) on both sides

Scoring is computed on the **adjusted window** (including overlap). Output is trimmed back to the **core interval** (`original_start`–`original_end`) to avoid duplication across adjacent chunks.

---

## Outputs

All outputs are prefixed with:

```
<out_prefix>_prot<PROTECTION>_lower<LOWER>_upper<UPPER>_maxdup<MAXDUP>
```

If randomization is used, an additional suffix is appended:

```
..._rand<MODE>
```

### Score tracks

#### bedGraph (combined)
If `--score-format bedgraph` or `both`, this is written:

- `<prefix>_combined_scores.bedGraph`

The exact column layout is produced by the shared PNS `write_bedgraph()` helper; it includes the tracks computed by this script (coverage/dyad/wps/wps_smoothed/mWPS/sm_mWPS) in a consistent per-base format.

#### `.wig.gz` (one file per track)
If `--score-format wiggz` or `both`, and `--score-tracks` is non-empty, one gzipped fixedStep WIG is written per requested track:

- `<prefix>_coverage.wig.gz`
- `<prefix>_dyad.wig.gz`
- `<prefix>_wps.wig.gz`
- `<prefix>_wps_smoothed.wig.gz`
- `<prefix>_mWPS.wig.gz`
- `<prefix>_sm_mWPS.wig.gz`

(Only the tracks listed in `--score-tracks` are written.)

### Peak calls

Two BED8 files are written (Kircher-style fields):

- `<prefix>_nucleosome_regions.bed`
- `<prefix>_breakpoint_peaks.bed`

Columns:
1. chrom (e.g. `chr12`)
2. start (0-based)
3. end (end-exclusive)
4. name (`<chrom_nochr>:<start>-<end>` in 1-based display style)
5. score (integer; max score in reported window)
6. strand (always `.`)
7. thickStart (midpoint-1)
8. thickEnd (midpoint)

Example:
```
chr12   51748712   51748745   12:51748713-51748745   6   .   51748728   51748729
```

### Fragment summaries

Two files are always written:

- `<prefix>_fragment_summary.txt`
  - `total_fragments_filtered_all`
  - `total_fragments_used_in_range` *(only those within frag range used for WPS)*
  - `unique_bases_covered_by_used_fragments`
- `<prefix>_fragment_length_counts.tsv`
  - two columns: `fragment_length`, `count`

---

## Conda / Mamba environment setup

Minimal environment:

```bash
mamba create -n wps_env -c conda-forge python=3.11 numpy scipy pysam tqdm
conda activate wps_env
```

---

## Notes

- Input BAMs must be indexed (`.bai` in the same directory).
- `-c contig:start-end` uses **0-based start** and **end-exclusive** semantics.
