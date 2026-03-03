# PNS with Nucleosome + Breakpoint Peak Calling (Fragmentomics scoring + peak calling)

This repository/script computes **PNS (Probabilistic Nucleosome Score)** tracks from paired-end BAM files, plus coverage and dyad tracks, optionally randomizes fragments, smooths the PNS signal, and calls:

- **Positive peaks** → **nucleosome regions**
- **Negative peaks (troughs)** → **breakpoint peaks** (called by flipping the sign and re-using the same peak caller)

---

## Table of contents

- [Overview](#overview)
- [Requirements](#requirements)
- [Concept: PNS scoring](#concept-pns-scoring)
- [Fragment extraction and filtering](#fragment-extraction-and-filtering)
- [Fragment randomization modes](#fragment-randomization-modes)
- [Smoothing](#smoothing)
- [Peak calling](#peak-calling)
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

1. Reads **paired-end fragments** from one or more BAMs across contigs or user-specified regions.
2. Filters and pairs reads (generator logic in-script), collapses duplicate fragments by coordinate tuple to `--max-duplicates`, and optionally subsamples (`--subsample`).
3. Optionally randomizes fragment positions within each processed window (see [Fragment randomization modes](#fragment-randomization-modes)).
4. For each fragment length in `[--frag-lower, --frag-upper]`:
   - Adds a **precomputed, length-specific PNS kernel** across the fragment footprint.
   - Computes a **coverage** track (+1 for each covered base).
   - Computes a **dyad** track (+1 at fragment center).
5. Smooths the raw PNS track using **Savitzky–Golay** (**window=21, polyorder=2**).
6. Calls peaks/regions on the smoothed track:
   - **positive peaks** → nucleosome regions
   - **negative troughs** (by sign flip) → breakpoint peaks
7. Writes (configurable):
   - **bedGraph** (combined multi-track bedGraph)
   - **wig.gz** (one file per track; `fixedStep`)
   - peak calls in either **BED8** or **rich**
   - fragment summary outputs

---

## Requirements

### Python
- Python 3

### Python libraries
- `numpy`
- `scipy`
- `pysam`
- `tqdm`

---

## Concept: PNS scoring

The script implements a **PNS scoring** approach:

- For each fragment length in `[--frag-lower, --frag-upper]`, a **kernel** is precomputed.
- The kernel is formed by:
  - building a **triangle distribution** across the first `--mode-length` bases,
  - mirroring it from the fragment end,
  - summing start + end triangles,
  - and **mean-centering** so the kernel has ~0 mean (prevent baseline drift).
- For fragments shorter than `--mode-length`, the kernel is expanded symmetrically so the “mode window” still fits.

This produces:
- positive signal across nucleosome-protected regions
- negative flanking contributions associated with fragment ends

---

## Fragment extraction and filtering

Paired-end fragment extraction uses `pysam.fetch(..., multiple_iterators=True)` and pairs reads by `query_name`.

Filters applied (both mates defensively):

- skip unmapped or mate-unmapped
- skip `is_duplicate` or `is_qcfail`
- skip reads with **S/H/P** CIGAR ops (soft clip / hard clip / pad)
- skip same-strand pairs
- fragment coordinates are `(contig, min(start), max(end))`

Duplicate handling:

- duplicate fragments are tracked by `(contig, frag_start, frag_end)`
- a fragment coordinate is allowed while `frag_counts[key] <= --max-duplicates`
  - `--max-duplicates 0` keeps **1** instance
  - `--max-duplicates 1` keeps **up to 2**, etc.

Optional subsampling:

- if `--subsample p` is set, each fragment is kept with probability `p`

---

## Fragment randomization modes

Randomization is applied **after filtering/subsampling/dup logic**, and **before scoring**, within each processed window.

`--randomize-mode` options:

- `none` (default): no randomization
- `uniform`: uniformly randomize fragment start positions inside the window, preserving fragment lengths
- `dinuc_anchor` *(requires `--fasta`)*:
  - for each fragment, choose to anchor on its **start** or **end** dinucleotide (probability set by `--anchor-prob-start`)
  - pick a random occurrence of that dinucleotide in the reference sequence for the window
  - place the fragment so its anchored boundary matches that occurrence
  - requires a reference FASTA with `.fai`

Additional controls:

- `--fasta <ref.fa>` (required for `dinuc_anchor`)
- `--anchor-prob-start` (default `0.5`)
- `--max-anchor-tries` (default `30`)
- `--randomize-fallback {uniform,keep,skip}` (default `uniform`)

---

## Smoothing

The script smooths the raw PNS track with Savitzky–Golay:

- window size: `21`
- polynomial order: `2`

---

## Peak calling

Peak calling is performed using a simple positive-region strategy on the pns_smoothed track.

Process:

1. Build runs where PNS value > 0.
2. Permit up to 5 consecutive bases with value ≤ 0 before terminating the run.
3. Retain runs that are at least 50 bp long.

### What gets called
- **Nucleosome regions**: positive runs in pns_smoothed
- **Breakpoint peaks**: positive runs in -1 * pns_smoothed

---

## Script usage

### Basic examples

#### Whole genome (all contigs in BAM header)
```bash
python3 PNS_with_nucleosome_peak_calling.py   -b sample1.bam
```

#### Restrict to one contig
```bash
python3 PNS_with_nucleosome_peak_calling.py   -b sample1.bam   -c 2
```

#### Restrict to a genomic interval
Use `contig:start-end` (**0-based start**, **end-exclusive**).
```bash
python3 PNS_with_nucleosome_peak_calling.py   -b sample1.bam   -c 2:100000-200000
```

#### Multiple BAMs (pooled fragments)
```bash
python3 PNS_with_nucleosome_peak_calling.py   -b sample1.bam sample2.bam   -c 1
```

#### Duplicate filtering + subsampling
```bash
python3 PNS_with_nucleosome_peak_calling.py   -b sample1.bam   --max-duplicates 1   --subsample 0.25
```

#### Write `.wig.gz` tracks (default in this version)
```bash
python3 PNS_with_nucleosome_peak_calling.py   -b sample1.bam   -c 12:52621135-52641135   --score-format wiggz
```

#### Write bedGraph + wig.gz
```bash
python3 PNS_with_nucleosome_peak_calling.py   -b sample1.bam   -c 12:52621135-52641135   --score-format both
```

#### Disable per-base score tracks (peaks + summaries only)
```bash
python3 PNS_with_nucleosome_peak_calling.py   -b sample1.bam   -c 12:52621135-52641135   --score-format none
```

#### Select which tracks to write
```bash
python3 PNS_with_nucleosome_peak_calling.py   -b sample1.bam   -c 12:52621135-52641135   --score-format wiggz   --score-tracks pns_smoothed coverage
```

#### Randomize fragments (uniform)
```bash
python3 PNS_with_nucleosome_peak_calling.py   -b sample1.bam   -c 12:52621135-52641135   --randomize-mode uniform   --seed 1
```

#### Randomize fragments (dinucleotide anchored)
```bash
python3 PNS_with_nucleosome_peak_calling.py   -b sample1.bam   -c 12:52621135-52641135   --randomize-mode dinuc_anchor   --fasta /path/to/ref.fa   --seed 1
```

---

## Arguments

### Core inputs
- `-b/--bamfiles` **(required)**: one or more paired-end BAMs
- `-o/--out_prefix`: output prefix (default derived from BAM basenames + region info)
- `-c/--contigs`: contig(s) and optional interval(s), e.g. `2` or `2:100000-200000`

### Kernel / fragment range
- `--mode-length` (default `167`): mode fragment length used to define kernel geometry
- `--frag-lower` (default `127`): lower fragment length (inclusive)
- `--frag-upper` (default `207`): upper fragment length (inclusive)

### Filtering / sampling / chunking
- `--max-duplicates` (default `0`): max allowed duplicate fragments with same `(contig,start,end)`
- `--subsample` (default `None`): keep each fragment with probability `p`
- `--chunk-bp` (default `100000`): chunk size for windowing
- `--overlap-bp` (default `1000`): overlap padding on each side of chunk
- `--seed` (default `None`): random seed

### Randomization
- `--randomize-mode {none,uniform,dinuc_anchor}` (default `none`)
- `--fasta <ref.fa>`: required for `dinuc_anchor` (must have `.fai`)
- `--anchor-prob-start` (default `0.5`)
- `--max-anchor-tries` (default `30`)
- `--randomize-fallback {uniform,keep,skip}` (default `uniform`)

### Score-track output controls
- `--score-format {bedgraph,wiggz,both,none}` (default `wiggz`)
  - `bedgraph`: write combined multi-track bedGraph
  - `wiggz`: write one `<prefix>_<track>.wig.gz` per track (`fixedStep`)
  - `both`: write both
  - `none`: write no per-base score tracks
- `--score-tracks` (default: `coverage pns_smoothed pns dyad`)
  - Valid: `coverage pns_smoothed pns dyad`
  - Use `--score-tracks none` to disable track writing (empty list)

### Peak output controls
- `--peak-format {rich,bed8}` (default `bed8`)
- `--peak-score-scale` (default `1.0`)
  - only used when `--peak-format bed8`
  - BED8 score is `round(prominence * scale)` as an integer

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
<out_prefix>_mode<MODE>_lower<LOWER>_upper<UPPER>
```

### Score tracks

#### 1) Combined per-base scores bedGraph (legacy)
If `--score-format bedgraph` or `both`:

- `<prefix>_combined_scores.bedGraph`

Rows are 1 bp wide (end = start + 1). Columns:

1. chrom
2. start
3. end
4. coverage *(int)*
5. pns_smoothed *(float)*
6. pns *(float)*
7. dyad *(int)*

#### 2) `.wig.gz` per-track outputs (NEW)
If `--score-format wiggz` or `both`, one gzipped **fixedStep** WIG is written per selected track:

- `<prefix>_coverage.wig.gz`
- `<prefix>_pns_smoothed.wig.gz`
- `<prefix>_pns.wig.gz`
- `<prefix>_dyad.wig.gz`

Notes:
- WIG `fixedStep start=` is **1-based**.
- Only the core (non-overlap) region is written per chunk.

### Peak calls

Two peak call files are written:

- `<prefix>_nucleosome_regions.bed`
- `<prefix>_breakpoint_peaks.bed`

The contents depend on `--peak-format`:

#### `--peak-format bed8` (default)
A BED8-like format:

1. chrom  
2. start *(region start, 0-based)*  
3. end *(region end, end-exclusive)*  
4. name *(e.g. `chr12:52630000_nuc` or `..._brk`)*  
5. score *(int; prominence × scale)*  
6. strand *(always `.`)*  
7. thickStart *(peak position)*  
8. thickEnd *(peak position + 1)*  

### Fragment summaries

Two files are always written:

- `<prefix>_fragment_summary.txt`
  - `total_fragments_filtered_all`
  - `total_fragments_used_in_range`
  - `unique_bases_covered_by_used_fragments`
- `<prefix>_fragment_length_counts.tsv`
  - columns: `fragment_length`, `count`
  - counts only fragments within `[--frag-lower, --frag-upper]`

---

## Notes

- Input BAMs must be indexed: each BAM requires a `.bai` index in the same directory as the BAM.
- `-c contig:start-end` uses **0-based start** and **end-exclusive** semantics.
- Coverage/dyad counting is only applied to fragments:
  - fully contained within the current scoring window, **and**
  - within the allowed fragment-length range.
- `.wig.gz` outputs are **fixedStep** and therefore assume 1 bp step with contiguous positions per chunk.
