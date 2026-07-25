# PNS with Nucleosome and Breakpoint Peak Calling

**PNS (Probabilistic Nucleosome Scoring)** is a fragmentomics pipeline
for generating high-resolution nucleosome protection and breakpoint maps
from paired-end sequencing data.

The pipeline reads one or more coordinate-sorted BAM files, filters
paired-end fragments, computes probabilistic nucleosome scores (PNS),
coverage, dyad and fragment-end tracks, and identifies
nucleosome protection peaks and breakpoint peaks. Additional modules
allow fragment randomisation, sequence-based WW/SS classification, and
aligned dinucleotide profiling.

Outputs are written directly as **BigWig** files by default, with
optional **WIG** output and BED-formatted peak calls.

------------------------------------------------------------------------

# Features

## Core scoring

-   Calculates **Probabilistic Nucleosome Score (PNS)** tracks from
    paired-end fragments
-   Generates raw and Savitzky--Golay smoothed PNS tracks
-   Produces:
    -   coverage
    -   dyad density
    -   fragment-end density
    -   left-end density
    -   right-end density

## Peak calling

Calls two complementary classes of peaks from the smoothed PNS signal:

-   **Nucleosome protection peaks** (positive PNS regions)
-   **Breakpoint peaks** (negative PNS regions)

Both are written as BED files suitable for downstream analysis or genome
browsers.

## Fragment processing

Supports:

-   multiple BAM inputs
-   coordinate-based duplicate filtering
-   deduplication across all BAMs or within individual BAMs
-   fragment length filtering
-   optional random subsampling

## Fragment randomisation

  -----------------------------------------------------------------------
  Mode                     Description
  ------------------------ ----------------------------------------------
  none                     No randomisation

  uniform                  Randomises fragment positions while preserving
                           fragment length

  dinuc_anchor             Preserves fragment boundary dinucleotides by
                           relocating fragments to matching reference
                           sequence positions
  -----------------------------------------------------------------------

The dinucleotide-anchor mode requires a reference FASTA.

## WW/SS sequence classification

Fragments can optionally be classified into the four canonical WW/SS
nucleosome sequence classes using a centred 147 bp reference sequence.

When enabled the pipeline automatically produces independent outputs
for:

-   Type 1
-   Type 2
-   Type 3
-   Type 4

while also retaining the standard all-fragment outputs.

## Dinucleotide profiling

Optionally calculates observed dinucleotide frequencies aligned relative
to fragment dyads.

Outputs include:

-   all 16 dinucleotides
-   combined WW frequencies
-   combined SS frequencies

Values may be written as percentages or fractions.

## Flexible genome selection

Examples include:

``` text
chr1
chr1,chr2,chr3
chr1-22
autosomes
all
chr12:100000-200000
```

Large chromosomes are processed automatically in overlapping windows to
minimise memory usage.

------------------------------------------------------------------------

# Requirements

Python ≥3.9

Required packages:

``` text
numpy
scipy
pysam
pyBigWig
```

------------------------------------------------------------------------

# Quick Start

## Whole genome

``` bash
python PNS_with_nucleosome_peak_calling.py \
    -b sample.bam
```

## Selected chromosomes

``` bash
python PNS_with_nucleosome_peak_calling.py \
    -b sample.bam \
    -c chr1-22,chrX
```

## Pool multiple BAMs

``` bash
python PNS_with_nucleosome_peak_calling.py \
    -b sample1.bam sample2.bam sample3.bam
```

## Generate WW type outputs

``` bash
python PNS_with_nucleosome_peak_calling.py \
    -b sample.bam \
    --split-ww-types \
    --fasta hg38.fa
```

## Generate aligned dinucleotide profiles

``` bash
python PNS_with_nucleosome_peak_calling.py \
    -b sample.bam \
    --dinuc-profile \
    --fasta hg38.fa
```

------------------------------------------------------------------------

# Main Options

## Input

``` text
-b / --bamfiles
-c / --contigs
-o / --out_prefix
```

## Fragment filtering

``` text
--frag-lower
--frag-upper
--max-duplicates
--dedup-scope
--subsample
```

## PNS

``` text
--mode-length
--pns-mode
```

## Randomisation

``` text
--randomize-mode
--fasta
--anchor-prob-start
--randomize-fallback
```

## Sequence analysis

``` text
--split-ww-types
--dinuc-profile
--dinuc-fraction
```

## Output

``` text
--pns-format
--other-format
--pns-tracks
--other-tracks
```

## Peak calling

``` text
--peak-format
--min-region-length
--max-neg-run
```

------------------------------------------------------------------------

# Output Files

## Score tracks

Direct **BigWig** output (default):

``` text
*_pns_smoothed.bw
*_pns.bw
*_posPNS.bw
*_coverage.bw
*_dyad.bw
*_fragment_ends.bw
*_fragment_left_ends.bw
*_fragment_right_ends.bw
```

Optional compressed WIG output:

``` text
*.wig.gz
```

## Peak calls

``` text
*_nucleosome_regions.bed
*_breakpoint_peaks.bed
```

## Sequence analysis

``` text
*_dinuc_profile.tsv
*_ww_type_summary.tsv
```

## Fragment summaries

``` text
*_fragment_summary.txt
*_fragment_length_counts.tsv
```

When `--split-ww-types` is enabled, each WW/SS class receives its own
complete set of outputs in addition to the combined all-fragment
results.

------------------------------------------------------------------------

# Pipeline Overview

1.  Read paired-end fragments from one or more BAM files.
2.  Apply filtering, deduplication and optional subsampling.
3.  Optionally randomise fragment positions.
4.  Optionally classify fragments into WW/SS Types 1--4.
5.  Compute PNS, coverage, dyad and fragment-end tracks.
6.  Optionally generate aligned dinucleotide profiles.
7.  Smooth the PNS signal using a Savitzky--Golay filter.
8.  Identify nucleosome protection peaks and breakpoint peaks.
9.  Write all requested tracks and summary files.
