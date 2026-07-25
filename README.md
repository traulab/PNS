# PNS with Nucleosome and Breakpoint Peak Calling

**PNS (Probabilistic Nucleosome Scoring)** is a fragmentomics pipeline for generating high-resolution nucleosome protection and breakpoint maps from paired-end sequencing data.

The pipeline reads one or more coordinate-sorted BAM files, filters paired-end fragments, computes probabilistic nucleosome scores (PNS), identifies nucleosome protection and breakpoint peaks, and outputs coverage, dyad and fragment-end tracks. Additional modules support fragment randomisation, WW/SS sequence classification, and aligned dinucleotide profiling.

Outputs are written directly as **BigWig** files by default, with optional compressed WIG output and BED-formatted peak calls.

---

## Features

### Core scoring

The pipeline can generate:

- raw PNS
- Savitzky–Golay-smoothed PNS
- positive-only PNS
- fragment coverage
- dyad density
- combined fragment-end density
- left fragment-end density
- right fragment-end density

### Peak calling

Two complementary classes of peaks are called from the smoothed PNS signal:

- **Nucleosome protection peaks** from positive PNS regions
- **Breakpoint peaks** from negative PNS regions

Both peak files are written in **BED8** format with the following columns:

| Column | Name | Description |
|:------:|------|-------------|
| 1 | **chrom** | Chromosome or contig name. |
| 2 | **start** | Start coordinate of the called region (0-based, inclusive). |
| 3 | **end** | End coordinate of the called region (0-based, end-exclusive). |
| 4 | **name** | Peak identifier, indicating the genomic position and peak type (nucleosome or breakpoint). |
| 5 | **score** | Maximum smoothed PNS signal within the called region. |
| 6 | **strand** | Strand field (`.`); peaks are not strand-specific. |
| 7 | **thickStart** | Genomic coordinate of the peak centre. |
| 8 | **thickEnd** | Peak centre coordinate + 1. |

### Fragment processing

The pipeline supports:

- one or more paired-end BAM files
- pooling all BAMs in a directory with a shell wildcard such as `*.bam`
- fragment-length filtering
- coordinate-based duplicate filtering
- deduplication within each BAM or across all BAMs
- optional random subsampling
- chromosome, contig-range, interval, autosome, or whole-genome analysis

### Fragment randomisation

Two randomisation modes are available:

- **`uniform`**  
  Randomise fragment positions within each processed window while preserving fragment length.

- **`dinuc_anchor`**  
  Relocate fragments to reference positions with matching fragment-boundary dinucleotides. This mode preserves the selected start- or end-boundary dinucleotide and requires an indexed reference FASTA.

### WW/SS sequence classification

Fragments can optionally be classified into four WW/SS nucleosome sequence classes using a centred 147 bp sequence around the fragment dyad:

- Type 1
- Type 2
- Type 3
- Type 4

When enabled, the pipeline produces separate outputs for each class while retaining the combined all-fragment outputs.

### Dinucleotide profiling

The pipeline can calculate aligned dinucleotide frequencies relative to fragment dyads.

Outputs include:

- all 16 dinucleotides
- combined WW frequencies
- combined SS frequencies

Values can be written as percentages or fractions.

### Flexible contig selection

Examples of valid contig selections include:

```text
chr1
chr1,chr2,chr3
chr1-22
chr1-22,chrX,chrY
autosomes
all
chr12:100000-200000
```

Large contigs are processed automatically in overlapping windows to reduce memory usage.

---

## Requirements

- Python 3
- `numpy`
- `scipy`
- `pysam`
- `pyBigWig`

Install the Python dependencies with:

```bash
python -m pip install numpy scipy pysam pyBigWig
```

Input BAM files must be coordinate sorted and indexed.

---

## Quick start

### Single BAM

```bash
python PNS_with_nucleosome_peak_calling.py \
    -b sample.bam
```

### All BAMs in the current directory

The shell expands `*.bam` into all matching BAM filenames before running the script:

```bash
python PNS_with_nucleosome_peak_calling.py \
    -b *.bam
```

All BAMs in another directory can be supplied in the same way:

```bash
python PNS_with_nucleosome_peak_calling.py \
    -b /path/to/bams/*.bam
```

Do not place quotes around `*.bam`, because the shell must expand the wildcard.

### Multiple explicitly named BAMs

```bash
python PNS_with_nucleosome_peak_calling.py \
    -b sample1.bam sample2.bam sample3.bam
```

### Selected chromosomes

```bash
python PNS_with_nucleosome_peak_calling.py \
    -b sample.bam \
    -c chr1-22,chrX
```

### A genomic interval

Coordinates use a 0-based start and end-exclusive end:

```bash
python PNS_with_nucleosome_peak_calling.py \
    -b sample.bam \
    -c chr12:52621135-52641135
```

### Custom fragment-length range

```bash
python PNS_with_nucleosome_peak_calling.py \
    -b sample.bam \
    --frag-lower 147 \
    --frag-upper 187
```

### Pool BAMs and deduplicate across all inputs

```bash
python PNS_with_nucleosome_peak_calling.py \
    -b /path/to/bams/*.bam \
    --dedup-scope all_bams \
    --max-duplicates 1
```

### Generate WW/SS type outputs

```bash
python PNS_with_nucleosome_peak_calling.py \
    -b sample.bam \
    --split-ww-types \
    --fasta hg38.fa
```

### Generate aligned dinucleotide profiles

```bash
python PNS_with_nucleosome_peak_calling.py \
    -b sample.bam \
    --dinuc-profile \
    --fasta hg38.fa
```

### Generate both WW/SS types and dinucleotide profiles

```bash
python PNS_with_nucleosome_peak_calling.py \
    -b sample.bam \
    --split-ww-types \
    --dinuc-profile \
    --fasta hg38.fa
```

### Disable PNS scoring

This can be useful when only coverage, dyad, fragment-end, WW/SS, or dinucleotide outputs are required:

```bash
python PNS_with_nucleosome_peak_calling.py \
    -b sample.bam \
    --pns-mode off
```

---

## Main options

### Input and regions

```text
-b, --bamfiles       One or more paired-end BAM files
-o, --out-prefix     Output prefix
-c, --contigs        Contigs, ranges, aliases, or genomic intervals
```

### Fragment filtering

```text
--frag-lower         Minimum fragment length
--frag-upper         Maximum fragment length
--max-duplicates     Number of additional identical coordinate copies allowed
--dedup-scope        Apply coordinate deduplication within each BAM or across all BAMs
--subsample          Randomly retain fragments with probability p
```

`--max-duplicates 0` disables coordinate-based deduplication and retains all fragments.  
`--max-duplicates 1` retains one copy of each identical fragment coordinate.  
`--max-duplicates 2` retains up to two copies, and so on.

### PNS scoring

```text
--mode-length        Fragment length used to define the PNS kernel geometry
--pns-mode           Enable or disable PNS scoring
```

### Fragment randomisation

```text
--randomize-mode     none, uniform, or dinuc_anchor
--fasta              Indexed reference FASTA
--anchor-prob-start  Probability of preserving the start-boundary dinucleotide
--max-anchor-tries   Maximum placement attempts for dinucleotide anchoring
--randomize-fallback Action when no valid anchored placement is found
--seed               Random seed
```

### Sequence analysis

```text
--split-ww-types     Produce separate Type 1-4 WW/SS outputs
--dinuc-profile      Calculate dyad-aligned dinucleotide profiles
--dinuc-fraction     Write dinucleotide values as fractions rather than percentages
--fasta              Indexed reference FASTA required for sequence-based analyses
```

### Output controls

```text
--pns-format         Output format for PNS tracks
--other-format       Output format for non-PNS tracks
--pns-tracks         Select which PNS tracks to write
--other-tracks       Select coverage, dyad, and fragment-end tracks
```

BigWig is the default output format. Compressed WIG output can also be requested where supported.

### Peak calling

```text
--peak-format        Peak output format
--min-region-length  Minimum length of a called positive region
--max-neg-run        Number of non-positive bases tolerated within a region
--peak-score-scale   Scale factor used for BED scores
```

---

## Output files

Output names include the selected PNS mode length and fragment-length range.

A typical prefix has the form:

```text
<out_prefix>_PNS_mode<MODE>_lower<LOWER>_upper<UPPER>
```

The exact files produced depend on the selected tracks and analysis options.

### PNS tracks

```text
*_pns.bw
*_pns_smoothed.bw
*_posPNS.bw
```

### Fragment tracks

```text
*_coverage.bw
*_dyad.bw
*_fragment_ends.bw
*_fragment_left_ends.bw
*_fragment_right_ends.bw
```

Optional compressed WIG outputs use:

```text
*.wig.gz
```

### Peak calls

```text
*_nucleosome_regions.bed
*_breakpoint_peaks.bed
```

### Sequence-analysis outputs

```text
*_dinuc_profile.tsv
*_ww_type_summary.tsv
```

When `--split-ww-types` is enabled, Type 1-4 outputs are written separately in addition to the combined all-fragment outputs.

### Fragment summaries

```text
*_fragment_summary.txt
*_fragment_length_counts.tsv
```

The summary files report fragment totals, fragment-length counts, and other run-level statistics.

---

## Pipeline overview

1. Read paired-end alignments from one or more indexed BAM files.
2. Reconstruct fragments and apply alignment and fragment-length filters.
3. Apply coordinate deduplication and optional subsampling.
4. Optionally randomise fragment positions.
5. Optionally classify fragments into WW/SS Types 1-4.
6. Calculate the requested PNS and fragment-based tracks.
7. Optionally calculate dyad-aligned dinucleotide profiles.
8. Smooth the raw PNS signal.
9. Call nucleosome protection and breakpoint peaks.
10. Write BigWig or WIG tracks, BED peak calls, and summary files.

---

## Notes

- BAM files must be coordinate sorted and indexed.
- A wildcard such as `*.bam` is expanded by the shell, not by Python.
- Region coordinates use a 0-based start and an end-exclusive end.
- Reference-based modes require a FASTA file with a corresponding `.fai` index.
- Large contigs are processed in overlapping chunks, but only the non-overlapping core of each chunk is written.
