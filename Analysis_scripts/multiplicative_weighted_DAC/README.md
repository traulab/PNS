# Streaming DAC from bigWig

This script calculates Distance Autocorrelation (DAC) values directly from bigWig signal tracks over genomic regions defined in a chromatin-state BED file.

---

## Requirements

### Python packages

Install required Python packages:

```bash
pip install numpy pandas tqdm matplotlib
```

### UCSC tools

This script requires:

- `bigWigToBedGraph`

Install UCSC utilities and ensure they are available in your `PATH`.

Example (Linux):

```bash
mamba install -c bioconda ucsc-bigwigtobedgraph
```

---

## Input Files

### 1. bigWig files

One or more bigWig signal tracks.

Example:

```bash
CH01_chr1.bw
CH01_chr2.bw
```

Wildcards are supported.

### 2. Chromatin BED file

BED file containing at least 4 columns:

```text
chromosome    start    end    state
```

Example:

```text
chr1    100000    105000    Active
chr1    200000    206000    Repressed
```

---

## Usage

Basic example:

```bash
python multiplicative_weighted_DAC.py \
    --bigwig "CH01_chr*.bw" \
    --chromatin_bed chromatin_states.bed
```

---

## Arguments

| Argument | Description |
|---|---|
| `--bigwig` | Space-separated bigWig patterns |
| `--chromatin_bed` | BED file containing chromatin states |
| `--dmax` | Maximum distance for DAC calculation (default: 1500) |
| `--value_limit` | Optional cap for signal values |
| `--min_region_length` | Minimum interval size to include (default: 2000) |
| `--convert_to_euchromatin` | Collapse selected ChromHMM states into `Euchromatin` |
| `--normalize_dac` | Apply opportunity-based normalization |

---

## Example Commands

### Basic DAC calculation

```bash
python multiplicative_weighted_DAC.py \
    --bigwig "CH01_chr*.bw" \
    --chromatin_bed states.bed
```

### DAC calculation with normalization

```bash
python multiplicative_weighted_DAC.py \
    --bigwig "CH01_chr*.bw" \
    --chromatin_bed states.bed \
    --normalize_dac
```

### DAC calculation with signal clipping

```bash
python multiplicative_weighted_DAC.py \
    --bigwig "CH01_chr*.bw" \
    --chromatin_bed states.bed \
    --value_limit 100
```

---

## Output

The script outputs one TSV file per chromatin state.

Example:

```text
CH01_Active_bigwig_streaming_DAC_values_normalized.tsv
CH01_Repressed_bigwig_streaming_DAC_values_normalized.tsv
```

Each file contains:

| Column | Description |
|---|---|
| `Distance` | Base-pair separation |
| `DAC Value` | Raw DAC score |
| `DAC Value Percent` | Percentage contribution to total DAC |

---

## DAC Calculation

For each pair of positions separated by distance `d`:

```text
DAC[d] += value1 × value2
```

The script processes intervals independently and accumulates DAC values across all intervals belonging to the same chromatin state.

---

## Notes

- Regions shorter than `--min_region_length` are skipped.
- Chromosome names are automatically standardized to `chr` format where possible.
- Normalization rescales DAC values according to the number of positional opportunities available at each distance.

---