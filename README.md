# Rupture

A Python pipeline for processing [Droplet Hi-C](https://www.nature.com/articles/s41587-024-02447-1) data. Takes three-read FASTQ input (R1 DNA, R2 barcode, R3 DNA mate) and produces multi-resolution contact matrices in mcool format.

## Overview

Rupture wraps a 13-step Snakemake workflow into an installable Python package with a single `rupture` CLI. The pipeline:

1. Aligns barcode reads (R2) to a whitelist and stamps barcodes into R1/R3 read headers
2. Trims and aligns DNA reads to a reference genome
3. Processes contacts with pairtools (parse, sort, per-cell dedup, split)
4. Generates QC summaries and multi-resolution `.mcool` contact matrices

## Installation

### With conda (recommended)

This installs both the external bioinformatics tools and the Python package:

```bash
git clone https://github.com/Xieeeee/Rupture.git
cd Rupture

conda env create -f environment.yaml
conda activate rupture
```

### With pip only

If you already have the required external tools available (bowtie, bwa, samtools, htslib, trim-galore, pairtools, pairix, cooler):

```bash
pip install .
```

### Verify installation

```bash
rupture --version
# rupture 0.1.0
```

## Quick start

```bash
# 1. Set up a working directory
rupture init --dir my_project
cd my_project

# 2. Edit config.yaml with your sample names and reference paths
#    (see Configuration below)

# 3. Place your FASTQ files
#    01.rawdata/sample1_R1.fq.gz
#    01.rawdata/sample1_R2.fq.gz
#    01.rawdata/sample1_R3.fq.gz

# 4. Dry run to verify the DAG
rupture run --configfile config.yaml --dryrun

# 5. Run the pipeline
rupture run --configfile config.yaml --cores 16
```

## Configuration

`rupture init` generates a `config.yaml` template. All parameters are required:

```yaml
sample_id:
  - sample1
  - sample2

bin_size: 5000       # base resolution for cooler (bp)
genome: "hg38"       # assembly label (used in output filenames)

bc_ref: "/path/to/barcode/bowtie_index_prefix"
ref: "/path/to/genome.fa"
chrom_size: "/path/to/chrom.sizes"
```

| Parameter | Description |
|-----------|-------------|
| `sample_id` | List of sample names matching FASTQ prefixes in `01.rawdata/` |
| `bin_size` | Base resolution in bp for cooler (default: 5000) |
| `genome` | Genome assembly label (e.g. `hg38`, `mm10`) |
| `bc_ref` | Bowtie v1 index prefix for barcode whitelist (e.g. 10x ARC or ATAC) |
| `ref` | Reference genome FASTA (must have a corresponding bwa index) |
| `chrom_size` | Tab-separated chromosome sizes file |

### Input files

Each sample requires three gzipped FASTQ files in `01.rawdata/`:

```
01.rawdata/{sample}_R1.fq.gz   # DNA read 1
01.rawdata/{sample}_R2.fq.gz   # Barcode read
01.rawdata/{sample}_R3.fq.gz   # DNA read 2 (mate)
```

If your data is split across multiple sequencing lanes, use `rupture merge-fastq` to combine them first (see [CLI reference](#cli-reference) below).

### Incremental merging

To merge new sequencing data with a previous run, place the prior BAM in `07.prev_mapping/`:

```
07.prev_mapping/{sample}_{genome}.bam
```

The pipeline will automatically detect and merge it before contact calling.

## Output

```
00.qc/                                       # Barcode stamping QC logs
03.mapping/{sample}_{genome}.bam             # Aligned BAM with CB tags
03.mapping/{sample}_{genome}.sc.pairs.gz     # Deduplicated contact pairs
03.mapping/{sample}_{genome}.sc.pairs.gz.px2 # Pairix index
03.mapping/{sample}_{genome}.sc.pairdedup.summary.txt   # QC summary
03.mapping/{sample}_{genome}.PairCount.stat.csv         # Per-cell contact counts
04.matrices/{sample}_{genome}.mcool          # Multi-resolution contact matrix
```

The `.mcool` file contains contact matrices at 9 resolutions: 5kb, 10kb, 25kb, 50kb, 100kb, 250kb, 500kb, 1Mb, and 2.5Mb.

## CLI reference

### `rupture run`

Run the full pipeline via Snakemake.

```bash
rupture run --configfile config.yaml --cores 16
rupture run --configfile config.yaml --dryrun
rupture run --configfile config.yaml --cores 8 --snakemake-args --forceall
```

### `rupture init`

Initialize a working directory with a config template and `01.rawdata/` folder.

```bash
rupture init --dir my_project
rupture init --dir my_project --force   # overwrite existing config.yaml
```

### `rupture stamp`

Filter R1/R3 FASTQ reads by mapped R2 barcodes and stamp barcodes into read headers.

```bash
rupture stamp \
    --bam R2_mapped.bam \
    --r1 sample_R1.fq.gz --r3 sample_R3.fq.gz \
    --out-r1 sample_R1_stamped.fq.gz --out-r3 sample_R3_stamped.fq.gz \
    --bc-tag BC --min-mapq 40
```

### `rupture add-tag`

Extract a field from each read name and write it as a BAM tag.

```bash
rupture add-tag --input raw.bam --output tagged.bam --field -1 --tag CB --sep ":"
```

### `rupture count-pairs`

Count per-cell contacts from a `.pairs` file.

```bash
rupture count-pairs --input sample.sc.pairs --output sample.PairCount
# writes sample.PairCount.stat.csv
```

### `rupture summarize-qc`

Compute QC summary metrics from pairtools dedup stats.

```bash
rupture summarize-qc --input sample.pairdedup.txt --output sample.summary.txt
```

### `rupture merge-fastq`

Merge FASTQ files across sequencing lanes for a library.

```bash
rupture merge-fastq --prefix CL010 --input-dir raw_lanes/ --output-dir 01.rawdata/
# Detects read types (R1, R2, R3, etc.) automatically
```

## Pipeline steps

| Step | Rule | Description |
|------|------|-------------|
| 1 | `align_bc` | Align R2 barcode reads to whitelist index (bowtie) |
| 2 | `stamp_hic` | Stamp matched barcodes into R1/R3 read names |
| 3 | `trim_dna` | Quality trim DNA reads (trim_galore, Q20) |
| 4 | `bwa` | Align to reference genome (bwa mem -SP5M) |
| 5 | `add_field` | Extract barcode from read name into CB BAM tag |
| 6 | `prepare_work_bam` | Merge with previous BAMs if available |
| 7 | `pairtools_parse` | Parse BAM to pairs format with CB columns |
| 8 | `pairtools_sort` | Sort pairs |
| 9 | `pairtools_dedup` | Per-cell deduplication using CB tags |
| 10 | `pairtools_split` | Split into .pairs + .bam |
| 11 | `QC` | Generate QC summary and per-cell contact counts |
| 12 | `bgzip` | Compress and index pairs (bgzip + pairix) |
| 13 | `pair2cool` | Build multi-resolution .mcool contact matrix |

## Dependencies

External tools (installed via conda):

- **Aligners**: bowtie (v1), bwa
- **BAM/SAM**: samtools, htslib (bgzip)
- **Trimming**: trim-galore
- **Hi-C**: pairtools, pairix, cooler
- **Workflow**: snakemake >=7.0

Python libraries (installed via pip):

- pysam >=0.20
- pandas >=1.4
- numpy >=1.21
- snakemake >=7.0
- pyyaml >=6.0

