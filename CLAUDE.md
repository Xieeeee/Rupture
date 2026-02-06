# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

**Rupture** — an installable Python package for single-cell Hi-C (scHiC) processing of droplet-based Hi-C data. Converts raw FASTQ reads into multi-resolution contact matrices (cooler/mcool format) for 3D genome structure analysis. Authored by Yang Xie.

## Installation

```bash
# Create conda environment with external tools
conda env create -f environment.yaml
conda activate rupture

# Or install just the Python package (requires external tools already available)
pip install -e .
```

## Running the Pipeline

```bash
# Initialize a working directory
rupture init --dir /path/to/workdir

# Edit config.yaml with your sample IDs and reference paths
# Place FASTQs in 01.rawdata/{sample}_R1.fq.gz, _R2.fq.gz, _R3.fq.gz

# Run the pipeline
rupture run --configfile config.yaml --cores 16

# Dry run to check DAG
rupture run --configfile config.yaml --dryrun
```

## CLI Subcommands

| Command | Description |
|---------|-------------|
| `rupture run` | Invoke the bundled Snakefile via Snakemake |
| `rupture init` | Generate config template in a working directory |
| `rupture stamp` | Filter R1/R3 by mapped R2 BAM, stamp barcodes into headers |
| `rupture add-tag` | Extract field from read name into a BAM tag |
| `rupture count-pairs` | Count per-cell contacts from .pairs file |
| `rupture summarize-qc` | Compute QC summary from pairtools dedup stats |
| `rupture merge-fastq` | Merge FASTQ files across lanes |

## Configuration

Edit `config.yaml` with these parameters:

| Parameter | Required | Description |
|-----------|----------|-------------|
| `sample_id` | yes | List of sample names matching FASTQ prefixes |
| `bin_size` | yes | Base resolution in bp (default: 5000) |
| `genome` | yes | Genome assembly label (used in output filenames) |
| `bc_ref` | yes | Barcode bowtie v1 index prefix |
| `ref` | yes | Genome FASTA path (must have bwa index) |
| `chrom_size` | yes | Chromosome sizes file |

## Pipeline Architecture

Three-read input: R1 (DNA), R2 (barcode), R3 (DNA mate). The pipeline:

1. **Barcode extraction** (steps 1-2): Aligns R2 to barcode whitelist via bowtie, stamps barcodes into R1/R3 read headers
2. **Alignment** (steps 3-5): Trims reads (trim_galore), aligns to genome (bwa mem -SP5M), adds CB BAM tag
3. **Pairs processing** (steps 6-10): Optionally merges with prior BAMs from `07.prev_mapping/`, then pairtools parse → sort → dedup (using CB tags for single-cell dedup) → split into .pairs + .bam
4. **Output generation** (steps 11-13): QC stats, bgzip/pairix indexing, cooler cload + zoomify to produce .mcool at 9 resolutions (5kb–2.5Mb)

Key design: barcodes are stamped into read names early (step 2) with `:` separator and carried through as CB tags for per-cell deduplication in pairtools.

## Package Structure

```
Rupture/
├── pyproject.toml                # Package metadata & entry points
├── environment.yaml              # Conda env for external tools
├── src/rupture/
│   ├── __init__.py               # __version__
│   ├── cli/                      # CLI subcommands (argparse-based)
│   │   ├── main.py               # Entry point dispatching subcommands
│   │   ├── cmd_run.py            # rupture run
│   │   ├── cmd_init.py           # rupture init
│   │   ├── cmd_stamp.py          # rupture stamp
│   │   ├── cmd_addtag.py         # rupture add-tag
│   │   ├── cmd_count.py          # rupture count-pairs
│   │   ├── cmd_summarize.py      # rupture summarize-qc
│   │   └── cmd_merge_fastq.py    # rupture merge-fastq
│   ├── core/                     # Reusable logic (importable)
│   │   ├── stamp.py              # Barcode stamping
│   │   ├── addtag.py             # BAM tag from read name field
│   │   ├── count_pairs.py        # Per-cell pair counting
│   │   ├── summarize_qc.py       # QC summary (Python rewrite of R script)
│   │   └── merge_fastq.py        # Lane merging
│   ├── workflow/
│   │   ├── __init__.py           # get_snakefile_path(), get_config_template_path()
│   │   └── Snakefile             # Adapted pipeline (no hardcoded paths)
│   └── data/
│       └── config_template.yaml  # Template for rupture init
└── tests/
```

## Directory Layout (working directory)

```
00.qc/          - Barcode stamping QC logs
01.rawdata/     - Input FASTQs + barcode-stamped FASTQs
02.trimmed/     - Trimmed reads (temp, auto-cleaned)
03.mapping/     - BAMs, pairs, dedup stats, QC summaries
04.matrices/    - Final .cool and .mcool contact matrices
07.prev_mapping/ - (optional) Previous BAMs for incremental merging
```

## Dependencies

Conda environment `rupture` (defined in `environment.yaml`):
- **Workflow**: snakemake >=7.0
- **Alignment**: bowtie (v1), bwa, samtools
- **Trimming**: trim-galore
- **Hi-C tools**: pairtools, pairix, cooler
- **Python**: >=3.9 (numpy, pandas, pysam, pyyaml)
