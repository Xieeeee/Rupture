### Droplet Hi-C Processing Pipeline
### Author: Yang Xie (yangxie@nygenome.org)
###
### Usage (run from your working directory, no need to copy this file):
###   snakemake -s /gpfs/commons/home/yangxie/scripts/Paired-HiC/scHiC_snakemake/Snakefile_fixed.py \
###             --configfile config.yaml -j 16
###
### Required directory layout in working directory:
###   01.rawdata/{sample}_R1.fq.gz
###   01.rawdata/{sample}_R2.fq.gz
###   01.rawdata/{sample}_R3.fq.gz
###   config.yaml

import os

# ── Resolve paths relative to this Snakefile (not the working directory) ──
PIPELINE_DIR = workflow.basedir
PAIRED_HIC_DIR = os.path.dirname(PIPELINE_DIR)  # .../Paired-HiC/
SCRIPTS_DIR = os.path.dirname(PAIRED_HIC_DIR)   # .../scripts/

STAMP_HIC_PY      = os.path.join(PIPELINE_DIR, "bin", "stamp_hic.py")
HEADER_TO_FIELD_PY = os.path.join(SCRIPTS_DIR, "scifi", "scifi.header_to_field2.py")
SUMMARIZE_PAIRS_R  = os.path.join(PAIRED_HIC_DIR, "phc.summarize_pairs_lec.R")
COUNT_PAIRS_PY     = os.path.join(PAIRED_HIC_DIR, "phc.count_pairs_sc.py")

# ── Configuration ──
bin_size = int(config["bin_size"])

# Barcode reference (override in config with "bc_ref" key)
if "bc_ref" in config:
    bc_ref = config["bc_ref"]
elif config["mode"] == "arc":
    bc_ref = "/gpfs/commons/home/yangxie/genome/others/737K-crarc-whitelist/737K-arc-v1"
elif config["mode"] == "atac":
    bc_ref = "/gpfs/commons/home/yangxie/genome/others/737K-cratac-whitelist/737K-cratac-v1"
else:
    raise ValueError(f"Unknown mode: {config['mode']}. Use 'arc', 'atac', or set 'bc_ref' in config.")

# Genome reference (override in config with "ref" and "chrom_size" keys)
if "ref" in config and "chrom_size" in config:
    ref = config["ref"]
    chrom_size = config["chrom_size"]
elif config["genome"] == "mm10":
    ref = "/gpfs/commons/home/yangxie/genome/mm10/bwa/mm10.fa"
    chrom_size = "/gpfs/commons/home/yangxie/genome/mm10/mm10.main.chrom.sizes"
elif config["genome"] == "hg38":
    ref = "/gpfs/commons/home/yangxie/genome/hg38/bwa/hg38.fa"
    chrom_size = "/gpfs/commons/home/yangxie/genome/hg38/hg38.main.chrom.sizes"
elif config["genome"] == "mix":
    ref = "/projects/ps-renlab/y2xie/projects/genome_ref/GRCh38_and_mm10/fasta/genome.fa"
    chrom_size = "/projects/ps-renlab/y2xie/projects/genome_ref/GRCh38_and_mm10/star/chrNameLength.txt"
else:
    raise ValueError(f"Unknown genome: {config['genome']}. Use 'mm10', 'hg38', 'mix', or set 'ref'+'chrom_size' in config.")

# ── Helper: runtime check for previous mapping BAMs ──
def get_merge_inputs(wildcards):
    """If a previous BAM exists in 07.prev_mapping/, include it for merging."""
    new_bam = f"03.mapping/{wildcards.sample}_{wildcards.genome}.bam"
    old_bam = f"07.prev_mapping/{wildcards.sample}_{wildcards.genome}.bam"
    if os.path.isfile(old_bam):
        return [new_bam, old_bam]
    return [new_bam]

###
### Rules
###

wildcard_constraints:
    genome="[^._/]+",
    sample="[^._/]+"

rule all:
    input:
        expand("03.mapping/{sample}_{genome}.sc.pairs.gz", sample=config["sample_id"], genome=config["genome"]),
        expand("03.mapping/{sample}_{genome}.sc.pairs.gz.px2", sample=config["sample_id"], genome=config["genome"]),
        expand("04.matrices/{sample}_{genome}.mcool", sample=config["sample_id"], genome=config["genome"]),
        expand("03.mapping/{sample}_{genome}.sc.pairdedup.summary.txt", sample=config["sample_id"], genome=config["genome"]),
        expand("03.mapping/{sample}_{genome}.PairCount.stat.csv", sample=config["sample_id"], genome=config["genome"]),
        expand("03.mapping/{sample}_{genome}.sc.pairdedup.txt", sample=config["sample_id"], genome=config["genome"]),
        expand("03.mapping/{sample}_{genome}.bam", sample=config["sample_id"], genome=config["genome"]),
        expand("00.qc/{sample}_qc.log", sample=config["sample_id"]),
        expand("01.rawdata/{sample}_R1_BC_cov.fq.gz", sample=config["sample_id"]),
        expand("01.rawdata/{sample}_R3_BC_cov.fq.gz", sample=config["sample_id"])

### Step 1: align R2 (barcode reads) to barcode index
rule align_bc:
    input:
        R2_fastq="01.rawdata/{sample}_R2.fq.gz"
    output:
        R2_bam=temp("01.rawdata/{sample}_R2_BC.bam")
    threads: 16
    shell:
        """
        bowtie -p {threads} -S {bc_ref} {input.R2_fastq} | \
            samtools view -b -F 4 -@ {threads} -o {output.R2_bam} -
        """

### Step 2: stamp barcodes onto R1/R3 read headers
rule stamp_hic:
    input:
        R1_fastq="01.rawdata/{sample}_R1.fq.gz",
        R3_fastq="01.rawdata/{sample}_R3.fq.gz",
        R2_bam="01.rawdata/{sample}_R2_BC.bam"
    output:
        R1_bc_fastq="01.rawdata/{sample}_R1_BC_cov.fq.gz",
        R3_bc_fastq="01.rawdata/{sample}_R3_BC_cov.fq.gz",
        stats="00.qc/{sample}_qc.log"
    shell:
        """
        python {STAMP_HIC_PY} \
            --bam {input.R2_bam} \
            --r1 {input.R1_fastq} \
            --r3 {input.R3_fastq} \
            --out-r1 {output.R1_bc_fastq} \
            --out-r3 {output.R3_bc_fastq} \
            --bc-tag BC --min-mapq 40 \
            2>{output.stats}
        """

### Step 3: quality trimming
rule trim_dna:
    input:
        R1_cov="01.rawdata/{sample}_R1_BC_cov.fq.gz",
        R3_cov="01.rawdata/{sample}_R3_BC_cov.fq.gz"
    output:
        R1_trimmed=temp("02.trimmed/{sample}_R1_BC_cov_val_1.fq.gz"),
        R3_trimmed=temp("02.trimmed/{sample}_R3_BC_cov_val_2.fq.gz"),
        R1_stats=temp("02.trimmed/{sample}_R1_BC_cov.fq.gz_trimming_report.txt"),
        R3_stats=temp("02.trimmed/{sample}_R3_BC_cov.fq.gz_trimming_report.txt")
    threads: 16
    shell:
        "trim_galore --cores {threads} -q 20 --paired {input.R1_cov} {input.R3_cov} -o 02.trimmed/"

### Step 4: align to reference genome
rule bwa:
    input:
        R1_trimmed="02.trimmed/{sample}_R1_BC_cov_val_1.fq.gz",
        R3_trimmed="02.trimmed/{sample}_R3_BC_cov_val_2.fq.gz"
    output:
        bam=temp("03.mapping/{sample}_{genome}.raw.bam"),
        stats=temp("03.mapping/{sample}_{genome}.bwa.log")
    threads: 16
    shell:
        "(bwa mem -SP5M -T0 -t{threads} {ref} {input.R1_trimmed} {input.R3_trimmed} | samtools view -bhS - > {output.bam}) 2>{output.stats}"

### Step 5: extract barcode from read name into CB BAM tag
rule add_field:
    input:
        "03.mapping/{sample}_{genome}.raw.bam"
    output:
        "03.mapping/{sample}_{genome}.bam"
    shell:
        """
        python {HEADER_TO_FIELD_PY} --in {input} --header -1 --bam CB
        mv {input}.addfield.bam {output}
        """

### Step 6: merge with previous mapping if available, otherwise just link
rule prepare_work_bam:
    input:
        get_merge_inputs
    output:
        temp("03.mapping/{sample}_work_{genome}.bam")
    threads: 16
    run:
        if len(input) > 1:
            shell("samtools merge -@ {threads} {output} {input}")
        else:
            shell("ln -sf $(realpath {input[0]}) {output[0]}")

### Step 7: parse BAM to pairs format
rule pairtools_parse:
    input:
        "03.mapping/{sample}_work_{genome}.bam"
    output:
        pairsam=temp("03.mapping/{sample}_{genome}.pairsam"),
        stats="03.mapping/{sample}_{genome}.pairparse.txt"
    threads: 16
    shell:
        '''
        samtools view -h {input} | \
        pairtools parse --min-mapq 30 --walks-policy all \
            --nproc-in {threads} --nproc-out {threads} \
            --max-inter-align-gap 30 \
            --chroms-path {chrom_size} \
            --assembly {wildcards.genome} \
            --output-stats {output.stats} \
            --add-columns CB \
            -o {output.pairsam}
        '''

### Step 8: sort pairs
rule pairtools_sort:
    input:
        "03.mapping/{sample}_{genome}.pairsam"
    output:
        temp("03.mapping/{sample}_{genome}_sorted.pairsam")
    threads: 16
    shell:
        "pairtools sort --nproc {threads} --memory 16G --tmpdir=03.mapping -o {output} {input}"

### Step 9: deduplicate
rule pairtools_dedup:
    input:
        "03.mapping/{sample}_{genome}_sorted.pairsam"
    output:
        dedup_pairsam=temp("03.mapping/{sample}_{genome}_dedup.pairsam"),
        stats="03.mapping/{sample}_{genome}.sc.pairdedup.txt"
    threads: 16
    shell:
        '''
        pairtools dedup --nproc-in {threads} --nproc-out {threads} \
            --extra-col-pair "CB1" "CB2" \
            --mark-dups \
            -o {output.dedup_pairsam} \
            --output-stats {output.stats} \
            {input}
        '''

### Step 10: split deduped pairsam into pairs + BAM
rule pairtools_split:
    input:
        "03.mapping/{sample}_{genome}_dedup.pairsam"
    output:
        pairs=temp("03.mapping/{sample}_{genome}.sc.pairs"),
        pairbam="03.mapping/{sample}_{genome}.sc.pairtools.bam"
    threads: 16
    shell:
        '''
        pairtools split --nproc-in {threads} --nproc-out {threads} \
            --output-pairs {output.pairs} \
            --output-sam - {input} | \
            samtools view -bS -@ {threads} | \
            samtools sort -T 03.mapping/ -@ {threads} -o {output.pairbam}
        '''

### Step 11: QC statistics
rule QC:
    input:
        stats="03.mapping/{sample}_{genome}.sc.pairdedup.txt",
        pairs="03.mapping/{sample}_{genome}.sc.pairs"
    output:
        summary="03.mapping/{sample}_{genome}.sc.pairdedup.summary.txt",
        count_file="03.mapping/{sample}_{genome}.PairCount.stat.csv"
    params:
        count_prefix=lambda w: f"03.mapping/{w.sample}_{w.genome}.PairCount"
    shell:
        '''
        Rscript {SUMMARIZE_PAIRS_R} {input.stats}
        python {COUNT_PAIRS_PY} --input {input.pairs} --output {params.count_prefix}
        '''

### Step 12: compress and index pairs
rule bgzip:
    input:
        "03.mapping/{sample}_{genome}.sc.pairs"
    output:
        pairs="03.mapping/{sample}_{genome}.sc.pairs.gz",
        index="03.mapping/{sample}_{genome}.sc.pairs.gz.px2"
    shell:
        '''
        bgzip -c {input} > {output.pairs}
        pairix -f {output.pairs}
        '''

### Step 13: generate multi-resolution contact matrices
rule pair2cool:
    input:
        "03.mapping/{sample}_{genome}.sc.pairs.gz"
    output:
        cool=temp(f"04.matrices/{{sample}}_{{genome}}_{bin_size}.cool"),
        mcool="04.matrices/{sample}_{genome}.mcool"
    threads: 16
    shell:
        '''
        cooler cload pairix {chrom_size}:{bin_size} {input} {output.cool}
        cooler zoomify --balance --balance-args '--convergence-policy store_nan' \
            -p {threads} -o {output.mcool} \
            -r 5000,10000,25000,50000,100000,250000,500000,1000000,2500000 \
            {output.cool}
        '''
