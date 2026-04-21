# Chunked Pipeline Implementation Plan

> **For Claude:** REQUIRED SUB-SKILL: Use superpowers:executing-plans to implement this plan task-by-task.

**Goal:** Split large FASTQs into chunks, process steps 1-5 per-chunk in parallel, merge chunk BAMs, then continue with the single-stream pairtools pipeline.

**Architecture:** At Snakefile load time, compute chunk counts per sample from FASTQ file sizes and a configurable `chunk_size` (GB). A new `split_fastq` rule splits R1/R2/R3 via `seqkit split2`. Existing rules (align_bc through add_field) gain a `{chunk}` wildcard. A new `merge_chunk_bams` rule combines chunk BAMs with `samtools merge`. The existing `prepare_work_bam` through `pair2cool` rules remain unchanged.

**Tech Stack:** seqkit (splitting), samtools (BAM merge), Snakemake wildcards

---

### Task 1: Add `chunk_size` to config template and `seqkit` to environment.yaml

**Files:**
- Modify: `src/rupture/data/config_template.yaml`
- Modify: `environment.yaml`

**Step 1: Update config template**

Add to `src/rupture/data/config_template.yaml` after the `chrom_size` line:

```yaml

# ── Optional: chunked processing for large FASTQs ──
chunk_size: 0            # max FASTQ size in GB per chunk (0 = no splitting)
```

**Step 2: Add seqkit to environment.yaml**

Add `- seqkit` under the `# BAM/SAM utilities` section in `environment.yaml`.

**Step 3: Commit**

```bash
git add src/rupture/data/config_template.yaml environment.yaml
git commit -m "feat: add chunk_size config option and seqkit dependency"
```

---

### Task 2: Rewrite Snakefile with chunked pipeline

**Files:**
- Modify: `src/rupture/workflow/Snakefile`

This is the core change. The Snakefile gets:
1. A chunk computation block at the top
2. A new `split_fastq` rule
3. Per-chunk versions of rules 1-5 (with `{chunk}` wildcard)
4. A new `merge_chunk_bams` rule
5. Rules 6-13 unchanged (operating on merged BAM)

**Step 1: Add chunk computation after config section**

Insert after the `chrom_size = config["chrom_size"]` line:

```python
import math

# ── Chunked processing ──
chunk_size_gb = float(config.get("chunk_size", 0))

def _get_chunk_ids(sample):
    """Compute chunk IDs for a sample based on R1 FASTQ size."""
    r1 = f"01.rawdata/{sample}_R1.fq.gz"
    if chunk_size_gb <= 0 or not os.path.isfile(r1):
        return ["part_001"]
    size_gb = os.path.getsize(r1) / (1024 ** 3)
    n = max(1, math.ceil(size_gb / chunk_size_gb))
    return [f"part_{i:03d}" for i in range(1, n + 1)]

# Pre-compute chunk IDs for all samples
CHUNKS = {s: _get_chunk_ids(s) for s in config["sample_id"]}
```

**Step 2: Add `split_fastq` rule**

This rule splits R1, R2, R3 into N chunks using `seqkit split2`. For samples with 1 chunk, it symlinks instead.

```python
def get_n_chunks(wildcards):
    return len(CHUNKS[wildcards.sample])

### Step 0: split FASTQs into chunks
rule split_fastq:
    input:
        R1="01.rawdata/{sample}_R1.fq.gz",
        R2="01.rawdata/{sample}_R2.fq.gz",
        R3="01.rawdata/{sample}_R3.fq.gz"
    output:
        R1_chunks=temp(expand("01.rawdata/chunks/{{sample}}/{{sample}}_R1.{chunk}.fq.gz",
                              chunk=lambda w: CHUNKS[w.sample])),
        R2_chunks=temp(expand("01.rawdata/chunks/{{sample}}/{{sample}}_R2.{chunk}.fq.gz",
                              chunk=lambda w: CHUNKS[w.sample])),
        R3_chunks=temp(expand("01.rawdata/chunks/{{sample}}/{{sample}}_R3.{chunk}.fq.gz",
                              chunk=lambda w: CHUNKS[w.sample]))
    params:
        n_chunks=get_n_chunks,
        outdir="01.rawdata/chunks/{sample}"
    threads: 4
    run:
        import os, glob
        os.makedirs(params.outdir, exist_ok=True)
        if params.n_chunks == 1:
            # Single chunk: symlink
            for read in ["R1", "R2", "R3"]:
                src = os.path.realpath(f"01.rawdata/{wildcards.sample}_{read}.fq.gz")
                dst = f"{params.outdir}/{wildcards.sample}_{read}.part_001.fq.gz"
                os.symlink(src, dst)
        else:
            # Split each read file
            for read in ["R1", "R2", "R3"]:
                fq = f"01.rawdata/{wildcards.sample}_{read}.fq.gz"
                shell(
                    "seqkit split2 -p {params.n_chunks} -j {threads} "
                    "-O {params.outdir} -1 {fq}"
                )
        # Rename seqkit output to match expected Snakemake output names
        for f in glob.glob(f"{params.outdir}/*.part_*.fq.gz"):
            # seqkit produces: sample_R1.part_001.fq.gz -> we want sample_R1.part_001.fq.gz
            # naming already matches, nothing to do
            pass
```

**NOTE:** The output naming needs careful handling. seqkit produces `{basename}.part_NNN.{ext}` where basename is `sample_R1.fq` and ext is `.gz`. The actual output will be `sample_R1.fq.part_001.gz`. We need to account for this. The plan will use a rename step in the `run:` block to normalize names to `{sample}_{read}.part_NNN.fq.gz`.

Actually, let me simplify this. Rather than fighting seqkit's naming, we'll use `--by-part-prefix` or just rename after splitting. The cleaner approach: use the `run:` directive to split then rename.

**REVISED Step 2: Simplified split rule using input functions**

Since Snakemake `expand()` in output doesn't support lambdas, we use an input function approach instead. Each per-chunk rule will request its chunk file, and a single split rule produces all chunks:

```python
def get_chunk_inputs_R1(wildcards):
    return expand("01.rawdata/chunks/{sample}/{sample}_R1.{chunk}.fq.gz",
                  sample=wildcards.sample, chunk=CHUNKS[wildcards.sample])

def get_chunk_inputs_R2(wildcards):
    return expand("01.rawdata/chunks/{sample}/{sample}_R2.{chunk}.fq.gz",
                  sample=wildcards.sample, chunk=CHUNKS[wildcards.sample])

def get_chunk_inputs_R3(wildcards):
    return expand("01.rawdata/chunks/{sample}/{sample}_R3.{chunk}.fq.gz",
                  sample=wildcards.sample, chunk=CHUNKS[wildcards.sample])
```

Given the complexity of dynamic outputs with seqkit naming, the cleanest Snakemake pattern is:

1. Use a **checkpoint** rule for splitting (even though we pre-compute counts, this handles the file creation cleanly)
2. OR use a simpler approach: a `split_fastq` rule per sample that creates a **directory**, and downstream rules reference files in that directory.

**FINAL APPROACH — directory-based split:**

```python
rule split_fastq:
    input:
        R1="01.rawdata/{sample}_R1.fq.gz",
        R2="01.rawdata/{sample}_R2.fq.gz",
        R3="01.rawdata/{sample}_R3.fq.gz"
    output:
        directory("01.rawdata/chunks/{sample}")
    params:
        n_chunks=get_n_chunks
    threads: 4
    run:
        ...split and rename logic...
```

Then per-chunk rules reference `01.rawdata/chunks/{sample}/{sample}_{read}.{chunk}.fq.gz`.

**Step 3: Modify rules 1-5 to add `{chunk}` wildcard**

All per-chunk rules get a `{chunk}` wildcard and operate in `01.rawdata/chunks/` and temp chunk directories:

- `align_bc_chunk`: input `chunks/{sample}/{sample}_R2.{chunk}.fq.gz` -> output `chunks/{sample}/{chunk}_R2_BC.bam`
- `stamp_hic_chunk`: input chunk R1, R3, R2 BAM -> output chunk stamped R1, R3
- `trim_dna_chunk`: input chunk stamped -> output `02.trimmed/chunks/{sample}/{chunk}_*`
- `bwa_chunk`: input chunk trimmed -> output `03.mapping/chunks/{sample}/{chunk}_{genome}.raw.bam`
- `add_field_chunk`: input chunk raw BAM -> output `03.mapping/chunks/{sample}/{chunk}_{genome}.bam`

**Step 4: Add `merge_chunk_bams` rule**

```python
rule merge_chunk_bams:
    input:
        lambda w: expand("03.mapping/chunks/{sample}/{chunk}_{genome}.bam",
                         sample=w.sample, chunk=CHUNKS[w.sample], genome=w.genome)
    output:
        "03.mapping/{sample}_{genome}.bam"
    threads: 16
    run:
        if len(input) == 1:
            shell("ln -sf $(realpath {input[0]}) {output[0]}")
        else:
            shell("samtools merge -@ {threads} {output} {input}")
```

**Step 5: Rules 6-13 remain unchanged**

`prepare_work_bam` takes `03.mapping/{sample}_{genome}.bam` as before. The rest of the pipeline is identical.

**Step 6: Update `rule all`**

Remove `01.rawdata/{sample}_R1_BC_cov.fq.gz` and `01.rawdata/{sample}_R3_BC_cov.fq.gz` from `rule all` since stamped FASTQs now live in chunk directories and are temp. The final outputs (pairs, mcool, QC) remain the same.

**Step 7: Verify with dry run**

```bash
rupture run --configfile config.yaml --dryrun
```

**Step 8: Commit**

```bash
git add src/rupture/workflow/Snakefile
git commit -m "feat: chunked pipeline - split large FASTQs and merge BAMs before pairtools"
```

---

### Task 3: Update CLAUDE.md and version

**Files:**
- Modify: `CLAUDE.md`
- Modify: `src/rupture/__init__.py`
- Modify: `pyproject.toml`

**Step 1: Add chunked pipeline docs to CLAUDE.md**

Add a section describing the `chunk_size` config option and the chunked pipeline flow.

**Step 2: Commit and push**

```bash
git add CLAUDE.md src/rupture/__init__.py pyproject.toml
git commit -m "docs: update CLAUDE.md with chunked pipeline documentation"
git push origin main
```
