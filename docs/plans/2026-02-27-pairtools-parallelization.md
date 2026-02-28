# pairtools Parallelization Implementation Plan

> **For Claude:** REQUIRED SUB-SKILL: Use superpowers:executing-plans to implement this plan task-by-task.

**Goal:** Parallelize the pairtools parse+sort bottleneck by running it per chunk, then merge sorted pairsams before a single-stream dedup.

**Architecture:** Remove the intermediate uncompressed `.pairsam` file by piping parse into sort. Keep chunks separate through parse+sort (instead of merging BAMs first), then use `pairtools merge` to merge-sort all chunk outputs before dedup. When `chunk_size: 0` (single chunk), the pipeline degenerates to a simple pipe with no merge overhead.

**Tech Stack:** Snakemake, pairtools, samtools, bgzip (already in environment via htslib)

---

## Overview of Snakefile changes

| Old rule | Action | New rule |
|---|---|---|
| `pairtools_parse` | remove | — |
| `pairtools_sort` | remove | — |
| `prepare_work_bam` | remove | — |
| `merge_chunk_bams` | keep, but no longer temp | unchanged |
| *(new)* | add | `pairtools_parse_sort_chunk` |
| *(new)* | add | `pairtools_parse_sort_prev` |
| *(new)* | add | `pairtools_merge_sorted_chunks` |
| `pairtools_dedup` | update input path | unchanged logic |

New intermediate file: `03.mapping/chunks/{sample}/{sample}.{chunk}_{genome}_sorted.pairsam.gz`
(temp, bgzip-compressed sorted pairsam per chunk)

---

### Task 1: Add helper functions for chunk pairsam collection

**Files:**
- Modify: `src/rupture/workflow/Snakefile` — add two helpers after the existing `get_chunk_qc_logs` function (around line 75)

**Step 1: Add `get_chunk_sorted_pairsams` and `get_all_sorted_pairsams_for_merge` helpers**

Insert after the `get_chunk_qc_logs` function:

```python
def get_chunk_sorted_pairsams(wildcards):
    """Collect all per-chunk sorted pairsams for merging."""
    chunks = _get_chunks_for_sample(wildcards)
    return expand(
        "03.mapping/chunks/{sample}/{sample}.{chunk}_{genome}_sorted.pairsam.gz",
        sample=wildcards.sample, chunk=chunks, genome=wildcards.genome,
    )

def get_all_sorted_pairsams_for_merge(wildcards):
    """Collect chunk pairsams + optional prev_mapping pairsam."""
    pairsams = get_chunk_sorted_pairsams(wildcards)
    prev_bam = f"07.prev_mapping/{wildcards.sample}_{wildcards.genome}.bam"
    if os.path.isfile(prev_bam):
        pairsams = list(pairsams) + [
            f"03.mapping/{wildcards.sample}_{wildcards.genome}.prev_sorted.pairsam.gz"
        ]
    return pairsams
```

**Step 2: Verify syntax**

```bash
cd /gpfs/commons/groups/ren_lab/yangxie/projects/git/Rupture
python -c "import ast; ast.parse(open('src/rupture/workflow/Snakefile').read()); print('OK')"
```

Expected: `OK`

---

### Task 2: Add `pairtools_parse_sort_chunk` rule

**Files:**
- Modify: `src/rupture/workflow/Snakefile` — add after the `### Step 5: add_field_chunk` rule (around line 203), before the merge stage section

**Step 1: Add the rule**

```python
### Step 5c: parse + sort per chunk — streams BAM → sorted compressed pairsam
rule pairtools_parse_sort_chunk:
    input:
        bam="03.mapping/chunks/{sample}/{sample}.{chunk}_{genome}.bam"
    output:
        sorted_pairsam=temp("03.mapping/chunks/{sample}/{sample}.{chunk}_{genome}_sorted.pairsam.gz"),
        stats="03.mapping/chunks/{sample}/{sample}.{chunk}_{genome}.pairparse.txt"
    threads: 16
    shell:
        '''
        samtools view -h {input.bam} | \
        pairtools parse --min-mapq 30 --walks-policy all \
            --nproc-in {threads} --nproc-out {threads} \
            --max-inter-align-gap 30 \
            --chroms-path {chrom_size} \
            --assembly {wildcards.genome} \
            --output-stats {output.stats} \
            --add-columns CB | \
        pairtools sort --nproc {threads} --memory 16G \
            --tmpdir 03.mapping/chunks/{wildcards.sample} \
            --compress-prog bgzip \
            -o {output.sorted_pairsam}
        '''
```

**Step 2: Verify syntax**

```bash
python -c "import ast; ast.parse(open('src/rupture/workflow/Snakefile').read()); print('OK')"
```

Expected: `OK`

---

### Task 3: Add `pairtools_parse_sort_prev` rule

**Files:**
- Modify: `src/rupture/workflow/Snakefile` — add immediately after `pairtools_parse_sort_chunk`

**Step 1: Add the rule**

```python
### Optional: parse + sort a previous-mapping BAM from 07.prev_mapping/
rule pairtools_parse_sort_prev:
    input:
        bam="07.prev_mapping/{sample}_{genome}.bam"
    output:
        temp("03.mapping/{sample}_{genome}.prev_sorted.pairsam.gz")
    threads: 16
    shell:
        '''
        samtools view -h {input.bam} | \
        pairtools parse --min-mapq 30 --walks-policy all \
            --nproc-in {threads} --nproc-out {threads} \
            --max-inter-align-gap 30 \
            --chroms-path {chrom_size} \
            --assembly {wildcards.genome} \
            --add-columns CB | \
        pairtools sort --nproc {threads} --memory 16G \
            --tmpdir 03.mapping \
            --compress-prog bgzip \
            -o {output}
        '''
```

**Step 2: Verify syntax**

```bash
python -c "import ast; ast.parse(open('src/rupture/workflow/Snakefile').read()); print('OK')"
```

Expected: `OK`

---

### Task 4: Add `pairtools_merge_sorted_chunks` rule

**Files:**
- Modify: `src/rupture/workflow/Snakefile` — add in the `### Merge stage` section, after `merge_qc_logs`

**Step 1: Add the rule**

```python
### Merge sorted per-chunk pairsams into one sorted pairsam
rule pairtools_merge_sorted_chunks:
    input:
        get_all_sorted_pairsams_for_merge
    output:
        temp("03.mapping/{sample}_{genome}_sorted.pairsam.gz")
    threads: 16
    run:
        if len(input) == 1:
            shell("ln -sf $(realpath {input[0]}) {output[0]}")
        else:
            shell(
                "pairtools merge --nproc-in {threads} --nproc-out {threads} "
                "--memory 16G --tmpdir 03.mapping "
                "--compress-prog bgzip "
                "-o {output} {input}"
            )
```

**Step 2: Verify syntax**

```bash
python -c "import ast; ast.parse(open('src/rupture/workflow/Snakefile').read()); print('OK')"
```

Expected: `OK`

---

### Task 5: Update `pairtools_dedup` input and remove old rules

**Files:**
- Modify: `src/rupture/workflow/Snakefile`

**Step 1: Update `pairtools_dedup` to consume the new merged sorted pairsam**

The rule currently reads:
```python
rule pairtools_dedup:
    input:
        "03.mapping/{sample}_{genome}_sorted.pairsam"
```

Change the input path to:
```python
rule pairtools_dedup:
    input:
        "03.mapping/{sample}_{genome}_sorted.pairsam.gz"
```

(Only the input filename changes — the `.gz` extension. Everything else stays identical.)

**Step 2: Remove `prepare_work_bam` rule entirely**

Delete the entire `rule prepare_work_bam:` block (currently around lines 236–246):

```python
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
```

**Step 3: Remove `pairtools_parse` rule entirely**

Delete the entire `rule pairtools_parse:` block (currently around lines 249–267).

**Step 4: Remove `pairtools_sort` rule entirely**

Delete the entire `rule pairtools_sort:` block (currently around lines 270–278).

**Step 5: Verify syntax**

```bash
python -c "import ast; ast.parse(open('src/rupture/workflow/Snakefile').read()); print('OK')"
```

Expected: `OK`

---

### Task 6: Clean up now-unused helpers and update `merge_chunk_bams`

**Files:**
- Modify: `src/rupture/workflow/Snakefile`

**Step 1: Remove `get_merge_inputs` helper** (used only by the now-deleted `prepare_work_bam`)

Delete the `get_merge_inputs` function (around lines 43–49):

```python
def get_merge_inputs(wildcards):
    """If a previous BAM exists in 07.prev_mapping/, include it for merging."""
    new_bam = f"03.mapping/{wildcards.sample}_{wildcards.genome}.bam"
    old_bam = f"07.prev_mapping/{wildcards.sample}_{wildcards.genome}.bam"
    if os.path.isfile(old_bam):
        return [new_bam, old_bam]
    return [new_bam]
```

**Step 2: Remove `temp()` wrapper from `merge_chunk_bams` output**

The merged BAM is a final pipeline output, not a temp file. It was already non-temp (check the current rule — `output: "03.mapping/{sample}_{genome}.bam"` with no `temp()`), so this step may be a no-op. Confirm it is not wrapped in `temp()`.

**Step 3: Verify syntax**

```bash
python -c "import ast; ast.parse(open('src/rupture/workflow/Snakefile').read()); print('OK')"
```

Expected: `OK`

---

### Task 7: Dry-run test

**Step 1: Set up a minimal test config**

```bash
cd /gpfs/commons/groups/ren_lab/yangxie/projects/colab/Doege_obesity_DHC
```

**Step 2: Dry-run with existing config**

```bash
rupture run --configfile config.yaml --dryrun -j 32 2>&1 | head -80
```

Expected output: Snakemake prints a DAG summary showing `pairtools_parse_sort_chunk` jobs (one per chunk per sample), `pairtools_merge_sorted_chunks`, and `pairtools_dedup`. No mention of old `pairtools_parse`, `pairtools_sort`, or `prepare_work_bam` rules. No errors.

**Step 3: If there are errors**

Common issues and fixes:
- `KeyError` on a wildcard → check wildcard_constraints covers new file patterns
- `MissingInputException` → check the helper function returns the right paths
- Circular dependency → check no rule still references the old `{sample}_work_{genome}.bam` temp file

---

### Task 8: Commit

```bash
cd /gpfs/commons/groups/ren_lab/yangxie/projects/git/Rupture
git add src/rupture/workflow/Snakefile
git commit -m "feat: parallelize pairtools parse+sort per chunk, merge before dedup

- Pipe parse→sort per chunk (no uncompressed .pairsam intermediate)
- bgzip-compress sorted per-chunk pairsams (~4-5x smaller)
- pairtools_merge_sorted_chunks merge-sorts all chunks before dedup
- Remove prepare_work_bam; handle 07.prev_mapping at pairsam level
- Remove pairtools_parse and pairtools_sort single-stream rules

Co-Authored-By: Claude Sonnet 4.6 <noreply@anthropic.com>"
```

---

## Notes for the implementer

- **`pairtools merge` behavior:** It performs a merge-sort on already-sorted inputs — this is O(N) not O(N log N), so it's fast even for large files.
- **Single-chunk case:** When `chunk_size: 0` or the FASTQ is small enough to be one chunk, `get_all_sorted_pairsams_for_merge` returns a list of length 1. The `pairtools_merge_sorted_chunks` rule symlinks instead of running pairtools merge — zero overhead.
- **bgzip threading:** `pairtools sort --compress-prog bgzip` uses bgzip for output. bgzip is single-threaded by default but the `--nproc-out` flag in pairtools controls the number of output decompressor processes used when reading the result back in downstream steps.
- **Existing large files:** The current `CL010_hg38_sorted.pairsam` (2.9T uncompressed) and `CL011_hg38.pairsam` (4.5T) in the working directory are stale intermediates. After implementing this plan and re-running, Snakemake will recompute them. The stale files can be deleted manually to free disk space before re-running.
