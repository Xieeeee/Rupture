# Design: pairtools Parallelization

Date: 2026-02-27

## Problem

For large datasets (400–750 GB merged BAMs), the single-stream `pairtools parse → sort → dedup`
block is the dominant bottleneck. Intermediate uncompressed `.pairsam` files balloon to multi-TB
on disk (4.5 TB observed for a 742 GB BAM), and sorting them sequentially exceeds 96-hour
SLURM wall-time limits. The per-chunk parallelism that already exists for alignment (steps 0–5)
stops at BAM merge and does not extend into pairtools.

## Approach

Two complementary changes, implemented together:

### A: Pipe parse→sort + lz4 compression

Collapse `pairtools_parse` and `pairtools_sort` into a single rule that pipes directly from
parse into sort, with no uncompressed intermediate written to disk. Sorted output is written
with `--compress-prog lz4` (fast compression, ~4-5x ratio at near-disk speed, natively
supported by pairtools).

- Eliminates the multi-TB uncompressed `.pairsam` intermediate
- Reduces sorted pairsam disk footprint by ~4-5x
- No structural pipeline change

### C: Extend chunk pipeline through pairtools parse+sort

Instead of merging chunk BAMs into one sample BAM and then running pairtools on it, keep
chunks separate through parse+sort, then merge all sorted chunk pairsams before dedup.

**New pipeline flow:**

```
per chunk BAM  →  pairtools parse | sort  →  sorted_chunk.pairsam.lz4
                                                       ↓ (all chunks)
                                             pairtools merge  →  merged sorted pairsam
                                                       ↓
                                             pairtools dedup  →  dedup.pairsam
```

- `merge_chunk_bams` (step 5b) is removed as a pipeline step; the merged BAM is no longer
  a required intermediate. The final `{sample}_{genome}.bam` can still be produced as a
  separate output rule if needed for downstream use.
- `prepare_work_bam` (optional `07.prev_mapping/` merging) is reworked: a previous BAM
  is also parsed+sorted as its own chunk and merged at the pairsam level.
- When `chunk_size: 0` (single chunk), the pipeline degenerates correctly to Approach A:
  one chunk → one parse+sort → straight to dedup with no merge overhead.

## Rules Changed / Added

| Old rule | New rule(s) | Notes |
|---|---|---|
| `pairtools_parse` | removed | collapsed into chunk rule |
| `pairtools_sort` | removed | collapsed into chunk rule |
| `merge_chunk_bams` | removed (or made optional output) | no longer a pipeline intermediate |
| `prepare_work_bam` | reworked | merges at pairsam level, not BAM level |
| *(new)* `pairtools_parse_sort_chunk` | per-chunk parallel parse+sort | streams BAM chunk → sorted lz4 pairsam |
| *(new)* `pairtools_merge_chunks` | merge-sorts all chunk pairsams | uses `pairtools merge` |
| `pairtools_dedup` | unchanged | receives merged sorted pairsam |

## Config

No new config parameters. `chunk_size` controls chunk count as before.

## Expected Impact

- Eliminates multi-TB uncompressed intermediates
- parse+sort steps parallelize across N chunks (N = ceil(R1_size / chunk_size))
- Dedup remains single-stream (unavoidable for CB-based single-cell dedup)
- Estimated wall-time reduction: 3–5x for large datasets
