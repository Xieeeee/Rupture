"""Barcode stamping: filter R1/R3 by mapped R2 BAM and append barcode to headers.

Uses a streaming sort-merge join to avoid loading all read IDs into memory.
Delegates heavy sorting to GNU sort (external merge sort), keeping Python
memory at O(1). Trade-off: uses temp disk instead of RAM.

Performance optimisations
-------------------------
- FASTQ linearisation via awk (avoids per-record Python overhead).
- gzip output via pigz when available (multi-threaded compression).
- R1 and R3 are processed in parallel (ThreadPoolExecutor).
"""

import gzip
import os
import shutil
import subprocess
import sys
import tempfile
from concurrent.futures import ThreadPoolExecutor

import pysam


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def core_name(n: str) -> str:
    """Normalize read IDs: strip whitespace; remove trailing /1, /2, /3."""
    n = n.split()[0]
    if len(n) >= 2 and n[-2] == "/" and n[-1] in "123":
        n = n[:-2]
    return n


_PIGZ = shutil.which("pigz")


class _PigzWriter:
    """Context manager for writing gzip-compressed data via pigz subprocess."""

    def __init__(self, path, threads=4, compresslevel=5):
        self._outfile = open(path, "wb")
        self._proc = subprocess.Popen(
            [_PIGZ, "-p", str(threads), f"-{compresslevel}", "-c"],
            stdin=subprocess.PIPE,
            stdout=self._outfile,
            text=True,
        )

    def write(self, data):
        return self._proc.stdin.write(data)

    def __enter__(self):
        return self

    def __exit__(self, *exc):
        self._proc.stdin.close()
        self._proc.wait()
        self._outfile.close()
        return False


def open_fastq_writer(path: str, threads=4):
    """Open a FASTQ writer.  Uses pigz for .gz when available."""
    if path.endswith(".gz"):
        if _PIGZ:
            return _PigzWriter(path, threads=threads)
        return gzip.open(path, "wt", compresslevel=5)
    return open(path, "w")


# ---------------------------------------------------------------------------
# Phase 1: extract barcodes from BAM
# ---------------------------------------------------------------------------

def _extract_barcodes_to_tsv(
    bam_path, tsv_path, bc_tag="BC", min_mapq=0, primary_only=True, threads=None
):
    """Stream BAM and write rid\\tbarcode\\tmapq lines to a TSV file.

    Returns the number of records written.
    """
    if threads is None:
        threads = os.cpu_count()
    n = 0
    mode = "rb" if bam_path.endswith((".bam", ".cram")) else "r"
    with pysam.AlignmentFile(bam_path, mode) as aln, open(tsv_path, "w") as out:
        try:
            aln.set_threads(max(1, int(threads)))
        except Exception:
            pass
        for rec in aln.fetch(until_eof=True):
            if rec.is_unmapped:
                continue
            if primary_only and (rec.is_secondary or rec.is_supplementary):
                continue
            if rec.mapping_quality < min_mapq:
                continue

            rid = core_name(rec.query_name)

            if bc_tag and rec.has_tag(bc_tag):
                bc = str(rec.get_tag(bc_tag))
            else:
                bc = rec.reference_name or "NA"

            mq = int(rec.mapping_quality)
            out.write(f"{rid}\t{bc}\t{mq}\n")
            n += 1
    return n


# ---------------------------------------------------------------------------
# Phase 2: sort and dedup barcodes
# ---------------------------------------------------------------------------

def _run_sort(input_path, output_path, sort_key, tmp_dir,
              parallel=None, sort_mem="2G"):
    """Sort a TSV file using GNU sort with LC_ALL=C.

    Parameters
    ----------
    sort_key : str
        Sort key spec for GNU sort, e.g. '-k1,1 -k3,3nr'.
    sort_mem : str
        Memory buffer for GNU sort (``-S`` flag), e.g. '2G'.
    """
    if parallel is None:
        parallel = os.cpu_count()
    cmd = [
        "sort",
        "--parallel", str(parallel),
        "-T", tmp_dir,
        "-S", sort_mem,
    ] + sort_key.split() + [
        "-o", output_path,
        input_path,
    ]
    env = os.environ.copy()
    env["LC_ALL"] = "C"
    subprocess.run(cmd, env=env, check=True)


def _dedup_sorted_barcodes(sorted_tsv, deduped_tsv):
    """Streaming single-pass dedup of barcode TSV sorted by rid asc, mapq desc.

    Input:  rid\\tbarcode\\tmapq  (sorted by -k1,1 -k3,3nr)
    Output: rid\\tbarcode         (one line per unique rid, highest mapq kept)

    Returns the number of unique read IDs.
    """
    n = 0
    prev_rid = None
    with open(sorted_tsv) as inp, open(deduped_tsv, "w") as out:
        for line in inp:
            rid, bc, _ = line.split("\t", 2)
            if rid != prev_rid:
                out.write(f"{rid}\t{bc}\n")
                prev_rid = rid
                n += 1
    return n


# ---------------------------------------------------------------------------
# Phase 3: linearise FASTQ via awk
# ---------------------------------------------------------------------------

# awk script: read 4-line FASTQ records, normalise read ID (strip @, strip
# trailing /1 /2 /3), and emit rid\tseq\tqual.  Print record count to stderr.
_AWK_LINEARIZE = (
    'NR%4==1{name=$1; sub(/^@/,"",name); '
    'if(name~/\\/[123]$/)name=substr(name,1,length(name)-2)} '
    'NR%4==2{seq=$0} '
    'NR%4==0{n++; print name"\\t"seq"\\t"$0} '
    'END{print n+0 > "/dev/stderr"}'
)


def _linearize_fastq(fastq_path, tsv_path):
    """Linearize FASTQ -> rid\\tseq\\tqual TSV using awk for speed.

    Falls back to a pure-Python implementation if awk is not available.
    Returns the number of records written.
    """
    if fastq_path.endswith(".gz"):
        decomp_bin = _PIGZ or "gzip"
        decomp = subprocess.Popen(
            [decomp_bin, "-dc", fastq_path],
            stdout=subprocess.PIPE,
        )
        awk_stdin = decomp.stdout
    else:
        decomp = None
        awk_stdin = open(fastq_path, "rb")

    with open(tsv_path, "w") as outf:
        awk = subprocess.Popen(
            ["awk", _AWK_LINEARIZE],
            stdin=awk_stdin,
            stdout=outf,
            stderr=subprocess.PIPE,
        )

    # Close parent's copy so child processes get SIGPIPE properly.
    if decomp:
        decomp.stdout.close()
    else:
        awk_stdin.close()

    _, awk_stderr = awk.communicate()
    if decomp:
        decomp.wait()

    if awk.returncode != 0:
        raise subprocess.CalledProcessError(awk.returncode, "awk")

    try:
        n = int(awk_stderr.decode().strip())
    except (ValueError, AttributeError):
        n = 0
    return n


# ---------------------------------------------------------------------------
# Phase 5: merge-join and write stamped FASTQ
# ---------------------------------------------------------------------------

def _merge_join_and_write(bc_tsv, fq_tsv, out_fastq, read_suffix,
                          sep=":", pigz_threads=4):
    """Walk sorted barcode + sorted FASTQ in lockstep, write stamped FASTQ on match.

    Both files must be sorted by rid (LC_ALL=C byte order).

    Returns the number of records written.
    """
    n_out = 0
    with open(bc_tsv) as bc_f, \
         open(fq_tsv) as fq_f, \
         open_fastq_writer(out_fastq, threads=pigz_threads) as w:
        bc_line = bc_f.readline()
        fq_line = fq_f.readline()

        while bc_line and fq_line:
            bc_rid, bc_val = bc_line.rstrip("\n").split("\t", 1)
            parts = fq_line.rstrip("\n").split("\t", 2)
            fq_rid = parts[0]
            fq_seq = parts[1]
            fq_qual = parts[2] if len(parts) > 2 else ""

            if bc_rid == fq_rid:
                new_name = f"{fq_rid}{sep}{bc_val}{read_suffix}"
                w.write(f"@{new_name}\n{fq_seq}\n+\n{fq_qual}\n")
                n_out += 1
                # Advance FASTQ; barcode stays in case of duplicate FASTQ rids
                fq_line = fq_f.readline()
            elif bc_rid < fq_rid:
                bc_line = bc_f.readline()
            else:
                fq_line = fq_f.readline()

    return n_out


# ---------------------------------------------------------------------------
# Orchestration for one FASTQ
# ---------------------------------------------------------------------------

def _stamp_one_fastq(fastq_in, fastq_out, bc_tsv, read_suffix, tmp_dir,
                     sep=":", threads=None, sort_mem="2G"):
    """Orchestrate linearise -> sort -> merge-join for a single FASTQ file.

    Returns (n_in, n_out).
    """
    if threads is None:
        threads = os.cpu_count()
    pigz_threads = max(1, threads // 2)

    suffix = read_suffix.strip("/")
    linear_tsv = os.path.join(tmp_dir, f"linearized_{suffix}.tsv")
    sorted_fq_tsv = os.path.join(tmp_dir, f"sorted_fq_{suffix}.tsv")

    # Phase 3: linearize FASTQ
    n_in = _linearize_fastq(fastq_in, linear_tsv)

    # Phase 4: sort linearized FASTQ by core_name
    _run_sort(linear_tsv, sorted_fq_tsv, sort_key="-k1,1",
              tmp_dir=tmp_dir, parallel=threads, sort_mem=sort_mem)
    os.unlink(linear_tsv)

    # Phase 5: merge-join and write
    n_out = _merge_join_and_write(
        bc_tsv, sorted_fq_tsv, fastq_out,
        read_suffix=read_suffix, sep=sep, pigz_threads=pigz_threads,
    )
    os.unlink(sorted_fq_tsv)

    return n_in, n_out


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------

def stamp(
    bam,
    r1,
    r3,
    out_r1,
    out_r3,
    bc_tag="BC",
    min_mapq=0,
    primary_only=True,
    threads=None,
    sep=":",
    tmp_dir=None,
):
    """Full stamping workflow using streaming sort-merge join.

    Keeps Python memory at O(1) by delegating sorting to GNU sort.
    R1 and R3 are processed in parallel to reduce wall-clock time.

    Peak memory: ~4 GB (two GNU sort processes each using -S 2G), plus
    minor overhead from pigz and awk.
    """
    if threads is None:
        threads = os.cpu_count()

    with tempfile.TemporaryDirectory(dir=tmp_dir, prefix="rupture_stamp_") as td:
        raw_bc_tsv = os.path.join(td, "raw_barcodes.tsv")
        sorted_bc_tsv = os.path.join(td, "sorted_barcodes.tsv")
        deduped_bc_tsv = os.path.join(td, "barcodes.tsv")

        # Phase 1: extract barcodes from BAM
        n_extracted = _extract_barcodes_to_tsv(
            bam, raw_bc_tsv,
            bc_tag=bc_tag, min_mapq=min_mapq,
            primary_only=primary_only, threads=threads,
        )
        print(
            f"[info] extracted {n_extracted:,} barcode records from BAM.",
            file=sys.stderr,
        )

        # Phase 2: sort by rid asc, mapq desc; then dedup
        _run_sort(
            raw_bc_tsv, sorted_bc_tsv,
            sort_key="-k1,1 -k3,3nr", tmp_dir=td, parallel=threads,
        )
        os.unlink(raw_bc_tsv)

        n_unique = _dedup_sorted_barcodes(sorted_bc_tsv, deduped_bc_tsv)
        os.unlink(sorted_bc_tsv)
        print(
            f"[info] collected barcodes for {n_unique:,} read IDs.",
            file=sys.stderr,
        )

        # Phases 3-5 for R1 and R3 in parallel
        per_task_threads = max(1, threads // 2)

        with ThreadPoolExecutor(max_workers=2) as pool:
            fut_r1 = pool.submit(
                _stamp_one_fastq,
                r1, out_r1, deduped_bc_tsv,
                read_suffix="/1", tmp_dir=td, sep=sep,
                threads=per_task_threads,
            )
            fut_r3 = pool.submit(
                _stamp_one_fastq,
                r3, out_r3, deduped_bc_tsv,
                read_suffix="/2", tmp_dir=td, sep=sep,
                threads=per_task_threads,
            )

            n1_in, n1_out = fut_r1.result()
            n3_in, n3_out = fut_r3.result()

        print(f"[R1] {n1_out:,}/{n1_in:,} kept.", file=sys.stderr)
        print(f"[R3] {n3_out:,}/{n3_in:,} kept.", file=sys.stderr)
