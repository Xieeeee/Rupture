"""Barcode stamping: filter R1/R3 by mapped R2 BAM and append barcode to headers."""

import gzip
import os
import sys

import pysam


def core_name(n: str) -> str:
    """Normalize read IDs: strip whitespace; remove trailing /1, /2, /3."""
    n = n.split()[0]
    if len(n) >= 2 and n[-2] == "/" and n[-1] in "123":
        n = n[:-2]
    return n


def open_fastq_writer(path: str):
    """Open a FASTQ writer, gzip-compressed if path ends with .gz."""
    if path.endswith(".gz"):
        return gzip.open(path, "wt", compresslevel=5)
    return open(path, "w")


def build_id_to_barcode(
    bam_path,
    bc_tag="BC",
    min_mapq=0,
    primary_only=True,
    threads=None,
):
    """Build a mapping of read ID -> barcode from a BAM file."""
    if threads is None:
        threads = os.cpu_count()
    id2bc = {}
    id2mq = {}
    mode = "rb" if bam_path.endswith((".bam", ".cram")) else "r"
    with pysam.AlignmentFile(bam_path, mode) as aln:
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

            # prefer tag (e.g., BC); fallback to reference name
            if bc_tag and rec.has_tag(bc_tag):
                bc = str(rec.get_tag(bc_tag))
            else:
                bc = rec.reference_name or "NA"

            # keep best MAPQ per read id
            mq = int(rec.mapping_quality)
            if rid not in id2mq or mq > id2mq[rid]:
                id2mq[rid] = mq
                id2bc[rid] = bc
    return id2bc


def stamp_and_filter_fastq(in_fastq, out_fastq, id2bc, read_suffix, sep=":"):
    """Filter FASTQ by barcode map and stamp barcode into read name.

    Parameters
    ----------
    sep : str
        Separator between read name and barcode (default ':').
    """
    n_in = n_out = 0
    with pysam.FastxFile(in_fastq) as fq, open_fastq_writer(out_fastq) as w:
        for rec in fq:
            n_in += 1
            rid = core_name(rec.name)
            bc = id2bc.get(rid)
            if bc is None:
                continue
            new_name = f"{rid}{sep}{bc}{read_suffix}"
            w.write(f"@{new_name}\n{rec.sequence}\n+\n{rec.quality or ''}\n")
            n_out += 1
    return n_in, n_out


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
):
    """Full stamping workflow: build barcode map then filter+stamp R1 and R3."""
    if threads is None:
        threads = os.cpu_count()

    id2bc = build_id_to_barcode(
        bam,
        bc_tag=bc_tag,
        min_mapq=min_mapq,
        primary_only=primary_only,
        threads=threads,
    )
    print(f"[info] collected barcodes for {len(id2bc):,} read IDs.", file=sys.stderr)

    n1_in, n1_out = stamp_and_filter_fastq(r1, out_r1, id2bc, read_suffix="/1", sep=sep)
    print(f"[R1] {n1_out:,}/{n1_in:,} kept.", file=sys.stderr)

    n3_in, n3_out = stamp_and_filter_fastq(r3, out_r3, id2bc, read_suffix="/2", sep=sep)
    print(f"[R3] {n3_out:,}/{n3_in:,} kept.", file=sys.stderr)
