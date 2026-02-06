"""Add n-th field from read name to a BAM tag."""

import os

import pysam


def nth_field(name: str, sep: str, idx: int) -> str:
    """Extract the n-th field from a string split by sep.

    Supports negative indexing (split fully, index from end).
    For positive indices, only splits up to the needed part for speed.
    """
    if idx < 0:
        parts = name.split(sep)
        return parts[idx] if len(parts) >= abs(idx) else ""
    parts = name.split(sep, idx + 1)
    return parts[idx] if len(parts) > idx else ""


def add_tag_to_bam(input_bam, output_bam, field, tag, sep=":", threads=None):
    """Read input BAM, extract field from read name, write as SAM tag.

    Parameters
    ----------
    input_bam : str
        Path to input BAM file.
    output_bam : str
        Path to output BAM file.
    field : int
        0-based field index (or negative for indexing from end).
    tag : str
        SAM tag to write (e.g. 'CB').
    sep : str
        Delimiter in read name (default ':').
    threads : int or None
        BGZF threads for read/write.
    """
    if threads is None:
        threads = os.cpu_count()
    with pysam.AlignmentFile(input_bam, "rb") as inb:
        try:
            inb.set_threads(max(1, threads))
        except Exception:
            pass
        with pysam.AlignmentFile(output_bam, "wb", template=inb) as outb:
            try:
                outb.set_threads(max(1, threads))
            except Exception:
                pass
            for read in inb.fetch(until_eof=True):
                val = nth_field(read.query_name, sep, field)
                read.set_tag(tag, val, value_type="Z")
                outb.write(read)
