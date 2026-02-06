"""Merge FASTQ files across lanes by read type for a given library prefix."""

import re
import shutil
import sys
from pathlib import Path


def detect_read_types(input_dir, prefix):
    """Auto-detect read types (I1, I2, R1, R2, R3, etc.) from FASTQ files.

    Expects filenames matching {prefix}_*_{ReadType}_*.fastq.gz
    """
    input_dir = Path(input_dir)
    pattern = re.compile(rf"^{re.escape(prefix)}_.*_([IR]\d)_\d+\.fastq\.gz$")
    read_types = set()
    for f in input_dir.iterdir():
        m = pattern.match(f.name)
        if m:
            read_types.add(m.group(1))
    return sorted(read_types)


def merge_lanes(prefix, input_dir, output_dir):
    """Merge FASTQ files across lanes for each detected read type.

    Parameters
    ----------
    prefix : str
        Library prefix (e.g. 'CL010').
    input_dir : str or Path
        Directory containing per-lane FASTQ files.
    output_dir : str or Path
        Directory for merged output files.
    """
    input_dir = Path(input_dir)
    output_dir = Path(output_dir)

    if not input_dir.is_dir():
        raise FileNotFoundError(f"Input directory does not exist: {input_dir}")

    output_dir.mkdir(parents=True, exist_ok=True)

    read_types = detect_read_types(input_dir, prefix)
    if not read_types:
        raise FileNotFoundError(
            f"No fastq.gz files found matching {prefix}_* in {input_dir}"
        )

    print(f"Library: {prefix}", file=sys.stderr)
    print(f"Input:   {input_dir}", file=sys.stderr)
    print(f"Output:  {output_dir}", file=sys.stderr)
    print(f"Read types detected: {' '.join(read_types)}", file=sys.stderr)

    for rt in read_types:
        outfile = output_dir / f"{prefix}_{rt}.fq.gz"
        files = sorted(input_dir.glob(f"{prefix}_*_{rt}_*.fastq.gz"))
        print(
            f"Merging {len(files)} files for {rt} -> {outfile.name}",
            file=sys.stderr,
        )
        with open(outfile, "wb") as out_fh:
            for f in files:
                with open(f, "rb") as in_fh:
                    shutil.copyfileobj(in_fh, out_fh)
