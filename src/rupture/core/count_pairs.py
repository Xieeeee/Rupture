"""Count contacts per cell from a .pairs file."""

import gzip

import numpy as np
import pandas as pd


def open_pairs_file(filename, mode="r"):
    """Open a .pairs or .pairs.gz file."""
    if filename.endswith(".gz"):
        return gzip.open(filename, mode)
    return open(filename, mode)


def count_contacts(input_path, output_prefix):
    """Count per-cell contact types from a .pairs file.

    Parameters
    ----------
    input_path : str
        Path to .pairs or .pairs.gz file.
    output_prefix : str
        Output prefix; writes {output_prefix}.stat.csv.
    """
    cid = {}
    with open_pairs_file(input_path, "r") as infile:
        for line in infile:
            try:
                dline = line.decode("utf8")
            except AttributeError:
                dline = line
            if dline[:1] == "#":
                continue
            fields = dline.strip("\n").split("\t")[:10]
            (
                readname,
                chrom1,
                pos1,
                chrom2,
                pos2,
                strand1,
                strand2,
                pair_type,
                cell1,
                cell2,
            ) = fields
            if cell1 in cid:
                cid[cell1]["total"] += 1
            else:
                cid[cell1] = {
                    "total": 1,
                    "mapped": 0,
                    "unmapped": 0,
                    "duplicate": 0,
                    "cis": 0,
                    "cis_1kb-": 0,
                    "cis_1kb+": 0,
                    "cis_10kb+": 0,
                    "trans": 0,
                }
            if pair_type.upper() in ["NN", "XX"]:
                cid[cell1]["unmapped"] += 1
            elif chrom1 != "!" and chrom2 != "!":
                cid[cell1]["mapped"] += 1
                if pair_type.upper() == "DD":
                    cid[cell1]["duplicate"] += 1
                elif chrom1 == chrom2:
                    cid[cell1]["cis"] += 1
                    dist = np.abs(int(pos2) - int(pos1))
                    if dist < 1000:
                        cid[cell1]["cis_1kb-"] += 1
                    if dist > 1000:
                        cid[cell1]["cis_1kb+"] += 1
                    if dist > 10000:
                        cid[cell1]["cis_10kb+"] += 1
                else:
                    cid[cell1]["trans"] += 1

    cidd = pd.DataFrame.from_dict(cid)
    cidd.T.to_csv(output_prefix + ".stat.csv", sep="\t")
