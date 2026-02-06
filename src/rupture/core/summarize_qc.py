"""QC summary from pairtools dedup stats — Python rewrite of phc.summarize_pairs_lec.R."""

import csv
import os


def parse_pairtools_stats(stats_path):
    """Parse a pairtools stats file into a dict of {key: value}.

    The file is tab-separated with two columns: metric name and count.
    Lines starting with '#' are skipped.
    """
    data = {}
    with open(stats_path) as fh:
        for line in fh:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split("\t", 1)
            if len(parts) == 2:
                try:
                    data[parts[0]] = int(parts[1])
                except ValueError:
                    pass
    return data


def compute_summary(data):
    """Compute 9 derived QC metrics + ratios from pairtools stats.

    Returns list of (metric_name, count, ratio) tuples.
    """
    total_unmapped = data.get("total_unmapped", 0)
    total_single = data.get("total_single_sided_mapped", 0)
    total_dups = data.get("total_dups", 0)
    total_nodups = data.get("total_nodups", 0)
    cis = data.get("cis", 0)
    cis_1kb_plus = data.get("cis_1kb+", 0)
    trans = data.get("trans", 0)

    sequenced = total_unmapped + total_single + total_dups + total_nodups
    mapped = total_dups + total_nodups
    nonclonal = total_nodups
    intra = cis
    short_range = cis - cis_1kb_plus
    long_range = cis_1kb_plus
    inter = trans
    unmapped = total_unmapped + total_single
    duplicates = total_dups

    def safe_div(a, b):
        return a / b if b else 0.0

    metrics = [
        # (name, count, ratio)
        ("Sequenced_Read_Pairs", sequenced, safe_div(sequenced, sequenced)),
        ("Mapped_Read_Pairs", mapped, safe_div(mapped, sequenced)),
        ("Nonclonal_Read_Pairs", nonclonal, safe_div(nonclonal, mapped)),
        ("Intra_chromosomal", intra, safe_div(intra, nonclonal)),
        ("short_range_1k", short_range, safe_div(short_range, intra)),
        ("long_range_1k", long_range, safe_div(long_range, intra)),
        ("Inter_chromosomal", inter, safe_div(inter, nonclonal)),
        ("Unmapped", unmapped, safe_div(unmapped, sequenced)),
        ("Removed_Duplicates", duplicates, safe_div(duplicates, mapped)),
    ]
    return metrics


def summarize_qc(input_path, output_path=None):
    """Read pairtools stats, compute summary, write tab-separated output.

    Parameters
    ----------
    input_path : str
        Path to pairtools dedup stats file (e.g. *.pairdedup.txt).
    output_path : str or None
        Output path. If None, derives from input_path by replacing
        the last .txt with .summary.txt.
    """
    if output_path is None:
        base, _ = os.path.splitext(input_path)
        output_path = base + ".summary.txt"

    data = parse_pairtools_stats(input_path)
    metrics = compute_summary(data)

    with open(output_path, "w", newline="") as fh:
        writer = csv.writer(fh, delimiter="\t")
        for name, count, ratio in metrics:
            writer.writerow([name, count, f"{ratio:.6f}"])
