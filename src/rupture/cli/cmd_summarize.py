"""rupture summarize-qc — QC summary from pairtools stats."""


def add_parser_summarize(subparsers):
    p = subparsers.add_parser(
        "summarize-qc",
        help="Compute QC summary from pairtools dedup stats.",
    )
    p.add_argument(
        "--input", required=True, help="Pairtools stats file (e.g. .pairdedup.txt)."
    )
    p.add_argument(
        "--output",
        default=None,
        help="Output summary file (default: input with .summary.txt suffix).",
    )
    p.set_defaults(func=summarize_cmd)


def summarize_cmd(args):
    from rupture.core.summarize_qc import summarize_qc

    summarize_qc(args.input, args.output)
