"""rupture merge-fastq — merge FASTQ files across lanes."""


def add_parser_merge_fastq(subparsers):
    p = subparsers.add_parser(
        "merge-fastq",
        help="Merge FASTQ files across lanes by read type.",
    )
    p.add_argument("--prefix", required=True, help="Library prefix (e.g. CL010).")
    p.add_argument("--input-dir", required=True, help="Directory with per-lane FASTQs.")
    p.add_argument("--output-dir", required=True, help="Directory for merged output.")
    p.set_defaults(func=merge_cmd)


def merge_cmd(args):
    from rupture.core.merge_fastq import merge_lanes

    merge_lanes(args.prefix, args.input_dir, args.output_dir)
