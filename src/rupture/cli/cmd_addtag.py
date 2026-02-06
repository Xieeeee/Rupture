"""rupture add-tag — extract field from read name into a BAM tag."""

import os


def add_parser_addtag(subparsers):
    p = subparsers.add_parser(
        "add-tag",
        help="Add n-th field from read name to a BAM tag.",
    )
    p.add_argument("--input", required=True, help="Input BAM file.")
    p.add_argument("--output", required=True, help="Output BAM file.")
    p.add_argument(
        "--field",
        type=int,
        required=True,
        help="0-based field index in read name (negative for indexing from end).",
    )
    p.add_argument(
        "--tag", default="CB", help="SAM tag to write (default: CB)."
    )
    p.add_argument(
        "--sep", default=":", help="Delimiter in read name (default: ':')."
    )
    p.add_argument(
        "--threads",
        type=int,
        default=os.cpu_count(),
        help="BGZF threads for BAM read/write.",
    )
    p.set_defaults(func=addtag_cmd)


def addtag_cmd(args):
    from rupture.core.addtag import add_tag_to_bam

    add_tag_to_bam(
        input_bam=args.input,
        output_bam=args.output,
        field=args.field,
        tag=args.tag,
        sep=args.sep,
        threads=args.threads,
    )
