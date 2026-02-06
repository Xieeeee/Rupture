"""rupture count-pairs — per-cell pair counts from .pairs file."""


def add_parser_count(subparsers):
    p = subparsers.add_parser(
        "count-pairs",
        help="Count contacts per cell from a .pairs file.",
    )
    p.add_argument("--input", required=True, help=".pairs or .pairs.gz input.")
    p.add_argument("--output", required=True, help="Output prefix (writes {output}.stat.csv).")
    p.set_defaults(func=count_cmd)


def count_cmd(args):
    from rupture.core.count_pairs import count_contacts

    count_contacts(args.input, args.output)
