"""Rupture CLI entry point — dispatches subcommands."""

import argparse
import sys

from rupture import __version__


def main(argv=None):
    parser = argparse.ArgumentParser(
        prog="rupture",
        description="Single-cell Hi-C pipeline for droplet-based Hi-C data.",
    )
    parser.add_argument(
        "--version", action="version", version=f"%(prog)s {__version__}"
    )

    subparsers = parser.add_subparsers(dest="command")

    # Import and register each subcommand
    from rupture.cli.cmd_run import add_parser_run
    from rupture.cli.cmd_init import add_parser_init
    from rupture.cli.cmd_stamp import add_parser_stamp
    from rupture.cli.cmd_addtag import add_parser_addtag
    from rupture.cli.cmd_count import add_parser_count
    from rupture.cli.cmd_summarize import add_parser_summarize
    from rupture.cli.cmd_merge_fastq import add_parser_merge_fastq

    add_parser_run(subparsers)
    add_parser_init(subparsers)
    add_parser_stamp(subparsers)
    add_parser_addtag(subparsers)
    add_parser_count(subparsers)
    add_parser_summarize(subparsers)
    add_parser_merge_fastq(subparsers)

    args = parser.parse_args(argv)

    if args.command is None:
        parser.print_help()
        sys.exit(1)

    args.func(args)


if __name__ == "__main__":
    main()
