"""rupture run — invoke the bundled Snakefile via snakemake."""

import argparse
import subprocess
import sys


def add_parser_run(subparsers):
    p = subparsers.add_parser("run", help="Run the Hi-C pipeline via Snakemake.")
    p.add_argument(
        "-c", "--configfile", required=True, help="Path to config.yaml."
    )
    p.add_argument(
        "-j", "--cores", type=int, default=16, help="Number of cores (default: 16)."
    )
    p.add_argument(
        "-n", "--dryrun", action="store_true", help="Dry-run (do not execute)."
    )
    p.add_argument(
        "--snakemake-args",
        nargs=argparse.REMAINDER,
        default=[],
        help="Additional arguments passed directly to snakemake.",
    )
    p.set_defaults(func=run)


def run(args):
    from rupture.workflow import get_snakefile_path

    snakefile = get_snakefile_path()

    cmd = [
        "snakemake",
        "-s",
        str(snakefile),
        "--configfile",
        args.configfile,
        "-j",
        str(args.cores),
    ]
    if args.dryrun:
        cmd.append("-n")
    cmd.extend(args.snakemake_args)

    sys.exit(subprocess.call(cmd))
