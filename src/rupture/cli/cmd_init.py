"""rupture init — generate config template in a working directory."""

import shutil
import sys
from pathlib import Path


def add_parser_init(subparsers):
    p = subparsers.add_parser(
        "init", help="Initialize a working directory with config template."
    )
    p.add_argument(
        "-d",
        "--dir",
        default=".",
        help="Target directory (default: current directory).",
    )
    p.add_argument(
        "-f", "--force", action="store_true", help="Overwrite existing config.yaml."
    )
    p.set_defaults(func=init)


def init(args):
    from rupture.workflow import get_config_template_path

    target = Path(args.dir)
    target.mkdir(parents=True, exist_ok=True)

    config_dest = target / "config.yaml"
    if config_dest.exists() and not args.force:
        print(
            f"Error: {config_dest} already exists. Use --force to overwrite.",
            file=sys.stderr,
        )
        sys.exit(1)

    template = get_config_template_path()
    shutil.copy2(template, config_dest)
    print(f"Created {config_dest}", file=sys.stderr)

    rawdata = target / "01.rawdata"
    rawdata.mkdir(exist_ok=True)
    print(f"Created {rawdata}/", file=sys.stderr)
