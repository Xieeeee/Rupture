"""rupture stamp — barcode stamping onto FASTQ headers."""

import os


def add_parser_stamp(subparsers):
    p = subparsers.add_parser(
        "stamp",
        help="Filter R1/R3 by mapped R2 BAM and stamp barcodes into headers.",
    )
    p.add_argument("--bam", required=True, help="Mapped R2 BAM/SAM (barcode index).")
    p.add_argument("--r1", required=True, help="Input R1 FASTQ(.gz).")
    p.add_argument("--r3", required=True, help="Input R3 FASTQ(.gz).")
    p.add_argument("--out-r1", required=True, help="Output stamped R1 FASTQ(.gz).")
    p.add_argument("--out-r3", required=True, help="Output stamped R3 FASTQ(.gz).")
    p.add_argument(
        "--bc-tag",
        default="BC",
        help="SAM tag carrying barcode (default: BC).",
    )
    p.add_argument(
        "--min-mapq", type=int, default=0, help="Minimum MAPQ to keep (default: 0)."
    )
    p.add_argument(
        "--sep",
        default=":",
        help="Separator between read name and barcode (default: ':').",
    )
    p.add_argument(
        "--threads",
        type=int,
        default=os.cpu_count(),
        help="BGZF threads for BAM reading.",
    )
    p.add_argument(
        "--tmp-dir",
        default=None,
        help="Directory for temporary sort files (default: system tmpdir).",
    )
    p.set_defaults(func=stamp_cmd)


def stamp_cmd(args):
    from rupture.core.stamp import stamp

    stamp(
        bam=args.bam,
        r1=args.r1,
        r3=args.r3,
        out_r1=args.out_r1,
        out_r3=args.out_r3,
        bc_tag=args.bc_tag,
        min_mapq=args.min_mapq,
        threads=args.threads,
        sep=args.sep,
        tmp_dir=args.tmp_dir,
    )
