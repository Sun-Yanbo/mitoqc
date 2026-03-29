"""
MitoQC command-line interface.

Usage
-----
  # Run on a local GenBank flat file
  python -m mitoqc --input sequences.gb --output results.csv

  # Run on NCBI accessions (fetches via Entrez)
  python -m mitoqc --accessions NC_012920 NC_001643 --output results.csv

  # Run with BLAST verification enabled (slow; requires network)
  python -m mitoqc --input sequences.gb --output results.csv --blast

  # Specify NCBI API key for higher rate limits
  python -m mitoqc --input sequences.gb --output results.csv \\
                   --blast --api-key YOUR_NCBI_KEY

  # Use M2 v2 (ORF scan) and M4 v2 (BLAST interface)
  python -m mitoqc --input sequences.gb --output results.csv --v2-modules
"""

from __future__ import annotations

import argparse
import sys
import pandas as pd
from pathlib import Path

from mitoqc.utils.io import iter_genbank_file, fetch_genbank_records
from mitoqc.pipeline import run_pipeline


def _build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        prog="mitoqc",
        description="MitoQC: Systematic quality control for metazoan mitogenomes.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__,
    )
    src = p.add_mutually_exclusive_group(required=True)
    src.add_argument(
        "--input", "-i", metavar="FILE",
        help="Input GenBank flat file (.gb or .gb.gz).",
    )
    src.add_argument(
        "--accessions", "-a", nargs="+", metavar="ACC",
        help="One or more NCBI accession numbers or UIDs.",
    )
    p.add_argument(
        "--output", "-o", metavar="FILE", required=True,
        help="Output CSV file path.",
    )
    p.add_argument(
        "--email", default="mitoqc@example.com",
        help="Email address for NCBI Entrez (required by NCBI policy).",
    )
    p.add_argument(
        "--blast", action="store_true",
        help="Enable BLAST species verification in M4 (slow; requires network).",
    )
    p.add_argument(
        "--api-key", metavar="KEY",
        help="NCBI API key for higher rate limits (10 req/s vs 3 req/s).",
    )
    p.add_argument(
        "--v2-modules", action="store_true",
        help="Use M2 v2 (ORF scan) and M4 v2 (BLAST interface) instead of v1.",
    )
    p.add_argument(
        "--batch-size", type=int, default=50,
        help="Records per Entrez fetch batch (default: 50).",
    )
    p.add_argument(
        "--verbose", "-v", action="store_true",
        help="Print progress to stderr.",
    )
    return p


def main(argv: list[str] | None = None) -> int:
    args = _build_parser().parse_args(argv)

    # Build config
    config: dict = {}
    if args.blast:
        config["m4"] = {"run_blast": True}
        if args.api_key:
            config["m4"]["api_key"] = args.api_key

    # Optionally swap in v2 modules
    if args.v2_modules:
        from mitoqc.modules import m2_annotation_v2, m4_species_v2
        from mitoqc import pipeline as _pl
        # Monkey-patch the module registry
        _pl._MODULES = [
            ("m1", __import__("mitoqc.modules.m1_integrity", fromlist=["m1_integrity"])),
            ("m2", m2_annotation_v2),
            ("m3", __import__("mitoqc.modules.m3_chimera",   fromlist=["m3_chimera"])),
            ("m4", m4_species_v2),
            ("m5", __import__("mitoqc.modules.m5_pcg",       fromlist=["m5_pcg"])),
        ]

    # Source records
    if args.input:
        records = iter_genbank_file(args.input)
    else:
        records = fetch_genbank_records(
            args.accessions,
            email=args.email,
            batch_size=args.batch_size,
        )

    # Process
    rows = []
    for i, rec in enumerate(records):
        row = run_pipeline(rec, config=config)
        rows.append(row)
        if args.verbose and (i + 1) % 100 == 0:
            print(f"  Processed {i+1} records...", file=sys.stderr)

    if not rows:
        print("ERROR: No records processed.", file=sys.stderr)
        return 1

    df = pd.DataFrame(rows)
    df.to_csv(args.output, index=False)

    if args.verbose:
        print(f"\nDone. {len(df)} records → {args.output}", file=sys.stderr)
        print(df["quality_grade"].value_counts().to_string(), file=sys.stderr)

    return 0


if __name__ == "__main__":
    sys.exit(main())
