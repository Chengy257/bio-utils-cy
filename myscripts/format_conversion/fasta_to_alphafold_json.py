#!/usr/bin/env python3
#########################################################################
# File Name: fasta_to_alphafold_json.py
# Author: ChengYu
# Description: Convert FASTA protein sequences to AlphaFold server
#              JSON input format.
# Created Time: 2026
#########################################################################
"""Convert FASTA protein sequences to AlphaFold server JSON format.

Generates JSON input compatible with the AlphaFold server API, with
each FASTA entry becoming a separate protein chain entry.
"""

import argparse
import json
import logging
import sys
from pathlib import Path

__version__ = "1.0.0"


def parse_fasta(file_path: str) -> list:
    """Parse a FASTA file preserving full headers.

    Args:
        file_path: Path to FASTA file.

    Returns:
        List of dicts with 'header' and 'sequence' keys.
    """
    entries = []
    header = None
    seq_parts = []

    with open(file_path, "r") as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if header is not None:
                    entries.append({"header": header, "sequence": "".join(seq_parts)})
                header = line[1:].strip()
                seq_parts = []
            else:
                seq_parts.append(line)
        if header is not None:
            entries.append({"header": header, "sequence": "".join(seq_parts)})

    return entries


def build_alphafold_json(entries: list, chain_count: int = 1) -> list:
    """Build AlphaFold server JSON structure.

    Args:
        entries: List of {'header': str, 'sequence': str} dicts.
        chain_count: Number of copies per chain.

    Returns:
        List of AlphaFold-formatted JSON dicts.
    """
    json_output = []
    for entry in entries:
        json_entry = {
            "name": entry["header"],
            "modelSeeds": [],
            "sequences": [
                {
                    "proteinChain": {
                        "sequence": entry["sequence"],
                        "glycans": [],
                        "modifications": [],
                        "count": chain_count,
                    }
                }
            ],
        }
        json_output.append(json_entry)
    return json_output


def build_parser() -> argparse.ArgumentParser:
    """Build argument parser."""
    parser = argparse.ArgumentParser(
        description="Convert FASTA protein sequences to AlphaFold server JSON format.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""\
examples:
  python fasta_to_alphafold_json.py -i proteins.fa -o alphafold_input.json
  python fasta_to_alphafold_json.py -i proteins.fa -o input.json --chain-count 2
""",
    )
    parser.add_argument("-i", "--input", type=str, required=True, help="Input protein FASTA file.")
    parser.add_argument("-o", "--output", type=str, required=True, help="Output JSON file.")
    parser.add_argument(
        "--chain-count", type=int, default=1,
        help="Number of copies per protein chain (default: 1).",
    )
    parser.add_argument(
        "--log-level", type=str, default="INFO",
        choices=["DEBUG", "INFO", "WARNING", "ERROR"],
        help="Logging level (default: INFO).",
    )
    parser.add_argument("--version", action="version", version=f"%(prog)s {__version__}")
    return parser


def main() -> None:
    """Entry point."""
    parser = build_parser()
    args = parser.parse_args()

    logging.basicConfig(
        format="%(asctime)s - %(levelname)s - %(message)s",
        level=getattr(logging, args.log_level),
    )

    if not Path(args.input).is_file():
        logging.error("Input file not found: %s", args.input)
        sys.exit(1)

    entries = parse_fasta(args.input)
    if not entries:
        logging.error("No sequences found in %s", args.input)
        sys.exit(1)

    logging.info("Parsed %d sequences from %s", len(entries), args.input)

    json_output = build_alphafold_json(entries, args.chain_count)

    with open(args.output, "w") as fh:
        json.dump(json_output, fh, indent=4, ensure_ascii=False)

    logging.info("Wrote %d entries -> %s", len(json_output), args.output)


if __name__ == "__main__":
    main()
