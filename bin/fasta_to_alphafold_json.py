#!/usr/bin/env python3
#########################################################################
# File Name: fasta_to_alphafold_json.py
# Author: ChengYu
# Description: Convert FASTA protein sequences to AlphaFold server
#              JSON input format.
# Created Time: 2026
#
# Changelog:
#   v1.1.0  2026-09-25
#   - FIX: sequences are validated before writing — empty entries and
#     entries with non-amino-acid characters are skipped with a warning
#     (previously silent; the JSON only failed at submission time).
#     Sequences are uppercased; allowed letters: 20 standard AA + X.
#   - FIX: --chain-count below 1 fails with a clean error.
#   - DOC: one JSON job is produced PER FASTA record (single-chain jobs);
#     multi-chain complexes cannot be expressed from a multi-record
#     FASTA — noted in the help.
#########################################################################
"""Convert FASTA protein sequences to AlphaFold server JSON format.

Generates JSON input compatible with the AlphaFold server API. Each
FASTA record becomes one single-chain JSON job (multi-chain complexes
are not expressible from a multi-record FASTA).
"""

import argparse
import json
import logging
import re
import sys
from pathlib import Path

__version__ = "1.1.0"

# Standard amino acids + X (unspecified); anything else makes the entry invalid
VALID_AA = set("ACDEFGHIKLMNPQRSTVWYX")


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

    with open(file_path, "r", encoding="utf-8") as fh:
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


def validate_entries(entries: list) -> list:
    """Filter entries: drop empty sequences and non-amino-acid characters.

    Sequences are uppercased. Entries that are empty or contain characters
    outside the standard amino acids (plus X) are dropped with a warning.

    Returns the cleaned entry list.
    """
    cleaned = []
    for entry in entries:
        seq = entry["sequence"].upper()
        if not seq:
            logging.warning("Skipping '%s': empty sequence.", entry["header"])
            continue
        bad = sorted(set(seq) - VALID_AA)
        if bad:
            logging.warning(
                "Skipping '%s': non-amino-acid character(s): %s",
                entry["header"], ", ".join(bad),
            )
            continue
        cleaned.append({"header": entry["header"], "sequence": seq})
    return cleaned


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

notes:
  Each FASTA record becomes one single-chain JSON job. Sequences must
  contain only standard amino-acid letters (plus X) and are uppercased;
  invalid or empty records are skipped with a warning.
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

    if args.chain_count < 1:
        parser.error("--chain-count must be >= 1")

    entries = parse_fasta(args.input)
    if not entries:
        logging.error("No sequences found in %s", args.input)
        sys.exit(1)

    logging.info("Parsed %d sequences from %s", len(entries), args.input)

    entries = validate_entries(entries)
    if not entries:
        logging.error("No valid sequences remained after validation.")
        sys.exit(1)

    json_output = build_alphafold_json(entries, args.chain_count)

    with open(args.output, "w") as fh:
        json.dump(json_output, fh, indent=4, ensure_ascii=False)

    logging.info("Wrote %d entries -> %s", len(json_output), args.output)


if __name__ == "__main__":
    main()
