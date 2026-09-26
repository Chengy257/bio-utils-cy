#!/usr/bin/env python3
#########################################################################
# File Name: maf_extract_regions.py
# Author: ChengYu
# Description: Extract genomic regions from MAF alignment files based
#              on BED12 coordinates and translate to protein.
# Created Time: 2026
#
# Changelog:
#   v1.1.0  2026-09-26
#   - REWRITE: hand-rolled MAF parser (no bx-python). Component
#     coordinates now come from the s-line start/size fields; the old
#     code parsed a ":start-end" suffix out of the src string, which
#     standard MAF files do not carry, so every component was silently
#     skipped.
#   - FIX: reference gap columns misaligned extraction. Bases are now
#     mapped by walking the component text (gaps do not advance the
#     genomic position) instead of slicing the text with genomic
#     offsets.
#   - FIX: regions spanning multiple alignment blocks were truncated to
#     the first matching block; contributions from all overlapping
#     blocks are now concatenated in BED12 block order. Genomic
#     positions not covered by any block are filled with 'N'.
#   - FIX: reference component strand ('-' in the s line) is now
#     honored: the text is reverse-complemented to reconstruct the
#     plus-strand genome sequence.
#   - FIX: component matching is an exact species/chrom match (the old
#     src prefix match collided, e.g. chr1 matching chr10).
#   - Duplicate BED names, malformed lines and blockCount mismatches
#     are warned about; exit 1 when no region is found in the MAF.
#########################################################################
"""Extract genomic regions from MAF alignment files using BED12 coordinates.

For each region in a BED12 file, extracts the reference species' genomic
sequence from a MAF file, optionally reverse-complementing for negative
strand. Outputs both nucleotide and translated protein FASTA.

Positions of a region not covered by any alignment block are filled
with 'N'; regions never covered by the MAF are omitted from the output
(and reported). The MAF is read twice-free: a single streaming pass.

Requires: biopython (for translation)
"""

import argparse
import logging
import sys
from pathlib import Path
from typing import Dict, List, Optional, Tuple

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

__version__ = "1.1.0"

_COMP_TRANS = str.maketrans("ACGTNacgtn-", "TGCANtgcan-")


def _revcomp(text: str) -> str:
    """Reverse complement keeping '-' gap columns in place."""
    return text.translate(_COMP_TRANS)[::-1]


def read_bed12(bed_file: str) -> Dict[str, dict]:
    """Parse a BED12 file into region dictionaries.

    Args:
        bed_file: Path to BED12 file.

    Returns:
        Dict mapping region name -> {chrom, strand, blocks}, blocks a
        list of ascending (start, end) genomic intervals.
    """
    regions = {}
    n_malformed = 0
    with open(bed_file) as fh:
        for line in fh:
            line = line.rstrip("\n")
            if not line.strip() or line.startswith(("#", "track", "browser")):
                continue
            fields = line.split("\t")
            if len(fields) < 12:
                logging.warning("Skipping incomplete BED12 line: %s", line)
                n_malformed += 1
                continue
            chrom = fields[0]
            try:
                start = int(fields[1])
                block_sizes = [int(x) for x in fields[10].rstrip(",").split(",") if x]
                block_starts = [int(x) for x in fields[11].rstrip(",").split(",") if x]
            except ValueError:
                logging.warning("Skipping malformed BED12 line: %s", line)
                n_malformed += 1
                continue
            strand = fields[5]
            name = fields[3]
            if strand not in ("+", "-"):
                logging.warning("Skipping %s: unrecognized strand %r.", name, strand)
                n_malformed += 1
                continue
            if not block_sizes or len(block_sizes) != len(block_starts):
                logging.warning("Skipping %s: inconsistent block arrays.", name)
                n_malformed += 1
                continue
            if fields[9].isdigit() and int(fields[9]) != len(block_sizes):
                logging.warning("Skipping %s: blockCount != len(blockSizes).", name)
                n_malformed += 1
                continue
            if name in regions:
                logging.warning("Duplicate region name %r; keeping the last occurrence.", name)
            blocks = sorted((start + bs, start + bs + sz) for bs, sz in zip(block_starts, block_sizes))
            regions[name] = {"chrom": chrom, "strand": strand, "blocks": blocks}

    if n_malformed:
        logging.warning("Skipped %d malformed BED12 lines.", n_malformed)
    logging.info("Loaded %d regions from %s", len(regions), bed_file)
    return regions


def _parse_s_line(line: str) -> Optional[dict]:
    """Parse an MAF 's' line into a component dict, or None if invalid."""
    fields = line.split()
    if len(fields) < 7:
        return None
    src = fields[1]
    try:
        start, size = int(fields[2]), int(fields[3])
    except ValueError:
        return None
    strand = fields[4]
    text = fields[6]
    if strand not in ("+", "-"):
        return None
    if len(text.replace("-", "")) != size:
        logging.warning("MAF s-line text length != size for %s; skipping component.", src)
        return None
    species, _, chrom = src.partition(".")
    chrom = chrom.split(":")[0]  # tolerate optional ':start-end' suffix
    return {
        "species": species,
        "chrom": chrom,
        "start": start,
        "size": size,
        "strand": strand,
        "text": text,
    }


def _iter_maf_components(maf_file: str):
    """Yield reference-matching component groups as (component, None) per
    alignment block. Yields every component of every block; callers
    filter by species/chrom."""
    block: List[dict] = []
    with open(maf_file) as fh:
        for line in fh:
            if line.startswith("s"):
                comp = _parse_s_line(line)
                if comp is not None:
                    block.append(comp)
            elif not line.strip() and block:
                yield block
                block = []
    if block:
        yield block


def _component_pos_map(comp: dict) -> Dict[int, str]:
    """Map genomic (plus-strand) positions to bases for one component.

    Minus-strand components are reverse-complemented first so that the
    map always yields plus-strand genome sequence; gap columns do not
    advance the position.
    """
    text = comp["text"] if comp["strand"] == "+" else _revcomp(comp["text"])
    pos_map = {}
    pos = comp["start"]
    for ch in text:
        if ch != "-":
            pos_map[pos] = ch
            pos += 1
    return pos_map


def extract_sequences(
    maf_file: str,
    regions: Dict[str, dict],
    reference_species: str,
) -> Dict[str, str]:
    """Extract sequences from MAF for each BED12 region.

    Args:
        maf_file: Path to MAF alignment file.
        regions: BED12 region dicts from read_bed12().
        reference_species: Reference species prefix (e.g., 'hg38').

    Returns:
        Dict mapping region name -> assembled nucleotide sequence, in
        transcript orientation (reverse-complemented for '-' regions).
    """
    regions_by_chrom: Dict[str, List[Tuple[str, dict]]] = {}
    for name, info in regions.items():
        regions_by_chrom.setdefault(info["chrom"], []).append((name, info))

    # (region_name, interval_index) -> {genome_pos: base}, first match wins.
    collected: Dict[Tuple[str, int], Dict[int, str]] = {}
    hit_regions = set()

    for block in _iter_maf_components(maf_file):
        for comp in block:
            if comp["species"] != reference_species:
                continue
            candidates = regions_by_chrom.get(comp["chrom"])
            if not candidates:
                continue
            pos_map = _component_pos_map(comp)
            comp_end = comp["start"] + comp["size"]
            for name, info in candidates:
                for idx, (bs, be) in enumerate(info["blocks"]):
                    if be <= comp["start"] or bs >= comp_end:
                        continue
                    hit_regions.add(name)
                    target = collected.setdefault((name, idx), {})
                    for pos in range(max(bs, comp["start"]), min(be, comp_end)):
                        base = pos_map.get(pos)
                        if base is not None:
                            target.setdefault(pos, base)

    sequences = {}
    for name, info in regions.items():
        if name not in hit_regions:
            continue
        plus = "".join(
            "".join(collected.get((name, idx), {}).get(p, "N") for p in range(bs, be))
            for idx, (bs, be) in enumerate(info["blocks"])
        )
        sequences[name] = (
            str(Seq(plus).reverse_complement()) if info["strand"] == "-" else plus
        )

    logging.info("Extracted %d/%d regions.", len(sequences), len(regions))
    return sequences


def write_fasta(sequences: Dict[str, str], output_path: str, translate: bool = False) -> None:
    """Write sequences to FASTA file.

    Args:
        sequences: Dict of name -> sequence.
        output_path: Output FASTA path.
        translate: If True, translate nucleotide to protein.
    """
    records = []
    for name, seq in sequences.items():
        if translate:
            seq = str(Seq(seq).translate(to_stop=True))
        records.append(SeqRecord(Seq(seq), id=name, description=""))
    SeqIO.write(records, output_path, "fasta")
    logging.info("Wrote %d sequences -> %s", len(records), output_path)


def build_parser() -> argparse.ArgumentParser:
    """Build argument parser."""
    parser = argparse.ArgumentParser(
        description="Extract genomic regions from MAF using BED12 coordinates.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""\
examples:
  python maf_extract_regions.py --bed cds.bed12 --maf alignment.maf --ref hg38 --cds output.fa --protein protein.fa

notes:
  The reference component is matched by exact species prefix and
  chromosome (src 'hg38.chr1' with --ref hg38). Minus-strand components
  and minus-strand BED regions are reverse-complemented so output is in
  transcript orientation. Uncovered positions are filled with 'N'.
""",
    )
    parser.add_argument("--bed", type=str, required=True, help="Input BED12 file with regions to extract.")
    parser.add_argument("--maf", type=str, required=True, help="Input MAF alignment file.")
    parser.add_argument("--ref", type=str, required=True, help="Reference species name (e.g., hg38, mm10).")
    parser.add_argument("--cds", type=str, required=True, help="Output nucleotide FASTA file.")
    parser.add_argument("--protein", type=str, required=True, help="Output protein FASTA file.")
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

    for label, path in [("BED12", args.bed), ("MAF", args.maf)]:
        if not Path(path).is_file():
            logging.error("%s file not found: %s", label, path)
            sys.exit(1)

    regions = read_bed12(args.bed)
    if not regions:
        logging.error("No valid regions in BED12 file: %s", args.bed)
        sys.exit(1)

    sequences = extract_sequences(args.maf, regions, args.ref)
    if not sequences:
        logging.error(
            "No regions found in the MAF; check --ref (%s) and chromosome names.",
            args.ref,
        )
        sys.exit(1)

    write_fasta(sequences, args.cds, translate=False)
    write_fasta(sequences, args.protein, translate=True)


if __name__ == "__main__":
    main()
