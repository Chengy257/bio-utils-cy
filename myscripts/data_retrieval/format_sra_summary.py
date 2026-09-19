#!/usr/bin/env python3
"""
File Name: format_sra_summary.py
Author: ChengYu
Description: Format comprehensive SRA metadata into a standardized summary TSV
             with four-layer SeqType classification, field mapping, text cleaning,
             and synonym resolution.
Created Time: 2026
"""

from __future__ import annotations

import argparse
import csv
import logging
import re
import sys
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

__version__ = "1.0.0"

logger = logging.getLogger(__name__)

# ---------------------------------------------------------------------------
# Synonym mapping: normalize common terms to canonical names
# ---------------------------------------------------------------------------

SYNONYM_MAP: Dict[str, str] = {
    # Instrument synonyms
    "illumina hiseq 2500": "Illumina HiSeq 2500",
    "illumina hiseq 2000": "Illumina HiSeq 2000",
    "illumina hiseq 1500": "Illumina HiSeq 1500",
    "illumina hiseq 4000": "Illumina HiSeq 4000",
    "illumina hiseq x ten": "Illumina HiSeq X Ten",
    "hiseq 2500": "Illumina HiSeq 2500",
    "hiseq 2000": "Illumina HiSeq 2000",
    "hiseq 4000": "Illumina HiSeq 4000",
    "illumina miseq": "Illumina MiSeq",
    "miseq": "Illumina MiSeq",
    "illumina novaseq 6000": "Illumina NovaSeq 6000",
    "novaseq 6000": "Illumina NovaSeq 6000",
    "illumina novaseq x plus": "Illumina NovaSeq X Plus",
    "illumina genome analyzer": "Illumina Genome Analyzer",
    "illumina genome analyzer ii": "Illumina Genome Analyzer II",
    "illumina genome analyzer iix": "Illumina Genome Analyzer IIx",
    "nextseq 500": "Illumina NextSeq 500",
    "nextseq 550": "Illumina NextSeq 550",
    "illumina nextseq 500": "Illumina NextSeq 500",
    "illumina nextseq 550": "Illumina NextSeq 550",
    "illumina nextseq 1000": "Illumina NextSeq 1000",
    "illumina nextseq 2000": "Illumina NextSeq 2000",
    "bgiseq-500": "BGISeq-500",
    "bgiseq-50": "BGISeq-50",
    "ion torrent pgm": "Ion Torrent PGM",
    "ion torrent s5": "Ion Torrent S5",
    "ion torrent proton": "Ion Torrent Proton",
    "pacbio rs ii": "PacBio RS II",
    "pacbio sequel": "PacBio Sequel",
    "pacbio sequel ii": "PacBio Sequel II",
    "pacbio sequel iie": "PacBio Sequel IIe",
    "ont minion": "ONT MinION",
    "ont gridion": "ONT GridION",
    "ont promethion": "ONT PromethION",
    "minion": "ONT MinION",
    "gridion": "ONT GridION",
    "promethion": "ONT PromethION",
    "oxford nanopore minion": "ONT MinION",
    "oxford nanopore gridion": "ONT GridION",
    "ls454 titanium": "LS454 Titanium",
    "454 gs flx titanium": "LS454 Titanium",
    "ab 5500 genetic analyzer": "AB 5500 Genetic Analyzer",
    # Library strategy synonyms
    "rna-seq": "RNA-Seq",
    "rnaseq": "RNA-Seq",
    "rna seq": "RNA-Seq",
    "wgs": "WGS",
    "whole genome sequencing": "WGS",
    "whole-genome shotgun": "WGS",
    "wes": "WES",
    "whole exome sequencing": "WES",
    "exome sequencing": "WES",
    "amplicon": "Amplicon",
    "16s": "16S",
    "16s rrna": "16S",
    "16s rrna sequencing": "16S",
    "mirna-seq": "miRNA-Seq",
    "mirna seq": "miRNA-Seq",
    "chip-seq": "ChIP-Seq",
    "chip seq": "ChIP-Seq",
    "atac-seq": "ATAC-Seq",
    "atac seq": "ATAC-Seq",
    "bisulfite-seq": "Bisulfite-Seq",
    "bisulfite seq": "Bisulfite-Seq",
    "methyl-seq": "Methyl-Seq",
}

# ---------------------------------------------------------------------------
# Field mapping: input column names -> output standard names
# ---------------------------------------------------------------------------

FIELD_MAP: Dict[str, str] = {
    "Run": "RunAccession",
    "run_accession": "RunAccession",
    "run": "RunAccession",
    "Experiment": "ExperimentAccession",
    "experiment_accession": "ExperimentAccession",
    "experiment": "ExperimentAccession",
    "BioSample": "BioSample",
    "sample_accession": "BioSample",
    "biosample": "BioSample",
    "SampleName": "SampleName",
    "bio_sample": "BioSample",
    "BioProject": "BioProject",
    "study_accession": "BioProject",
    "Study": "BioProject",
    "bioproject": "BioProject",
    "Organism": "Organism",
    "organism": "Organism",
    "ScientificName": "Organism",
    "Instrument": "Instrument",
    "instrument_model": "Instrument",
    "InstrumentModel": "Instrument",
    "LibraryStrategy": "LibraryStrategy",
    "library_strategy": "LibraryStrategy",
    "LibrarySource": "LibrarySource",
    "library_source": "LibrarySource",
    "LibrarySelection": "LibrarySelection",
    "library_selection": "LibrarySelection",
    "LibraryLayout": "LibraryLayout",
    "library_layout": "LibraryLayout",
    "LibraryName": "LibraryName",
    "library_name": "LibraryName",
    "RunInfo": "RunInfo",
    "InsertSize": "InsertSize",
    "insert_size": "InsertSize",
    "InsertDev": "InsertDev",
    "read_count": "ReadCount",
    "ReadCount": "ReadCount",
    "base_count": "BaseCount",
    "BaseCount": "BaseCount",
    "Bases": "BaseCount",
    "Spots": "ReadCount",
    "AvgSpotLen": "AvgSpotLen",
    "avg_spot_len": "AvgSpotLen",
    "fastq_ftp": "FastQ_FTP",
    "FastqFTP": "FastQ_FTP",
    "fastq_bytes": "FastQ_Bytes",
    "fastq_md5": "FastQ_MD5",
    "ReleaseDate": "ReleaseDate",
    "release_date": "ReleaseDate",
    "LoadDate": "LoadDate",
    "load_date": "LoadDate",
    "ExperimentTitle": "ExperimentTitle",
    "experiment_title": "ExperimentTitle",
    "experiment_title": "ExperimentTitle",
    "StudyTitle": "StudyTitle",
    "study_title": "StudyTitle",
    "SampleTitle": "SampleTitle",
    "sample_title": "SampleTitle",
    "biosample_description_title": "BioSampleTitle",
    "biosample_description": "BioSampleDescription",
    "description_title": "BioSampleTitle",
}


# ---------------------------------------------------------------------------
# Text cleaning
# ---------------------------------------------------------------------------

def clean_text(value: str) -> str:
    """Clean a text value: strip whitespace, collapse internal spaces, normalize quotes.

    Args:
        value: Raw text string.

    Returns:
        Cleaned string.
    """
    if not value:
        return ""
    value = value.strip()
    value = re.sub(r"\s+", " ", value)
    value = value.replace("\u201c", '"').replace("\u201d", '"')
    value = value.replace("\u2018", "'").replace("\u2019", "'")
    return value


def resolve_synonym(value: str) -> str:
    """Resolve a value to its canonical form using the synonym map.

    Args:
        value: Raw value string.

    Returns:
        Canonical form, or the cleaned original if no synonym found.
    """
    cleaned = clean_text(value)
    if not cleaned:
        return cleaned
    lookup = cleaned.lower()
    return SYNONYM_MAP.get(lookup, cleaned)


# ---------------------------------------------------------------------------
# Four-layer SeqType classification
# ---------------------------------------------------------------------------

# Layer 1: Direct LibraryStrategy match
_STRATEGY_CLASS: Dict[str, str] = {
    "WGS": "WGS",
    "Whole Genome Sequencing": "WGS",
    "Whole-Genome Shotgun": "WGS",
    "WES": "WES",
    "Whole Exome Sequencing": "WES",
    "Exome Sequencing": "WES",
    "RNA-Seq": "RNA-Seq",
    "RNA Seq": "RNA-Seq",
    "Transcriptome": "RNA-Seq",
    "EST": "RNA-Seq",
    "FL-cDNA": "RNA-Seq",
    "ChIP-Seq": "ChIP-Seq",
    "ChIP Seq": "ChIP-Seq",
    "ATAC-Seq": "ATAC-Seq",
    "ATAC Seq": "ATAC-Seq",
    "miRNA-Seq": "miRNA-Seq",
    "miRNA Seq": "miRNA-Seq",
    "small RNA-Seq": "smallRNA-Seq",
    "smallRNA-Seq": "smallRNA-Seq",
    "ncRNA-Seq": "ncRNA-Seq",
    "ncRNA Seq": "ncRNA-Seq",
    "lncRNA-Seq": "lncRNA-Seq",
    "Bisulfite-Seq": "Bisulfite-Seq",
    "Bisulfite Seq": "Bisulfite-Seq",
    "Methyl-Seq": "Bisulfite-Seq",
    "Methylation": "Bisulfite-Seq",
    "MBD-Seq": "Bisulfite-Seq",
    "MeDIP-Seq": "Bisulfite-Seq",
    "Amplicon": "Amplicon",
    "Targeted Capture": "Amplicon",
    "16S": "16S",
    "16S rRNA": "16S",
    "16S rRNA Sequencing": "16S",
    "AMR": "Amplicon",
    "Other": "Other",
    "Synthetic-Long-Read": "Long-Read",
    "Synthetic Long Read": "Long-Read",
    "Tethered sequencing": "Long-Read",
    "CLONE": "Clone",
    "CLONEEND": "Clone-End",
    "Clone End": "Clone-End",
    "POOLCLONE": "Clone",
    "ChIA-PET": "ChIA-PET",
    "Hi-C": "Hi-C",
    "MNase-Seq": "MNase-Seq",
    "DNase-Seq": "DNase-Seq",
    "FAIRE-seq": "FAIRE-Seq",
    "SELEX": "SELEX",
    "RIP-Seq": "RIP-Seq",
    "RIP Seq": "RIP-Seq",
    "CLIP-Seq": "CLIP-Seq",
    "CLIP Seq": "CLIP-Seq",
    "ReDDi-Seq": "Other",
    "RAD-Seq": "RAD-Seq",
    "GBS": "RAD-Seq",
    "Genotyping by Sequencing": "RAD-Seq",
    "Restriction Site Associated DNA Sequencing": "RAD-Seq",
    "WCS": "WCS",
    "WGA": "WGA",
    "Pseudo-RNA-Seq": "Other",
    "Pseudo-WGS": "Other",
}

# Layer 2: LibrarySource-based heuristics
_SOURCE_CLASS: Dict[str, str] = {
    "GENOMIC": "WGS",
    "GENOMIC SINGLE CELL": "scDNA-Seq",
    "TRANSCRIPTOMIC": "RNA-Seq",
    "TRANSCRIPTOMIC SINGLE CELL": "scRNA-Seq",
    "METAGENOMIC": "Metagenomics",
    "METATRANSCRIPTOMIC": "MetaTranscriptomics",
    "SYNTHETIC": "Other",
    "VIRAL RNA": "RNA-Seq",
}

# Layer 3: Instrument-based heuristics
_LONG_READ_INSTRUMENTS = {
    "PacBio RS II",
    "PacBio Sequel",
    "PacBio Sequel II",
    "PacBio Sequel IIe",
    "PacBio Revio",
    "ONT MinION",
    "ONT GridION",
    "ONT PromethION",
    "Oxford Nanopore MinION",
    "Oxford Nanopore GridION",
    "Oxford Nanopore PromethION",
}

# Layer 4: Title/description keyword heuristics
_TITLE_KEYWORDS: List[Tuple[str, str]] = [
    (r"\b16[sS]\b", "16S"),
    (r"\bmetagenom", "Metagenomics"),
    (r"\bmeta-?genom", "Metagenomics"),
    (r"\bsingle[\s_-]?cell\b", "scRNA-Seq"),
    (r"\bscRNA\b", "scRNA-Seq"),
    (r"\bscDNA\b", "scDNA-Seq"),
    (r"\bchi[ap][\s_-]?seq\b", "ChIA-PET"),
    (r"\bhi[\s_-]?c\b", "Hi-C"),
    (r"\bchip[\s_-]?seq\b", "ChIP-Seq"),
    (r"\batac[\s_-]?seq\b", "ATAC-Seq"),
    (r"\brna[\s_-]?seq\b", "RNA-Seq"),
    (r"\bwhole[\s_-]genome\b", "WGS"),
    (r"\bwgs\b", "WGS"),
    (r"\bexome\b", "WES"),
    (r"\bwes\b", "WES"),
    (r"\bamplicon\b", "Amplicon"),
    (r"\bbisulfite\b", "Bisulfite-Seq"),
    (r"\bmethyl", "Bisulfite-Seq"),
    (r"\bmirna\b", "miRNA-Seq"),
    (r"\blong[\s_-]?read\b", "Long-Read"),
]


def classify_seqtype(
    row: Dict[str, str],
    override: Optional[str] = None,
) -> str:
    """Classify the sequencing type for a row using four-layer logic.

    Layer 1: Direct LibraryStrategy match.
    Layer 2: LibrarySource-based heuristics.
    Layer 3: Instrument-based long-read detection.
    Layer 4: Title/description keyword heuristics.

    Args:
        row: A metadata row dict.
        override: If non-empty, forces this value as the SeqType.

    Returns:
        Classified SeqType string.
    """
    if override:
        return override

    # Gather relevant fields
    strategy = ""
    for key in ("LibraryStrategy", "library_strategy", "LIBRARY_STRATEGY"):
        if key in row and row[key]:
            strategy = row[key].strip()
            break
    source = ""
    for key in ("LibrarySource", "library_source", "LIBRARY_SOURCE"):
        if key in row and row[key]:
            source = row[key].strip()
            break
    instrument = ""
    for key in ("Instrument", "instrument_model", "InstrumentModel", "instrument"):
        if key in row and row[key]:
            instrument = row[key].strip()
            break
    title = ""
    for key in ("ExperimentTitle", "experiment_title", "StudyTitle", "study_title"):
        if key in row and row[key]:
            title = row[key].strip()
            break

    # Layer 1: LibraryStrategy direct match
    if strategy:
        norm = _STRATEGY_CLASS.get(strategy)
        if norm:
            return norm
        # Try case-insensitive
        for k, v in _STRATEGY_CLASS.items():
            if k.lower() == strategy.lower():
                return v

    # Layer 2: LibrarySource heuristics
    if source:
        src_norm = _SOURCE_CLASS.get(source)
        if src_norm:
            return src_norm
        for k, v in _SOURCE_CLASS.items():
            if k.lower() == source.lower():
                return v

    # Layer 3: Instrument-based long-read
    if instrument:
        canon_instr = resolve_synonym(instrument)
        if canon_instr in _LONG_READ_INSTRUMENTS:
            return "Long-Read"

    # Layer 4: Title/description keyword heuristics
    title_lower = title.lower()
    if title_lower:
        for pattern, label in _TITLE_KEYWORDS:
            if re.search(pattern, title_lower):
                return label

    return "Other"


# ---------------------------------------------------------------------------
# Field remapping
# ---------------------------------------------------------------------------

def remap_row(row: Dict[str, str]) -> Dict[str, str]:
    """Remap input column names to standardized output names.

    Also applies text cleaning and synonym resolution to appropriate fields.

    Args:
        row: Raw metadata row.

    Returns:
        Row with standardized column names and cleaned values.
    """
    out: Dict[str, str] = {}
    for key, val in row.items():
        std_name = FIELD_MAP.get(key, key)
        cleaned = clean_text(val)
        # Apply synonym resolution to instrument and strategy fields
        if std_name in ("Instrument", "LibraryStrategy"):
            cleaned = resolve_synonym(cleaned)
        if std_name not in out:
            out[std_name] = cleaned
        else:
            # Append if duplicate mapping
            if cleaned and cleaned != out[std_name]:
                out[std_name] = f"{out[std_name]};{cleaned}"
    return out


# ---------------------------------------------------------------------------
# Output field ordering
# ---------------------------------------------------------------------------

OUTPUT_COLUMNS: List[str] = [
    "RunAccession",
    "ExperimentAccession",
    "BioSample",
    "SampleName",
    "BioProject",
    "Organism",
    "Instrument",
    "LibraryStrategy",
    "LibrarySource",
    "LibrarySelection",
    "LibraryLayout",
    "LibraryName",
    "InsertSize",
    "ReadCount",
    "BaseCount",
    "AvgSpotLen",
    "SeqType",
    "ReleaseDate",
    "LoadDate",
    "ExperimentTitle",
    "StudyTitle",
    "SampleTitle",
    "FastQ_FTP",
    "FastQ_Bytes",
    "FastQ_MD5",
]


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main() -> None:
    """Main entry point."""
    parser = argparse.ArgumentParser(
        description="Format comprehensive SRA metadata into standardized summary TSV.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""\
examples:
  %(prog)s --input metadata.tsv --output summary.tsv
  %(prog)s --input metadata.tsv --output summary.tsv --seq-type-override RNA-Seq
  %(prog)s --input metadata.tsv --output summary.tsv --quiet
  %(prog)s --input metadata.tsv --output - --log-level DEBUG
""",
    )
    parser.add_argument(
        "--input", "-i", required=True, help="Input comprehensive metadata TSV."
    )
    parser.add_argument(
        "--output", "-o", default="-", help="Output summary TSV path (default: stdout)."
    )
    parser.add_argument(
        "--seq-type-override",
        default="",
        help="Force all rows to this SeqType value (e.g. 'RNA-Seq').",
    )
    parser.add_argument(
        "--quiet", action="store_true", help="Suppress progress messages."
    )
    parser.add_argument(
        "--log-level",
        default="INFO",
        choices=["DEBUG", "INFO", "WARNING", "ERROR", "CRITICAL"],
        help="Set logging level (default: INFO).",
    )
    parser.add_argument("--version", action="version", version=f"%(prog)s {__version__}")
    args = parser.parse_args()

    logging.basicConfig(
        level=getattr(logging, args.log_level),
        format="%(asctime)s [%(levelname)s] %(name)s: %(message)s",
    )

    in_path = Path(args.input)
    if not in_path.exists():
        logger.error("Input file not found: %s", args.input)
        sys.exit(1)

    # Read input TSV
    rows: List[Dict[str, str]] = []
    with open(in_path, "r", encoding="utf-8") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            rows.append(dict(row))
    logger.info("Read %d rows from %s", len(rows), args.input)

    if not rows:
        logger.error("No data rows in input file.")
        sys.exit(1)

    # Process rows
    processed: List[Dict[str, str]] = []
    extra_cols: Dict[str, None] = {}
    for row in rows:
        remapped = remap_row(row)
        seqtype = classify_seqtype(remapped, override=args.seq_type_override or None)
        remapped["SeqType"] = seqtype
        processed.append(remapped)
        for k in remapped:
            if k not in OUTPUT_COLUMNS:
                extra_cols[k] = None

    if not args.quiet:
        seqtype_counts: Dict[str, int] = {}
        for r in processed:
            st = r.get("SeqType", "Other")
            seqtype_counts[st] = seqtype_counts.get(st, 0) + 1
        logger.info("SeqType distribution:")
        for st, cnt in sorted(seqtype_counts.items(), key=lambda x: -x[1]):
            logger.info("  %s: %d", st, cnt)

    # Determine final column order: OUTPUT_COLUMNS first, then extras sorted
    all_cols = [c for c in OUTPUT_COLUMNS if any(c in r for r in processed)]
    extras_sorted = sorted(extra_cols.keys())
    all_cols.extend(extras_sorted)

    # Write output
    if args.output == "-":
        writer = csv.DictWriter(
            sys.stdout, fieldnames=all_cols, delimiter="\t", extrasaction="ignore"
        )
        writer.writeheader()
        writer.writerows(processed)
    else:
        out = Path(args.output)
        out.parent.mkdir(parents=True, exist_ok=True)
        with open(out, "w", newline="", encoding="utf-8") as fh:
            writer = csv.DictWriter(
                fh, fieldnames=all_cols, delimiter="\t", extrasaction="ignore"
            )
            writer.writeheader()
            writer.writerows(processed)
        logger.info("Wrote %d rows to %s", len(processed), args.output)


if __name__ == "__main__":
    main()
