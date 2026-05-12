#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#########################################################################
# File Name: search_riboseq_bioproject.py
# Author: ChengYu
# Description: Multi-layer Ribo-seq BioProject retrieval from NCBI
# Created Time: 2026
#########################################################################
"""
Multi-layer Ribo-seq BioProject retrieval from NCBI.

Searches the NCBI BioProject database using a configurable multi-tier
keyword strategy (strict / medium / broad), resolves species names to
TaxIDs, fetches BioProject metadata via XML, links each project to SRA,
and writes a consolidated CSV report.

Features:
  - Configurable keyword tiers (strict, medium, broad)
  - Species name to TaxID auto-resolution
  - BioProject -> SRA linking
  - Robust XML parsing for NCBI / ENA / DDBJ variants
  - CSV output
"""

import argparse
import csv
import logging
import os
import sys
import time
from typing import Dict, List, Optional
from urllib.error import HTTPError, URLError
from xml.etree import ElementTree as ET

from Bio import Entrez

__version__ = "1.0.0"

logger = logging.getLogger(__name__)

DEFAULT_SLEEP = 0.34  # seconds between Entrez calls

# ============================================================================
# DEFAULT KEYWORD DICTIONARIES
# ============================================================================

DEFAULT_STRICT_KEYWORDS = [
    "ribosome profiling",
    "ribo-seq", "riboseq", "ribo seq",
    "ribosome footprinting", "ribosome footprints",
    "ribosome protected fragments", "ribosome-protected fragments",
    "RPF", "RPF-seq",
]

DEFAULT_MEDIUM_KEYWORDS = [
    "translatome", "translatomics",
    "translatomic profiling", "translatome sequencing",
    "translation-level profiling",
    "actively translated mRNA",
    "actively translating",
    "polysome profiling", "polysome-seq", "polysome sequencing",
    "polyribosome profiling",
    "polysome fractionation",
    "ribosome-associated mRNA",
    "ribosome bound mRNA",
    "ribosome-associated transcriptome",
    "TE profiling",
    "translation-level analysis",
    "ribosome occupancy",
    "ribosome release",
    "ribosome pausing",
    "ribosome stalling",
    "translation dynamics",
]

DEFAULT_BROAD_PATTERNS = [
    "(ribosome[All Fields] AND (profiling[All Fields] OR footprinting[All Fields]))",
    "(translation[All Fields] AND (profiling[All Fields] OR footprint[All Fields]))",
]


# ============================================================================
# TAXONOMY
# ============================================================================

def species_to_taxid(species: str, sleep: float = DEFAULT_SLEEP) -> Optional[int]:
    """Resolve a species name to an NCBI TaxID."""
    try:
        logger.info("Resolving TaxID for: %s", species)
        handle = Entrez.esearch(db="taxonomy", term=f'"{species}"[Scientific Name]')
        data = Entrez.read(handle)
        handle.close()
        time.sleep(sleep)
        ids = data.get("IdList", [])
        if ids:
            taxid = int(ids[0])
            logger.info("Resolved TaxID: %d", taxid)
            return taxid
        logger.warning("Species not found: %s", species)
    except Exception as exc:
        logger.warning("Taxonomy lookup failed for %s: %s", species, exc)
    return None


# ============================================================================
# QUERY BUILDING
# ============================================================================

def build_queries(
    taxid: Optional[int],
    strict_keywords: Optional[List[str]] = None,
    medium_keywords: Optional[List[str]] = None,
    broad_patterns: Optional[List[str]] = None,
) -> List[tuple]:
    """
    Build (tier_name, query_string) tuples for each keyword tier.

    Args:
        taxid: Optional NCBI TaxID to restrict results.
        strict_keywords: Override strict keyword list.
        medium_keywords: Override medium keyword list.
        broad_patterns: Override broad pattern list.

    Returns:
        List of ``(tier, query)`` pairs.
    """
    strict = strict_keywords or DEFAULT_STRICT_KEYWORDS
    medium = medium_keywords or DEFAULT_MEDIUM_KEYWORDS
    broad = broad_patterns or DEFAULT_BROAD_PATTERNS

    suffix = f" AND txid{taxid}[Organism]" if taxid else ""

    strict_q = " OR ".join(f'"{k}"[All Fields]' for k in strict)
    medium_q = " OR ".join(f'"{k}"[All Fields]' for k in medium)
    broad_q = " OR ".join(broad)

    return [
        ("STRICT", f"({strict_q}){suffix}"),
        ("MEDIUM", f"({medium_q}){suffix}"),
        ("BROAD",  f"({broad_q}){suffix}"),
    ]


# ============================================================================
# ESEARCH (paginated)
# ============================================================================

def esearch_all(query: str, retmax: int = 500, sleep: float = DEFAULT_SLEEP) -> List[str]:
    """Paginate through ``esearch`` on the bioproject database."""
    all_ids: List[str] = []
    retstart = 0
    logger.debug("Query: %s", query)

    while True:
        try:
            handle = Entrez.esearch(
                db="bioproject",
                term=query,
                retmax=retmax,
                retstart=retstart,
                usehistory="n",
            )
            data = Entrez.read(handle)
            handle.close()
            time.sleep(sleep)
        except (HTTPError, URLError) as exc:
            logger.error("esearch network error: %s", exc)
            break
        except Exception as exc:
            logger.error("esearch error: %s", exc)
            break

        batch = data.get("IdList", [])
        if not batch:
            break
        all_ids.extend(batch)
        if len(batch) < retmax:
            break
        retstart += retmax

    return all_ids


# ============================================================================
# EFETCH BIOPROJECT DETAILS
# ============================================================================

def fetch_bioproject(ids: List[str], sleep: float = DEFAULT_SLEEP) -> List[Dict]:
    """
    Fetch and parse BioProject XML for a list of BioProject UIDs.

    Robustly handles NCBI, ENA, and DDBJ XML variants.
    """
    if not ids:
        return []

    results: List[Dict] = []
    chunk_size = 50

    for i in range(0, len(ids), chunk_size):
        sub = ids[i : i + chunk_size]
        try:
            handle = Entrez.efetch(db="bioproject", id=",".join(sub), rettype="xml")
            xml_bytes = handle.read()
            handle.close()
            time.sleep(sleep)
        except Exception as exc:
            logger.warning("efetch error at batch %d: %s", i, exc)
            continue

        try:
            root = ET.fromstring(xml_bytes)
        except ET.ParseError as exc:
            logger.warning("XML parse error at batch %d: %s", i, exc)
            continue

        for project in root.findall(".//Project"):
            rec = _parse_project(project)
            if rec:
                results.append(rec)

    return results


def _parse_project(project: ET.Element) -> Optional[Dict]:
    """Extract fields from a single ``<Project>`` XML element."""
    try:
        rec: Dict = {}

        # Accession (always present)
        arch = project.find(".//ArchiveID")
        rec["Accession"] = arch.attrib.get("accession", "") if arch is not None else ""

        # Numeric BioProject ID
        pid = project.find(".//ProjectID")
        if pid is not None and pid.text:
            rec["BioProjectID"] = pid.text.strip()
        else:
            rec["BioProjectID"] = rec["Accession"]

        # Title
        rec["Title"] = _first_text(project, [
            ".//Title",
            ".//ProjectDescr/Title",
            ".//Description/Title",
            ".//Name",
        ])

        # Organism
        rec["Organism"] = _first_text(project, [
            ".//Organism/OrganismName",
            ".//Organism/Name",
            ".//Organism",
        ])

        # Description
        rec["Description"] = _first_text(project, [
            ".//ProjectDescr/Description",
            ".//ProjectDescription",
            ".//Description/ProjectDescription",
            ".//Description",
        ])

        # Release date
        rec["ReleaseDate"] = _first_text(project, [
            ".//ReleaseDate",
            ".//Submission/ReleaseDate",
            ".//ProjectReleaseDate",
        ])

        return rec
    except Exception as exc:
        logger.warning("Failed to parse Project element: %s", exc)
        return None


def _first_text(element: ET.Element, paths: List[str]) -> str:
    """Return stripped text from the first matching XPath."""
    for path in paths:
        node = element.find(path)
        if node is not None and node.text:
            return node.text.strip()
    return ""


# ============================================================================
# BIOPROJECT -> SRA LINKING
# ============================================================================

def link_sra(bioproject_id: str, sleep: float = DEFAULT_SLEEP) -> List[str]:
    """Link a BioProject UID to SRA UIDs via ``elink``."""
    sra_ids: List[str] = []
    try:
        handle = Entrez.elink(dbfrom="bioproject", db="sra", id=bioproject_id)
        data = Entrez.read(handle)
        handle.close()
        time.sleep(sleep)
        for block in data:
            for linkset in block.get("LinkSetDb", []):
                if linkset.get("DbTo") == "sra":
                    for link in linkset.get("Link", []):
                        sra_ids.append(link.get("Id"))
    except Exception as exc:
        logger.debug("elink failed for %s: %s", bioproject_id, exc)
    return sra_ids


# ============================================================================
# CSV OUTPUT
# ============================================================================

CSV_FIELDS = [
    "QuerySpecies", "BioProjectID", "Accession",
    "Title", "Organism", "ReleaseDate", "Description", "LinkedSRA",
]


def save_csv(records: List[Dict], outfile: str) -> None:
    """Write records to a CSV file."""
    with open(outfile, "w", newline="", encoding="utf-8") as fh:
        writer = csv.DictWriter(fh, CSV_FIELDS)
        writer.writeheader()
        for rec in records:
            rec["LinkedSRA"] = ";".join(rec.get("LinkedSRA", []))
            writer.writerow(rec)
    logger.info("Saved %d records to %s", len(records), outfile)


# ============================================================================
# CLI
# ============================================================================

def build_parser() -> "argparse.ArgumentParser":
    """Construct and return the argument parser."""
    parser = argparse.ArgumentParser(
        prog="search_riboseq_bioproject.py",
        description="Multi-layer Ribo-seq BioProject retrieval from NCBI",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""\
Examples:
  # Search a single species
  %(prog)s -s "Arabidopsis thaliana" -o results.csv

  # Search multiple species / TaxIDs
  %(prog)s -s "Oryza sativa" "Zea mays" 3702 -o multi_results.csv

  # All organisms (no species filter)
  %(prog)s -o all_riboseq.csv

  # Custom Entrez email
  %(prog)s -s "Homo sapiens" --email user@example.com -o human.csv

  # Verbose logging
  %(prog)s -s "Homo sapiens" -o human.csv --log-level DEBUG
""",
    )

    parser.add_argument("--version", action="version", version=f"%(prog)s {__version__}")
    parser.add_argument(
        "-s", "--species", nargs="+",
        help="Species name(s) or TaxID(s). Omit to search all organisms.",
        default=[],
    )
    parser.add_argument(
        "-o", "--output", default="ribo_bioproject_recall.csv",
        help="Output CSV file path (default: ribo_bioproject_recall.csv)",
    )
    parser.add_argument(
        "--email", type=str, default=None,
        help="NCBI Entrez email (default: $NCBI_EMAIL env var)",
    )
    parser.add_argument(
        "--api-key", type=str, default=None,
        help="NCBI API key (default: $NCBI_API_KEY env var)",
    )
    parser.add_argument(
        "--api-delay", type=float, default=DEFAULT_SLEEP,
        help="Delay between Entrez calls in seconds (default: 0.34)",
    )
    parser.add_argument(
        "--tiers", type=str,
        choices=["strict", "medium", "broad", "all"], default="all",
        help="Which keyword tiers to use (default: all)",
    )
    parser.add_argument(
        "--log-level", type=str,
        choices=["DEBUG", "INFO", "WARNING", "ERROR", "CRITICAL"],
        default="INFO",
        help="Logging level (default: INFO)",
    )
    return parser


def main() -> None:
    """Entry point."""
    parser = build_parser()
    args = parser.parse_args()

    # Logging
    logging.basicConfig(
        level=getattr(logging, args.log_level),
        format="[%(levelname)s] %(message)s",
    )

    # Entrez credentials
    email = args.email or os.getenv("NCBI_EMAIL")
    if not email:
        parser.error("NCBI email required. Use --email or set NCBI_EMAIL env var.")
    Entrez.email = email

    api_key = args.api_key or os.getenv("NCBI_API_KEY")
    if api_key:
        Entrez.api_key = api_key

    sleep = args.api_delay

    # Determine which tiers to run
    tier_filter = None if args.tiers == "all" else {args.tiers}

    species_list = args.species if args.species else [None]
    final_records: List[Dict] = []

    for sp in species_list:
        # Resolve TaxID
        taxid: Optional[int] = None
        label: str = ""

        if sp is not None:
            try:
                taxid = int(sp)
                label = f"taxid:{taxid}"
            except ValueError:
                taxid = species_to_taxid(sp, sleep=sleep)
                label = f"{sp}(txid={taxid})"
        else:
            label = "ALL_SPECIES"

        logger.info("=" * 50)
        logger.info("Species: %s", label)

        # Build queries
        queries = build_queries(taxid)
        if tier_filter:
            queries = [(t, q) for t, q in queries if t.lower() in tier_filter]

        # Search
        id_set: set = set()
        for tier_name, query in queries:
            logger.info("Running %s tier ...", tier_name)
            ids = esearch_all(query, sleep=sleep)
            logger.info("%s: found %d IDs", tier_name, len(ids))
            id_set.update(ids)

        logger.info("Unique BioProject IDs: %d", len(id_set))

        # Fetch details
        details = fetch_bioproject(list(id_set), sleep=sleep)

        # Link SRA
        for rec in details:
            bid = rec.get("BioProjectID", "")
            rec["LinkedSRA"] = link_sra(bid, sleep=sleep)
            rec["QuerySpecies"] = label
            final_records.append(rec)

    # Write output
    save_csv(final_records, args.output)
    logger.info("Done. %d records written to %s", len(final_records), args.output)


if __name__ == "__main__":
    main()
