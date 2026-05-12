#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#########################################################################
# File Name: search_sra_riboseq.py
# Author: ChengYu
# Description: Search NCBI SRA for Ribo-seq data by species name or TaxID
# Created Time: 2026
#########################################################################
"""
Search NCBI for Ribo-seq data by species name or TaxID.

Uses a two-tier keyword strategy (STRICT + MEDIUM) to search BioProjects,
extracts linked SRR accessions, and generates summary reports.

Features:
  - Search by species name or TaxID
  - Two-tier keyword strategy (STRICT + MEDIUM)
  - Automatic BioProject to SRA linking
  - Custom keyword support
  - Detailed search reports in TSV and plain text
"""

import argparse
import csv
import logging
import os
import re
import sys
import time
from datetime import datetime
from pathlib import Path
from typing import Dict, List, Optional, Set
from xml.etree import ElementTree as ET

from Bio import Entrez
from tqdm import tqdm

__version__ = "1.0.0"

logger = logging.getLogger(__name__)

# ============================================================================
# BUILT-IN KEYWORD CONFIGURATION
# ============================================================================

# Strict keywords -- high specificity for ribosome profiling
DEFAULT_STRICT_KEYWORDS = [
    "ribosome profiling",
    "ribo-seq", "riboseq", "ribo seq",
    "ribosome footprinting", "ribosome footprints",
    "ribosome protected fragments", "ribosome-protected fragments",
    "RPF", "RPF-seq",
]

# Medium keywords -- broader translation-related terms
DEFAULT_MEDIUM_KEYWORDS = [
    "translatome", "translatomics", "translatomic profiling",
    "translatome sequencing", "translation-level profiling",
    "actively translated mRNA",
    "polysome profiling", "polysome-seq", "polysome sequencing",
    "polyribosome profiling", "polysome fractionation",
    "ribosome-associated mRNA", "ribosome bound mRNA",
    "ribosome-associated transcriptome",
    "TE profiling", "translation efficiency",
    "ribosome occupancy", "ribosome release", "ribosome pausing",
    "ribosome stalling", "translation dynamics",
]


# ============================================================================
# CORE SEARCHER CLASS
# ============================================================================

class NCBISearcher:
    """Search NCBI for Ribo-seq data."""

    def __init__(
        self,
        email: str,
        api_key: Optional[str] = None,
        api_delay: float = 0.34,
        custom_keywords: Optional[Dict[str, List[str]]] = None,
    ):
        """
        Initialize NCBI Searcher.

        Args:
            email: Email for NCBI Entrez (required by NCBI).
            api_key: Optional NCBI API key for higher rate limits.
            api_delay: Seconds to sleep between API calls.
            custom_keywords: Dict with ``'strict'`` and ``'medium'`` keys
                overriding built-in keyword lists.
        """
        Entrez.email = email
        if api_key:
            Entrez.api_key = api_key

        self.api_delay = api_delay

        if custom_keywords:
            self.strict_keywords = custom_keywords.get("strict", DEFAULT_STRICT_KEYWORDS)
            self.medium_keywords = custom_keywords.get("medium", DEFAULT_MEDIUM_KEYWORDS)
        else:
            self.strict_keywords = DEFAULT_STRICT_KEYWORDS
            self.medium_keywords = DEFAULT_MEDIUM_KEYWORDS

    # ------------------------------------------------------------------
    # Taxonomy helpers
    # ------------------------------------------------------------------

    def species_to_taxid(self, species: str) -> Optional[int]:
        """Convert a species name to an NCBI TaxID."""
        try:
            logger.info("Resolving TaxID for species: %s", species)
            handle = Entrez.esearch(
                db="taxonomy",
                term=f'"{species}"[Scientific Name]',
                retmax=1,
            )
            record = Entrez.read(handle)
            handle.close()
            time.sleep(self.api_delay)

            if record["IdList"]:
                taxid = int(record["IdList"][0])
                logger.info("Resolved TaxID: %d", taxid)
                return taxid

            logger.warning("Species not found in taxonomy database: %s", species)
            return None
        except Exception as exc:
            logger.error("Failed to convert species to TaxID: %s", exc)
            return None

    def detect_input_type(self, input_value: str) -> Dict[str, Optional[int]]:
        """Return dict describing whether *input_value* is a TaxID or species name."""
        try:
            taxid = int(input_value)
            return {
                "type": "taxid",
                "value": input_value,
                "taxid": taxid,
                "label": f"taxid:{taxid}",
            }
        except ValueError:
            pass

        taxid = self.species_to_taxid(input_value)
        if taxid:
            return {
                "type": "species",
                "value": input_value,
                "taxid": taxid,
                "label": f"{input_value}(txid={taxid})",
            }
        return {
            "type": "species",
            "value": input_value,
            "taxid": None,
            "label": input_value,
        }

    # ------------------------------------------------------------------
    # Query building
    # ------------------------------------------------------------------

    def build_queries(self, taxid: Optional[int]) -> List[tuple]:
        """
        Build (label, query) tuples for STRICT and MEDIUM tiers.

        Returns:
            List of ``(tier_name, query_string)`` tuples.
        """
        org_filter = f" AND txid{taxid}[Organism]" if taxid else ""

        strict_terms = " OR ".join(f'"{k}"[All Fields]' for k in self.strict_keywords)
        medium_terms = " OR ".join(f'"{k}"[All Fields]' for k in self.medium_keywords)

        return [
            ("STRICT", f"({strict_terms}){org_filter}"),
            ("MEDIUM", f"({medium_terms}){org_filter}"),
        ]

    # ------------------------------------------------------------------
    # Entrez helpers
    # ------------------------------------------------------------------

    def esearch_all(self, database: str, query: str, retmax: int = 10000) -> List[str]:
        """Paginate through ``esearch`` to collect every matching ID."""
        all_ids: List[str] = []
        retstart = 0

        while True:
            try:
                handle = Entrez.esearch(
                    db=database,
                    term=query,
                    retmax=retmax,
                    retstart=retstart,
                )
                data = Entrez.read(handle)
                handle.close()
                time.sleep(self.api_delay)
            except Exception as exc:
                logger.warning("esearch failed at retstart=%d: %s", retstart, exc)
                break

            batch = data.get("IdList", [])
            if not batch:
                break
            all_ids.extend(batch)
            if len(batch) < retmax:
                break
            retstart += retmax

        return all_ids

    # ------------------------------------------------------------------
    # BioProject details
    # ------------------------------------------------------------------

    def fetch_bioproject_details(self, ids: List[str]) -> List[Dict]:
        """Fetch and parse XML details for a list of BioProject IDs."""
        bioprojects: List[Dict] = []
        batch_size = 50

        logger.info("Fetching details for %d BioProjects ...", len(ids))

        for i in tqdm(range(0, len(ids), batch_size), desc="BioProject details"):
            batch = ids[i : i + batch_size]
            try:
                handle = Entrez.efetch(
                    db="bioproject",
                    id=",".join(batch),
                    rettype="xml",
                    retmode="xml",
                )
                xml_data = handle.read()
                handle.close()
                time.sleep(self.api_delay)

                if isinstance(xml_data, bytes):
                    xml_data = xml_data.decode("utf-8")

                root = ET.fromstring(xml_data)
                for project in root.findall(".//Project"):
                    record = self._parse_bioproject_xml(project)
                    if record:
                        bioprojects.append(record)
            except Exception as exc:
                logger.warning("Failed to fetch batch starting at %d: %s", i, exc)

        return bioprojects

    @staticmethod
    def _parse_bioproject_xml(project: ET.Element) -> Optional[Dict]:
        """Parse a single BioProject XML ``<Project>`` element."""
        record: Dict = {}

        try:
            # Accession
            acc_elem = project.find(".//ArchiveID")
            record["Accession"] = (
                acc_elem.attrib.get("accession", "") if acc_elem is not None else ""
            )

            # BioProject ID
            pid_elem = project.find(".//ProjectID")
            if pid_elem is not None and pid_elem.text:
                record["BioProjectID"] = pid_elem.text.strip()
            else:
                record["BioProjectID"] = record["Accession"]

            # Title (try multiple paths)
            record["Title"] = _find_first_text(project, [
                ".//Title",
                ".//ProjectDescr/Title",
                ".//Description/Title",
                ".//Name",
            ])

            # Organism
            record["Organism"] = _find_first_text(project, [
                ".//Organism/OrganismName",
                ".//Organism",
                ".//Name",
            ])

            # Description
            record["Description"] = _find_first_text(project, [
                ".//ProjectDescr/Description",
                ".//Description",
                ".//Comment",
            ])

            # Release date
            record["ReleaseDate"] = _find_first_text(project, [
                ".//ProjectReleaseDate",
                ".//ReleaseDate",
                ".//Date",
            ])

            return record
        except Exception as exc:
            logger.warning("Failed to parse BioProject XML: %s", exc)
            return None

    # ------------------------------------------------------------------
    # BioProject -> SRA linking
    # ------------------------------------------------------------------

    def link_sra(self, bioproject_id: str) -> List[str]:
        """
        Retrieve SRR accessions linked to a BioProject.

        Uses ``esearch`` + ``esummary`` against the SRA database.
        """
        srr_ids: Set[str] = set()

        try:
            query = f"{bioproject_id}[BioProject]"
            handle = Entrez.esearch(db="sra", term=query, retmax=10000)
            search_data = Entrez.read(handle)
            handle.close()
            time.sleep(self.api_delay)

            uid_list = search_data.get("IdList", [])
            if not uid_list:
                return []

            for i in range(0, len(uid_list), 100):
                batch_uids = uid_list[i : i + 100]
                try:
                    handle = Entrez.esummary(db="sra", id=",".join(batch_uids))
                    summary_data = Entrez.read(handle)
                    handle.close()
                    time.sleep(self.api_delay)

                    for item in summary_data:
                        runs_field = str(item.get("Runs", ""))
                        srr_ids.update(re.findall(r'acc="(SRR\d+)"', runs_field))
                except Exception as exc:
                    logger.warning("esummary batch %d failed: %s", i, exc)

        except Exception as exc:
            logger.warning("Failed to link SRA for %s: %s", bioproject_id, exc)

        return sorted(srr_ids)

    # ------------------------------------------------------------------
    # Main search workflow
    # ------------------------------------------------------------------

    def search(self, species_or_taxid: str, output_dir: Path) -> Dict:
        """
        Execute the full search pipeline.

        Returns:
            Dict with keys ``bioprojects``, ``srr_ids``, ``statistics``,
            ``output_files``.
        """
        logger.info("=" * 60)
        logger.info("NCBI SRA RIBO-SEQ SEARCH")
        logger.info("=" * 60)

        # 1. Detect input type
        input_info = self.detect_input_type(species_or_taxid)
        logger.info("Input type: %s  |  Query: %s", input_info["type"], input_info["label"])

        # 2. Build queries
        queries = self.build_queries(input_info["taxid"])
        for tier_name, _ in queries:
            kw_count = (
                len(self.strict_keywords) if tier_name == "STRICT"
                else len(self.medium_keywords)
            )
            logger.info("%s tier: %d keywords", tier_name, kw_count)

        # 3. Search BioProjects
        all_bioproject_ids: Set[str] = set()
        for tier_name, query in queries:
            logger.info("Running %s query ...", tier_name)
            ids = self.esearch_all("bioproject", query)
            logger.info("%s: found %d BioProjects", tier_name, len(ids))
            all_bioproject_ids.update(ids)

        logger.info("Total unique BioProjects: %d", len(all_bioproject_ids))
        if not all_bioproject_ids:
            logger.error("No BioProjects found!")
            return {"error": "No BioProjects found"}

        # 4. Fetch details
        bioprojects = self.fetch_bioproject_details(list(all_bioproject_ids))
        logger.info("Retrieved details for %d BioProjects", len(bioprojects))

        # 5. Link to SRA
        logger.info("Linking BioProjects to SRA ...")
        all_srr_ids: List[str] = []
        for bp in tqdm(bioprojects, desc="Linking SRA"):
            bp["QuerySpecies"] = input_info["label"]
            srr_ids = self.link_sra(bp["BioProjectID"])
            bp["LinkedSRA"] = srr_ids
            bp["SRR_Count"] = len(srr_ids)
            all_srr_ids.extend(srr_ids)

        logger.info("Total SRR IDs extracted: %d", len(all_srr_ids))

        # 6. Save results
        output_dir = Path(output_dir)
        output_dir.mkdir(parents=True, exist_ok=True)
        output_files = self._save_results(bioprojects, all_srr_ids, input_info["label"], output_dir)

        # 7. Statistics report
        stats = self._generate_statistics(bioprojects, all_srr_ids, input_info["label"])
        self._print_summary(stats, output_files)

        return {
            "bioprojects": bioprojects,
            "srr_ids": all_srr_ids,
            "statistics": stats,
            "output_files": output_files,
        }

    # ------------------------------------------------------------------
    # Output helpers
    # ------------------------------------------------------------------

    @staticmethod
    def _save_results(
        bioprojects: List[Dict],
        srr_ids: List[str],
        label: str,
        output_dir: Path,
    ) -> Dict[str, str]:
        """Write SRR IDs, BioProject summary TSV, and a text report."""
        output_files: Dict[str, str] = {}

        # SRR ID list
        srr_file = output_dir / "srr_ids.txt"
        with open(srr_file, "w") as fh:
            for srr_id in srr_ids:
                fh.write(f"{srr_id}\n")
        output_files["srr_ids"] = str(srr_file)
        logger.info("Saved SRR IDs to %s", srr_file)

        # BioProject summary TSV
        bp_file = output_dir / "bioprojects_summary.tsv"
        fieldnames = ["BioProjectID", "Accession", "Title", "Organism", "SRR_Count", "ReleaseDate"]
        with open(bp_file, "w", newline="", encoding="utf-8") as fh:
            writer = csv.DictWriter(fh, fieldnames=fieldnames, delimiter="\t")
            writer.writeheader()
            for bp in sorted(bioprojects, key=lambda x: x["SRR_Count"], reverse=True):
                writer.writerow({
                    "BioProjectID": bp.get("BioProjectID", ""),
                    "Accession": bp.get("Accession", ""),
                    "Title": bp.get("Title", ""),
                    "Organism": bp.get("Organism", ""),
                    "SRR_Count": bp["SRR_Count"],
                    "ReleaseDate": bp.get("ReleaseDate", ""),
                })
        output_files["bioprojects"] = str(bp_file)
        logger.info("Saved BioProject summary to %s", bp_file)

        # Text report
        report_file = output_dir / "search_report.txt"
        stats = NCBISearcher._generate_statistics(bioprojects, srr_ids, label)
        with open(report_file, "w", encoding="utf-8") as fh:
            fh.write(stats)
        output_files["report"] = str(report_file)
        logger.info("Saved report to %s", report_file)

        return output_files

    @staticmethod
    def _generate_statistics(bioprojects: List[Dict], srr_ids: List[str], label: str) -> str:
        """Build a plain-text statistics report."""
        lines = [
            "=" * 60,
            "SRA RIBO-SEQ SEARCH REPORT",
            "=" * 60,
            "",
            f"Query Species: {label}",
            f"Search Time:  {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}",
            "",
            "Results:",
            f"  Total BioProjects found:    {len(bioprojects)}",
            f"  Total SRR IDs extracted:     {len(srr_ids)}",
            f"  BioProjects with SRA data:   "
            f"{sum(1 for bp in bioprojects if bp['SRR_Count'] > 0)}",
            f"  BioProjects without SRA:     "
            f"{sum(1 for bp in bioprojects if bp['SRR_Count'] == 0)}",
            "",
            "Top BioProjects (by SRR count):",
        ]
        sorted_bps = sorted(bioprojects, key=lambda x: x["SRR_Count"], reverse=True)[:10]
        for idx, bp in enumerate(sorted_bps, 1):
            title_preview = bp.get("Title", "")[:50]
            lines.append(f"  {idx}. {bp['BioProjectID']}: {bp['SRR_Count']} SRRs - {title_preview}...")
        lines.append("")
        lines.append("=" * 60)
        return "\n".join(lines)

    @staticmethod
    def _print_summary(stats: str, output_files: Dict[str, str]) -> None:
        """Print a short summary to stdout."""
        print("\n" + "=" * 60)
        print("SUMMARY")
        print("=" * 60)
        for line in stats.splitlines():
            if line.startswith("Query Species:") or line.startswith("Total "):
                print(f"  {line.strip()}")
        print("\nOutput Files:")
        for ftype, fpath in output_files.items():
            print(f"  - {ftype}: {fpath}")


# ============================================================================
# MODULE-LEVEL HELPERS
# ============================================================================

def _find_first_text(element: ET.Element, paths: List[str]) -> str:
    """Return stripped text from the first matching XPath, or empty string."""
    for path in paths:
        node = element.find(path)
        if node is not None and node.text:
            return node.text.strip()
    return ""


def parse_keywords(keywords_str: str) -> Dict[str, List[str]]:
    """
    Split a comma-separated keyword string into strict/medium buckets.

    Keywords containing ``ribo``, ``rpf``, ``footprint``, or ``protected``
    are treated as strict; everything else as medium.
    """
    strict_indicators = ("ribo", "rpf", "footprint", "protected")
    strict: List[str] = []
    medium: List[str] = []

    for kw in keywords_str.split(","):
        kw = kw.strip()
        if not kw:
            continue
        if any(tok in kw.lower() for tok in strict_indicators):
            strict.append(kw)
        else:
            medium.append(kw)

    return {"strict": strict, "medium": medium}


# ============================================================================
# CLI
# ============================================================================

def build_parser() -> "argparse.ArgumentParser":
    """Construct and return the argument parser."""
    parser = argparse.ArgumentParser(
        prog="search_sra_riboseq.py",
        description="Search NCBI for Ribo-seq data by species name or TaxID",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""\
Examples:
  # Search by species name
  %(prog)s --species "Arabidopsis thaliana" --output-dir results/

  # Search by TaxID
  %(prog)s --species 3702 --output-dir results/

  # Use custom keywords
  %(prog)s --species "Oryza sativa" \\
      --keywords "ribosome profiling,translatome,polysome" \\
      --output-dir results/

  # Restrict to strict tier only
  %(prog)s --species "Zea mays" --query-type strict --output-dir results/

  # Verbose logging to console
  %(prog)s --species "Homo sapiens" --output-dir results/ --log-level DEBUG
""",
    )

    parser.add_argument("--version", action="version", version=f"%(prog)s {__version__}")
    parser.add_argument(
        "--species", type=str, required=True,
        help='Species name (e.g., "Arabidopsis thaliana") or TaxID (e.g., 3702)',
    )
    parser.add_argument(
        "--output-dir", type=str, required=True,
        help="Output directory for results",
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
        "--keywords", type=str, default=None,
        help="Comma-separated custom keywords (overrides built-in lists)",
    )
    parser.add_argument(
        "--query-type", type=str,
        choices=["strict", "medium", "both"], default="both",
        help="Which keyword tier(s) to use (default: both)",
    )
    parser.add_argument(
        "--api-delay", type=float, default=0.34,
        help="Delay between NCBI API calls in seconds (default: 0.34)",
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

    # Logging setup
    logging.basicConfig(
        level=getattr(logging, args.log_level),
        format="[%(levelname)s] %(message)s",
    )

    # NCBI email (required)
    email = args.email or os.getenv("NCBI_EMAIL")
    if not email:
        parser.error(
            "NCBI email is required. Use --email or set the NCBI_EMAIL env var."
        )

    api_key = args.api_key or os.getenv("NCBI_API_KEY")

    # Custom keywords
    custom_keywords = None
    if args.keywords:
        custom_keywords = parse_keywords(args.keywords)
        logger.info("Custom strict keywords: %s", custom_keywords["strict"])
        logger.info("Custom medium keywords: %s", custom_keywords["medium"])

    searcher = NCBISearcher(
        email=email,
        api_key=api_key,
        api_delay=args.api_delay,
        custom_keywords=custom_keywords,
    )

    # Optionally filter tiers
    if args.query_type != "both":
        if args.query_type == "strict":
            searcher.medium_keywords = []
        else:
            searcher.strict_keywords = []

    try:
        searcher.search(args.species, Path(args.output_dir))
        logger.info("Done.")
    except Exception as exc:
        logger.critical("Search failed: %s", exc, exc_info=True)
        sys.exit(1)


if __name__ == "__main__":
    main()
