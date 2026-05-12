#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
NCBI SRA Ribo-seq Searcher - Standalone Script
================================================
Search NCBI for Ribo-seq data by species name or TaxID

This script searches NCBI BioProject database for Ribo-seq studies,
extracts linked SRR accessions, and generates summary reports.

Features:
- Search by species name or TaxID
- Two-tier keyword strategy (STRICT + MEDIUM)
- Automatic BioProject to SRA linking
- Custom keyword support
- Detailed search reports

Author: Auto-generated from ribo_metadata_workflow
Date: 2025
"""

import os
import re
import sys
import time
import csv
import argparse
from pathlib import Path
from typing import List, Dict, Optional, Set
from datetime import datetime
from xml.etree import ElementTree as ET

from tqdm import tqdm
from Bio import Entrez


# ============================================================================
# BUILT-IN KEYWORD CONFIGURATION
# ============================================================================

# Default Ribo-seq Keywords (STRICT - High Specificity)
DEFAULT_STRICT_KEYWORDS = [
    'ribosome profiling',
    'ribo-seq', 'riboseq', 'ribo seq',
    'ribosome footprinting', 'ribosome footprints',
    'ribosome protected fragments', 'ribosome-protected fragments',
    'RPF', 'RPF-seq'
]

# Default Translation-related Keywords (MEDIUM - Broader)
DEFAULT_MEDIUM_KEYWORDS = [
    'translatome', 'translatomics', 'translatomic profiling',
    'translatome sequencing', 'translation-level profiling',
    'actively translated mRNA',
    'polysome profiling', 'polysome-seq', 'polysome sequencing',
    'polyribosome profiling', 'polysome fractionation',
    'ribosome-associated mRNA', 'ribosome bound mRNA',
    'ribosome-associated transcriptome',
    'TE profiling', 'translation efficiency',
    'ribosome occupancy', 'ribosome release', 'ribosome pausing',
    'ribosome stalling', 'translation dynamics'
]


# ============================================================================
# CORE SEARCHER CLASS
# ============================================================================

class NCBISearcher:
    """Search NCBI for Ribo-seq data"""

    def __init__(
        self,
        email: str,
        api_key: Optional[str] = None,
        api_delay: float = 0.34,
        custom_keywords: Optional[Dict[str, List[str]]] = None,
        verbose: bool = True
    ):
        """
        Initialize NCBI Searcher

        Args:
            email: Email for NCBI API (required)
            api_key: NCBI API Key (optional)
            api_delay: Delay between API calls
            custom_keywords: Custom keywords dict with 'strict' and 'medium' keys
            verbose: Print progress information
        """
        Entrez.email = email
        if api_key:
            Entrez.api_key = api_key

        self.api_delay = api_delay
        self.verbose = verbose

        # Use custom keywords or default
        if custom_keywords:
            self.strict_keywords = custom_keywords.get('strict', DEFAULT_STRICT_KEYWORDS)
            self.medium_keywords = custom_keywords.get('medium', DEFAULT_MEDIUM_KEYWORDS)
        else:
            self.strict_keywords = DEFAULT_STRICT_KEYWORDS
            self.medium_keywords = DEFAULT_MEDIUM_KEYWORDS

    def species_to_taxid(self, species: str) -> Optional[int]:
        """
        Convert species name to NCBI TaxID

        Args:
            species: Species name (e.g., "Arabidopsis thaliana")

        Returns:
            TaxID as integer, or None if not found
        """
        try:
            if self.verbose:
                print(f"[INFO] Converting species name to TaxID: {species}")

            handle = Entrez.esearch(
                db="taxonomy",
                term=f'"{species}"[Scientific Name]',
                retmax=1
            )
            record = Entrez.read(handle)
            handle.close()
            time.sleep(self.api_delay)

            if record["IdList"]:
                taxid = int(record["IdList"][0])
                if self.verbose:
                    print(f"[INFO] Found TaxID: {taxid}")
                return taxid
            else:
                print(f"[WARN] Species not found: {species}")
                return None

        except Exception as e:
            print(f"[ERROR] Failed to convert species to TaxID: {e}")
            return None

    def detect_input_type(self, input_value: str) -> Dict[str, any]:
        """
        Detect if input is species name or TaxID

        Args:
            input_value: Species name or TaxID

        Returns:
            Dict with 'type' ('species' or 'taxid'), 'value', and 'taxid' (if applicable)
        """
        # Try to parse as integer (TaxID)
        try:
            taxid = int(input_value)
            return {
                'type': 'taxid',
                'value': input_value,
                'taxid': taxid,
                'label': f'taxid:{taxid}'
            }
        except ValueError:
            # It's a species name
            taxid = self.species_to_taxid(input_value)
            if taxid:
                return {
                    'type': 'species',
                    'value': input_value,
                    'taxid': taxid,
                    'label': f'{input_value}(txid={taxid})'
                }
            else:
                # Species not found, still proceed with label
                return {
                    'type': 'species',
                    'value': input_value,
                    'taxid': None,
                    'label': input_value
                }

    def build_queries(self, taxid: Optional[int]) -> List[str]:
        """
        Build search queries for NCBI

        Args:
            taxid: TaxID (optional)

        Returns:
            List of query strings
        """
        queries = []

        # Add organism filter if TaxID is available
        org_filter = f" AND txid{taxid}[Organism]" if taxid else ""

        # STRICT query (high specificity)
        strict_terms = " OR ".join([f'"{k}"[All Fields]' for k in self.strict_keywords])
        strict_query = f"({strict_terms}){org_filter}"
        queries.append(('STRICT', strict_query))

        # MEDIUM query (translation-related)
        medium_terms = " OR ".join([f'"{k}"[All Fields]' for k in self.medium_keywords])
        medium_query = f"({medium_terms}){org_filter}"
        queries.append(('MEDIUM', medium_query))

        return queries

    def esearch_all(
        self,
        database: str,
        query: str,
        retmax: int = 10000
    ) -> List[str]:
        """
        Perform esearch with pagination to get all results

        Args:
            database: NCBI database name
            query: Search query
            retmax: Maximum results per batch

        Returns:
            List of IDs
        """
        all_ids = []
        retstart = 0

        while True:
            try:
                handle = Entrez.esearch(
                    db=database,
                    term=query,
                    retmax=retmax,
                    retstart=retstart
                )
                data = Entrez.read(handle)
                handle.close()
                time.sleep(self.api_delay)

                batch_ids = data.get("IdList", [])

                if not batch_ids:
                    break

                all_ids.extend(batch_ids)

                if len(batch_ids) < retmax:
                    # Last batch
                    break

                retstart += retmax

            except Exception as e:
                print(f"[WARN] Search failed at retstart={retstart}: {e}")
                break

        return all_ids

    def fetch_bioproject_details(self, ids: List[str]) -> List[Dict]:
        """
        Fetch detailed information for BioProjects

        Args:
            ids: List of BioProject IDs

        Returns:
            List of BioProject records
        """
        bioprojects = []
        batch_size = 50

        if self.verbose:
            print(f"[INFO] Fetching details for {len(ids)} BioProjects...")

        for i in tqdm(range(0, len(ids), batch_size), desc="BioProject details"):
            batch = ids[i:i + batch_size]

            try:
                handle = Entrez.efetch(
                    db="bioproject",
                    id=",".join(batch),
                    rettype="xml",
                    retmode="xml"
                )
                xml_data = handle.read()
                handle.close()
                time.sleep(self.api_delay)

                # Parse XML
                if isinstance(xml_data, bytes):
                    xml_data = xml_data.decode('utf-8')

                root = ET.fromstring(xml_data)

                # Extract projects
                for project in root.findall(".//Project"):
                    record = self._parse_bioproject_xml(project)
                    if record:
                        bioprojects.append(record)

            except Exception as e:
                print(f"[WARN] Failed to fetch batch {i}: {e}")
                continue

        return bioprojects

    def _parse_bioproject_xml(self, project: ET.Element) -> Optional[Dict]:
        """
        Parse BioProject XML element

        Args:
            project: XML Element

        Returns:
            Dict with BioProject information
        """
        record = {}

        try:
            # Accession
            acc_elem = project.find(".//ArchiveID")
            if acc_elem is not None:
                record["Accession"] = acc_elem.attrib.get("accession", "")
            else:
                record["Accession"] = ""

            # BioProject ID
            pid_elem = project.find(".//ProjectID")
            if pid_elem is not None and pid_elem.text:
                record["BioProjectID"] = pid_elem.text
            else:
                record["BioProjectID"] = record["Accession"]

            # Title (try multiple paths)
            record["Title"] = self._find_first_text(project, [
                ".//Title",
                ".//ProjectDescr/Title",
                ".//Description/Title",
                ".//Name"
            ])

            # Organism
            record["Organism"] = self._find_first_text(project, [
                ".//Organism/OrganismName",
                ".//Organism",
                ".//Name"
            ])

            # Description
            record["Description"] = self._find_first_text(project, [
                ".//ProjectDescr/Description",
                ".//Description",
                ".//Comment"
            ])

            # ReleaseDate
            record["ReleaseDate"] = self._find_first_text(project, [
                ".//ProjectReleaseDate",
                ".//ReleaseDate",
                ".//Date"
            ])

            return record

        except Exception as e:
            print(f"[WARN] Failed to parse BioProject XML: {e}")
            return None

    def _find_first_text(self, element: ET.Element, paths: List[str]) -> str:
        """
        Find first non-empty text from multiple XML paths

        Args:
            element: XML Element
            paths: List of XPath expressions

        Returns:
            Text content or empty string
        """
        for path in paths:
            elem = element.find(path)
            if elem is not None and elem.text:
                return elem.text.strip()
        return ""

    def link_sra(self, bioproject_id: str) -> List[str]:
        """
        Link BioProject to SRA and extract SRR accessions

        Args:
            bioproject_id: BioProject ID (e.g., "PRJNA123456")

        Returns:
            Sorted list of SRR accessions
        """
        srr_ids = set()

        try:
            # Search in SRA database
            query = f"{bioproject_id}[BioProject]"
            handle = Entrez.esearch(db="sra", term=query, retmax=10000)
            search_data = Entrez.read(handle)
            handle.close()
            time.sleep(self.api_delay)

            uid_list = search_data.get("IdList", [])

            if not uid_list:
                return []

            # Fetch summaries in batches
            for i in range(0, len(uid_list), 100):
                batch_uids = uid_list[i:i+100]

                try:
                    handle = Entrez.esummary(db="sra", id=",".join(batch_uids))
                    summary_data = Entrez.read(handle)
                    handle.close()
                    time.sleep(self.api_delay)

                    # Extract SRR IDs using regex
                    for item in summary_data:
                        runs_field = str(item.get("Runs", ""))
                        srr_matches = re.findall(r'acc="(SRR\d+)"', runs_field)
                        srr_ids.update(srr_matches)

                except Exception as e:
                    print(f"[WARN] Failed to fetch summary batch {i}: {e}")
                    continue

        except Exception as e:
            print(f"[WARN] Failed to link SRA for {bioproject_id}: {e}")

        return sorted(list(srr_ids))

    def search(self, species_or_taxid: str, output_dir: Path) -> Dict:
        """
        Main search workflow

        Args:
            species_or_taxid: Species name or TaxID
            output_dir: Output directory path

        Returns:
            Dict with results and statistics
        """
        # Step 1: Detect input type
        print("\n" + "=" * 60)
        print("NCBI SRA RIBO-SEQ SEARCH")
        print("=" * 60)

        input_info = self.detect_input_type(species_or_taxid)

        if self.verbose:
            print(f"\n[INFO] Input type: {input_info['type']}")
            print(f"[INFO] Query: {input_info['label']}")

        # Step 2: Build queries
        print("\n" + "-" * 60)
        print("[INFO] Building search queries...")
        queries = self.build_queries(input_info['taxid'])

        for query_type, query in queries:
            print(f"  - {query_type}: {len(self.strict_keywords if query_type == 'STRICT' else self.medium_keywords)} keywords")

        # Step 3: Search BioProjects
        print("\n" + "-" * 60)
        print("[INFO] Searching BioProjects...")

        all_bioproject_ids = set()

        for query_type, query in queries:
            print(f"\n[INFO] Running {query_type} query...")
            ids = self.esearch_all("bioproject", query)
            print(f"[INFO] {query_type}: Found {len(ids)} BioProjects")
            all_bioproject_ids.update(ids)

        print(f"\n[INFO] Total unique BioProjects: {len(all_bioproject_ids)}")

        if not all_bioproject_ids:
            print("[ERROR] No BioProjects found!")
            return {'error': 'No BioProjects found'}

        # Step 4: Fetch BioProject details
        print("\n" + "-" * 60)
        bioprojects = self.fetch_bioproject_details(list(all_bioproject_ids))
        print(f"[INFO] Retrieved details for {len(bioprojects)} BioProjects")

        # Step 5: Link to SRA
        print("\n" + "-" * 60)
        print("[INFO] Linking to SRA database...")

        all_srr_ids = []

        for bp in tqdm(bioprojects, desc="Linking SRA"):
            bp["QuerySpecies"] = input_info['label']
            srr_ids = self.link_sra(bp["BioProjectID"])
            bp["LinkedSRA"] = srr_ids
            bp["SRR_Count"] = len(srr_ids)
            all_srr_ids.extend(srr_ids)

        print(f"\n[INFO] Total SRR IDs extracted: {len(all_srr_ids)}")

        # Step 6: Save results
        print("\n" + "-" * 60)
        output_dir = Path(output_dir)
        output_dir.mkdir(parents=True, exist_ok=True)

        output_files = self._save_results(
            bioprojects,
            all_srr_ids,
            input_info['label'],
            output_dir
        )

        # Step 7: Generate statistics
        stats = self._generate_statistics(
            bioprojects,
            all_srr_ids,
            input_info['label']
        )

        # Step 8: Print summary
        self._print_summary(stats, output_files)

        return {
            'bioprojects': bioprojects,
            'srr_ids': all_srr_ids,
            'statistics': stats,
            'output_files': output_files
        }

    def _save_results(
        self,
        bioprojects: List[Dict],
        srr_ids: List[str],
        label: str,
        output_dir: Path
    ) -> Dict[str, str]:
        """
        Save search results to files

        Args:
            bioprojects: List of BioProject records
            srr_ids: List of SRR IDs
            label: Query label
            output_dir: Output directory

        Returns:
            Dict mapping file type to file path
        """
        output_files = {}
        safe_label = label.replace("(", "_").replace(")", "_").replace(" ", "_").replace("/", "_")

        # 1. Save SRR IDs
        srr_file = output_dir / "srr_ids.txt"
        with open(srr_file, 'w') as f:
            for srr_id in srr_ids:
                f.write(f"{srr_id}\n")
        output_files['srr_ids'] = str(srr_file)
        print(f"[INFO] Saved: {srr_file}")

        # 2. Save BioProject summary
        bp_file = output_dir / "bioprojects_summary.tsv"
        with open(bp_file, 'w', newline='', encoding='utf-8') as f:
            fieldnames = ['BioProjectID', 'Accession', 'Title', 'Organism', 'SRR_Count', 'ReleaseDate']
            writer = csv.DictWriter(f, fieldnames=fieldnames, delimiter='\t')
            writer.writeheader()

            for bp in sorted(bioprojects, key=lambda x: x['SRR_Count'], reverse=True):
                writer.writerow({
                    'BioProjectID': bp.get('BioProjectID', ''),
                    'Accession': bp.get('Accession', ''),
                    'Title': bp.get('Title', ''),
                    'Organism': bp.get('Organism', ''),
                    'SRR_Count': bp['SRR_Count'],
                    'ReleaseDate': bp.get('ReleaseDate', '')
                })
        output_files['bioprojects'] = str(bp_file)
        print(f"[INFO] Saved: {bp_file}")

        # 3. Save detailed report
        report_file = output_dir / "search_report.txt"
        stats = self._generate_statistics(bioprojects, srr_ids, label)

        with open(report_file, 'w', encoding='utf-8') as f:
            f.write(stats)
        output_files['report'] = str(report_file)
        print(f"[INFO] Saved: {report_file}")

        return output_files

    def _generate_statistics(
        self,
        bioprojects: List[Dict],
        srr_ids: List[str],
        label: str
    ) -> str:
        """
        Generate search statistics report

        Args:
            bioprojects: List of BioProject records
            srr_ids: List of SRR IDs
            label: Query label

        Returns:
            Statistics report as string
        """
        lines = []
        lines.append("=" * 60)
        lines.append("SRA RIBO-SEQ SEARCH REPORT")
        lines.append("=" * 60)
        lines.append("")
        lines.append(f"Query Species: {label}")
        lines.append(f"Search Time: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
        lines.append("")

        lines.append("Query Strategy:")
        lines.append(f"  - STRICT: {len(self.strict_keywords)} keywords (high specificity)")
        lines.append(f"  - MEDIUM: {len(self.medium_keywords)} keywords (translation-related)")
        lines.append("")

        lines.append("Results:")
        lines.append(f"  - Total BioProjects found: {len(bioprojects)}")
        lines.append(f"  - Total SRR IDs extracted: {len(srr_ids)}")
        lines.append(f"  - BioProjects with SRA data: {sum(1 for bp in bioprojects if bp['SRR_Count'] > 0)}")
        lines.append(f"  - BioProjects without SRA: {sum(1 for bp in bioprojects if bp['SRR_Count'] == 0)}")
        lines.append("")

        # Top BioProjects
        lines.append("Top BioProjects (by SRR count):")
        sorted_bps = sorted(bioprojects, key=lambda x: x['SRR_Count'], reverse=True)[:10]
        for i, bp in enumerate(sorted_bps, 1):
            lines.append(f"  {i}. {bp['BioProjectID']}: {bp['SRR_Count']} SRRs - {bp['Title'][:50]}...")
        lines.append("")

        lines.append("=" * 60)

        return "\n".join(lines)

    def _print_summary(self, stats: str, output_files: Dict[str, str]):
        """Print search summary"""
        print("\n" + "=" * 60)
        print("SUMMARY")
        print("=" * 60)

        # Parse statistics string
        newline = '\n'
        query_species = stats.split('Query Species: ')[1].split(newline)[0]
        total_bp = stats.split('Total BioProjects found: ')[1].split(newline)[0]
        total_srr = stats.split('Total SRR IDs extracted: ')[1].split(newline)[0]

        print(f"Query Species: {query_species}")
        print(f"Total BioProjects: {total_bp}")
        print(f"Total SRR IDs: {total_srr}")
        print(f"\nOutput Files:")
        for file_type, file_path in output_files.items():
            print(f"  - {file_type}: {file_path}")


# ============================================================================
# COMMAND LINE INTERFACE
# ============================================================================

def parse_keywords(keywords_str: str) -> Dict[str, List[str]]:
    """
    Parse custom keywords from command line string

    Args:
        keywords_str: Comma-separated keywords

    Returns:
        Dict with 'strict' and 'medium' keys
    """
    keywords = [k.strip() for k in keywords_str.split(',')]

    # Split into strict and medium (simple heuristic)
    # If keyword contains "ribo", "rpf", "footprint", "profiling" → strict
    # Otherwise → medium
    strict = []
    medium = []

    for kw in keywords:
        kw_lower = kw.lower()
        if any(x in kw_lower for x in ['ribo', 'rpf', 'footprint', 'protected']):
            strict.append(kw)
        else:
            medium.append(kw)

    return {'strict': strict, 'medium': medium}


def main():
    """Main function"""
    parser = argparse.ArgumentParser(
        description="Search NCBI for Ribo-seq data by species name or TaxID",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Search by species name
  python search_sra_riboseq.py --species "Arabidopsis thaliana" --output-dir results/

  # Search by TaxID
  python search_sra_riboseq.py --species 3702 --output-dir results/

  # Use custom keywords
  python search_sra_riboseq.py \\
    --species "Oryza sativa" \\
    --keywords "ribosome profiling,translatome,polysome" \\
    --output-dir results/

  # Specify query type
  python search_sra_riboseq.py \\
    --species "Zea mays" \\
    --query-type strict \\
    --output-dir results/

Workflow:
  1. Search for SRR IDs
     python search_sra_riboseq.py --species "Arabidopsis" --output-dir search/

  2. Fetch metadata
     python batch_fetch_sra_metadata.py \\
       --input-file search/srr_ids.txt \\
       --output-dir fetch/

  3. Format metadata
     python format_sra_metadata.py \\
       --input fetch/metadata_comprehensive.tsv \\
       --output fetch/samples_summary.csv
        """
    )

    parser.add_argument(
        '--species',
        type=str,
        required=True,
        help='Species name (e.g., "Arabidopsis thaliana") or TaxID (e.g., 3702)'
    )

    parser.add_argument(
        '--output-dir',
        type=str,
        required=True,
        help='Output directory'
    )

    parser.add_argument(
        '--email',
        type=str,
        default=None,
        help='NCBI API email (default: from NCBI_EMAIL env var)'
    )

    parser.add_argument(
        '--api-key',
        type=str,
        default=None,
        help='NCBI API key (default: from NCBI_API_KEY env var)'
    )

    parser.add_argument(
        '--keywords',
        type=str,
        default=None,
        help='Custom keywords (comma-separated, overrides built-in)'
    )

    parser.add_argument(
        '--query-type',
        type=str,
        choices=['strict', 'medium', 'both'],
        default='both',
        help='Query type (default: both)'
    )

    parser.add_argument(
        '--api-delay',
        type=float,
        default=0.34,
        help='Delay between API calls in seconds (default: 0.34)'
    )

    args = parser.parse_args()

    # Get email with priority: CLI > env var > error
    email = args.email or os.getenv('NCBI_EMAIL')
    if not email:
        print("[ERROR] NCBI email is required!")
        print("[INFO] Set NCBI_EMAIL environment variable or use --email argument")
        sys.exit(1)

    # Get API key
    api_key = args.api_key or os.getenv('NCBI_API_KEY')

    # Parse custom keywords if provided
    custom_keywords = None
    if args.keywords:
        custom_keywords = parse_keywords(args.keywords)
        print(f"[INFO] Using custom keywords:")
        print(f"  - Strict: {custom_keywords['strict']}")
        print(f"  - Medium: {custom_keywords['medium']}")

    # Initialize searcher
    searcher = NCBISearcher(
        email=email,
        api_key=api_key,
        api_delay=args.api_delay,
        custom_keywords=custom_keywords
    )

    # Run search
    try:
        result = searcher.search(
            species_or_taxid=args.species,
            output_dir=Path(args.output_dir)
        )
        print("\n[INFO] Done!")
        sys.exit(0)

    except Exception as e:
        print(f"\n[ERROR] Search failed: {e}")
        import traceback
        traceback.print_exc()
        sys.exit(1)


if __name__ == "__main__":
    main()
