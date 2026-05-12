#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
SRA Metadata Fetcher - Standalone Script
==========================================
Fetch SRA metadata for batch SRA Run IDs

This script fetches three types of metadata from NCBI:
1. RunInfo basic information (Run, Experiment, BioSample, BioProject IDs, Library strategy, etc.)
2. BioSample detailed attributes (tissue, treatment, developmental stage, cultivar, etc.)
3. Experiment detailed information (library design, platform information, etc.)

Output format: TSV (tab-separated values)

Author: Auto-generated from ribo_metadata_workflow
Date: 2025
"""

import os
import sys
import time
import csv
import argparse
import logging
from pathlib import Path
from typing import List, Dict, Optional
from io import StringIO

import pandas as pd
from tqdm import tqdm
from Bio import Entrez
from xml.etree import ElementTree as ET


class SRAMetadataFetcher:
    """Independent SRA Metadata Fetcher class"""

    def __init__(
        self,
        email: str,
        api_key: Optional[str] = None,
        api_delay: float = 0.34,
        batch_size: int = 500
    ):
        """
        Initialize SRA Metadata Fetcher

        Args:
            email: Email for NCBI API (required)
            api_key: NCBI API Key (optional, increases rate limits)
            api_delay: Delay between API calls in seconds
            batch_size: Number of accessions per batch
        """
        Entrez.email = email
        if api_key:
            Entrez.api_key = api_key

        self.api_delay = api_delay
        self.batch_size = batch_size
        self.failed_srr_ids = []

    def fetch_runinfo(self, srr_ids: List[str]) -> pd.DataFrame:
        """
        Batch fetch SRA RunInfo for SRR accessions

        Args:
            srr_ids: List of SRR accessions (e.g., ['SRR123456', 'SRR789012'])

        Returns:
            DataFrame with RunInfo data
        """
        if not srr_ids:
            return pd.DataFrame()

        all_runinfo = []

        print(f"[INFO] Fetching RunInfo for {len(srr_ids)} SRA runs...")

        for i in tqdm(range(0, len(srr_ids), self.batch_size), desc="RunInfo"):
            batch = srr_ids[i:i + self.batch_size]

            try:
                # STEP 1: Convert SRR to UIDs using esearch
                search_terms = " OR ".join([f"{acc}[Accession]" for acc in batch])
                handle = Entrez.esearch(db="sra", term=search_terms, retmax=len(batch))
                search_data = Entrez.read(handle)
                handle.close()
                time.sleep(self.api_delay)

                uids = search_data.get("IdList", [])
                if not uids:
                    print(f"[WARN] No UIDs found for batch {i}-{i+len(batch)}")
                    self.failed_srr_ids.extend(batch)
                    continue

                # STEP 2: Fetch runinfo using efetch (CSV format)
                handle = Entrez.efetch(
                    db="sra",
                    id=",".join(uids),
                    rettype="runinfo",
                    retmode="text"
                )
                raw_data = handle.read()
                handle.close()
                time.sleep(self.api_delay)

                # STEP 3: Parse CSV
                runinfo_text = raw_data.decode('utf-8') if isinstance(raw_data, bytes) else raw_data

                if runinfo_text.strip():
                    reader = csv.DictReader(StringIO(runinfo_text))
                    batch_rows = list(reader)
                    for row in batch_rows:
                        all_runinfo.append(dict(row))

            except Exception as e:
                print(f"[WARN] Batch {i}-{i+len(batch)} failed: {e}")
                self.failed_srr_ids.extend(batch)
                continue

        if all_runinfo:
            df = pd.DataFrame(all_runinfo)
            print(f"[INFO] Total runs retrieved: {len(df)}")
            return df

        return pd.DataFrame()

    def fetch_biosample(self, biosample_ids: List[str]) -> Dict[str, Dict]:
        """
        Fetch BioSample metadata with robust error handling

        Args:
            biosample_ids: List of BioSample accessions (SAMN*)

        Returns:
            Dict mapping biosample_id -> {attribute_name: value}
        """
        metadata = {}

        print(f"[INFO] Fetching BioSample metadata for {len(biosample_ids)} samples...")

        for i in tqdm(range(0, len(biosample_ids), self.batch_size), desc="BioSample"):
            batch = biosample_ids[i:i + self.batch_size]

            try:
                handle = Entrez.efetch(
                    db="biosample",
                    id=",".join(batch),
                    rettype="xml",
                    retmode="xml"
                )
                xml_data = handle.read()
                handle.close()
                time.sleep(self.api_delay)

            except Exception as e:
                print(f"[WARN] BioSample fetch error for batch {i}: {e}")
                continue

            # Parse XML
            try:
                if isinstance(xml_data, bytes):
                    xml_data = xml_data.decode('utf-8')

                root = ET.fromstring(xml_data)
                samples = root.findall(".//BioSample") or root.findall(".//Sample") or []

                if not samples:
                    print(f"[WARN] No BioSample elements found in batch {i}")
                    continue

                for sample in samples:
                    # Get accession - try multiple strategies
                    bio_id = sample.attrib.get("accession", "")

                    if not bio_id:
                        acc = sample.find(".//Accession") or sample.find(".//ID") or sample.find(".//primary_id")
                        if acc is not None and acc.text:
                            bio_id = acc.text.strip()

                    if not bio_id:
                        for id_elem in sample.findall(".//Id"):
                            if id_elem.attrib.get("is_primary") == "1" and id_elem.text:
                                bio_id = id_elem.text.strip()
                                break

                    if not bio_id:
                        continue

                    # Extract ALL attributes
                    attrs = {}

                    # Strategy 1: Attributes with attribute_name attribute
                    for attr in sample.findall(".//Attribute"):
                        attr_name = attr.attrib.get("attribute_name", "")
                        if not attr_name:
                            attr_name = attr.attrib.get("harmonized_name", "")
                        attr_val = attr.text if attr.text else ""
                        if attr_name:
                            attrs[attr_name] = attr_val

                    # Strategy 2: Description/Title fields
                    title = sample.find(".//Title")
                    if title is not None and title.text:
                        attrs['title'] = title.text

                    organism = sample.find(".//OrganismName")
                    if organism is not None and organism.text:
                        attrs['organism'] = organism.text

                    if attrs:
                        metadata[bio_id] = attrs

            except ET.ParseError as e:
                print(f"[WARN] XML parse error in batch {i}: {e}")
                continue
            except Exception as e:
                print(f"[WARN] Unexpected error parsing batch {i}: {e}")
                continue

        print(f"[INFO] Successfully fetched metadata for {len(metadata)}/{len(biosample_ids)} BioSamples")
        return metadata

    def fetch_experiment(self, experiment_ids: List[str]) -> Dict[str, Dict]:
        """
        Fetch SRA Experiment metadata with robust error handling

        Args:
            experiment_ids: List of Experiment accessions (SRX*)

        Returns:
            Dict mapping experiment_id -> {attribute_name: value}
        """
        metadata = {}

        print(f"[INFO] Fetching metadata for {len(experiment_ids)} experiments...")

        for i in tqdm(range(0, len(experiment_ids), self.batch_size), desc="Experiment"):
            batch = experiment_ids[i:i + self.batch_size]

            try:
                # STEP 1: Convert SRX to UIDs
                search_terms = " OR ".join([f"{acc}[Accession]" for acc in batch])
                handle = Entrez.esearch(db="sra", term=search_terms, retmax=len(batch))
                search_data = Entrez.read(handle)
                handle.close()
                time.sleep(self.api_delay)

                uids = search_data.get("IdList", [])
                if not uids:
                    print(f"[WARN] No UIDs found for experiment batch {i}-{i+len(batch)}")
                    continue

                # STEP 2: Fetch metadata via esummary
                handle = Entrez.esummary(db="sra", id=",".join(uids))
                summary_data = Entrez.read(handle)
                handle.close()
                time.sleep(self.api_delay)

                # Parse each experiment summary
                for exp_summary in summary_data:
                    exp_id = None

                    # Strategy 1: Try direct 'experiment' field
                    if 'experiment' in exp_summary:
                        exp_id = exp_summary['experiment']

                    # Strategy 2: Extract from runs array
                    if not exp_id and 'runs' in exp_summary:
                        runs = exp_summary['runs']
                        if runs and len(runs) > 0:
                            if isinstance(runs[0], dict):
                                exp_id = runs[0].get('experiment', '')

                    if not exp_id:
                        continue

                    attrs = {}

                    # Extract experiment title
                    if 'title' in exp_summary:
                        attrs['experiment_title'] = exp_summary['title']

                    # Extract library design information
                    if 'library' in exp_summary:
                        lib = exp_summary['library']
                        if isinstance(lib, dict):
                            if 'strategy' in lib:
                                attrs['library_strategy'] = lib['strategy']
                            if 'source' in lib:
                                attrs['library_source'] = lib['source']
                            if 'selection' in lib:
                                attrs['library_selection'] = lib['selection']
                            if 'layout' in lib:
                                attrs['library_layout'] = lib['layout']

                    # Extract platform information
                    if 'platform' in exp_summary:
                        plat = exp_summary['platform']
                        if isinstance(plat, dict):
                            if 'name' in plat:
                                attrs['platform'] = plat['name']
                            if 'instrument' in plat:
                                attrs['instrument_model'] = plat['instrument']

                    # Extract sample link
                    if 'sample' in exp_summary:
                        attrs['linked_sample'] = exp_summary['sample']

                    if exp_id and attrs:
                        metadata[exp_id] = attrs

            except Exception as e:
                print(f"[WARN] Failed to fetch experiment metadata for batch {i}: {e}")
                continue

        print(f"[INFO] Successfully fetched metadata for {len(metadata)}/{len(experiment_ids)} experiments")
        return metadata

    def merge_metadata(
        self,
        runinfo_df: pd.DataFrame,
        biosample_metadata: Dict[str, Dict],
        experiment_metadata: Dict[str, Dict]
    ) -> pd.DataFrame:
        """
        Merge BioSample and Experiment metadata into RunInfo DataFrame

        Args:
            runinfo_df: RunInfo DataFrame
            biosample_metadata: Dict mapping BioSample ID -> attributes
            experiment_metadata: Dict mapping Experiment ID -> attributes

        Returns:
            Updated RunInfo DataFrame
        """
        print(f"[INFO] Merging all metadata into RunInfo...")

        # 1. Merge BioSample metadata
        if biosample_metadata and 'BioSample' in runinfo_df.columns:
            for bio_id, attrs in biosample_metadata.items():
                mask = runinfo_df['BioSample'] == bio_id
                if not mask.any():
                    continue

                for attr_name, attr_val in attrs.items():
                    col_name = 'biosample_' + attr_name.replace(' ', '_').replace('-', '_').replace('/', '_')
                    col_name = ''.join(c if c.isalnum() or c == '_' else '_' for c in col_name)

                    if col_name not in runinfo_df.columns:
                        runinfo_df[col_name] = pd.Series(dtype=object)

                    runinfo_df.loc[mask, col_name] = attr_val

            print(f"[INFO] BioSample coverage: {runinfo_df['BioSample'].nunique()} unique samples")

        # 2. Merge Experiment metadata
        if experiment_metadata and 'Experiment' in runinfo_df.columns:
            for exp_id, attrs in experiment_metadata.items():
                mask = runinfo_df['Experiment'] == exp_id
                if not mask.any():
                    continue

                for attr_name, attr_val in attrs.items():
                    col_name = 'experiment_' + attr_name.replace(' ', '_').replace('-', '_').replace('/', '_')
                    col_name = ''.join(c if c.isalnum() or c == '_' else '_' for c in col_name)

                    if col_name not in runinfo_df.columns:
                        runinfo_df[col_name] = pd.Series(dtype=object)

                    runinfo_df.loc[mask, col_name] = attr_val

            print(f"[INFO] Experiment coverage: {runinfo_df['Experiment'].nunique()} unique experiments")

        return runinfo_df


def read_srr_ids_from_file(file_path: str) -> List[str]:
    """
    Read SRR IDs from a text file

    Args:
        file_path: Path to text file (one SRR ID per line)

    Returns:
        List of SRR IDs
    """
    srr_ids = []
    with open(file_path, 'r') as f:
        for line in f:
            line = line.strip()
            if line and not line.startswith('#'):
                srr_ids.append(line)
    return srr_ids


def save_metadata_tsv(
    output_dir: Path,
    runinfo_df: pd.DataFrame,
    biosample_metadata: Dict[str, Dict],
    experiment_metadata: Dict[str, Dict]
) -> Dict[str, str]:
    """
    Save comprehensive metadata to a single TSV file

    Args:
        output_dir: Output directory path
        runinfo_df: RunInfo DataFrame (already merged with BioSample and Experiment metadata)
        biosample_metadata: BioSample metadata dict (unused, kept for compatibility)
        experiment_metadata: Experiment metadata dict (unused, kept for compatibility)

    Returns:
        Dict mapping file type to file path
    """
    output_dir.mkdir(parents=True, exist_ok=True)

    output_files = {}

    # Save comprehensive metadata file (all information in one file)
    comprehensive_file = output_dir / "metadata_comprehensive.tsv"

    # Replace NaN with "NA" for better readability
    runinfo_df_filled = runinfo_df.fillna('NA')

    # Save to TSV
    runinfo_df_filled.to_csv(comprehensive_file, sep='\t', index=False, encoding='utf-8')
    output_files['comprehensive'] = str(comprehensive_file)
    print(f"[INFO] Saved comprehensive metadata file: {comprehensive_file}")

    return output_files


def main():
    """Main function"""
    parser = argparse.ArgumentParser(
        description="Fetch SRA metadata for batch SRA Run IDs",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # From file
  python fetch_sra_metadata.py --input-file srr_list.txt --output-dir results/

  # From command line
  python fetch_sra_metadata.py --srr-ids SRR123456,SRR789012 --output-dir results/

  # With custom email and API key
  python fetch_sra_metadata.py --input-file srr_list.txt --email user@example.com --api-key YOUR_KEY

  # Set environment variables (recommended)
  export NCBI_EMAIL="your-email@example.com"
  export NCBI_API_KEY="your-api-key"
  python fetch_sra_metadata.py --input-file srr_list.txt
        """
    )

    # Input arguments
    input_group = parser.add_mutually_exclusive_group(required=True)
    input_group.add_argument('--input-file', type=str, help='Text file with SRR IDs (one per line)')
    input_group.add_argument('--srr-ids', type=str, help='Comma-separated SRR IDs')

    # Output arguments
    parser.add_argument('--output-dir', type=str, default='.', help='Output directory (default: current directory)')

    # NCBI API configuration
    parser.add_argument('--email', type=str, default=None,
                        help='NCBI API email (default: from NCBI_EMAIL env var or your-email@example.com)')
    parser.add_argument('--api-key', type=str, default=None,
                        help='NCBI API key (default: from NCBI_API_KEY env var)')

    # Performance arguments
    parser.add_argument('--batch-size', type=int, default=500,
                        help='Number of accessions per batch (default: 500)')
    parser.add_argument('--api-delay', type=float, default=0.34,
                        help='Delay between API calls in seconds (default: 0.34)')

    args = parser.parse_args()

    # Get email with priority: CLI > env var > default
    email = args.email or os.getenv('NCBI_EMAIL') or 'your-email@example.com'
    if email == 'your-email@example.com':
        print("[WARN] Using default email. Please set NCBI_EMAIL environment variable or use --email argument.")

    # Get API key with priority: CLI > env var > None
    api_key = args.api_key or os.getenv('NCBI_API_KEY')

    # Get SRR ID list
    if args.input_file:
        if not Path(args.input_file).exists():
            print(f"[ERROR] Input file not found: {args.input_file}")
            sys.exit(1)
        srr_ids = read_srr_ids_from_file(args.input_file)
        print(f"[INFO] Loaded {len(srr_ids)} SRR IDs from {args.input_file}")
    else:
        srr_ids = [s.strip() for s in args.srr_ids.split(',')]
        print(f"[INFO] Processing {len(srr_ids)} SRR IDs from command line")

    if not srr_ids:
        print("[ERROR] No SRR IDs provided")
        sys.exit(1)

    # Initialize fetcher
    print(f"[INFO] Initializing SRA Metadata Fetcher...")
    print(f"[INFO] Email: {email}")
    if api_key:
        print(f"[INFO] API Key: {api_key[:10]}...")
    fetcher = SRAMetadataFetcher(
        email=email,
        api_key=api_key,
        api_delay=args.api_delay,
        batch_size=args.batch_size
    )

    # Step 1: Fetch RunInfo
    print("\n" + "="*60)
    print("STEP 1: Fetching RunInfo...")
    print("="*60)
    runinfo_df = fetcher.fetch_runinfo(srr_ids)

    if runinfo_df.empty:
        print("[ERROR] No RunInfo data retrieved. Please check your SRR IDs.")
        sys.exit(1)

    # Step 2: Extract BioSample and Experiment IDs
    biosample_ids = []
    if 'BioSample' in runinfo_df.columns:
        biosample_ids = runinfo_df['BioSample'].dropna().unique().tolist()
        biosample_ids = [bid for bid in biosample_ids if bid and str(bid).startswith('SAMN')]

    experiment_ids = []
    if 'Experiment' in runinfo_df.columns:
        experiment_ids = runinfo_df['Experiment'].dropna().unique().tolist()
        experiment_ids = [eid for eid in experiment_ids if eid and str(eid).startswith('SRX')]

    print(f"\n[INFO] Found {len(biosample_ids)} unique BioSamples")
    print(f"[INFO] Found {len(experiment_ids)} unique Experiments")

    # Step 3: Fetch BioSample metadata
    biosample_metadata = {}
    if biosample_ids:
        print("\n" + "="*60)
        print("STEP 2: Fetching BioSample metadata...")
        print("="*60)
        biosample_metadata = fetcher.fetch_biosample(biosample_ids)

    # Step 4: Fetch Experiment metadata
    experiment_metadata = {}
    if experiment_ids:
        print("\n" + "="*60)
        print("STEP 3: Fetching Experiment metadata...")
        print("="*60)
        experiment_metadata = fetcher.fetch_experiment(experiment_ids)

    # Step 5: Merge metadata
    print("\n" + "="*60)
    print("STEP 4: Merging metadata...")
    print("="*60)
    runinfo_df = fetcher.merge_metadata(runinfo_df, biosample_metadata, experiment_metadata)

    # Step 6: Save results
    print("\n" + "="*60)
    print("STEP 5: Saving results...")
    print("="*60)
    output_dir = Path(args.output_dir)
    output_files = save_metadata_tsv(output_dir, runinfo_df, biosample_metadata, experiment_metadata)

    # Print summary
    print("\n" + "="*60)
    print("SUMMARY")
    print("="*60)
    print(f"Total SRR IDs processed: {len(srr_ids)}")
    print(f"Total runs retrieved: {len(runinfo_df)}")
    print(f"BioSamples with metadata: {len(biosample_metadata)}/{len(biosample_ids)}")
    print(f"Experiments with metadata: {len(experiment_metadata)}/{len(experiment_ids)}")
    print(f"\nOutput file:")
    for file_type, file_path in output_files.items():
        if file_type == 'comprehensive':
            print(f"  ✓ Comprehensive metadata: {file_path}")
            print(f"    (Contains all RunInfo, BioSample, and Experiment metadata)")
        else:
            print(f"  - {file_type}: {file_path}")

    if fetcher.failed_srr_ids:
        print(f"\n[WARN] Failed to retrieve {len(fetcher.failed_srr_ids)} SRR IDs")
        print(f"[INFO] Check log for details")

    print("\n[INFO] Done!")


if __name__ == "__main__":
    main()
