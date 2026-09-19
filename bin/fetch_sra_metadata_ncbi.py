#!/usr/bin/env python3
# ==============================================================================
# File Name:    fetch_sra_metadata_ncbi.py
# Author:       ChengYu
# Description:  Fetch SRA metadata from NCBI by batch-querying SRR accession
#               IDs via the Entrez E-Utilities API, parsing the returned
#               EXPERIMENT_PACKAGE XML, and extracting sample attributes into
#               a structured TSV/CSV report.
# Created Time: 2026
# ==============================================================================

from __future__ import annotations

import argparse
import csv
import logging
import os
import sys
import time
import xml.etree.ElementTree as ET
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path
from typing import Dict, List, Optional, Set, TextIO
from urllib.error import HTTPError, URLError
from urllib.parse import urlencode
from urllib.request import Request, urlopen

__version__ = "2.0.0"

logger = logging.getLogger("fetch_sra_metadata_ncbi")

# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------
DEFAULT_BATCH_SIZE = 50
DEFAULT_THREADS = 1
DEFAULT_RATE_LIMIT = 0.5  # seconds between API calls
NCBI_EFETCH_URL = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
NCBI_ESRCH_URL = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi"

# Fields that are always extracted
CORE_FIELDS = [
    "accession",
    "experiment_accession",
    "sample_accession",
    "study_accession",
    "run_alias",
    "sample_alias",
    "experiment_alias",
    "instrument_model",
    "library_strategy",
    "library_source",
    "library_selection",
    "library_layout",
    "spots",
    "bases",
    "avg_length",
    "biosample",
    "bioproject",
    "center_name",
    "tax_id",
    "scientific_name",
]

# Default metadata output fields
DEFAULT_OUTPUT_FIELDS = CORE_FIELDS + ["sample_attributes"]


# ---------------------------------------------------------------------------
# NCBI API helpers
# ---------------------------------------------------------------------------
def esearch_sra(
    term: str,
    email: str,
    retmax: int = 10000,
    api_key: Optional[str] = None,
    timeout: int = 30,
) -> List[str]:
    """Search NCBI SRA via ESearch and return a list of accession IDs.

    Parameters
    ----------
    term : str
        Entrez search term, e.g. a study accession or free text.
    email : str
        Email address (required by NCBI policy).
    retmax : int
        Maximum number of UIDs to return.
    api_key : str or None
        NCBI API key for higher rate limits.
    timeout : int
        HTTP request timeout in seconds.

    Returns
    -------
    list[str]
        List of SRA Run accession strings.

    Raises
    ------
    RuntimeError
        If the ESearch response cannot be parsed.
    """
    params: Dict[str, str] = {
        "db": "sra",
        "term": term,
        "retmax": str(retmax),
        "usehistory": "n",
        "retmode": "json",
        "email": email,
        "tool": "fetch_sra_metadata_ncbi",
    }
    if api_key:
        params["api_key"] = api_key

    url = f"{NCBI_ESRCH_URL}?{urlencode(params)}"
    logger.info("ESearch: %s", term)

    import json as _json

    req = Request(url)
    req.add_header("User-Agent", "fetch_sra_metadata_ncbi/2.0")
    try:
        with urlopen(req, timeout=timeout) as resp:
            data = _json.loads(resp.read().decode())
    except (URLError, _json.JSONDecodeError) as exc:
        raise RuntimeError(f"ESearch failed: {exc}") from exc

    ids = data.get("esearchresult", {}).get("idlist", [])
    logger.info("ESearch returned %d UID(s)", len(ids))
    return ids


def efetch_sra_xml(
    run_ids: List[str],
    email: str,
    api_key: Optional[str] = None,
    timeout: int = 120,
) -> str:
    """Fetch SRA run metadata as XML via EFetch.

    Parameters
    ----------
    run_ids : list[str]
        List of SRA Run accessions or UIDs.
    email : str
        Email for NCBI policy.
    api_key : str or None
        Optional NCBI API key.
    timeout : int
        HTTP request timeout in seconds.

    Returns
    -------
    str
        Raw XML response string.

    Raises
    ------
    RuntimeError
        If the request fails.
    """
    params: Dict[str, str] = {
        "db": "sra",
        "id": ",".join(run_ids),
        "rettype": "full",
        "retmode": "xml",
        "email": email,
        "tool": "fetch_sra_metadata_ncbi",
    }
    if api_key:
        params["api_key"] = api_key

    url = f"{NCBI_EFETCH_URL}?{urlencode(params)}"
    logger.debug("EFetch URL: %s", url)

    req = Request(url)
    req.add_header("User-Agent", "fetch_sra_metadata_ncbi/2.0")
    try:
        with urlopen(req, timeout=timeout) as resp:
            return resp.read().decode()
    except (URLError, HTTPError) as exc:
        raise RuntimeError(f"EFetch failed for batch: {exc}") from exc


# ---------------------------------------------------------------------------
# XML parsing
# ---------------------------------------------------------------------------
def _text(elem: ET.Element, tag: str) -> str:
    """Extract text from a child element, returning empty string if absent."""
    child = elem.find(tag)
    if child is not None and child.text:
        return child.text.strip()
    return ""


def _findtext_recursive(elem: ET.Element, tag: str) -> str:
    """Find the first descendant matching *tag* and return its text."""
    found = elem.find(f".//{tag}")
    if found is not None and found.text:
        return found.text.strip()
    return ""


def parse_experiment_package(pkg: ET.Element) -> Dict[str, str]:
    """Parse a single EXPERIMENT_PACKAGE element into a flat dict.

    Parameters
    ----------
    pkg : xml.etree.ElementTree.Element
        An ``<EXPERIMENT_PACKAGE>`` element.

    Returns
    -------
    dict[str, str]
        Flat dictionary with standardised field names.
    """
    row: Dict[str, str] = {}

    # --- EXPERIMENT ---
    exp = pkg.find("EXPERIMENT")
    if exp is not None:
        row["experiment_accession"] = exp.get("accession", "")
        row["experiment_alias"] = exp.get("alias", "")
        row["center_name"] = exp.get("center_name", "")

        # Descriptor
        desc = exp.find("DESIGN/DESIGN_DESCRIPTOR")
        if desc is not None:
            row["library_strategy"] = _text(desc, "LIBRARY_STRATEGY")
            row["library_source"] = _text(desc, "LIBRARY_SOURCE")
            row["library_selection"] = _text(desc, "LIBRARY_SELECTION")

        # Library layout
        layout = exp.find("DESIGN/LIBRARY_DESCRIPTOR/LIBRARY_LAYOUT")
        if layout is not None:
            if layout.find("SINGLE") is not None:
                row["library_layout"] = "SINGLE"
            elif layout.find("PAIRED") is not None:
                row["library_layout"] = "PAIRED"
            else:
                row["library_layout"] = ""

        # Instrument
        plat = exp.find("PLATFORM")
        if plat is not None:
            for child in plat:
                row["instrument_model"] = _text(child, "INSTRUMENT_MODEL")
                if row["instrument_model"]:
                    break

    # --- RUN_SET ---
    run_set = pkg.find("RUN_SET")
    if run_set is not None:
        run = run_set.find("RUN")
        if run is not None:
            row["accession"] = run.get("accession", "")
            row["run_alias"] = run.get("alias", "")
            row["spots"] = run.get("total_spots", "")
            row["bases"] = run.get("total_bases", "")
            row["avg_length"] = run.get("avg_length", "")

            # Experiment ref
            exp_ref = run.find("EXPERIMENT_REF")
            if exp_ref is not None:
                if not row.get("experiment_accession"):
                    row["experiment_accession"] = exp_ref.get("accession", "")

    # --- SAMPLE ---
    sample = pkg.find("SAMPLE")
    if sample is not None:
        row["sample_accession"] = sample.get("accession", "")
        row["sample_alias"] = sample.get("alias", "")
        row["tax_id"] = _text(sample, "SAMPLE_NAME/TAXON_ID")
        row["scientific_name"] = _text(sample, "SAMPLE_NAME/SCIENTIFIC_NAME")

        # Sample attributes -> semicolon-separated key=value pairs
        attrs = []
        sa_block = sample.find("SAMPLE_ATTRIBUTES")
        if sa_block is not None:
            for sa in sa_block.findall("SAMPLE_ATTRIBUTE"):
                k = _text(sa, "TAG")
                v = _text(sa, "VALUE")
                if k:
                    attrs.append(f"{k}={v}")
        row["sample_attributes"] = ";".join(attrs)

    # --- STUDY ---
    study = pkg.find("STUDY")
    if study is not None:
        ext_ids = study.find("IDENTIFIERS")
        if ext_ids is not None:
            for eid in ext_ids.findall("EXTERNAL_ID"):
                ns = eid.get("namespace", "")
                if ns.upper() == "BIOPROJECT":
                    row["bioproject"] = (eid.text or "").strip()
        if not row.get("bioproject"):
            row["bioproject"] = ""

        # Also try viaDescriptor
        for rid in study.iter("PRIMARY_ID"):
            if rid.text and rid.text.startswith("PRJ"):
                row["study_accession"] = rid.text.strip()
                break
        for rid in study.iter("EXTERNAL_ID"):
            if rid.text and rid.text.startswith("SRP"):
                row["study_accession"] = rid.text.strip()
                break
        if not row.get("study_accession"):
            row["study_accession"] = study.get("accession", "")

    # --- Pool / BioSample ---
    if not row.get("biosample"):
        pool = pkg.find("Pool")
        if pool is not None:
            for member in pool.findall("Member"):
                bs = member.get("biosample_accession", "")
                if bs:
                    row["biosample"] = bs
                    break

    return row


def parse_sra_xml(xml_text: str) -> List[Dict[str, str]]:
    """Parse a full EFetch XML response into a list of record dicts.

    Parameters
    ----------
    xml_text : str
        Raw XML from NCBI EFetch.

    Returns
    -------
    list[dict[str, str]]
    """
    root = ET.fromstring(xml_text)
    results: List[Dict[str, str]] = []
    for pkg in root.iter("EXPERIMENT_PACKAGE"):
        try:
            record = parse_experiment_package(pkg)
            results.append(record)
        except Exception as exc:
            logger.warning("Failed to parse an EXPERIMENT_PACKAGE: %s", exc)
    return results


# ---------------------------------------------------------------------------
# Batch processing
# ---------------------------------------------------------------------------
def fetch_metadata_batched(
    accessions: List[str],
    email: str,
    api_key: Optional[str],
    batch_size: int,
    timeout: int,
    rate_limit: float,
) -> List[Dict[str, str]]:
    """Fetch metadata for accessions in batches.

    Parameters
    ----------
    accessions : list[str]
        SRR/ERR/DRR accession identifiers.
    email : str
        Email for NCBI.
    api_key : str or None
        Optional NCBI API key.
    batch_size : int
        Number of accessions per EFetch request.
    timeout : int
        HTTP timeout per request.
    rate_limit : float
        Minimum seconds between requests.

    Returns
    -------
    list[dict[str, str]]
    """
    all_records: List[Dict[str, str]] = []
    total = len(accessions)
    failed_batches: List[List[str]] = []

    for i in range(0, total, batch_size):
        batch = accessions[i : i + batch_size]
        batch_num = i // batch_size + 1
        total_batches = (total + batch_size - 1) // batch_size
        logger.info(
            "Fetching batch %d/%d (%d accessions)",
            batch_num,
            total_batches,
            len(batch),
        )

        try:
            xml_text = efetch_sra_xml(batch, email, api_key, timeout)
            records = parse_sra_xml(xml_text)
            logger.info("  Parsed %d record(s)", len(records))
            all_records.extend(records)
        except RuntimeError as exc:
            logger.error("  Batch %d failed: %s", batch_num, exc)
            failed_batches.append(batch)

        # Rate limiting
        if i + batch_size < total:
            time.sleep(rate_limit)

    if failed_batches:
        logger.warning("%d batch(es) failed.", len(failed_batches))
        # Retry failed batches once
        logger.info("Retrying failed batches ...")
        for batch in failed_batches:
            try:
                xml_text = efetch_sra_xml(batch, email, api_key, timeout)
                records = parse_sra_xml(xml_text)
                logger.info("  Retry parsed %d record(s)", len(records))
                all_records.extend(records)
            except RuntimeError as exc:
                logger.error("  Retry failed: %s", exc)

    return all_records


# ---------------------------------------------------------------------------
# Output writing
# ---------------------------------------------------------------------------
def write_output(
    records: List[Dict[str, str]],
    output_path: str,
    fields: List[str],
) -> None:
    """Write metadata records to a TSV file.

    Parameters
    ----------
    records : list[dict]
        Parsed metadata records.
    output_path : str
        Path to the output file.
    fields : list[str]
        Column names to write.
    """
    outdir = os.path.dirname(os.path.abspath(output_path))
    if outdir:
        os.makedirs(outdir, exist_ok=True)

    # Determine delimiter from file extension
    delim = "\t"
    if output_path.endswith(".csv"):
        delim = ","

    with open(output_path, "w", newline="", encoding="utf-8") as fh:
        writer = csv.DictWriter(
            fh, fieldnames=fields, delimiter=delim, extrasaction="ignore"
        )
        writer.writeheader()
        for rec in records:
            # Ensure all fields are present
            for f in fields:
                rec.setdefault(f, "")
            writer.writerow(rec)

    logger.info("Wrote %d records to %s", len(records), output_path)


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------
def build_parser() -> argparse.ArgumentParser:
    """Construct the argument parser."""
    parser = argparse.ArgumentParser(
        prog="fetch_sra_metadata_ncbi.py",
        description=(
            "Fetch SRA metadata from NCBI by querying SRR accession IDs via "
            "Entrez E-Utilities. Parses EXPERIMENT_PACKAGE XML and outputs "
            "a tab-delimited (or CSV) table of sample attributes."
        ),
        epilog=(
            "Examples:\n"
            "  # Basic usage with an accession list\n"
            "  %(prog)s -i srr_list.txt -o metadata.tsv -e your@email.com\n"
            "\n"
            "  # Batch size 100, 3 requests/sec with API key\n"
            "  %(prog)s -i srr_list.txt -o metadata.tsv -e your@email.com \\\n"
            "       -b 100 --api-key YOUR_NCBI_API_KEY\n"
            "\n"
            "  # Custom field list\n"
            "  %(prog)s -i srr_list.txt -o metadata.tsv -e your@email.com \\\n"
            "       --fields accession,sample_accession,scientific_name,sample_attributes\n"
            "\n"
            "  # Single accession from command line\n"
            "  %(prog)s -i SRR1234567 -o metadata.tsv -e your@email.com\n"
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument(
        "-i", "--input", required=True,
        help=(
            "Input file with one SRR/ERR/DRR accession per line, "
            "or a single accession / comma-separated list."
        ),
    )
    parser.add_argument(
        "-o", "--output", default="sra_metadata.tsv",
        help="Output file path (TSV or CSV based on extension, default: sra_metadata.tsv)",
    )
    parser.add_argument(
        "-e", "--email", required=True,
        help="Email address for NCBI Entrez (required by NCBI policy)",
    )
    parser.add_argument(
        "-b", "--batch-size", type=int, default=DEFAULT_BATCH_SIZE,
        help=f"Number of accessions per EFetch request (default: {DEFAULT_BATCH_SIZE})",
    )
    parser.add_argument(
        "-t", "--threads", type=int, default=DEFAULT_THREADS,
        help=f"Number of threads (default: {DEFAULT_THREADS})",
    )
    parser.add_argument(
        "--api-key", default=None,
        help="NCBI API key for higher rate limits",
    )
    parser.add_argument(
        "--rate-limit", type=float, default=DEFAULT_RATE_LIMIT,
        help=f"Seconds to wait between API calls (default: {DEFAULT_RATE_LIMIT})",
    )
    parser.add_argument(
        "--timeout", type=int, default=120,
        help="HTTP request timeout in seconds (default: 120)",
    )
    parser.add_argument(
        "--fields", default=None,
        help=(
            "Comma-separated list of output fields. Default: all standard fields "
            "plus sample_attributes."
        ),
    )
    parser.add_argument(
        "--log-level", default="INFO",
        choices=["DEBUG", "INFO", "WARNING", "ERROR", "CRITICAL"],
        help="Set logging verbosity (default: INFO)",
    )
    parser.add_argument(
        "--version", action="version",
        version=f"%(prog)s {__version__}",
    )
    return parser


def setup_logging(level: str) -> None:
    """Configure the root logger.

    Parameters
    ----------
    level : str
        Logging level string.
    """
    logging.basicConfig(
        level=getattr(logging, level),
        format="%(asctime)s [%(levelname)-8s] %(name)s: %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S",
    )


def read_accessions(input_arg: str) -> List[str]:
    """Read accession IDs from a file or treat input as a bare accession.

    Parameters
    ----------
    input_arg : str
        File path or accession identifier(s).

    Returns
    -------
    list[str]
    """
    path = Path(input_arg)
    if path.is_file():
        accs = [
            line.strip()
            for line in path.read_text().splitlines()
            if line.strip() and not line.strip().startswith("#")
        ]
        logger.info("Read %d accession(s) from %s", len(accs), input_arg)
        return accs
    # Single accession or comma-separated
    accs = [a.strip() for a in input_arg.split(",") if a.strip()]
    return accs


def parse_field_list(fields_str: Optional[str]) -> List[str]:
    """Parse a comma-separated field list, falling back to defaults.

    Parameters
    ----------
    fields_str : str or None
        Comma-separated field names, or None for defaults.

    Returns
    -------
    list[str]
    """
    if fields_str is None:
        return list(DEFAULT_OUTPUT_FIELDS)
    return [f.strip() for f in fields_str.split(",") if f.strip()]


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
def main() -> None:
    """Entry point for the NCBI SRA metadata fetcher."""
    parser = build_parser()
    args = parser.parse_args()
    setup_logging(args.log_level)

    logger.info("=== fetch_sra_metadata_ncbi.py %s ===", __version__)

    # Read accessions
    accessions = read_accessions(args.input)
    if not accessions:
        logger.error("No accessions to process.")
        sys.exit(1)

    logger.info("Total accessions: %d", len(accessions))
    logger.info("Batch size: %d", args.batch_size)
    logger.info("Rate limit: %.2f s", args.rate_limit)

    # Resolve output fields
    fields = parse_field_list(args.fields)
    logger.info("Output fields (%d): %s", len(fields), ", ".join(fields))

    # Fetch metadata
    try:
        records = fetch_metadata_batched(
            accessions=accessions,
            email=args.email,
            api_key=args.api_key,
            batch_size=args.batch_size,
            timeout=args.timeout,
            rate_limit=args.rate_limit,
        )
    except Exception as exc:
        logger.error("Metadata fetching failed: %s", exc)
        sys.exit(1)

    if not records:
        logger.warning("No metadata records retrieved.")
        sys.exit(0)

    # Write output
    write_output(records, args.output, fields)

    # Summary
    found_accessions: Set[str] = {r.get("accession", "") for r in records}
    missing = set(accessions) - found_accessions
    logger.info("=" * 50)
    logger.info("Records retrieved : %d", len(records))
    logger.info("Output file       : %s", args.output)
    if missing:
        logger.warning(
            "Accessions not found in results (%d): %s",
            len(missing),
            ", ".join(sorted(missing)[:10]),
        )
        miss_file = os.path.splitext(args.output)[0] + ".missing.txt"
        with open(miss_file, "w") as fh:
            for m in sorted(missing):
                fh.write(m + "\n")
        logger.warning("Missing accessions written to %s", miss_file)
    logger.info("Done.")


if __name__ == "__main__":
    main()
