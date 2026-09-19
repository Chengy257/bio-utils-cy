#!/usr/bin/env python3
"""
File Name: fetch_sra_metadata_comprehensive.py
Author: ChengYu
Description: Three-step comprehensive SRA metadata retrieval:
             RunInfo -> BioSample -> Experiment from NCBI.
             Fetches metadata in batches, merges results into a single TSV.
Created Time: 2026
"""

from __future__ import annotations

import argparse
import csv
import logging
import os
import sys
import time
import xml.etree.ElementTree as ET
from pathlib import Path
from typing import Any, Dict, List, Optional
from urllib.error import HTTPError, URLError
from urllib.parse import urlencode
from urllib.request import Request, urlopen

__version__ = "1.0.0"

NCBI_EUTILS_BASE = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
DEFAULT_BATCH_SIZE = 50
DEFAULT_TIMEOUT = 120
DEFAULT_RETRY_DELAY = 3.0
DEFAULT_MAX_RETRIES = 3

logger = logging.getLogger(__name__)


# ---------------------------------------------------------------------------
# Utility helpers
# ---------------------------------------------------------------------------

def _fetch_url(url: str, timeout: int = DEFAULT_TIMEOUT) -> str:
    """Fetch content from a URL.

    Args:
        url: The URL to request.
        timeout: Request timeout in seconds.

    Returns:
        Response body as a decoded string.

    Raises:
        RuntimeError: On any fetch failure.
    """
    logger.debug("Fetching: %s", url)
    req = Request(url)
    try:
        with urlopen(req, timeout=timeout) as resp:
            charset = resp.headers.get_content_charset() or "utf-8"
            return resp.read().decode(charset)
    except HTTPError as exc:
        raise RuntimeError(f"HTTP {exc.code} fetching {url}: {exc.reason}") from exc
    except URLError as exc:
        raise RuntimeError(f"URL error fetching {url}: {exc}") from exc


def _fetch_with_retry(
    url: str,
    timeout: int = DEFAULT_TIMEOUT,
    max_retries: int = DEFAULT_MAX_RETRIES,
    retry_delay: float = DEFAULT_RETRY_DELAY,
) -> str:
    """Fetch a URL with exponential-backoff retries.

    Args:
        url: URL to fetch.
        timeout: Per-request timeout in seconds.
        max_retries: Maximum number of attempts.
        retry_delay: Base delay in seconds between retries.

    Returns:
        Decoded response text.

    Raises:
        RuntimeError: If all retries are exhausted.
    """
    last_exc: Optional[Exception] = None
    for attempt in range(1, max_retries + 1):
        try:
            return _fetch_url(url, timeout=timeout)
        except RuntimeError as exc:
            last_exc = exc
            logger.warning("Attempt %d/%d failed: %s", attempt, max_retries, exc)
            if attempt < max_retries:
                sleep = retry_delay * (2 ** (attempt - 1))
                logger.info("Retrying in %.1f s ...", sleep)
                time.sleep(sleep)
    raise RuntimeError(f"All {max_retries} retries exhausted. Last: {last_exc}")


def _read_ids(path: str) -> List[str]:
    """Read IDs from a file, one per line, skipping blanks and comments.

    Args:
        path: Path to the text file.

    Returns:
        List of ID strings.
    """
    ids: List[str] = []
    with open(path, "r", encoding="utf-8") as fh:
        for line in fh:
            token = line.strip()
            if token and not token.startswith("#"):
                ids.append(token)
    return ids


# ---------------------------------------------------------------------------
# Step 1: RunInfo
# ---------------------------------------------------------------------------

def fetch_runinfo(
    srr_ids: List[str],
    batch_size: int = DEFAULT_BATCH_SIZE,
    timeout: int = DEFAULT_TIMEOUT,
    api_key: Optional[str] = None,
) -> List[Dict[str, str]]:
    """Fetch RunInfo CSV metadata from NCBI for given SRR IDs.

    Args:
        srr_ids: List of SRR/ERR/DRR accession strings.
        batch_size: Number of IDs per request.
        timeout: HTTP timeout per request.
        api_key: Optional NCBI API key for higher rate limits.

    Returns:
        List of dicts representing each run's RunInfo fields.
    """
    all_rows: List[Dict[str, str]] = []
    for start in range(0, len(srr_ids), batch_size):
        batch = srr_ids[start : start + batch_size]
        id_str = ",".join(batch)
        params = {
            "db": "sra",
            "id": id_str,
            "rettype": "runinfo",
            "retmode": "text",
        }
        if api_key:
            params["api_key"] = api_key
        url = f"{NCBI_EUTILS_BASE}/efetch.fcgi?{urlencode(params)}"
        logger.info(
            "Fetching RunInfo batch %d-%d of %d",
            start + 1,
            min(start + batch_size, len(srr_ids)),
            len(srr_ids),
        )
        csv_text = _fetch_with_retry(url, timeout=timeout)
        reader = csv.DictReader(csv_text.splitlines())
        for row in reader:
            all_rows.append(dict(row))
        if start + batch_size < len(srr_ids):
            time.sleep(0.4 if api_key else 0.6)
    logger.info("RunInfo: retrieved %d rows", len(all_rows))
    return all_rows


# ---------------------------------------------------------------------------
# Step 2: BioSample
# ---------------------------------------------------------------------------

def _extract_biosample_accessions(runinfo_rows: List[Dict[str, str]]) -> List[str]:
    """Extract unique BioSample accessions from RunInfo rows.

    Args:
        runinfo_rows: Rows from the RunInfo CSV.

    Returns:
        Deduplicated list of BioSample accessions (SAMN/SAME/ERS/...).
    """
    sample_set = set()
    for row in runinfo_rows:
        for key in ("BioSample", "SampleAccession", "bio_sample", "biosample"):
            val = row.get(key, "").strip()
            if val:
                sample_set.add(val)
                break
    return sorted(sample_set)


def fetch_biosample_xml(
    biosample_ids: List[str],
    batch_size: int = DEFAULT_BATCH_SIZE,
    timeout: int = DEFAULT_TIMEOUT,
    api_key: Optional[str] = None,
) -> Dict[str, Dict[str, str]]:
    """Fetch BioSample XML and extract sample attributes.

    Args:
        biosample_ids: List of BioSample accessions.
        batch_size: Number of IDs per request.
        timeout: HTTP timeout per request.
        api_key: Optional NCBI API key.

    Returns:
        Dict mapping BioSample accession -> attribute dict.
    """
    biosample_data: Dict[str, Dict[str, str]] = {}
    for start in range(0, len(biosample_ids), batch_size):
        batch = biosample_ids[start : start + batch_size]
        id_str = ",".join(batch)
        params = {
            "db": "biosample",
            "id": id_str,
            "rettype": "full",
            "retmode": "xml",
        }
        if api_key:
            params["api_key"] = api_key
        url = f"{NCBI_EUTILS_BASE}/efetch.fcgi?{urlencode(params)}"
        logger.info(
            "Fetching BioSample XML batch %d-%d of %d",
            start + 1,
            min(start + batch_size, len(biosample_ids)),
            len(biosample_ids),
        )
        xml_text = _fetch_with_retry(url, timeout=timeout)
        root = ET.fromstring(xml_text)
        for bs_elem in root.iter("BioSample"):
            acc = bs_elem.attrib.get("accession", bs_elem.attrib.get("id", ""))
            attrs: Dict[str, str] = {}
            for attr_elem in bs_elem.iter("Attribute"):
                attr_name = attr_elem.attrib.get(
                    "harmonized_name",
                    attr_elem.attrib.get("attribute_name", attr_elem.tag),
                )
                attrs[attr_name] = (attr_elem.text or "").strip()
            description = bs_elem.findtext("Description/Title", default="").strip()
            if description:
                attrs["description_title"] = description
            organism = bs_elem.findtext(
                "Description/Organism/OrganismName", default=""
            ).strip()
            if organism:
                attrs["organism"] = organism
            if acc:
                biosample_data[acc] = attrs
        if start + batch_size < len(biosample_ids):
            time.sleep(0.4 if api_key else 0.6)
    logger.info("BioSample: retrieved %d samples", len(biosample_data))
    return biosample_data


# ---------------------------------------------------------------------------
# Step 3: Experiment (SRA experiment XML)
# ---------------------------------------------------------------------------

def fetch_experiment_xml(
    srr_ids: List[str],
    batch_size: int = DEFAULT_BATCH_SIZE,
    timeout: int = DEFAULT_TIMEOUT,
    api_key: Optional[str] = None,
) -> Dict[str, Dict[str, str]]:
    """Fetch SRA experiment XML and extract experiment-level attributes.

    Args:
        srr_ids: List of SRR accessions.
        batch_size: Number of IDs per request.
        timeout: HTTP timeout per request.
        api_key: Optional NCBI API key.

    Returns:
        Dict mapping run accession -> experiment attribute dict.
    """
    experiment_data: Dict[str, Dict[str, str]] = {}
    for start in range(0, len(srr_ids), batch_size):
        batch = srr_ids[start : start + batch_size]
        id_str = ",".join(batch)
        params = {
            "db": "sra",
            "id": id_str,
            "rettype": "full",
            "retmode": "xml",
        }
        if api_key:
            params["api_key"] = api_key
        url = f"{NCBI_EUTILS_BASE}/efetch.fcgi?{urlencode(params)}"
        logger.info(
            "Fetching Experiment XML batch %d-%d of %d",
            start + 1,
            min(start + batch_size, len(srr_ids)),
            len(srr_ids),
        )
        xml_text = _fetch_with_retry(url, timeout=timeout)
        root = ET.fromstring(xml_text)
        for pkg in root.iter("EXPERIMENT_PACKAGE"):
            run_elem = pkg.find(".//RUN")
            if run_elem is None:
                continue
            run_acc = run_elem.attrib.get("accession", "")
            exp_info: Dict[str, str] = {}
            exp_elem = pkg.find(".//EXPERIMENT")
            if exp_elem is not None:
                exp_info["experiment_accession"] = exp_elem.attrib.get("accession", "")
                title_elem = exp_elem.find("TITLE")
                if title_elem is not None and title_elem.text:
                    exp_info["experiment_title"] = title_elem.text.strip()
            lib_elem = pkg.find(".//LIBRARY_DESCRIPTOR")
            if lib_elem is not None:
                for child in lib_elem:
                    tag = child.tag.replace("LIBRARY_", "")
                    if child.text:
                        exp_info[f"library_{tag.lower()}"] = child.text.strip()
            plat_elem = pkg.find(".//INSTRUMENT_MODEL")
            if plat_elem is not None and plat_elem.text:
                exp_info["instrument_model"] = plat_elem.text.strip()
            if run_acc:
                experiment_data[run_acc] = exp_info
        if start + batch_size < len(srr_ids):
            time.sleep(0.4 if api_key else 0.6)
    logger.info("Experiment: retrieved %d runs", len(experiment_data))
    return experiment_data


# ---------------------------------------------------------------------------
# Merge and write
# ---------------------------------------------------------------------------

def merge_metadata(
    runinfo_rows: List[Dict[str, str]],
    biosample_data: Dict[str, Dict[str, str]],
    experiment_data: Dict[str, Dict[str, str]],
) -> List[Dict[str, str]]:
    """Merge RunInfo, BioSample, and Experiment data into unified rows.

    Args:
        runinfo_rows: Rows from RunInfo CSV.
        biosample_data: BioSample accession -> attributes dict.
        experiment_data: Run accession -> experiment attributes dict.

    Returns:
        Merged list of flat dicts, one per run.
    """
    merged: List[Dict[str, str]] = []
    for row in runinfo_rows:
        combined = dict(row)
        run_acc = row.get("Run", row.get("run_accession", ""))
        bio_acc = row.get("BioSample", row.get("SampleAccession", ""))
        if bio_acc and bio_acc in biosample_data:
            for k, v in biosample_data[bio_acc].items():
                combined[f"biosample_{k}"] = v
        if run_acc and run_acc in experiment_data:
            for k, v in experiment_data[run_acc].items():
                combined[f"experiment_{k}"] = v
        merged.append(combined)
    return merged


def collect_fieldnames(rows: List[Dict[str, str]]) -> List[str]:
    """Collect all unique keys across rows, preserving first-seen order.

    Args:
        rows: List of dicts.

    Returns:
        Ordered list of all keys.
    """
    seen: Dict[str, None] = {}
    for row in rows:
        for key in row:
            if key not in seen:
                seen[key] = None
    return list(seen.keys())


def write_tsv(rows: List[Dict[str, str]], path: str) -> None:
    """Write merged rows to a TSV file.

    Args:
        rows: Merged data rows.
        path: Output file path.
    """
    fields = collect_fieldnames(rows)
    out_path = Path(path)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with open(out_path, "w", newline="", encoding="utf-8") as fh:
        writer = csv.DictWriter(fh, fieldnames=fields, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)
    logger.info("Wrote %d rows x %d columns to %s", len(rows), len(fields), path)


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def main() -> None:
    """Main entry point."""
    parser = argparse.ArgumentParser(
        description=(
            "Fetch comprehensive SRA metadata via NCBI E-utilities "
            "(RunInfo -> BioSample -> Experiment)."
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""\
examples:
  %(prog)s --srr-ids SRR123456 SRR789012 --output-dir ./meta --email you@example.com
  %(prog)s --input-file ids.txt --output-dir ./results --email you@example.com --api-key ABC123
  %(prog)s --input-file ids.txt -o out --email you@example.com --batch-size 20 --log-level DEBUG
""",
    )
    input_group = parser.add_mutually_exclusive_group(required=True)
    input_group.add_argument(
        "--input-file", help="Path to file with one SRR ID per line."
    )
    input_group.add_argument(
        "--srr-ids", nargs="+", help="One or more SRR accessions on the command line."
    )
    parser.add_argument(
        "--output-dir",
        "-o",
        default=".",
        help="Directory for output TSV files (default: current directory).",
    )
    parser.add_argument(
        "--email", required=True, help="Email address for NCBI E-utilities."
    )
    parser.add_argument(
        "--api-key",
        default=os.environ.get("NCBI_API_KEY", ""),
        help="NCBI API key (or set NCBI_API_KEY env var).",
    )
    parser.add_argument(
        "--batch-size",
        type=int,
        default=DEFAULT_BATCH_SIZE,
        help=f"Number of IDs per API request (default: {DEFAULT_BATCH_SIZE}).",
    )
    parser.add_argument(
        "--timeout",
        type=int,
        default=DEFAULT_TIMEOUT,
        help=f"HTTP timeout per request in seconds (default: {DEFAULT_TIMEOUT}).",
    )
    parser.add_argument(
        "--retries",
        type=int,
        default=DEFAULT_MAX_RETRIES,
        help=f"Retry attempts per request (default: {DEFAULT_MAX_RETRIES}).",
    )
    parser.add_argument(
        "--log-level",
        default="INFO",
        choices=["DEBUG", "INFO", "WARNING", "ERROR", "CRITICAL"],
        help="Set the logging level (default: INFO).",
    )
    parser.add_argument(
        "--version", action="version", version=f"%(prog)s {__version__}"
    )
    args = parser.parse_args()

    logging.basicConfig(
        level=getattr(logging, args.log_level),
        format="%(asctime)s [%(levelname)s] %(name)s: %(message)s",
    )

    if args.input_file:
        srr_ids = _read_ids(args.input_file)
    else:
        srr_ids = list(args.srr_ids)

    if not srr_ids:
        logger.error("No SRR IDs provided.")
        sys.exit(1)

    logger.info("Starting comprehensive metadata fetch for %d accession(s).", len(srr_ids))
    out_dir = Path(args.output_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    # Step 1: RunInfo
    logger.info("Step 1/3: Fetching RunInfo ...")
    runinfo_rows = fetch_runinfo(
        srr_ids,
        batch_size=args.batch_size,
        timeout=args.timeout,
        api_key=args.api_key or None,
    )
    if not runinfo_rows:
        logger.error("No RunInfo data retrieved. Exiting.")
        sys.exit(1)

    # Step 2: BioSample
    logger.info("Step 2/3: Fetching BioSample XML ...")
    biosample_ids = _extract_biosample_accessions(runinfo_rows)
    biosample_data: Dict[str, Dict[str, str]] = {}
    if biosample_ids:
        biosample_data = fetch_biosample_xml(
            biosample_ids,
            batch_size=args.batch_size,
            timeout=args.timeout,
            api_key=args.api_key or None,
        )
    else:
        logger.warning("No BioSample accessions found in RunInfo.")

    # Step 3: Experiment
    logger.info("Step 3/3: Fetching Experiment XML ...")
    experiment_data = fetch_experiment_xml(
        srr_ids,
        batch_size=args.batch_size,
        timeout=args.timeout,
        api_key=args.api_key or None,
    )

    # Merge
    merged = merge_metadata(runinfo_rows, biosample_data, experiment_data)
    output_path = str(out_dir / "sra_metadata_comprehensive.tsv")
    write_tsv(merged, output_path)
    logger.info("Done. Output: %s", output_path)


if __name__ == "__main__":
    main()
