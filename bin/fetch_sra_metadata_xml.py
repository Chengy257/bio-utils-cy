#!/usr/bin/env python3
"""
File Name: fetch_sra_metadata_xml.py
Author: ChengYu
Description: Fetch SRA metadata by combining RunInfo CSV with esummary XML.
             Extracts SAMPLE_ATTRIBUTE fields from XML, merges with RunInfo
             columns, and produces a unified TSV output.
Created Time: 2026
"""

from __future__ import annotations

import argparse
import csv
import io
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
DEFAULT_BATCH_SIZE = 100
DEFAULT_TIMEOUT = 120
DEFAULT_MAX_RETRIES = 3
DEFAULT_RETRY_DELAY = 3.0

logger = logging.getLogger(__name__)


# ---------------------------------------------------------------------------
# HTTP helpers
# ---------------------------------------------------------------------------

def _fetch(url: str, timeout: int = DEFAULT_TIMEOUT) -> str:
    """Perform a GET request and return decoded body text.

    Args:
        url: The URL to fetch.
        timeout: Timeout in seconds.

    Returns:
        Decoded response text.

    Raises:
        RuntimeError: On any HTTP/URL error.
    """
    logger.debug("GET %s", url)
    try:
        with urlopen(Request(url), timeout=timeout) as resp:
            charset = resp.headers.get_content_charset() or "utf-8"
            return resp.read().decode(charset)
    except HTTPError as exc:
        raise RuntimeError(f"HTTP {exc.code}: {exc.reason} ({url})") from exc
    except URLError as exc:
        raise RuntimeError(f"URL error: {exc} ({url})") from exc


def _fetch_retry(
    url: str,
    timeout: int = DEFAULT_TIMEOUT,
    retries: int = DEFAULT_MAX_RETRIES,
    delay: float = DEFAULT_RETRY_DELAY,
) -> str:
    """Fetch with exponential-backoff retries.

    Args:
        url: URL to fetch.
        timeout: Per-request timeout in seconds.
        retries: Maximum attempts.
        delay: Base delay between retries in seconds.

    Returns:
        Decoded response text.

    Raises:
        RuntimeError: After all retries exhausted.
    """
    last_exc: Optional[Exception] = None
    for attempt in range(1, retries + 1):
        try:
            return _fetch(url, timeout=timeout)
        except RuntimeError as exc:
            last_exc = exc
            logger.warning("Attempt %d/%d failed: %s", attempt, retries, exc)
            if attempt < retries:
                wait = delay * (2 ** (attempt - 1))
                logger.info("Sleeping %.1f s before retry ...", wait)
                time.sleep(wait)
    raise RuntimeError(f"All {retries} retries exhausted: {last_exc}")


# ---------------------------------------------------------------------------
# ID input
# ---------------------------------------------------------------------------

def _read_ids(path: str) -> List[str]:
    """Read accessions from a text file, one per line.

    Args:
        path: File path.

    Returns:
        List of accession strings.
    """
    ids: List[str] = []
    with open(path, "r", encoding="utf-8") as fh:
        for raw in fh:
            tok = raw.strip()
            if tok and not tok.startswith("#"):
                ids.append(tok)
    logger.info("Loaded %d IDs from %s", len(ids), path)
    return ids


# ---------------------------------------------------------------------------
# Step 1: esearch to convert accessions -> UIDs
# ---------------------------------------------------------------------------

def esearch_uids(
    accessions: List[str],
    batch_size: int = DEFAULT_BATCH_SIZE,
    timeout: int = DEFAULT_TIMEOUT,
    api_key: Optional[str] = None,
) -> List[str]:
    """Convert SRA accessions to NCBI UIDs via esearch.

    Args:
        accessions: List of SRA run accessions.
        batch_size: Number of accessions per request.
        timeout: HTTP timeout.
        api_key: Optional NCBI API key.

    Returns:
        List of UID strings.
    """
    all_uids: List[str] = []
    for start in range(0, len(accessions), batch_size):
        batch = accessions[start : start + batch_size]
        query = " OR ".join(batch)
        params: Dict[str, str] = {
            "db": "sra",
            "term": query,
            "retmax": str(len(batch)),
            "usehistory": "n",
            "retmode": "json",
        }
        if api_key:
            params["api_key"] = api_key
        url = f"{NCBI_EUTILS_BASE}/esearch.fcgi?{urlencode(params)}"
        logger.info("esearch batch %d-%d", start + 1, start + len(batch))
        text = _fetch_retry(url, timeout=timeout)
        import json
        data = json.loads(text)
        id_list = data.get("esearchresult", {}).get("idlist", [])
        all_uids.extend(id_list)
        time.sleep(0.4 if api_key else 0.6)
    logger.info("Resolved %d UIDs from %d accessions", len(all_uids), len(accessions))
    return all_uids


# ---------------------------------------------------------------------------
# Step 2: RunInfo CSV via efetch
# ---------------------------------------------------------------------------

def fetch_runinfo_csv(
    uids: List[str],
    batch_size: int = DEFAULT_BATCH_SIZE,
    timeout: int = DEFAULT_TIMEOUT,
    api_key: Optional[str] = None,
) -> List[Dict[str, str]]:
    """Fetch RunInfo CSV for given UIDs.

    Args:
        uids: NCBI SRA UIDs.
        batch_size: UIDs per request.
        timeout: HTTP timeout.
        api_key: Optional API key.

    Returns:
        List of dicts from CSV rows.
    """
    all_rows: List[Dict[str, str]] = []
    for start in range(0, len(uids), batch_size):
        batch = uids[start : start + batch_size]
        params: Dict[str, str] = {
            "db": "sra",
            "id": ",".join(batch),
            "rettype": "runinfo",
            "retmode": "text",
        }
        if api_key:
            params["api_key"] = api_key
        url = f"{NCBI_EUTILS_BASE}/efetch.fcgi?{urlencode(params)}"
        logger.info("RunInfo batch %d-%d", start + 1, start + len(batch))
        csv_text = _fetch_retry(url, timeout=timeout)
        reader = csv.DictReader(io.StringIO(csv_text))
        for row in reader:
            all_rows.append(dict(row))
        time.sleep(0.4 if api_key else 0.6)
    logger.info("RunInfo: %d rows", len(all_rows))
    return all_rows


# ---------------------------------------------------------------------------
# Step 3: esummary / efetch XML for SAMPLE_ATTRIBUTE
# ---------------------------------------------------------------------------

def fetch_sample_attributes(
    uids: List[str],
    batch_size: int = DEFAULT_BATCH_SIZE,
    timeout: int = DEFAULT_TIMEOUT,
    api_key: Optional[str] = None,
) -> Dict[str, Dict[str, str]]:
    """Fetch SRA XML and extract SAMPLE_ATTRIBUTE fields per run.

    Args:
        uids: NCBI SRA UIDs.
        batch_size: UIDs per request.
        timeout: HTTP timeout.
        api_key: Optional API key.

    Returns:
        Dict mapping run accession -> sample attribute dict.
    """
    sample_attrs: Dict[str, Dict[str, str]] = {}
    for start in range(0, len(uids), batch_size):
        batch = uids[start : start + batch_size]
        params: Dict[str, str] = {
            "db": "sra",
            "id": ",".join(batch),
            "rettype": "full",
            "retmode": "xml",
        }
        if api_key:
            params["api_key"] = api_key
        url = f"{NCBI_EUTILS_BASE}/efetch.fcgi?{urlencode(params)}"
        logger.info("XML batch %d-%d", start + 1, start + len(batch))
        xml_text = _fetch_retry(url, timeout=timeout)
        root = ET.fromstring(xml_text)
        for pkg in root.iter("EXPERIMENT_PACKAGE"):
            sample_elem = pkg.find(".//SAMPLE")
            run_elems = pkg.findall(".//RUN")
            if sample_elem is None or not run_elems:
                continue
            attrs: Dict[str, str] = {}
            for sa in sample_elem.iter("SAMPLE_ATTRIBUTE"):
                tag_elem = sa.find("TAG")
                val_elem = sa.find("VALUE")
                if tag_elem is not None and tag_elem.text:
                    key = tag_elem.text.strip()
                    val = val_elem.text.strip() if val_elem is not None and val_elem.text else ""
                    attrs[key] = val
            sample_acc = sample_elem.attrib.get("accession", "")
            if sample_acc:
                attrs["sample_accession"] = sample_acc
            for run_elem in run_elems:
                run_acc = run_elem.attrib.get("accession", "")
                if run_acc:
                    sample_attrs[run_acc] = dict(attrs)
        time.sleep(0.4 if api_key else 0.6)
    logger.info("Extracted sample attributes for %d runs", len(sample_attrs))
    return sample_attrs


# ---------------------------------------------------------------------------
# Merge & output
# ---------------------------------------------------------------------------

def merge(
    runinfo_rows: List[Dict[str, str]],
    sample_attrs: Dict[str, Dict[str, str]],
) -> List[Dict[str, str]]:
    """Merge RunInfo rows with per-run SAMPLE_ATTRIBUTE dicts.

    Args:
        runinfo_rows: List of RunInfo CSV dicts.
        sample_attrs: run accession -> sample attribute dict.

    Returns:
        Merged list of dicts.
    """
    merged: List[Dict[str, str]] = []
    for row in runinfo_rows:
        run_acc = row.get("Run", row.get("run_accession", ""))
        combined = dict(row)
        sa = sample_attrs.get(run_acc, {})
        for k, v in sa.items():
            col = f"sample_{k}" if not k.startswith("sample_") else k
            combined[col] = v
        merged.append(combined)
    return merged


def _all_keys(rows: List[Dict[str, str]]) -> List[str]:
    """Collect ordered unique keys from all rows.

    Args:
        rows: List of dicts.

    Returns:
        Ordered list of keys.
    """
    seen: Dict[str, None] = {}
    for row in rows:
        for k in row:
            if k not in seen:
                seen[k] = None
    return list(seen.keys())


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def main() -> None:
    """Main entry point."""
    parser = argparse.ArgumentParser(
        description="Fetch SRA metadata combining RunInfo CSV and esummary XML.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""\
examples:
  %(prog)s -i ids.txt -o result.tsv --email you@example.com
  %(prog)s -i ids.txt -o out.tsv --email you@example.com --api-key MYKEY --batch 50
  %(prog)s -i ids.txt -o out.tsv --email you@example.com --timeout 180 --log-level DEBUG
""",
    )
    parser.add_argument("-i", "--input", required=True, help="File with one SRA accession per line.")
    parser.add_argument("-o", "--output", default="-", help="Output TSV path (default: stdout).")
    parser.add_argument("--email", required=True, help="Email for NCBI E-utilities.")
    parser.add_argument(
        "--api-key",
        default=os.environ.get("NCBI_API_KEY", ""),
        help="NCBI API key (or set NCBI_API_KEY env var).",
    )
    parser.add_argument(
        "--batch",
        type=int,
        default=DEFAULT_BATCH_SIZE,
        help=f"Batch size per API request (default: {DEFAULT_BATCH_SIZE}).",
    )
    parser.add_argument(
        "--timeout",
        type=int,
        default=DEFAULT_TIMEOUT,
        help=f"HTTP timeout in seconds (default: {DEFAULT_TIMEOUT}).",
    )
    parser.add_argument(
        "--retries",
        type=int,
        default=DEFAULT_MAX_RETRIES,
        help=f"Retry attempts (default: {DEFAULT_MAX_RETRIES}).",
    )
    parser.add_argument(
        "--log-level",
        default="INFO",
        choices=["DEBUG", "INFO", "WARNING", "ERROR", "CRITICAL"],
        help="Logging level (default: INFO).",
    )
    parser.add_argument("--version", action="version", version=f"%(prog)s {__version__}")
    args = parser.parse_args()

    logging.basicConfig(
        level=getattr(logging, args.log_level),
        format="%(asctime)s [%(levelname)s] %(name)s: %(message)s",
    )

    accessions = _read_ids(args.input)
    if not accessions:
        logger.error("No accessions loaded from %s", args.input)
        sys.exit(1)

    api_key = args.api_key or None

    # Resolve UIDs
    logger.info("Resolving NCBI UIDs for %d accessions ...", len(accessions))
    uids = esearch_uids(accessions, batch_size=args.batch, timeout=args.timeout, api_key=api_key)
    if not uids:
        logger.error("No UIDs resolved. Check your accessions.")
        sys.exit(1)

    # RunInfo CSV
    logger.info("Fetching RunInfo CSV ...")
    runinfo_rows = fetch_runinfo_csv(uids, batch_size=args.batch, timeout=args.timeout, api_key=api_key)

    # Sample attributes from XML
    logger.info("Fetching sample attributes from XML ...")
    sample_attrs = fetch_sample_attributes(uids, batch_size=args.batch, timeout=args.timeout, api_key=api_key)

    # Merge
    merged = merge(runinfo_rows, sample_attrs)
    fields = _all_keys(merged)

    if args.output == "-":
        writer = csv.DictWriter(sys.stdout, fieldnames=fields, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(merged)
    else:
        out = Path(args.output)
        out.parent.mkdir(parents=True, exist_ok=True)
        with open(out, "w", newline="", encoding="utf-8") as fh:
            writer = csv.DictWriter(fh, fieldnames=fields, delimiter="\t", extrasaction="ignore")
            writer.writeheader()
            writer.writerows(merged)
        logger.info("Wrote %d rows to %s", len(merged), args.output)


if __name__ == "__main__":
    main()
