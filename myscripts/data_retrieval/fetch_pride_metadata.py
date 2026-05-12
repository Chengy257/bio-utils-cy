#!/usr/bin/env python3
"""
File Name: fetch_pride_metadata.py
Author: ChengYu
Description: Fetch PRIDE proteomics dataset metadata for PXD accessions.
             Queries the PRIDE API and outputs a merged TSV of project metadata.
Created Time: 2026
"""

from __future__ import annotations

import argparse
import csv
import json
import logging
import sys
import time
from pathlib import Path
from typing import Any, Dict, List, Optional
from urllib.error import HTTPError, URLError
from urllib.request import Request, urlopen

__version__ = "1.0.0"

PRIDE_PROJECT_API = "https://www.ebi.ac.uk/pride/ws/archive/v2/projects/{accession}"
PRIDE_FILE_API = "https://www.ebi.ac.uk/pride/ws/archive/v2/files/{accession}"
DEFAULT_TIMEOUT = 60
DEFAULT_MAX_RETRIES = 3
DEFAULT_RETRY_DELAY = 2.0

logger = logging.getLogger(__name__)


# ---------------------------------------------------------------------------
# HTTP helpers
# ---------------------------------------------------------------------------

def _fetch_json(url: str, timeout: int = DEFAULT_TIMEOUT) -> Any:
    """Fetch JSON content from a URL.

    Args:
        url: URL to fetch.
        timeout: Request timeout in seconds.

    Returns:
        Parsed JSON object (dict or list).

    Raises:
        RuntimeError: On fetch failure.
    """
    logger.debug("GET %s", url)
    req = Request(url, headers={"Accept": "application/json"})
    try:
        with urlopen(req, timeout=timeout) as resp:
            charset = resp.headers.get_content_charset() or "utf-8"
            body = resp.read().decode(charset)
            return json.loads(body)
    except HTTPError as exc:
        raise RuntimeError(f"HTTP {exc.code} for {url}: {exc.reason}") from exc
    except URLError as exc:
        raise RuntimeError(f"URL error for {url}: {exc}") from exc
    except json.JSONDecodeError as exc:
        raise RuntimeError(f"Invalid JSON from {url}: {exc}") from exc


def _fetch_with_retry(
    url: str,
    timeout: int = DEFAULT_TIMEOUT,
    retries: int = DEFAULT_MAX_RETRIES,
    delay: float = DEFAULT_RETRY_DELAY,
) -> Any:
    """Fetch JSON with retry and exponential backoff.

    Args:
        url: URL to fetch.
        timeout: Per-request timeout in seconds.
        retries: Maximum number of attempts.
        delay: Base delay in seconds between retries.

    Returns:
        Parsed JSON object.

    Raises:
        RuntimeError: After all retries exhausted.
    """
    last_exc: Optional[Exception] = None
    for attempt in range(1, retries + 1):
        try:
            return _fetch_json(url, timeout=timeout)
        except RuntimeError as exc:
            last_exc = exc
            logger.warning("Attempt %d/%d failed: %s", attempt, retries, exc)
            if attempt < retries:
                wait = delay * (2 ** (attempt - 1))
                logger.info("Retrying in %.1f s ...", wait)
                time.sleep(wait)
    raise RuntimeError(f"All {retries} retries exhausted: {last_exc}")


# ---------------------------------------------------------------------------
# ID loading
# ---------------------------------------------------------------------------

def _read_accessions(path: str) -> List[str]:
    """Read PXD accessions from a file, one per line.

    Args:
        path: Path to text file.

    Returns:
        List of accession strings.
    """
    accs: List[str] = []
    with open(path, "r", encoding="utf-8") as fh:
        for line in fh:
            tok = line.strip()
            if tok and not tok.startswith("#"):
                accs.append(tok)
    return accs


# ---------------------------------------------------------------------------
# PRIDE project metadata extraction
# ---------------------------------------------------------------------------

def _safe_get(data: Dict[str, Any], *keys: str, default: str = "") -> str:
    """Safely traverse nested dicts, returning default on any missing key.

    Args:
        data: Dict to traverse.
        *keys: Sequence of keys.
        default: Default value if any key is missing.

    Returns:
        Retrieved string value or default.
    """
    current: Any = data
    for key in keys:
        if isinstance(current, dict) and key in current:
            current = current[key]
        else:
            return default
    if current is None:
        return default
    return str(current).strip()


def extract_project_metadata(project_data: Dict[str, Any]) -> Dict[str, str]:
    """Extract key fields from a PRIDE project JSON response.

    Args:
        project_data: Parsed JSON dict from PRIDE project API.

    Returns:
        Flat dict of project-level metadata fields.
    """
    row: Dict[str, str] = {}

    row["accession"] = _safe_get(project_data, "accession")
    row["title"] = _safe_get(project_data, "title")
    row["project_description"] = _safe_get(project_data, "projectDescription")
    row["publication_date"] = _safe_get(project_data, "publicationDate")
    row["submission_date"] = _safe_get(project_data, "submissionDate")
    row["updated_date"] = _safe_get(project_data, "updatedDate")

    # Submitter information
    submitters = project_data.get("submitters", [])
    if submitters:
        first = submitters[0]
        row["submitter_name"] = _safe_get(first, "name")
        row["submitter_email"] = _safe_get(first, "email")
        row["submitter_affiliation"] = _safe_get(first, "affiliation")
    else:
        row["submitter_name"] = ""
        row["submitter_email"] = ""
        row["submitter_affiliation"] = ""

    # Lab PI
    lab_pis = project_data.get("labPIs", [])
    row["lab_pi"] = "; ".join(
        _safe_get(p, "name") for p in lab_pis if _safe_get(p, "name")
    )

    # Species / tissue / cell type / disease (CvParams)
    def _extract_cv_param_list(items: List[Dict[str, Any]]) -> str:
        names = []
        for item in items:
            name = item.get("name", item.get("value", ""))
            if name:
                names.append(str(name).strip())
        return "; ".join(names)

    row["species"] = _extract_cv_param_list(
        project_data.get("organisms", project_data.get("species", []))
    )
    row["tissue"] = _extract_cv_param_list(
        project_data.get("tissues", project_data.get("tissue", []))
    )
    row["cell_type"] = _extract_cv_param_list(
        project_data.get("cellTypes", project_data.get("cellType", []))
    )
    row["disease"] = _extract_cv_param_list(
        project_data.get("diseases", project_data.get("disease", []))
    )
    row["instrument"] = _extract_cv_param_list(
        project_data.get("instrumentNames", project_data.get("instruments", []))
    )

    # Quantification and modification methods
    row["quantification_method"] = _extract_cv_param_list(
        project_data.get("quantificationMethods", [])
    )
    row["modification"] = _extract_cv_param_list(
        project_data.get("ptmNames", project_data.get("modifications", []))
    )

    # Experiment type
    row["experiment_type"] = _extract_cv_param_list(
        project_data.get("experimentTypes", [])
    )

    # Keywords
    keywords = project_data.get("keywords", [])
    if isinstance(keywords, list):
        row["keywords"] = "; ".join(keywords)
    elif isinstance(keywords, str):
        row["keywords"] = keywords
    else:
        row["keywords"] = ""

    # Doi / publication
    row["doi"] = _safe_get(project_data, "doi")
    pubmed = project_data.get("pubmedIds", [])
    if isinstance(pubmed, list):
        row["pubmed_ids"] = "; ".join(str(p) for p in pubmed)
    else:
        row["pubmed_ids"] = ""

    # Links
    row["project_url"] = _safe_get(project_data, "_links", "self", "href")

    # File counts
    row["num_files"] = str(project_data.get("numFiles", ""))

    return row


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main() -> None:
    """Main entry point."""
    parser = argparse.ArgumentParser(
        description="Fetch PRIDE proteomics dataset metadata for PXD accessions.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""\
examples:
  %(prog)s -i pxd_accessions.txt -o pride_metadata.tsv
  %(prog)s -i ids.txt -o out.tsv --timeout 120 --retries 5
  %(prog)s -i ids.txt -o out.tsv --log-level DEBUG
""",
    )
    parser.add_argument(
        "-i", "--input-file", required=True,
        help="File with one PXD accession per line.",
    )
    parser.add_argument(
        "-o", "--output", default="-",
        help="Output TSV path (default: stdout).",
    )
    parser.add_argument(
        "--timeout", type=int, default=DEFAULT_TIMEOUT,
        help=f"HTTP timeout in seconds (default: {DEFAULT_TIMEOUT}).",
    )
    parser.add_argument(
        "--retries", type=int, default=DEFAULT_MAX_RETRIES,
        help=f"Retry attempts per request (default: {DEFAULT_MAX_RETRIES}).",
    )
    parser.add_argument(
        "--delay", type=float, default=1.0,
        help="Delay in seconds between requests (default: 1.0).",
    )
    parser.add_argument(
        "--log-level", default="INFO",
        choices=["DEBUG", "INFO", "WARNING", "ERROR", "CRITICAL"],
        help="Logging level (default: INFO).",
    )
    parser.add_argument(
        "--version", action="version", version=f"%(prog)s {__version__}"
    )
    args = parser.parse_args()

    logging.basicConfig(
        level=getattr(logging, args.log_level),
        format="%(asctime)s [%(levelname)s] %(name)s: %(message)s",
    )

    accessions = _read_accessions(args.input_file)
    if not accessions:
        logger.error("No accessions found in %s", args.input_file)
        sys.exit(1)

    logger.info("Fetching metadata for %d PXD accession(s).", len(accessions))

    all_rows: List[Dict[str, str]] = []
    errors: List[str] = []

    for idx, acc in enumerate(accessions, start=1):
        logger.info("[%d/%d] Fetching %s ...", idx, len(accessions), acc)
        url = PRIDE_PROJECT_API.format(accession=acc)
        try:
            data = _fetch_with_retry(url, timeout=args.timeout, retries=args.retries)
            if isinstance(data, dict):
                row = extract_project_metadata(data)
                all_rows.append(row)
                logger.info("  -> OK: %s", row.get("title", "")[:60])
            else:
                logger.warning("Unexpected response type for %s: %s", acc, type(data))
                errors.append(acc)
        except RuntimeError as exc:
            logger.error("  -> FAILED: %s", exc)
            errors.append(acc)

        if idx < len(accessions):
            time.sleep(args.delay)

    if not all_rows:
        logger.error("No metadata retrieved successfully.")
        sys.exit(1)

    # Collect fieldnames preserving order
    fieldnames_set: Dict[str, None] = {}
    for row in all_rows:
        for k in row:
            if k not in fieldnames_set:
                fieldnames_set[k] = None
    fieldnames = list(fieldnames_set.keys())

    if args.output == "-":
        writer = csv.DictWriter(
            sys.stdout, fieldnames=fieldnames, delimiter="\t", extrasaction="ignore"
        )
        writer.writeheader()
        writer.writerows(all_rows)
    else:
        out = Path(args.output)
        out.parent.mkdir(parents=True, exist_ok=True)
        with open(out, "w", newline="", encoding="utf-8") as fh:
            writer = csv.DictWriter(
                fh, fieldnames=fieldnames, delimiter="\t", extrasaction="ignore"
            )
            writer.writeheader()
            writer.writerows(all_rows)
        logger.info("Wrote %d rows to %s", len(all_rows), args.output)

    if errors:
        logger.warning(
            "%d accession(s) failed: %s", len(errors), ", ".join(errors)
        )


if __name__ == "__main__":
    main()
