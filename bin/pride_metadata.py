#!/usr/bin/env python3
"""
File Name: pride_metadata.py
Author: ChengYu
Description: PRIDE proteomics metadata tool — merged from fetch_pride_metadata.py
             and parse_pride_json.py (which it replaces).
             Subcommands:
               fetch  Query the PRIDE API for PXD accessions -> merged TSV
               parse  Extract metadata from local PRIDE JSON file(s) -> CSV/TSV/XLSX
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
from typing import Any, Dict, List, Optional, Tuple
from urllib.error import HTTPError, URLError
from urllib.request import Request, urlopen

__version__ = "2.0.0"

PRIDE_PROJECT_API = "https://www.ebi.ac.uk/pride/ws/archive/v2/projects/{accession}"
DEFAULT_TIMEOUT = 60
DEFAULT_MAX_RETRIES = 3
DEFAULT_RETRY_DELAY = 2.0

logger = logging.getLogger(__name__)


# ---------------------------------------------------------------------------
# Shared helpers
# ---------------------------------------------------------------------------

def _safe_get(data: Dict[str, Any], *keys: str, default: str = "") -> str:
    """Safely traverse nested dicts, returning default on any missing key."""
    current: Any = data
    for key in keys:
        if isinstance(current, dict) and key in current:
            current = current[key]
        else:
            return default
    if current is None:
        return default
    return str(current).strip()


def _join_cv_list(items: List[Dict[str, Any]], key: str = "name") -> str:
    """Join CV parameter values from a list of dicts (or plain strings)."""
    values: List[str] = []
    if not isinstance(items, list):
        return ""
    for item in items:
        if isinstance(item, dict):
            val = item.get(key, item.get("value", ""))
            if val:
                values.append(str(val).strip())
        elif isinstance(item, str):
            values.append(item.strip())
    return "; ".join(values)


def _extract_instruments(project_data: Dict[str, Any]) -> str:
    """Extract instrument names, checking multiple JSON locations."""
    names = project_data.get("instrumentNames", [])
    if names:
        if isinstance(names, list):
            return "; ".join(str(n).strip() for n in names if n)
        return str(names)

    instruments = project_data.get("instruments", [])
    if instruments:
        return _join_cv_list(instruments)

    # Fall back to additionalAttributes entries mentioning "instrument"
    attrs = project_data.get("additionalAttributes", [])
    for attr in attrs:
        if isinstance(attr, dict):
            cv = attr.get("cvParam", attr)
            name = cv.get("name", "")
            if name and "instrument" in name.lower():
                val = cv.get("value", "")
                if val:
                    return str(val)
    return ""


def _keywords_to_str(keywords: Any) -> str:
    """Normalize the keywords field (list or string) to one string."""
    if isinstance(keywords, list):
        return "; ".join(str(k) for k in keywords)
    if isinstance(keywords, str):
        return keywords
    return ""


def extract_project_metadata(project_data: Dict[str, Any]) -> Dict[str, str]:
    """Extract the unified metadata schema from a PRIDE project JSON dict.

    This is the single schema used by both the ``fetch`` and ``parse``
    subcommands.
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
    if submitters and isinstance(submitters, list):
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
    row["species"] = _join_cv_list(
        project_data.get("organisms", project_data.get("species", []))
    )
    row["tissue"] = _join_cv_list(
        project_data.get("tissues", project_data.get("tissue", []))
    )
    row["cell_type"] = _join_cv_list(
        project_data.get("cellTypes", project_data.get("cellType", []))
    )
    row["disease"] = _join_cv_list(
        project_data.get("diseases", project_data.get("disease", []))
    )
    row["instrument"] = _extract_instruments(project_data)

    # Quantification and modification methods
    row["quantification_method"] = _join_cv_list(
        project_data.get("quantificationMethods", [])
    )
    row["modification"] = _join_cv_list(
        project_data.get("ptmNames", project_data.get("modifications", []))
    )
    row["experiment_type"] = _join_cv_list(
        project_data.get("experimentTypes", [])
    )

    # Keywords / references / links / counts
    row["keywords"] = _keywords_to_str(project_data.get("keywords", []))
    row["doi"] = _safe_get(project_data, "doi")
    pubmed = project_data.get("pubmedIds", [])
    row["pubmed_ids"] = (
        "; ".join(str(p) for p in pubmed) if isinstance(pubmed, list) else ""
    )
    row["project_url"] = _safe_get(project_data, "_links", "self", "href")
    row["num_files"] = str(project_data.get("numFiles", ""))

    return row


def _collect_fieldnames(rows: List[Dict[str, str]]) -> List[str]:
    """Collect fieldnames across rows, preserving first-seen order."""
    fieldnames_set: Dict[str, None] = {}
    for row in rows:
        for k in row:
            if k not in fieldnames_set:
                fieldnames_set[k] = None
    return list(fieldnames_set.keys())


# ---------------------------------------------------------------------------
# fetch subcommand — PRIDE API
# ---------------------------------------------------------------------------

def _fetch_json(url: str, timeout: int = DEFAULT_TIMEOUT) -> Any:
    """Fetch JSON content from a URL."""
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
    """Fetch JSON with retry and exponential backoff."""
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


def _read_accessions(path: str) -> List[str]:
    """Read PXD accessions from a file, one per line (skips blanks/comments)."""
    accs: List[str] = []
    with open(path, "r", encoding="utf-8") as fh:
        for line in fh:
            tok = line.strip()
            if tok and not tok.startswith("#"):
                accs.append(tok)
    return accs


def cmd_fetch(args: argparse.Namespace) -> None:
    """Fetch metadata for PXD accessions from the PRIDE API."""
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

    _write_table(all_rows, args.output, fmt="tsv")

    if errors:
        logger.warning("%d accession(s) failed: %s", len(errors), ", ".join(errors))


# ---------------------------------------------------------------------------
# parse subcommand — local PRIDE JSON files
# ---------------------------------------------------------------------------

def find_json_files(paths: List[str], recursive: bool = True) -> List[Path]:
    """Find JSON files in the given files/directories."""
    json_files: List[Path] = []
    for p in paths:
        path = Path(p)
        if path.is_file():
            if path.suffix.lower() == ".json":
                json_files.append(path)
            else:
                logger.warning("Skipping non-JSON file: %s", path)
        elif path.is_dir():
            found = sorted(path.rglob("*.json")) if recursive else sorted(path.glob("*.json"))
            json_files.extend(found)
            logger.info("Found %d JSON file(s) in %s", len(found), path)
        else:
            logger.warning("Path not found: %s", p)
    return json_files


def parse_json_file(filepath: Path) -> Tuple[Optional[Dict[str, str]], Optional[str]]:
    """Parse a single PRIDE JSON file and extract metadata.

    Returns:
        Tuple of (extracted row dict or None, error message or None).
    """
    try:
        with open(filepath, "r", encoding="utf-8") as fh:
            data = json.load(fh)
    except json.JSONDecodeError as exc:
        return None, f"Invalid JSON in {filepath}: {exc}"
    except OSError as exc:
        return None, f"Cannot read {filepath}: {exc}"

    if not isinstance(data, dict):
        return None, f"Expected dict in {filepath}, got {type(data).__name__}"

    def _row_from(obj: Dict[str, Any]) -> Dict[str, str]:
        row = extract_project_metadata(obj)
        row["source_file"] = filepath.name
        return row

    # Project data at the root
    if "accession" in data:
        return _row_from(data), None

    # Common wrapper keys; lists of projects are merged into one row
    for wrapper_key in ("project", "projects", "data"):
        if wrapper_key in data:
            inner = data[wrapper_key]
            if isinstance(inner, dict):
                return _row_from(inner), None
            if isinstance(inner, list):
                rows_found = [
                    _row_from(item)
                    for item in inner
                    if isinstance(item, dict) and "accession" in item
                ]
                if rows_found:
                    if len(rows_found) == 1:
                        return rows_found[0], None
                    merged: Dict[str, str] = {}
                    for r in rows_found:
                        for k, v in r.items():
                            if k in merged and merged[k]:
                                if v and v != merged[k]:
                                    merged[k] = f"{merged[k]}; {v}"
                            elif v:
                                merged[k] = v
                            else:
                                merged.setdefault(k, "")
                    return merged, None

    # Fallback: extract whatever we can
    row = _row_from(data)
    if row.get("accession") or row.get("title"):
        return row, None

    return None, f"No recognizable PRIDE project data in {filepath}"


def write_xlsx(rows: List[Dict[str, str]], path: str) -> None:
    """Write rows to an XLSX file (falls back to CSV without openpyxl)."""
    try:
        from openpyxl import Workbook
    except ImportError:
        csv_path = path.rsplit(".", 1)[0] + ".csv"
        logger.warning("openpyxl not installed; falling back to CSV: %s", csv_path)
        _write_table(rows, csv_path, fmt="csv")
        return

    fieldnames = _collect_fieldnames(rows)
    wb = Workbook()
    ws = wb.active
    ws.title = "PRIDE Metadata"
    ws.append(fieldnames)
    for row in rows:
        ws.append([row.get(f, "") for f in fieldnames])

    out = Path(path)
    out.parent.mkdir(parents=True, exist_ok=True)
    wb.save(str(out))
    logger.info("Wrote %d rows to %s (xlsx)", len(rows), path)


def _write_table(rows: List[Dict[str, str]], output: str, fmt: str = "tsv") -> None:
    """Write rows as CSV/TSV to a file path or stdout ('-')."""
    if not rows:
        logger.warning("No rows to write.")
        return
    fieldnames = _collect_fieldnames(rows)
    delimiter = "\t" if fmt == "tsv" else ","
    if output == "-":
        writer = csv.DictWriter(
            sys.stdout, fieldnames=fieldnames, delimiter=delimiter, extrasaction="ignore"
        )
        writer.writeheader()
        writer.writerows(rows)
        return
    out = Path(output)
    out.parent.mkdir(parents=True, exist_ok=True)
    with open(out, "w", newline="", encoding="utf-8") as fh:
        writer = csv.DictWriter(
            fh, fieldnames=fieldnames, delimiter=delimiter, extrasaction="ignore"
        )
        writer.writeheader()
        writer.writerows(rows)
    logger.info("Wrote %d rows to %s", len(rows), output)


def cmd_parse(args: argparse.Namespace) -> None:
    """Extract metadata from local PRIDE JSON files."""
    json_files = find_json_files(args.inputs, recursive=not args.no_recursive)
    if not json_files:
        logger.error("No JSON files found in: %s", ", ".join(args.inputs))
        sys.exit(1)

    logger.info("Processing %d JSON file(s).", len(json_files))

    rows: List[Dict[str, str]] = []
    errors: List[str] = []
    for jf in json_files:
        row, error = parse_json_file(jf)
        if row:
            rows.append(row)
            logger.info("  OK: %s (%s)", jf.name, row.get("accession", "N/A"))
        if error:
            errors.append(error)
            logger.warning("  SKIP: %s", error)

    if not rows:
        logger.error("No valid PRIDE metadata extracted from any file.")
        sys.exit(1)

    if args.format == "xlsx":
        out_path = args.output if args.output != "-" else "pride_metadata.xlsx"
        write_xlsx(rows, out_path)
    else:
        _write_table(rows, args.output, fmt=args.format)

    if errors:
        logger.warning("%d file(s) had errors:\n  %s", len(errors), "\n  ".join(errors))
    logger.info("Done. %d/%d file(s) parsed successfully.", len(rows), len(json_files))


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def build_parser() -> argparse.ArgumentParser:
    """Build the argument parser with fetch/parse subcommands."""
    parser = argparse.ArgumentParser(
        description="PRIDE proteomics metadata tool: fetch from the PRIDE API "
                    "or parse local PRIDE JSON files.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""\
examples:
  %(prog)s fetch -i pxd_accessions.txt -o pride_metadata.tsv
  %(prog)s fetch -i ids.txt -o out.tsv --timeout 120 --retries 5
  %(prog)s parse ./pride_jsons/ -o results.csv
  %(prog)s parse file1.json file2.json -o output.xlsx -f xlsx
""",
    )
    parser.add_argument("--version", action="version", version=f"%(prog)s {__version__}")
    sub = parser.add_subparsers(dest="command", required=True)

    p_fetch = sub.add_parser(
        "fetch", help="Query the PRIDE API for PXD accessions.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    p_fetch.add_argument("-i", "--input-file", required=True,
                         help="File with one PXD accession per line.")
    p_fetch.add_argument("-o", "--output", default="-",
                         help="Output TSV path (default: stdout).")
    p_fetch.add_argument("--timeout", type=int, default=DEFAULT_TIMEOUT,
                         help=f"HTTP timeout in seconds (default: {DEFAULT_TIMEOUT}).")
    p_fetch.add_argument("--retries", type=int, default=DEFAULT_MAX_RETRIES,
                         help=f"Retry attempts per request (default: {DEFAULT_MAX_RETRIES}).")
    p_fetch.add_argument("--delay", type=float, default=1.0,
                         help="Delay in seconds between requests (default: 1.0).")
    p_fetch.add_argument("--log-level", default="INFO",
                         choices=["DEBUG", "INFO", "WARNING", "ERROR", "CRITICAL"])
    p_fetch.add_argument("--version", action="version", version=f"%(prog)s {__version__}")
    p_fetch.set_defaults(func=cmd_fetch)

    p_parse = sub.add_parser(
        "parse", help="Extract metadata from local PRIDE JSON files.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    p_parse.add_argument("inputs", nargs="+",
                         help="Input JSON file(s) or directory(ies).")
    p_parse.add_argument("-o", "--output", default="-",
                         help="Output file path (default: stdout as CSV).")
    p_parse.add_argument("-f", "--format", choices=["csv", "tsv", "xlsx"],
                         default="csv", help="Output format (default: csv).")
    p_parse.add_argument("--no-recursive", action="store_true",
                         help="Do not search directories recursively.")
    p_parse.add_argument("--log-level", default="INFO",
                         choices=["DEBUG", "INFO", "WARNING", "ERROR", "CRITICAL"])
    p_parse.add_argument("--version", action="version", version=f"%(prog)s {__version__}")
    p_parse.set_defaults(func=cmd_parse)

    return parser


def main() -> None:
    """Main entry point."""
    parser = build_parser()
    args = parser.parse_args()

    logging.basicConfig(
        level=getattr(logging, args.log_level),
        format="%(asctime)s [%(levelname)s] %(name)s: %(message)s",
    )
    args.func(args)


if __name__ == "__main__":
    main()
