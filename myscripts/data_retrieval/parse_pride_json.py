#!/usr/bin/env python3
"""
File Name: parse_pride_json.py
Author: ChengYu
Description: Parse PRIDE JSON metadata files. Extracts instruments, core info,
             and other metadata from PRIDE project JSON files. Supports recursive
             directory scanning and outputs to CSV or XLSX format.
Created Time: 2026
"""

from __future__ import annotations

import argparse
import csv
import json
import logging
import os
import sys
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

__version__ = "1.0.0"

logger = logging.getLogger(__name__)

# ---------------------------------------------------------------------------
# JSON value extraction helpers
# ---------------------------------------------------------------------------

def _safe_get(data: Any, *keys: str, default: str = "") -> str:
    """Safely traverse nested dicts/objects and return a string value.

    Args:
        data: Data structure to traverse.
        *keys: Sequence of keys to follow.
        default: Default if key path fails.

    Returns:
        Extracted string or default.
    """
    current = data
    for key in keys:
        if isinstance(current, dict) and key in current:
            current = current[key]
        else:
            return default
    if current is None:
        return default
    return str(current).strip()


def _join_cv_list(items: List[Dict[str, Any]], key: str = "name") -> str:
    """Join CV parameter values from a list of dicts.

    Args:
        items: List of CV parameter dicts.
        key: Key to extract from each dict (default: 'name').

    Returns:
        Semicolon-separated string of values.
    """
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
    """Extract instrument names from PRIDE project JSON.

    Checks multiple possible locations for instrument data.

    Args:
        project_data: Parsed PRIDE project JSON dict.

    Returns:
        Semicolon-separated instrument names.
    """
    # Try instrumentNames first
    names = project_data.get("instrumentNames", [])
    if names:
        if isinstance(names, list):
            return "; ".join(str(n).strip() for n in names if n)
        return str(names)

    # Try instruments as CV params
    instruments = project_data.get("instruments", [])
    if instruments:
        return _join_cv_list(instruments)

    # Try nested in additionalAttributes
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


def _extract_core_info(project_data: Dict[str, Any]) -> Dict[str, str]:
    """Extract core project information from PRIDE JSON.

    Args:
        project_data: Parsed PRIDE project JSON dict.

    Returns:
        Dict of standardized core metadata fields.
    """
    row: Dict[str, str] = {}

    row["accession"] = _safe_get(project_data, "accession")
    row["title"] = _safe_get(project_data, "title")
    row["description"] = _safe_get(project_data, "projectDescription")
    row["publication_date"] = _safe_get(project_data, "publicationDate")
    row["submission_date"] = _safe_get(project_data, "submissionDate")
    row["updated_date"] = _safe_get(project_data, "updatedDate")

    # Submitter
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

    # Lab PIs
    pis = project_data.get("labPIs", [])
    row["lab_pi"] = "; ".join(
        _safe_get(p, "name") for p in pis if _safe_get(p, "name")
    ) if pis else ""

    # Biological context
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

    # Instrument
    row["instrument"] = _extract_instruments(project_data)

    # Methods
    row["quantification_method"] = _join_cv_list(
        project_data.get("quantificationMethods", [])
    )
    row["modification"] = _join_cv_list(
        project_data.get("ptmNames", project_data.get("modifications", []))
    )
    row["experiment_type"] = _join_cv_list(
        project_data.get("experimentTypes", [])
    )

    # Keywords
    kw = project_data.get("keywords", [])
    if isinstance(kw, list):
        row["keywords"] = "; ".join(kw)
    elif isinstance(kw, str):
        row["keywords"] = kw
    else:
        row["keywords"] = ""

    # References
    row["doi"] = _safe_get(project_data, "doi")
    pubmed = project_data.get("pubmedIds", [])
    row["pubmed_ids"] = "; ".join(str(p) for p in pubmed) if isinstance(pubmed, list) else ""

    # File counts
    row["num_files"] = str(project_data.get("numFiles", ""))

    # Links
    row["project_url"] = _safe_get(project_data, "_links", "self", "href")

    return row


# ---------------------------------------------------------------------------
# File discovery
# ---------------------------------------------------------------------------

def find_json_files(
    paths: List[str],
    recursive: bool = True,
) -> List[Path]:
    """Find JSON files in the given paths.

    Args:
        paths: List of file or directory paths.
        recursive: Whether to search directories recursively.

    Returns:
        Sorted list of Path objects for discovered JSON files.
    """
    json_files: List[Path] = []
    for p in paths:
        path = Path(p)
        if path.is_file():
            if path.suffix.lower() == ".json":
                json_files.append(path)
            else:
                logger.warning("Skipping non-JSON file: %s", path)
        elif path.is_dir():
            if recursive:
                found = sorted(path.rglob("*.json"))
            else:
                found = sorted(path.glob("*.json"))
            json_files.extend(found)
            logger.info("Found %d JSON file(s) in %s", len(found), path)
        else:
            logger.warning("Path not found: %s", p)
    return json_files


# ---------------------------------------------------------------------------
# Parse a single JSON file
# ---------------------------------------------------------------------------

def parse_json_file(filepath: Path) -> Tuple[Optional[Dict[str, str]], Optional[str]]:
    """Parse a single PRIDE JSON file and extract metadata.

    Args:
        filepath: Path to the JSON file.

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

    # Handle PRIDE project-level JSON
    # Some files have the project data at the root, others wrap it
    if "accession" in data:
        row = _extract_core_info(data)
        row["source_file"] = filepath.name
        return row, None

    # Try common wrapper keys
    for wrapper_key in ("project", "projects", "data"):
        if wrapper_key in data:
            inner = data[wrapper_key]
            if isinstance(inner, dict):
                row = _extract_core_info(inner)
                row["source_file"] = filepath.name
                return row, None
            elif isinstance(inner, list) and inner:
                rows_found: List[Dict[str, str]] = []
                for item in inner:
                    if isinstance(item, dict) and "accession" in item:
                        r = _extract_core_info(item)
                        r["source_file"] = filepath.name
                        rows_found.append(r)
                if rows_found:
                    # Return first; store extras via side channel is complex,
                    # so we merge all into one row with semicolons for multi
                    if len(rows_found) == 1:
                        return rows_found[0], None
                    # Merge multiple projects into semicolon-separated values
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
                    merged["source_file"] = filepath.name
                    return merged, None

    # Fallback: try to extract whatever we can
    row = _extract_core_info(data)
    if row.get("accession") or row.get("title"):
        row["source_file"] = filepath.name
        return row, None

    return None, f"No recognizable PRIDE project data in {filepath}"


# ---------------------------------------------------------------------------
# Output writers
# ---------------------------------------------------------------------------

def write_csv(rows: List[Dict[str, str]], path: str) -> None:
    """Write rows to a CSV/TSV file.

    Args:
        rows: List of row dicts.
        path: Output file path.
    """
    if not rows:
        logger.warning("No rows to write.")
        return
    # Collect ordered fieldnames
    fieldnames_set: Dict[str, None] = {}
    for row in rows:
        for k in row:
            if k not in fieldnames_set:
                fieldnames_set[k] = None
    fieldnames = list(fieldnames_set.keys())

    out = Path(path)
    out.parent.mkdir(parents=True, exist_ok=True)
    with open(out, "w", newline="", encoding="utf-8") as fh:
        writer = csv.DictWriter(fh, fieldnames=fieldnames, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)
    logger.info("Wrote %d rows to %s", len(rows), path)


def write_xlsx(rows: List[Dict[str, str]], path: str) -> None:
    """Write rows to an XLSX file using openpyxl.

    Falls back to CSV if openpyxl is not available.

    Args:
        rows: List of row dicts.
        path: Output file path.
    """
    try:
        from openpyxl import Workbook
    except ImportError:
        csv_path = path.rsplit(".", 1)[0] + ".csv"
        logger.warning(
            "openpyxl not installed; falling back to CSV: %s", csv_path
        )
        write_csv(rows, csv_path)
        return

    if not rows:
        logger.warning("No rows to write.")
        return

    fieldnames_set: Dict[str, None] = {}
    for row in rows:
        for k in row:
            if k not in fieldnames_set:
                fieldnames_set[k] = None
    fieldnames = list(fieldnames_set.keys())

    wb = Workbook()
    ws = wb.active
    ws.title = "PRIDE Metadata"

    # Header
    ws.append(fieldnames)

    # Data
    for row in rows:
        ws.append([row.get(f, "") for f in fieldnames])

    out = Path(path)
    out.parent.mkdir(parents=True, exist_ok=True)
    wb.save(str(out))
    logger.info("Wrote %d rows to %s (xlsx)", len(rows), path)


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def main() -> None:
    """Main entry point."""
    parser = argparse.ArgumentParser(
        description="Parse PRIDE JSON metadata files and output structured tables.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""\
examples:
  %(prog)s ./pride_jsons/ -o results.csv
  %(prog)s file1.json file2.json -o output.xlsx -f xlsx
  %(prog)s ./data/ -o out.csv --no-recursive
  %(prog)s ./data/ -o out.csv --log-level DEBUG
""",
    )
    parser.add_argument(
        "inputs",
        nargs="+",
        help="Input JSON file(s) or directory(ies) to parse.",
    )
    parser.add_argument(
        "-o", "--output",
        default="-",
        help="Output file path (default: stdout as CSV).",
    )
    parser.add_argument(
        "-f", "--format",
        choices=["csv", "tsv", "xlsx"],
        default="csv",
        help="Output format (default: csv).",
    )
    parser.add_argument(
        "--no-recursive",
        action="store_true",
        help="Do not search directories recursively.",
    )
    parser.add_argument(
        "--log-level",
        default="INFO",
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

    # Discover JSON files
    json_files = find_json_files(args.inputs, recursive=not args.no_recursive)
    if not json_files:
        logger.error("No JSON files found in: %s", ", ".join(args.inputs))
        sys.exit(1)

    logger.info("Processing %d JSON file(s).", len(json_files))

    # Parse files
    rows: List[Dict[str, str]] = []
    errors: List[str] = []

    for jf in json_files:
        logger.debug("Parsing: %s", jf)
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

    # Determine output format
    if args.format == "xlsx":
        out_path = args.output if args.output != "-" else "pride_metadata.xlsx"
        write_xlsx(rows, out_path)
    elif args.format == "tsv":
        if args.output == "-":
            # stdout as TSV
            fieldnames_set: Dict[str, None] = {}
            for row in rows:
                for k in row:
                    if k not in fieldnames_set:
                        fieldnames_set[k] = None
            fieldnames = list(fieldnames_set.keys())
            writer = csv.DictWriter(
                sys.stdout, fieldnames=fieldnames, delimiter="\t", extrasaction="ignore"
            )
            writer.writeheader()
            writer.writerows(rows)
        else:
            # Write TSV to file by temporarily adjusting the csv writer
            if not args.output.endswith(".tsv"):
                out_path = args.output
            else:
                out_path = args.output
            out_p = Path(out_path)
            out_p.parent.mkdir(parents=True, exist_ok=True)
            fieldnames_set2: Dict[str, None] = {}
            for row in rows:
                for k in row:
                    if k not in fieldnames_set2:
                        fieldnames_set2[k] = None
            fieldnames2 = list(fieldnames_set2.keys())
            with open(out_p, "w", newline="", encoding="utf-8") as fh:
                writer = csv.DictWriter(
                    fh, fieldnames=fieldnames2, delimiter="\t", extrasaction="ignore"
                )
                writer.writeheader()
                writer.writerows(rows)
            logger.info("Wrote %d rows to %s (tsv)", len(rows), out_path)
    else:
        # csv format
        if args.output == "-":
            fieldnames_set3: Dict[str, None] = {}
            for row in rows:
                for k in row:
                    if k not in fieldnames_set3:
                        fieldnames_set3[k] = None
            fieldnames3 = list(fieldnames_set3.keys())
            writer = csv.DictWriter(
                sys.stdout, fieldnames=fieldnames3, extrasaction="ignore"
            )
            writer.writeheader()
            writer.writerows(rows)
        else:
            write_csv(rows, args.output)

    if errors:
        logger.warning(
            "%d file(s) had errors:\n  %s",
            len(errors),
            "\n  ".join(errors),
        )
    logger.info("Done. %d/%d file(s) parsed successfully.", len(rows), len(json_files))


if __name__ == "__main__":
    main()
