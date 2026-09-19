#!/usr/bin/env python3
"""
Visualize MS2 spectra with theoretical b/y ion annotations from mzID+mzML/MGF.

Author: ChengYu
Created Time: 2026
"""

__version__ = "1.0.0"

import argparse
import logging
import re
import sys
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

import matplotlib.pyplot as plt
import numpy as np
from lxml import etree
from pyteomics import mass, mgf

try:
    import pymzml
except ImportError:
    pymzml = None

logger = logging.getLogger(__name__)

NS_HINT = "http://psidev.info/psi/pi/mzIdentML/1.2"
ION_OFFSET = 0.1
COLORS = {"b": "#FF3030", "y": "#4876FF"}


# ---------------------------------------------------------------------------
# mzIdentML parsing helpers
# ---------------------------------------------------------------------------

def get_namespace(tree: etree._ElementTree) -> Dict[str, str]:
    """Detect the XML namespace used in the mzIdentML file."""
    for elem in tree.iter():
        if elem.tag.startswith("{"):
            ns = elem.tag.split("}")[0][1:]
            return {"ns": ns}
    return {"ns": NS_HINT}


def parse_mzid(
    mzid_path: str,
    target_protein_id: str,
) -> List[Dict[str, Any]]:
    """Parse an mzIdentML file and return PSMs matching *target_protein_id*.

    Parameters
    ----------
    mzid_path : str
        Path to the .mzid file.
    target_protein_id : str
        Protein accession to search for (fuzzy match on prefix).

    Returns
    -------
    list[dict]
        Each dict contains keys: scan, sequence, charge, mods, score.
    """
    tree = etree.parse(mzid_path)
    ns = get_namespace(tree)

    clean_target = re.split(r"[.\s_]", target_protein_id)[0].upper()
    logger.debug("Cleaned target protein ID: %s", clean_target)

    peptide_map: Dict[str, str] = {}
    for pep in tree.xpath("//ns:Peptide", namespaces=ns):
        pid = pep.get("id")
        seq_elem = pep.xpath("ns:PeptideSequence/text()", namespaces=ns)
        peptide_map[pid] = str(seq_elem[0]) if seq_elem else ""

    psm_list: List[Dict[str, Any]] = []

    for sir in tree.xpath("//ns:SpectrumIdentificationResult", namespaces=ns):
        spectrum_id = sir.get("spectrumID", "")
        scan_match = re.search(r"scan=(\d+)", spectrum_id)
        scan_num = int(scan_match.group(1)) if scan_match else None
        if scan_num is None:
            try:
                scan_num = int(spectrum_id)
            except ValueError:
                continue

        for sii in sir.xpath(".//ns:SpectrumIdentificationItem", namespaces=ns):
            protein_refs = sii.xpath(
                ".//ns:ProteinAccession/text() | "
                ".//ns:DBSequence/@accession | "
                ".//ns:ProteinDetectionHypothesis/@id",
                namespaces=ns,
            )

            matched = False
            for ref in protein_refs:
                clean_ref = re.split(r"[.\s_]", str(ref))[0].upper()
                if clean_ref == clean_target:
                    matched = True
                    break
            if not matched:
                continue

            pep_ref = sii.get("peptide_ref", "")
            sequence = peptide_map.get(pep_ref, "")
            if not sequence:
                seq_text = sii.xpath("ns:PeptideSequence/text()", namespaces=ns)
                sequence = str(seq_text[0]) if seq_text else ""
            sequence = re.sub(r"\[.*?\]", "", sequence).upper()

            charge = int(sii.get("chargeState", 1))
            score = float(sii.get("MS-GF:QValue", sii.get("MS-GF:EValue", 0)))

            # Parse modifications
            mods: List[Tuple[int, float]] = []
            for mod in sii.xpath(".//ns:Modification", namespaces=ns):
                location = int(mod.get("location", 0))
                mono = mod.get("monoisotopicMassDelta", "0")
                mods.append((location, float(mono)))

            psm_list.append({
                "scan": scan_num,
                "sequence": sequence,
                "charge": charge,
                "mods": mods,
                "score": score,
            })

    logger.debug("Found %d PSMs", len(psm_list))
    return psm_list


# ---------------------------------------------------------------------------
# Fragment ion calculation
# ---------------------------------------------------------------------------

def calculate_fragments(
    sequence: str,
    charge: int,
    modifications: Optional[List[Tuple[int, float]]] = None,
) -> Dict[str, List[float]]:
    """Calculate theoretical b/y ion m/z values.

    Parameters
    ----------
    sequence : str
        Peptide amino-acid sequence (one-letter codes).
    charge : int
        Precursor charge state (fragment ions reported at z=1).
    modifications : list[tuple[int, float]], optional
        Each tuple is (position, mass_delta).

    Returns
    -------
    dict[str, list[float]]
        Keys ``'b'`` and ``'y'`` with sorted m/z lists.
    """
    if modifications is None:
        modifications = []

    frags: Dict[str, List[float]] = {"b": [], "y": []}

    for i in range(1, len(sequence)):
        # b-ion
        b_mass = sum(mass.std_aa_mass.get(aa, 0) for aa in sequence[:i]) + mass.std_ion_comp["b"]
        for pos, delta in modifications:
            if pos < i:
                b_mass += delta
        frags["b"].append(b_mass)

        # y-ion
        y_mass = sum(mass.std_aa_mass.get(aa, 0) for aa in sequence[-i:]) + mass.std_ion_comp["y"]
        for pos, delta in modifications:
            if pos >= (len(sequence) - i):
                y_mass += delta
        frags["y"].append(y_mass)

    return frags


# ---------------------------------------------------------------------------
# Spectrum loading
# ---------------------------------------------------------------------------

def load_spectrum_mgf(
    file_path: str, scan_num: int
) -> Optional[Tuple[np.ndarray, np.ndarray]]:
    """Load a single spectrum from an MGF file by scan number."""
    for spec in mgf.read(file_path):
        scans = spec["params"].get("scans", "")
        if str(scan_num) in str(scans).split(","):
            return np.array(spec["m/z array"]), np.array(spec["intensity array"])
    return None


def load_spectrum_mzml(
    file_path: str, scan_num: int
) -> Optional[Tuple[np.ndarray, np.ndarray]]:
    """Load a single spectrum from an mzML file by scan number."""
    if pymzml is None:
        logger.error("pymzml is required for mzML files. Install with: pip install pymzml")
        return None
    run = pymzml.run.Reader(file_path)
    for spec in run:
        if spec.ms_level == 2:
            spec_scan = spec.get("MS:1000016", spec.ID)
            try:
                if int(spec_scan) == scan_num or int(spec.ID) == scan_num:
                    return np.array(spec.mz), np.array(spec.i)
            except (ValueError, TypeError):
                continue
    return None


# ---------------------------------------------------------------------------
# Plotting
# ---------------------------------------------------------------------------

def plot_ms2(
    mz: np.ndarray,
    intensity: np.ndarray,
    frags: Dict[str, List[float]],
    scan_num: int,
    sequence: str,
    protein_id: str,
    output_path: str,
    tolerance: float = 0.1,
) -> None:
    """Draw a single annotated MS2 spectrum and save to *output_path*."""
    fig, ax = plt.subplots(figsize=(12, 6))

    ax.stem(mz, intensity, linefmt="gray", markerfmt=" ", basefmt=" ")

    for ion_type in ("b", "y"):
        for i, theo_mz in enumerate(frags[ion_type]):
            nearest_idx = int(np.argmin(np.abs(mz - theo_mz)))
            if abs(mz[nearest_idx] - theo_mz) < tolerance:
                label = f"{ion_type}{i + 1}"
                ax.vlines(
                    theo_mz, 0, intensity[nearest_idx],
                    colors=COLORS[ion_type], linestyles="dashed", alpha=0.7,
                )
                ax.text(
                    theo_mz + ION_OFFSET,
                    intensity[nearest_idx],
                    label,
                    color=COLORS[ion_type],
                    fontsize=8,
                )

    ax.set_title(f"Scan {scan_num} | {sequence}\nProtein: {protein_id}", fontsize=12)
    ax.set_xlabel("m/z")
    ax.set_ylabel("Intensity")

    seq_text = " ".join(list(sequence))
    fig.text(
        0.5, 0.01, seq_text,
        ha="center", fontsize=10,
        bbox=dict(facecolor="lightgray", alpha=0.5),
    )

    fig.savefig(output_path, dpi=300, bbox_inches="tight")
    plt.close(fig)
    logger.info("Saved: %s", output_path)


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="plot_ms2_spectrum.py",
        description="Visualize MS2 spectra with theoretical b/y ion annotations "
                    "from mzIdentML + mzML/MGF files.",
        epilog=(
            "Examples:\n"
            "  %(prog)s --mzid result.mzid --spectra data.mzML --protein P12345 -o spec\n"
            "  %(prog)s --mzid result.mzid --spectra data.mgf  --protein sp|P12345 -o spec\n"
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument("--version", action="version", version=f"%(prog)s {__version__}")
    parser.add_argument(
        "--log-level", default="INFO",
        choices=["DEBUG", "INFO", "WARNING", "ERROR", "CRITICAL"],
        help="Set the logging level (default: INFO).",
    )
    parser.add_argument("--mzid", required=True, help="Path to the mzIdentML (.mzid) file.")
    parser.add_argument("--spectra", required=True, help="Path to the MS/MS file (.mzML or .mgf).")
    parser.add_argument("--protein", required=True, help="Target protein accession ID.")
    parser.add_argument(
        "-o", "--output", default="ms2_spectrum",
        help="Output prefix (figures saved as <prefix>_<N>.png).",
    )
    parser.add_argument(
        "--tolerance", type=float, default=0.1,
        help="Mass tolerance (Da) for matching theoretical fragments.",
    )
    return parser


def main(argv: Optional[List[str]] = None) -> None:
    parser = build_parser()
    args = parser.parse_args(argv)

    logging.basicConfig(
        level=getattr(logging, args.log_level),
        format="%(asctime)s [%(levelname)s] %(name)s: %(message)s",
    )

    # 1. Parse mzIdentML
    psms = parse_mzid(args.mzid, args.protein)
    if not psms:
        logger.error("No PSMs found for protein '%s' in %s", args.protein, args.mzid)
        sys.exit(1)

    # 2. Iterate over PSMs and plot
    spectra_path = Path(args.spectra)
    ext = spectra_path.suffix.lower()

    for idx, psm in enumerate(psms, 1):
        if ext == ".mgf":
            result = load_spectrum_mgf(str(spectra_path), psm["scan"])
        elif ext in (".mzml", ".mzxml"):
            result = load_spectrum_mzml(str(spectra_path), psm["scan"])
        else:
            logger.error("Unsupported spectrum format: %s", ext)
            sys.exit(1)

        if result is None:
            logger.warning("Scan %d not found in %s", psm["scan"], spectra_path)
            continue
        mz, intensity = result

        frags = calculate_fragments(psm["sequence"], psm["charge"], psm["mods"])
        out_path = f"{args.output}_{idx}.png"
        plot_ms2(
            mz, intensity, frags,
            psm["scan"], psm["sequence"], args.protein,
            out_path, tolerance=args.tolerance,
        )

    logger.info("Done — %d spectra plotted.", len(psms))


if __name__ == "__main__":
    main()
