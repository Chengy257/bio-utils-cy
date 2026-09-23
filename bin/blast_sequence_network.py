#!/usr/bin/env python3
#########################################################################
# File Name: blast_sequence_network.py
# Author: ChengYu
# Description: BLAST-based protein sequence clustering and network
#              visualization using MCL (Markov Clustering Algorithm).
# Created Time: 2026
#########################################################################
"""BLAST-based protein sequence clustering and network visualization.

Runs BLASTP all-vs-all, builds a similarity matrix, applies MCL clustering,
and outputs cluster assignments, a SIF network file, and a publication-quality
network visualization PDF.

Requires: blastp, makeblastdb (NCBI BLAST+), markov-clustering, networkx, matplotlib
"""

import argparse
import logging
import os
import shutil
import sys
from pathlib import Path
from typing import Dict, List, Tuple

__version__ = "1.0.0"


def resolve_tool(name: str) -> str:
    """Resolve an external tool: MYS_<NAME>_BIN (config/env.sh) first, then PATH."""
    env_path = os.environ.get(f"MYS_{name.upper()}_BIN", "")
    if env_path:
        if os.path.isfile(env_path) and os.access(env_path, os.X_OK):
            return env_path
        raise FileNotFoundError(f"MYS_{name.upper()}_BIN is set but not executable: {env_path}")
    found = shutil.which(name)
    if not found:
        raise FileNotFoundError(
            f"{name} not found in PATH; set MYS_{name.upper()}_BIN in config/env.local.sh"
        )
    return found


def _check_dependencies() -> None:
    """Import and validate all required dependencies."""
    try:
        import matplotlib  # noqa: F401
        import markov_clustering  # noqa: F401
        import networkx  # noqa: F401
        import numpy  # noqa: F401
    except ImportError as e:
        import sys
        print(f"ERROR: Missing dependency: {e.name}", file=sys.stderr)
        print("Install with: pip install markov-clustering networkx numpy matplotlib biopython", file=sys.stderr)
        sys.exit(1)


def read_sequences(fasta_file: str) -> list:
    """Read sequences from a FASTA file.

    Args:
        fasta_file: Path to FASTA file.

    Returns:
        List of Bio.SeqRecord objects.
    """
    records = list(SeqIO.parse(fasta_file, "fasta"))
    if not records:
        logging.error("No sequences found in %s", fasta_file)
        sys.exit(1)
    logging.info("Loaded %d sequences from %s", len(records), fasta_file)
    return records


def run_blast(
    fasta_file: str,
    blast_db: str = "protein_db",
    output_file: str = "blast_output.txt",
    evalue: float = 1e-5,
    threads: int = 1,
) -> str:
    """Run BLASTP all-vs-all and return output file path.

    Args:
        fasta_file: Input protein FASTA.
        blast_db: BLAST database name prefix.
        output_file: BLAST tabular output path.
        evalue: E-value threshold.
        threads: Number of BLAST threads.

    Returns:
        Path to BLAST output file.
    """
    if not os.path.exists(f"{blast_db}.pin"):
        logging.info("Creating BLAST database: %s", blast_db)
        makeblastdb = resolve_tool("makeblastdb")
        cmd = f"{makeblastdb} -in {fasta_file} -dbtype prot -out {blast_db}"
        run(cmd, shell=True, check=True)

    logging.info("Running BLASTP (evalue=%g, threads=%d)...", evalue, threads)
    blastp = resolve_tool("blastp")
    cmd = (
        f"{blastp} -query {fasta_file} -db {blast_db} "
        f"-evalue {evalue} -outfmt 6 -out {output_file} "
        f"-num_threads {threads}"
    )
    run(cmd, shell=True, check=True)
    logging.info("BLAST output: %s", output_file)
    return output_file


def parse_blast_results(
    blast_file: str,
    seq_ids: List[str],
    identity_threshold: float = 0.5,
) -> Tuple["np.ndarray", List[str]]:
    """Parse BLAST tabular output into a similarity matrix.

    Args:
        blast_file: Path to BLAST -outfmt 6 file.
        seq_ids: Ordered list of sequence IDs.
        identity_threshold: Minimum percent identity (0-1) to include.

    Returns:
        Tuple of (similarity_matrix, seq_ids).
    """
    n = len(seq_ids)
    matrix = np.zeros((n, n))
    id_to_idx = {sid: i for i, sid in enumerate(seq_ids)}
    threshold_pct = identity_threshold * 100

    with open(blast_file) as fh:
        for line in fh:
            parts = line.strip().split()
            if len(parts) < 3:
                continue
            qseqid, sseqid, pident = parts[0], parts[1], float(parts[2])
            if qseqid in id_to_idx and sseqid in id_to_idx:
                i, j = id_to_idx[qseqid], id_to_idx[sseqid]
                if pident >= threshold_pct:
                    matrix[i, j] = matrix[j, i] = pident

    n_edges = np.count_nonzero(matrix) // 2
    logging.info(
        "Similarity matrix: %d sequences, %d edges (threshold=%.0f%%).",
        n, n_edges, threshold_pct,
    )
    return matrix, seq_ids


def run_mcl(matrix: "np.ndarray", inflation: float = 2.0) -> list:
    """Apply Markov Clustering Algorithm.

    Args:
        matrix: Similarity matrix.
        inflation: MCL inflation parameter (higher = more clusters).

    Returns:
        List of clusters, each cluster is a list of node indices.
    """
    result = mc.run_mcl(matrix, inflation=inflation)
    clusters = mc.get_clusters(result)
    logging.info("MCL produced %d clusters (inflation=%.1f).", len(clusters), inflation)
    return clusters


def save_clusters(clusters: list, seq_ids: List[str], output_file: str) -> None:
    """Save cluster assignments to a text file.

    Args:
        clusters: List of clusters (each a list of indices).
        seq_ids: Sequence ID list.
        output_file: Output file path.
    """
    with open(output_file, "w") as fh:
        for i, cluster in enumerate(clusters, 1):
            members = ", ".join(seq_ids[idx] for idx in cluster)
            fh.write(f"Cluster {i} ({len(cluster)} members): {members}\n")
    logging.info("Clusters saved to %s (%d clusters)", output_file, len(clusters))


def build_network(clusters: list, seq_ids: List[str]) -> "nx.Graph":
    """Build a networkx graph from MCL clusters.

    Args:
        clusters: List of clusters.
        seq_ids: Sequence ID list.

    Returns:
        NetworkX Graph with 'cluster' node attribute.
    """
    G = nx.Graph()
    for i, cluster in enumerate(clusters, 1):
        for idx in cluster:
            G.add_node(seq_ids[idx], cluster=i)
        for j in range(len(cluster)):
            for k in range(j + 1, len(cluster)):
                G.add_edge(seq_ids[cluster[j]], seq_ids[cluster[k]])

    logging.info("Network: %d nodes, %d edges.", G.number_of_nodes(), G.number_of_edges())
    return G


def save_sif(G: "nx.Graph", output_path: str) -> None:
    """Save network in SIF format.

    Args:
        G: NetworkX graph.
        output_path: Output SIF file path.
    """
    with open(output_path, "w") as fh:
        for u, v in G.edges():
            fh.write(f"{u}\tpp\t{v}\n")
    logging.info("SIF file saved to %s", output_path)


def visualize_network(
    G: "nx.Graph",
    output_path: str,
    figsize: Tuple[int, int] = (6, 6),
    dpi: int = 300,
) -> None:
    """Visualize cluster network.

    Args:
        G: NetworkX graph with 'cluster' node attribute.
        output_path: Output PDF path.
        figsize: Figure size.
        dpi: Resolution.
    """
    if G.number_of_nodes() == 0:
        logging.warning("Empty network, skipping visualization.")
        return

    clusters_attr = nx.get_node_attributes(G, "cluster")
    unique_clusters = sorted(set(clusters_attr.values()))
    n_clusters = len(unique_clusters)
    cmap = cm.get_cmap("tab20", max(n_clusters, 1))

    cluster_colors = {cid: cmap(i / max(n_clusters, 1)) for i, cid in enumerate(unique_clusters)}
    node_colors = [cluster_colors[clusters_attr[n]] for n in G.nodes()]

    pos = nx.spring_layout(G, seed=42, k=0.15, iterations=50)

    plt.figure(figsize=figsize, dpi=dpi)
    nx.draw_networkx_nodes(
        G, pos, node_color=node_colors,
        node_size=20, alpha=0.6, edgecolors="white", linewidths=0,
    )
    nx.draw_networkx_edges(G, pos, edge_color="gray", alpha=0.2, width=0.3)

    plt.title("Protein Sequence Clustering Network", fontsize=12)
    plt.axis("off")
    plt.subplots_adjust(left=0, right=1, top=1, bottom=0)
    plt.savefig(output_path, format="pdf", dpi=dpi, bbox_inches="tight")
    plt.close()
    logging.info("Network plot saved to %s", output_path)


def build_parser() -> argparse.ArgumentParser:
    """Build argument parser."""
    parser = argparse.ArgumentParser(
        description="BLAST-based protein sequence clustering and network visualization using MCL.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""\
examples:
  python blast_sequence_network.py -i proteins.fa -o result
  python blast_sequence_network.py -i proteins.fa -o result --evalue 1e-10 --inflation 2.5 --threads 4
""",
    )
    parser.add_argument(
        "-i", "--input", type=str, required=True,
        help="Input protein FASTA file.",
    )
    parser.add_argument(
        "-o", "--output-prefix", type=str, required=True,
        help="Output prefix for cluster, SIF, and PDF files.",
    )
    parser.add_argument(
        "--evalue", type=float, default=1e-5,
        help="BLAST e-value threshold (default: 1e-5).",
    )
    parser.add_argument(
        "--identity-threshold", type=float, default=0.5,
        help="Minimum percent identity (0-1) for network edges (default: 0.5).",
    )
    parser.add_argument(
        "--inflation", type=float, default=2.0,
        help="MCL inflation parameter; higher = more clusters (default: 2.0).",
    )
    parser.add_argument(
        "-t", "--threads", type=int, default=1,
        help="Number of BLAST threads (default: 1).",
    )
    parser.add_argument(
        "--log-level", type=str, default="INFO",
        choices=["DEBUG", "INFO", "WARNING", "ERROR"],
        help="Logging level (default: INFO).",
    )
    parser.add_argument(
        "--version", action="version", version=f"%(prog)s {__version__}",
    )
    return parser


def main() -> None:
    """Entry point."""
    parser = build_parser()
    args = parser.parse_args()

    logging.basicConfig(
        format="%(asctime)s - %(levelname)s - %(message)s",
        level=getattr(logging, args.log_level),
    )

    _check_dependencies()

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.cm as cm
    import matplotlib.pyplot as plt
    import markov_clustering as mc
    import networkx as nx
    import numpy as np
    from Bio import SeqIO
    from subprocess import run, CalledProcessError

    input_path = Path(args.input)
    if not input_path.is_file():
        logging.error("Input file not found: %s", args.input)
        sys.exit(1)

    # Check BLAST+ availability (and resolve configured binaries)
    try:
        blastp_bin = resolve_tool("blastp")
        makeblastdb_bin = resolve_tool("makeblastdb")
        run([blastp_bin, "-version"], capture_output=True, check=True)
    except (CalledProcessError, FileNotFoundError) as exc:
        logging.error("BLAST+ not available: %s", exc)
        sys.exit(1)

    # Read sequences
    records = read_sequences(args.input)
    seq_ids = [rec.id for rec in records]

    # Run BLAST
    blast_output = f"{args.output_prefix}_blast.txt"
    run_blast(
        args.input,
        blast_db=f"{args.output_prefix}_db",
        output_file=blast_output,
        evalue=args.evalue,
        threads=args.threads,
    )

    # Parse and cluster
    matrix, _ = parse_blast_results(blast_output, seq_ids, args.identity_threshold)
    clusters = run_mcl(matrix, inflation=args.inflation)

    # Save outputs
    save_clusters(clusters, seq_ids, f"{args.output_prefix}_clusters.txt")

    G = build_network(clusters, seq_ids)
    save_sif(G, f"{args.output_prefix}_network.sif")
    visualize_network(G, f"{args.output_prefix}_network.pdf")

    # Clean up BLAST temp files
    for ext in (".phr", ".pin", ".pdb", ".psq", ".pog", ".pos", ".pot", ".pto"):
        db_file = f"{args.output_prefix}_db{ext}"
        if os.path.exists(db_file):
            os.unlink(db_file)


if __name__ == "__main__":
    main()
