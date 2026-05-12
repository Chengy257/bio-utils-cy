#!/usr/bin/env python3
#########################################################################
# File Name: sequence_similarity_network.py
# Author: ChengYu
# Description: Build and visualize sequence similarity networks using
#              k-mer Jaccard similarity and Louvain community detection.
# Created Time: 2026
#########################################################################
"""Build and visualize sequence similarity networks from FASTA files.

Uses k-mer based Jaccard similarity to construct a sparse similarity
matrix, then applies Louvain community detection for clustering. Outputs
a publication-quality network PDF and SIF format network file.

Requires: networkx, python-louvain, scikit-learn, scipy, biopython, matplotlib
"""

import argparse
import logging
import sys
from pathlib import Path
from typing import List, Tuple

__version__ = "1.0.0"


def _check_dependencies() -> None:
    """Import and validate all required dependencies."""
    try:
        import matplotlib
        matplotlib.use("Agg")
        import community.community_louvain  # noqa: F401
        import networkx  # noqa: F401
        import numpy  # noqa: F401
        import scipy  # noqa: F401
        import sklearn  # noqa: F401
    except ImportError as e:
        import sys
        print(f"ERROR: Missing dependency: {e.name}", file=sys.stderr)
        print("Install with: pip install networkx python-louvain scikit-learn scipy biopython matplotlib", file=sys.stderr)
        sys.exit(1)


def get_kmers(sequence: str, k: int) -> List[str]:
    """Convert a sequence into a list of k-mers.

    Args:
        sequence: Input sequence string.
        k: K-mer size.

    Returns:
        List of k-mer strings.
    """
    return [sequence[i:i + k] for i in range(len(sequence) - k + 1)]


def compute_similarity_matrix(
    sequences: List[str],
    k: int = 5,
    threshold: float = 0.5,
) -> "csr_matrix":
    """Compute sparse Jaccard similarity matrix from k-mer profiles.

    Args:
        sequences: List of sequence strings.
        k: K-mer size for profile construction.
        threshold: Minimum similarity to retain (values below are zeroed).

    Returns:
        Sparse similarity matrix.
    """
    if len(sequences) < 2:
        logging.error("Need at least 2 sequences for network construction.")
        sys.exit(1)

    vectorizer = CountVectorizer(analyzer=lambda seq: get_kmers(seq, k))
    kmer_matrix = vectorizer.fit_transform(sequences)

    jaccard_dist = pdist(kmer_matrix.toarray(), metric="jaccard")
    sim_matrix = squareform(1 - jaccard_dist)
    np.fill_diagonal(sim_matrix, 0)

    sim_matrix[sim_matrix < threshold] = 0
    sparse_matrix = csr_matrix(sim_matrix)

    n_edges = sparse_matrix.nnz // 2
    logging.info(
        "Similarity matrix: %d sequences, %d edges (threshold=%.2f, k=%d).",
        len(sequences), n_edges, threshold, k,
    )
    return sparse_matrix


def build_network(similarity_matrix: "csr_matrix", threshold: float) -> "nx.Graph":
    """Build a networkx graph from a sparse similarity matrix.

    Args:
        similarity_matrix: Sparse similarity matrix.
        threshold: Minimum edge weight.

    Returns:
        NetworkX Graph with 'weight' edge attribute.
    """
    G = nx.Graph()
    coo = similarity_matrix.tocoo()

    for i, j, weight in zip(coo.row, coo.col, coo.data):
        if weight >= threshold and i != j:
            G.add_edge(i, j, weight=float(weight))

    logging.info("Network: %d nodes, %d edges.", G.number_of_nodes(), G.number_of_edges())

    if G.number_of_nodes() == 0:
        logging.warning("Empty network — consider lowering the threshold.")

    return G


def detect_communities(G: "nx.Graph") -> dict:
    """Apply Louvain community detection.

    Args:
        G: NetworkX graph.

    Returns:
        Dict mapping node -> community id.
    """
    if G.number_of_nodes() == 0:
        return {}
    partition = community_louvain.best_partition(G, weight="weight")
    n_communities = len(set(partition.values()))
    logging.info("Detected %d communities.", n_communities)
    return partition


def save_sif(G: "nx.Graph", seq_ids: List[str], output_path: str) -> None:
    """Save network in SIF format.

    Args:
        G: NetworkX graph.
        seq_ids: Sequence ID list (node indices map to this).
        output_path: Output SIF file path.
    """
    with open(output_path, "w") as fh:
        for i, j, data in G.edges(data=True):
            fh.write(f"{seq_ids[i]}\tsimilarity\t{seq_ids[j]}\n")
    logging.info("SIF network saved to %s", output_path)


def visualize_network(
    G: "nx.Graph",
    seq_ids: List[str],
    partition: dict,
    output_prefix: str,
    figsize: Tuple[int, int] = (12, 10),
    dpi: int = 300,
) -> None:
    """Visualize network with community coloring.

    Args:
        G: NetworkX graph.
        seq_ids: Sequence ID list.
        partition: Community partition dict.
        output_prefix: Output file prefix (will append _network.pdf).
        figsize: Figure size in inches.
        dpi: Resolution.
    """
    if G.number_of_nodes() == 0:
        logging.warning("Empty network, skipping visualization.")
        return

    n_communities = len(set(partition.values()))
    cmap = cm.get_cmap("tab20", max(n_communities, 1))
    node_colors = [cmap(partition.get(n, 0) / max(n_communities, 1)) for n in G.nodes()]

    pos = nx.spring_layout(G, k=0.3, iterations=150, seed=42)

    plt.figure(figsize=figsize, dpi=dpi)
    nx.draw_networkx_nodes(
        G, pos, node_color=node_colors,
        node_size=500, alpha=0.9, edgecolors="black", linewidths=0.5,
    )
    nx.draw_networkx_edges(G, pos, edge_color="gray", alpha=0.3, width=0.8)

    # Label only first node of each community to avoid clutter
    community_first = {}
    for node, comm in partition.items():
        if comm not in community_first:
            community_first[comm] = node
    labels = {n: seq_ids[n] for n in community_first.values()}
    nx.draw_networkx_labels(G, pos, labels=labels, font_size=8, font_weight="bold")

    # Legend
    handles = [
        plt.Line2D([0], [0], marker="o", color="w",
                   label=f"Community {i+1}",
                   markerfacecolor=cmap(i / max(n_communities, 1)),
                   markersize=10)
        for i in range(n_communities)
    ]
    plt.legend(handles=handles, title="Communities", loc="best", fontsize=9)

    plt.title("Sequence Similarity Network", fontsize=14)
    plt.axis("off")
    plt.tight_layout()

    pdf_path = f"{output_prefix}_network.pdf"
    plt.savefig(pdf_path, format="pdf", dpi=dpi, bbox_inches="tight")
    plt.close()
    logging.info("Network plot saved to %s", pdf_path)


def build_parser() -> argparse.ArgumentParser:
    """Build argument parser."""
    parser = argparse.ArgumentParser(
        description="Build and visualize sequence similarity networks using k-mer Jaccard similarity.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""\
examples:
  python sequence_similarity_network.py -i proteins.fa -o result -t 0.5 -k 5
  python sequence_similarity_network.py -i proteins.fa -o result -t 0.3 -k 3
""",
    )
    parser.add_argument(
        "-i", "--input", type=str, required=True,
        help="Input FASTA file with sequences to compare.",
    )
    parser.add_argument(
        "-t", "--threshold", type=float, default=0.5,
        help="Similarity threshold for network edges (0-1, default: 0.5).",
    )
    parser.add_argument(
        "-k", "--kmer-size", type=int, default=5,
        help="K-mer size for similarity calculation (default: 5).",
    )
    parser.add_argument(
        "-o", "--output-prefix", type=str, required=True,
        help="Output prefix for PDF and SIF files.",
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
    import community.community_louvain as community_louvain
    import networkx as nx
    import numpy as np
    from Bio import SeqIO
    from scipy.sparse import csr_matrix
    from scipy.spatial.distance import pdist, squareform
    from sklearn.feature_extraction.text import CountVectorizer

    input_path = Path(args.input)
    if not input_path.is_file():
        logging.error("Input file not found: %s", args.input)
        sys.exit(1)

    # Parse sequences
    records = list(SeqIO.parse(args.input, "fasta"))
    if not records:
        logging.error("No sequences found in %s", args.input)
        sys.exit(1)

    sequences = [str(rec.seq) for rec in records]
    seq_ids = [rec.id for rec in records]
    logging.info("Loaded %d sequences from %s", len(sequences), args.input)

    # Compute similarity and build network
    sim_matrix = compute_similarity_matrix(sequences, k=args.kmer_size, threshold=args.threshold)
    G = build_network(sim_matrix, args.threshold)

    # Community detection
    partition = detect_communities(G)

    # Output
    visualize_network(G, seq_ids, partition, args.output_prefix)
    save_sif(G, seq_ids, f"{args.output_prefix}_network.sif")


if __name__ == "__main__":
    main()
