#!/usr/bin/env python3
"""
Cluster and visualize sequence similarities with community detection in a network graph.

Usage:
    python sequence_network.py -i input.fasta -t threshold -k kmer_size -o output_prefix

Options:
    -i, --input           Path to the input FASTA file.
    -t, --threshold       Similarity threshold (0-1) for edges in the network.
    -k, --kmer_size       Size of k-mers for approximate similarity (default=5).
    -o, --output_prefix   Prefix for output files (PDF and SIF).
    -h, --help            Show help message.
"""

import argparse
from Bio import SeqIO
import numpy as np
import networkx as nx
from scipy.spatial.distance import pdist, squareform
from sklearn.feature_extraction.text import CountVectorizer
from scipy.sparse import csr_matrix, coo_matrix
import matplotlib.pyplot as plt
import community.community_louvain as community_louvain  # 正确导入community_louvain
import matplotlib.cm as cm


def parse_fasta(fasta_file):
    """Parse sequences from a FASTA file."""
    return [str(record.seq) for record in SeqIO.parse(fasta_file, "fasta")]


def get_kmers(sequence, k):
    """Convert a sequence into a list of k-mers."""
    return [sequence[i:i + k] for i in range(len(sequence) - k + 1)]


def compute_sparse_similarity_matrix(sequences, k=5, threshold=0.5):
    """Compute a sparse similarity matrix based on K-mer Jaccard similarity using sklearn."""
    vectorizer = CountVectorizer(analyzer=lambda seq: get_kmers(seq, k))
    kmer_matrix = vectorizer.fit_transform(sequences)
    
    # Compute pairwise Jaccard similarity
    jaccard_similarities = 1 - pdist(kmer_matrix.toarray(), metric="jaccard")
    similarity_matrix = squareform(jaccard_similarities)

    # Threshold the similarity matrix to make it sparse
    similarity_matrix[similarity_matrix < threshold] = 0
    sparse_similarity_matrix = csr_matrix(similarity_matrix)
    
    return sparse_similarity_matrix


def build_network(similarity_matrix, threshold):
    """Build a network graph from a similarity matrix with edges above the threshold."""
    G = nx.Graph()
    coo_similarity_matrix = similarity_matrix.tocoo()

    # Add edges with weights above the threshold
    for i, j, weight in zip(coo_similarity_matrix.row, coo_similarity_matrix.col, coo_similarity_matrix.data):
        if weight >= threshold:
            G.add_edge(i, j, weight=weight)
    
    return G


# def plot_community_network(G, sequences, output_prefix):
#     """Detect communities in the graph and plot with different colors for each community."""
#     # Apply Louvain community detection
#     partition = community_louvain.best_partition(G, weight='weight')
    
#     # Generate colors for each community
#     num_communities = len(set(partition.values()))
#     cmap = cm.get_cmap('tab20', num_communities)
#     node_colors = [cmap(partition[node]) for node in G.nodes()]

#     # Draw network
#     pos = nx.spring_layout(G)
#     nx.draw(G, pos, node_color=node_colors, with_labels=True, labels={i: f"Seq{i+1}" for i in G.nodes()},
#             node_size=700, font_weight="bold", edge_color='gray', alpha=0.7)
#     plt.title("Sequence Similarity Network with Community Detection")
#     plt.savefig(f"{output_prefix}_network.pdf")  # Save to PDF
#     plt.show()

#     # Save SIF file
#     with open(f"{output_prefix}_network.sif", "w") as sif_file:
#         for i, j, w in G.edges(data="weight"):
#             sif_file.write(f"Seq{i+1} (similarity {w:.2f}) Seq{j+1}\n")

def plot_community_network(G, sequences, output_prefix):
    """Detect communities in the graph and plot with high-quality settings suitable for publication."""
    # Apply Louvain community detection
    partition = community_louvain.best_partition(G, weight='weight')
    
    # Generate a unique color map for each community
    num_communities = len(set(partition.values()))
    cmap = cm.get_cmap('tab20', num_communities)
    node_colors = [cmap(partition[node]) for node in G.nodes()]
    
    # Use a high spread for the layout to minimize overlap
    pos = nx.spring_layout(G, k=0.3, iterations=150, seed=42)  # Higher k and iterations for clearer layout

    # Determine node sizes based on community membership, larger for central nodes
    sizes = [800 if partition[node] == partition[list(G.nodes())[0]] else 500 for node in G.nodes()]

    # Show only the central node label for clarity
    center_node = list(G.nodes())[0]  # Assuming the first node as the center
    labels = {center_node: f"Seq{center_node+1}"}

    # Plot the network
    plt.figure(figsize=(12, 10), dpi=300)  # High resolution for publication
    nx.draw_networkx_nodes(G, pos, node_color=node_colors, node_size=sizes, edgecolors='black', linewidths=0.5, alpha=0.9)
    nx.draw_networkx_edges(G, pos, edge_color='gray', alpha=0.3, width=0.8)
    nx.draw_networkx_labels(G, pos, labels=labels, font_size=12, font_weight="bold")

    # Add a legend for community colors
    handles = [plt.Line2D([0], [0], marker='o', color='w', label=f'Community {i+1}',
                          markerfacecolor=cmap(i/num_communities), markersize=10) for i in range(num_communities)]
    plt.legend(handles=handles, title="Communities", loc="best", fontsize=10, title_fontsize=12, frameon=True)

    # Title and layout adjustments for clarity
    plt.title("Sequence Similarity Network with Community Detection", fontsize=16)
    plt.axis("off")  # Hide axes for a cleaner look
    plt.tight_layout()
    plt.savefig(f"{output_prefix}_network_high_res.pdf", format="pdf", dpi=300, bbox_inches="tight")
    #plt.show()

    # Save SIF file
    with open(f"{output_prefix}_network.sif", "w") as sif_file:
        for i, j, w in G.edges(data="weight"):
            sif_file.write(f"Seq{i+1} (similarity {w:.2f}) Seq{j+1}\n")


def main():
    parser = argparse.ArgumentParser(description="Cluster and visualize sequence similarities with community detection in a network graph.")
    parser.add_argument("-i", "--input", required=True, help="Path to the input FASTA file.")
    parser.add_argument("-t", "--threshold", type=float, default=0.7, help="Similarity threshold for network edges.")
    parser.add_argument("-k", "--kmer_size", type=int, default=5, help="Size of k-mers for approximate similarity.")
    parser.add_argument("-o", "--output_prefix", required=True, help="Prefix for output files (PDF and SIF).")
    args = parser.parse_args()

    # Step 1: Parse sequences
    sequences = parse_fasta(args.input)

    # Step 2: Compute sparse similarity matrix
    similarity_matrix = compute_sparse_similarity_matrix(sequences, args.kmer_size, args.threshold)

    # Step 3: Build network from similarity matrix
    G = build_network(similarity_matrix, args.threshold)

    # Step 4: Detect communities and plot network with colored clusters
    plot_community_network(G, sequences, args.output_prefix)


if __name__ == "__main__":
    main()
