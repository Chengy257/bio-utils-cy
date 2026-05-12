import os
import argparse
import networkx as nx
import matplotlib.pyplot as plt
from matplotlib import cm
from Bio import SeqIO
from subprocess import run
import numpy as np
import markov_clustering as mc

def read_fasta_sequences(fasta_file):
    """读取FASTA文件中的蛋白质序列"""
    sequences = list(SeqIO.parse(fasta_file, "fasta"))
    return sequences

def run_blast(fasta_file, blast_db="protein_db", output_file="blast_output.txt", evalue=1e-5, num_threads=1):
    """运行BLAST，创建相似性输出文件，使用subprocess模块"""
    if not os.path.exists(f"{blast_db}.pin"):  # 检查是否已有BLAST数据库
        makeblastdb_cmd = f"makeblastdb -in {fasta_file} -dbtype prot -out {blast_db}"
        run(makeblastdb_cmd, shell=True, check=True)
    blastp_cmd = f"blastp -query {fasta_file} -db {blast_db} -evalue {evalue} -outfmt 6 -out {output_file} -num_threads {num_threads}"
    run(blastp_cmd, shell=True, check=True)
    return output_file

def parse_blast_results(output_file, sequences, threshold=0.5):
    """将BLAST结果解析为相似性矩阵"""
    seq_ids = [seq.id for seq in sequences]
    matrix = np.zeros((len(seq_ids), len(seq_ids)))

    with open(output_file) as f:
        for line in f:
            qseqid, sseqid, pident, *_ = line.split()
            i, j = seq_ids.index(qseqid), seq_ids.index(sseqid)
            if float(pident) >= threshold * 100:
                matrix[i, j] = matrix[j, i] = float(pident)
    return matrix, seq_ids

def run_mcl_clustering(matrix, inflation=2.0):
    """使用markov_clustering包的MCL算法进行聚类"""
    result = mc.run_mcl(matrix, inflation=inflation)
    clusters = mc.get_clusters(result)
    return clusters

def save_clusters_to_file(clusters, seq_ids, output_file="clusters.txt"):
    """将聚类结果保存到文件中"""
    with open(output_file, 'w') as f:
        for i, cluster in enumerate(clusters, start=1):
            f.write(f"Cluster {i}: {', '.join([seq_ids[idx] for idx in cluster])}\n")

def generate_network_graph(clusters, seq_ids):
    """根据聚类结果生成网络图"""
    G = nx.Graph()
    for i, cluster in enumerate(clusters, start=1):
        for node in cluster:
            G.add_node(seq_ids[node], cluster=i)  # 为每个节点标记cluster ID
        for j in range(len(cluster)):
            for k in range(j + 1, len(cluster)):
                G.add_edge(seq_ids[cluster[j]], seq_ids[cluster[k]])
    return G

def save_sif_format(G, output_sif="network.sif"):
    """保存网络图为SIF格式文件"""
    with open(output_sif, 'w') as f:
        for edge in G.edges():
            f.write(f"{edge[0]} pp {edge[1]}\n")


def visualize_network(G, clusters):
    """优化的网络图可视化，适合学术展示，无图例"""
    plt.figure(figsize=(3, 3), dpi=600)  # 提高清晰度

    # 获取所有簇并设置颜色
    unique_clusters = set(nx.get_node_attributes(G, 'cluster').values())
    num_clusters = len(unique_clusters)
    
    # 使用colormap生成颜色，以确保色调适中且区分度高
    color_map = cm.get_cmap('tab20', num_clusters)
    cluster_colors = {cluster_id: color_map(i / num_clusters) for i, cluster_id in enumerate(unique_clusters)}
    
    # 设置节点位置和样式
    pos = nx.spring_layout(G, seed=42, k=0.15, iterations=50)  # 调整布局参数，确保节点均匀分布
    for cluster_id in unique_clusters:
        nodes_in_cluster = [node for node, attr in G.nodes(data=True) if attr['cluster'] == cluster_id]
        nx.draw_networkx_nodes(
            G, pos, nodelist=nodes_in_cluster,
            node_color=[cluster_colors[cluster_id]], 
            node_size=10, alpha=0.5,  # 更小的节点
            edgecolors='white', linewidths=0  # 边框更细
        )

    # 绘制边线并设置透明度和宽度
    nx.draw_networkx_edges(G, pos, edge_color="gray", alpha=0.2, width=0.2)
    
    # 移除图例，优化布局
    plt.title("Protein Sequence Clustering Network", fontsize=12)
    plt.axis('off')
    plt.subplots_adjust(left=0, right=1, top=1, bottom=0)  # 移除边距确保紧凑显示
    #plt.show()
    plt.savefig(f"Protein_Sequence_Clustering_Network.pdf", format="pdf", dpi=300, bbox_inches="tight")

def main():
    # 配置命令行参数
    parser = argparse.ArgumentParser(description="蛋白质序列聚类和网络图生成工具")
    parser.add_argument("fasta_file", type=str, help="输入蛋白质FASTA序列文件路径")
    parser.add_argument("--method", type=str, choices=["blast"], default="blast", help="相似性计算方法，目前支持blast")
    parser.add_argument("--output_sif", type=str, default="network.sif", help="输出SIF文件路径")
    parser.add_argument("--output_clusters", type=str, default="clusters.txt", help="保存聚类结果的文件路径")
    parser.add_argument("--evalue", type=float, default=1e-5, help="BLAST相似性e值阈值")
    parser.add_argument("--threshold", type=float, default=0.5, help="BLAST结果解析时的相似性百分比阈值")
    parser.add_argument("--num_threads", type=int, default=1, help="BLAST计算时使用的线程数")
    parser.add_argument("--inflation", type=float, default=2.0, help="MCL聚类的膨胀参数")
    args = parser.parse_args()

    # 读取FASTA序列
    sequences = read_fasta_sequences(args.fasta_file)

    # 运行BLAST并解析相似性矩阵
    blast_output = run_blast(args.fasta_file, evalue=args.evalue, num_threads=args.num_threads)
    matrix, seq_ids = parse_blast_results(blast_output, sequences, threshold=args.threshold)

    # 使用MCL进行聚类
    clusters = run_mcl_clustering(matrix, inflation=args.inflation)

    # 保存聚类结果到文件
    save_clusters_to_file(clusters, seq_ids, output_file=args.output_clusters)

    # 创建网络图
    G = generate_network_graph(clusters, seq_ids)
    
    # 可视化网络图
    visualize_network(G, clusters)

    save_sif_format(G, output_sif=args.output_sif)

if __name__ == "__main__":
    main()
