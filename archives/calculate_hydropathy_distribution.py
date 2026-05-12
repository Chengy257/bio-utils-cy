#!/usr/bin/
#########################################################################
# File Name: /home/chengyu/myscripts/calculate_hydropathy_distribution.py
# Author: ChengYu
# Description: 
# Created Time: Tue 22 Oct 2024 08:45:33 PM CST
#########################################################################
# file: calculate_hydropathy_distribution.py

from Bio import SeqIO
from Bio.SeqUtils.ProtParam import ProteinAnalysis
import csv
import argparse
import numpy as np
import matplotlib.pyplot as plt

def calculate_hydropathy_distribution(sequence, window_size=10):
    """Calculate hydropathy distribution using a sliding window."""
    hydropathy_scores = []

    # Define hydropathy scale (Kyte-Doolittle scale)
    hydropathy_scale = {
        'A': 1.8, 'C': 2.5, 'D': -3.5, 'E': -3.5,
        'F': 2.8, 'G': -0.4, 'H': -3.2, 'I': 1.8,
        'K': -3.9, 'L': 1.8, 'M': 1.9, 'N': -3.5,
        'P': -1.6, 'Q': -3.5, 'R': -4.5, 'S': -0.8,
        'T': -0.7, 'V': 4.2, 'W': -0.9, 'Y': -1.3
    }

    # Use a sliding window to calculate mean hydropathy for each window
    for i in range(len(sequence) - window_size + 1):
        window = sequence[i:i + window_size]
        score = np.mean([hydropathy_scale.get(amino, 0) for amino in window])
        hydropathy_scores.append(score)

    return hydropathy_scores

def process_fasta(fasta_file, output_csv, window_size):
    results = []

    # Read FASTA file
    try:
        sequences = [(record.id, str(record.seq)) for record in SeqIO.parse(fasta_file, "fasta")]
    except FileNotFoundError:
        print(f"Error: Unable to find file {fasta_file}. Please check the path.")
        return
    except Exception as e:
        print(f"An error occurred: {e}")
        return

    # Calculate hydropathy distribution for each sequence
    for record_id, sequence in sequences:
        scores = calculate_hydropathy_distribution(sequence, window_size)
        results.append({
            "id": record_id,
            "hydropathy_distribution": scores
        })
        # plot_hydropathy_distribution(record_id, scores)

    # Save results to CSV file
    try:
        with open(output_csv, mode='w', newline='') as file:
            writer = csv.writer(file)
            writer.writerow(["id", "hydropathy_distribution"])
            for result in results:
                writer.writerow([result["id"], ','.join(map(str, result["hydropathy_distribution"]))])
        print(f"Hydropathy distribution saved to {output_csv}")
    except Exception as e:
        print(f"Error writing to CSV file: {e}")

def plot_hydropathy_distribution(record_id, scores):
    """Plot the hydropathy distribution."""
    plt.figure(figsize=(10, 5))
    plt.plot(scores, marker='o', linestyle='-', label=f'Hydropathy Distribution: {record_id}')
    plt.title(f'Hydropathy Distribution for {record_id}')
    plt.xlabel('Position in Peptide (Window Size)')
    plt.ylabel('Hydropathy Score')
    plt.axhline(0, color='red', linestyle='--', label='Hydropathy Score = 0')
    plt.legend()
    plt.grid()
    plt.savefig(f'hydropathy_distribution_{record_id}.pdf')
    plt.close()
    print(f"Plot saved as hydropathy_distribution_{record_id}.png")

def main():
    parser = argparse.ArgumentParser(description='Calculate hydropathy distribution of peptides from a FASTA file.')
    parser.add_argument('fasta_file', type=str, help='Input FASTA file path')
    parser.add_argument('output_csv', type=str, help='Output CSV file path')
    parser.add_argument('--window_size', type=int, default=10, help='Window size for hydropathy calculation')

    args = parser.parse_args()
    
    process_fasta(args.fasta_file, args.output_csv, args.window_size)

if __name__ == "__main__":
    main()
