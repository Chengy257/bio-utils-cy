#!/usr/bin/
#########################################################################
# File Name: /home/chengyu/myscripts/calculate_peptide_properties.py
# Author: ChengYu
# Description: calculate_peptide_properties
# Created Time: Tue Oct 22 14:09:44 CST 2024
#########################################################################
# file: calculate_peptide_properties.py

from Bio import SeqIO
from Bio.SeqUtils.ProtParam import ProteinAnalysis
import csv
import argparse
from concurrent.futures import ThreadPoolExecutor, as_completed
from tqdm import tqdm
import matplotlib.pyplot as plt

def analyze_sequence(sequence, record_id):
    """Analyze a single peptide sequence and return its properties."""
    if not sequence or not all(c in 'ACDEFGHIKLMNPQRSTVWY' for c in sequence):
        print(f"Invalid sequence: {record_id}, skipping.")
        return None

    analysis = ProteinAnalysis(sequence)
    
    return {
        "id": record_id,
        "molecular_weight": analysis.molecular_weight(),
        "net_charge": analysis.charge_at_pH(7.0),
        "isoelectric_point": analysis.isoelectric_point(),
        "average_hydropathy": analysis.gravy(),
        "aromaticity": analysis.aromaticity(),
        # "flexibility": analysis.flexibility(),
        "instability_index": analysis.instability_index(),
        # "molar_extinction_coefficient": analysis.molar_extinction_coefficient(),
        "secondary_structure_fraction": analysis.secondary_structure_fraction(),
        "composition": analysis.get_amino_acids_percent(),
    }

def calculate_properties(fasta_file, output_csv, visualize):
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

    # Analyze sequences using multithreading and display a progress bar
    with ThreadPoolExecutor() as executor:
        futures = {executor.submit(analyze_sequence, seq, record_id): record_id for record_id, seq in sequences}
        for future in tqdm(as_completed(futures), total=len(futures), desc="Processing Progress"):
            result = future.result()
            if result:
                results.append(result)

    # Save results to CSV file
    try:
        with open(output_csv, mode='w', newline='') as file:
            writer = csv.DictWriter(file, fieldnames=results[0].keys(),delimiter='\t')
            writer.writeheader()
            writer.writerows(results)
        print(f"Calculation results saved to {output_csv}")
    except Exception as e:
        print(f"Error writing to CSV file: {e}")

    # Visualize properties distribution (if user chooses)
    if visualize:
        visualize_properties(results)



def visualize_properties(results):
    """Plot the distribution of various properties."""
    molecular_weights = [result["molecular_weight"]/1000 for result in results]
    net_charges = [result["net_charge"] for result in results]
    isoelectric_points = [result["isoelectric_point"] for result in results]
    instability_indexs = [result["instability_index"] for result in results]
    average_hydropathys = [result["average_hydropathy"] for result in results]
    aromaticitys = [result["aromaticity"] for result in results]

    # Amino acid composition
    amino_acid_composition = {}
    for result in results:
        composition = result.get("composition", {})
        for aa, percentage in composition.items():
            amino_acid_composition[aa] = amino_acid_composition.get(aa, 0) + percentage

    # Secondary structure fractions
    # Secondary structure fractions
    helix_fractions = []
    turn_fractions = []
    sheet_fractions = []

    for result in results:
        if isinstance(result.get("secondary_structure_fraction"), tuple):
            helix, turn, sheet = result["secondary_structure_fraction"]
            helix_fractions.append(helix)
            turn_fractions.append(turn)
            sheet_fractions.append(sheet)

    plt.figure(figsize=(12, 9))

    # Molecular weight distribution
    plt.subplot(3, 3, 1)
    plt.hist(molecular_weights, bins=30, color='skyblue')
    plt.title('Molecular Weight Distribution')
    plt.xlabel('Molecular Weight (kDa)')
    plt.ylabel('Frequency')
    plt.grid(axis='y', alpha=0.75)

    # Net charge distribution
    plt.subplot(3, 3, 2)
    plt.hist(net_charges, bins=30, color='lightgreen')
    plt.title('Net Charge Distribution')
    plt.xlabel('Net Charge')
    plt.ylabel('Frequency')
    plt.grid(axis='y', alpha=0.75)

    # Isoelectric point distribution
    plt.subplot(3, 3, 3)
    plt.hist(isoelectric_points, bins=30, color='salmon')
    plt.title('Isoelectric Point Distribution')
    plt.xlabel('Isoelectric Point')
    plt.ylabel('Frequency')
    plt.grid(axis='y', alpha=0.75)

    # Instability index distribution
    plt.subplot(3, 3, 4)
    plt.hist(instability_indexs, bins=30, color='grey')
    plt.title('Instability Indexs Distribution')
    plt.xlabel('Instability Indexs')
    plt.ylabel('Frequency')
    plt.grid(axis='y', alpha=0.75)

    # Plot amino acid composition
    plt.subplot(3, 3, 5)
    plt.bar(amino_acid_composition.keys(), amino_acid_composition.values(), color='lightblue')
    plt.title('Amino Acid Composition')
    plt.xlabel('Amino Acids')
    plt.ylabel('Percentage (%)')
    # plt.xticks(rotation=45)
    plt.grid(axis='y', alpha=0.75)

    # Plot secondary structure fractions
    plt.subplot(3, 3, 6)
    plt.bar(['Helix', 'Turn', 'Sheet'], 
            [sum(helix_fractions), sum(turn_fractions), sum(sheet_fractions)],
            color=['gold', 'lightgreen', 'lightcoral'])
    plt.title('Secondary Structure Fractions')
    plt.xlabel('Structure Type')
    plt.ylabel('Fraction (%)')
    plt.grid(axis='y', alpha=0.75)

    # Plot average hydropathys
    plt.subplot(3, 3, 7)
    plt.hist(average_hydropathys, bins=30, color='gold')
    plt.title('Average Hydropathys Distribution')
    plt.xlabel('Average Hydropathys')
    plt.ylabel('Frequency')
    plt.grid(axis='y', alpha=0.75)

    # Plot aromaticitys
    plt.subplot(3, 3, 8)
    plt.hist(aromaticitys, bins=30, color='lightcoral')
    plt.title('Aromaticitys Distribution')
    plt.xlabel('Aromaticity')
    plt.ylabel('Frequency')
    plt.grid(axis='y', alpha=0.75)

    plt.tight_layout()
    plt.savefig('properties_distribution.pdf')
    # plt.show()

def main():
    parser = argparse.ArgumentParser(description='Calculate various properties of peptides and save to a TSV file.')
    parser.add_argument('fasta_file', type=str, help='Input FASTA file path')
    parser.add_argument('output_csv', type=str, help='Output TSV file path')
    parser.add_argument('--visualize', action='store_true', help='Generate properties distribution plots')
    args = parser.parse_args()
    calculate_properties(args.fasta_file, args.output_csv, args.visualize)

# Call function example
if __name__ == "__main__":
    main()
