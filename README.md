# myscripts

Personal bioinformatics toolkit for next-generation sequencing (NGS) data analysis.

This repository contains standalone utility scripts covering RNA-seq, ChIP-seq, CUT&Tag, ATAC-seq, sRNA-seq, ribosome profiling, lncRNA analysis, proteomics, and general data processing.

- **`myscripts/`** — Rewritten, production-quality scripts (Python / R / Bash), optimized with robust error handling, full CLI help, and consistent code style.
- **`myscripts_raw/`** — Original scripts preserved as archive.

> The rewritten scripts in `myscripts/` were refactored and improved with the assistance of AI (Claude), adding proper argument parsing, logging, input validation, and comprehensive documentation while preserving the original algorithmic logic.

---

## Directory Structure

```
myscripts/
├── sequence_analysis/      # Sequence alignment, similarity networks, reverse complement
├── codon_analysis/         # Codon/AA frequency, CAI, dN/dS, Kozak score
├── format_conversion/      # GFF/GTF/BED/BAM format converters
├── data_processing/        # File joining, transposition, general utilities
├── gene_annotation/        # UTR extraction, intron filtering, MAF region extraction
├── rnaseq/                 # DEG analysis, expression quantification, count merging
├── data_retrieval/         # SRA/ENA/PRIDE downloading and metadata fetching
├── enrichment/             # GO/KEGG enrichment analysis
├── proteomics/             # MS/MS analysis, peptide properties, structure visualization
├── visualization/          # Heatmaps, volcano, Manhattan, UpSet, gene structure plots
└── ml_stats/               # Machine learning, statistics, clustering
```

## Quick Start

```bash
# Python scripts — use argparse
python myscripts/sequence_analysis/blast_align_analysis.py -h

# R scripts — use getopt
Rscript myscripts/rnaseq/deseq2_multigroup.R -h

# Bash scripts — use getopts
bash myscripts/format_conversion/bam_to_bigwig.sh -h
```

---

## Scripts by Category

### Sequence Analysis

> Sequence alignment, similarity networks, and sequence manipulation.

| Script | Language | Description |
|--------|----------|-------------|
| `blast_align_analysis.py` | Python | BLAST-based inter-species similarity metrics |
| `aligned_fasta_similarity.py` | Python | Direct pairwise similarity from aligned FASTA |
| `sequence_similarity_network.py` | Python | K-mer Jaccard similarity + Louvain clustering network |
| `blast_sequence_network.py` | Python | BLAST + MCL clustering network |
| `reverse_complement.py` | Python | Reverse complement / complement / reverse of FASTA or string |
| `sequence_complexity.py` | Python | Sequence complexity (block entropy, conditional entropy) |

<details>
<summary>Usage examples</summary>

```bash
# BLAST similarity between species
python blast_align_analysis.py -i alignments/ -r human -o results.csv -t 8

# Direct alignment similarity
python aligned_fasta_similarity.py -i aligned/ -r human -o similarity.csv

# K-mer similarity network
python sequence_similarity_network.py -i proteins.fa -o result -k 5 -t 0.5

# BLAST + MCL clustering
python blast_sequence_network.py -i proteins.fa -o result --evalue 1e-10 --inflation 2.5

# Reverse complement
python reverse_complement.py -i genes.fa -m revcomp -o output.fa
python reverse_complement.py -s ATCGATCG -m comp

# Sequence complexity
python sequence_complexity.py -i sequences.fa -o complexity.tsv -w 4 -b 2
```
</details>

---

### Codon Analysis

> Codon usage, amino acid composition, CAI, Kozak score, dN/dS, and nucleotide diversity.

| Script | Language | Description |
|--------|----------|-------------|
| `codon_frequency.py` | Python | Codon usage frequency table from CDS FASTA |
| `amino_acid_frequency.py` | Python | Amino acid composition from protein FASTA |
| `calculate_cai.py` | Python | Codon Adaptation Index (CAI) calculation |
| `kozak_similarity_score.py` | Python | Kozak consensus sequence similarity score |
| `calculate_dnds.py` | Python | Pairwise dN/dS ratios from aligned DNA FASTA |
| `calculate_nucleotide_diversity.py` | Python | Nucleotide diversity (Pi) via PLINK + VCFTools |

<details>
<summary>Usage examples</summary>

```bash
# Codon frequency
python codon_frequency.py -i cds.fa -o codon_freq.tsv

# Amino acid frequency
python amino_acid_frequency.py -i proteins.fa -o aa_freq.tsv

# CAI calculation
python calculate_cai.py -i target_cds.fa -r ribosomal_genes.fa -o cai_results.tsv

# Kozak score
python kozak_similarity_score.py -i tis_sequences.fa -o scores.tsv
python kozak_similarity_score.py -s GCCACCATGGCGATCGATCGATC

# dN/dS
python calculate_dnds.py -i gene_alignments/ -o dnds_results.tsv

# Nucleotide diversity
python calculate_nucleotide_diversity.py -b regions.bed -p sample_data -o pi.tsv -t 8
```
</details>

---

### Format Conversion

> Genome annotation format converters and FASTA utilities.

| Script | Language | Description |
|--------|----------|-------------|
| `genome_format_converter.sh` | Bash | UCSC tools wrapper: GFF↔GTF↔BED↔genePred |
| `gtf_standardize.sh` | Bash | GTF attribute standardization |
| `bam_to_bigwig.sh` | Bash | BAM → bigWig with multiple normalization methods |
| `bed12_effective_length.py` | Python | Non-redundant genomic length from BED12 |
| `fasta_to_alphafold_json.py` | Python | FASTA → AlphaFold server JSON format |

<details>
<summary>Usage examples</summary>

```bash
# GFF to GTF
bash genome_format_converter.sh -i annotation.gff3 -m gff2gtf -o annotation.gtf

# GTF to BED
bash genome_format_converter.sh -i annotation.gtf -m gtf2bed -o annotation.bed

# GTF standardization
bash gtf_standardize.sh -i input.gtf -o output.gtf

# BAM to bigWig (RPKM/CPM/BPM/RPGC/None)
bash bam_to_bigwig.sh -i bam_list.txt -n RPKM -t 8

# BED12 effective length
python bed12_effective_length.py -i regions.bed -o merged.bed

# FASTA to AlphaFold JSON
python fasta_to_alphafold_json.py -i proteins.fa -o alphafold_input.json
```
</details>

---

### Data Processing

> General-purpose file manipulation and data transformation utilities.

| Script | Language | Description |
|--------|----------|-------------|
| `multi_file_join.py` | Python | Join multiple TSV files by a common key column |
| `transpose_table.sh` | Bash | Transpose a tab-delimited table (swap rows/columns) |
| `number_duplicates.R` | R | Append sequential numbers to duplicate values |
| `tissue_specificity_tau.py` | Python | Tissue specificity index (tau) from expression matrix |

<details>
<summary>Usage examples</summary>

```bash
# Multi-file join
python multi_file_join.py -i file_list.txt -c 1,3 -o merged.tsv --how inner

# Transpose table
bash transpose_table.sh -i matrix.tsv -o transposed.tsv

# Number duplicates
Rscript number_duplicates.R input.tsv output.tsv

# Tissue specificity
python tissue_specificity_tau.py -i expression.tsv -o tau.tsv
```
</details>

---

### Gene Annotation

> Genomic coordinate handling, UTR extraction, and MAF region extraction.

| Script | Language | Description |
|--------|----------|-------------|
| `extract_utr.py` | Python | Extract 5'/3' UTR sequences from GTF + genome FASTA |
| `genepred_utr_to_bed12.py` | Python | UTR regions from genePred → BED12 |
| `filter_long_introns.py` | Python | Filter transcripts with abnormally long introns |
| `maf_extract_regions.py` | Python | Extract genomic regions from MAF using BED12 coordinates |

<details>
<summary>Usage examples</summary>

```bash
# Extract UTR
python extract_utr.py --gtf annotation.gtf --genome genome.fa -f utr.fa -c utr_lengths.csv

# genePred UTR to BED12
python genepred_utr_to_bed12.py -i annotation.gp -o utr.bed12 --utr 5UTR

# Filter long introns
python filter_long_introns.py -i annotation.gtf -o filtered.gtf -m 20000

# MAF region extraction
python maf_extract_regions.py --bed cds.bed12 --maf alignment.maf --ref hg38 --cds output.fa
```
</details>

---

### RNA-seq

> Differential expression, expression quantification, and count processing.

| Script | Language | Description |
|--------|----------|-------------|
| `deseq2_multigroup.R` | R | Full DESeq2 multi-group differential expression pipeline |
| `deg_wilcox_test.R` | R | Wilcoxon rank-sum DEG test with TMM normalization |
| `featurecounts_pipeline.sh` | Bash | Automated featureCounts pipeline with strandedness detection |
| `merge_featurecounts.py` | Python | Merge featureCounts results from multiple samples |
| `calculate_fpkm_tpm.R` | R | FPKM/TPM from STAR gene-level counts |
| `calculate_effective_length.R` | R | Non-redundant exon length per gene from GTF |
| `calculate_tpm.R` | R | Expression unit conversion functions (library) |
| `filter_expression.py` | Python | Filter low-expression and low-variability genes |

<details>
<summary>Usage examples</summary>

```bash
# DESeq2 multi-group
Rscript deseq2_multigroup.R -c counts.tsv -s samples.csv -o result -p 0.05 -f 2

# Wilcoxon DEG
Rscript deg_wilcox_test.R -c counts.tsv -t conditions.tsv -o result

# featureCounts pipeline
bash featurecounts_pipeline.sh -b bam_dir/ -g annotation.gtf -o output/ -r regions.bed -t 8

# Merge featureCounts
python merge_featurecounts.py -i featurecounts_output/ -o merged/

# FPKM/TPM from STAR
Rscript calculate_fpkm_tpm.R -g annotation.gtf -c ReadsPerGene.out.tab -s 0

# Effective length
Rscript calculate_effective_length.R -g annotation.gtf -o efflen.tsv

# Filter expression matrix
python filter_expression.py -i expr.tsv -o filtered.tsv --min_expr 5 --min_samples 3
```
</details>

---

### Data Retrieval

> Download and fetch metadata from NCBI SRA, ENA, and PRIDE databases.

| Script | Language | Description |
|--------|----------|-------------|
| `sra_prefetch_batch.sh` | Bash | Batch download FASTQ from NCBI SRA |
| `ena_ascp_download.py` | Python | Download FASTQ from ENA using Aspera ascp |
| `ena_ascp_download_batch.sh` | Bash | Batch ENA Aspera download by accession |
| `fetch_sra_metadata_ncbi.py` | Python | SRA metadata from NCBI Entrez E-Utilities |
| `fetch_sra_metadata_ena.py` | Python | SRA metadata from ENA filereport API |
| `fetch_sra_metadata_comprehensive.py` | Python | Comprehensive SRA metadata via NCBI |
| `fetch_sra_metadata_xml.py` | Python | SRA metadata combining RunInfo CSV + XML |
| `format_sra_summary.py` | Python | Format SRA metadata into standardized summary |
| `fetch_bioproject_metadata.sh` | Bash | BioProject metadata from ENA |
| `fetch_pride_metadata.py` | Python | PRIDE proteomics dataset metadata |
| `parse_pride_json.py` | Python | Parse PRIDE JSON metadata into tables |
| `search_sra_riboseq.py` | Python | Search NCBI for Ribo-seq data by species |
| `search_riboseq_bioproject.py` | Python | Multi-layer Ribo-seq BioProject retrieval |
| `extract_species_rna_rfam.sh` | Bash | Extract species RNA from Rfam database |

<details>
<summary>Usage examples</summary>

```bash
# SRA batch download
bash sra_prefetch_batch.sh -i srr_list.txt -d ./fastq -t 4

# ENA Aspera download
python ena_ascp_download.py -i SRR1234567 -o ./fastq -t 8
bash ena_ascp_download_batch.sh -i accession_list.txt -o ./download

# NCBI SRA metadata
python fetch_sra_metadata_ncbi.py -i srr_list.txt -o metadata.tsv -e user@email.com

# ENA SRA metadata
python fetch_sra_metadata_ena.py -i accessions.txt -o metadata.tsv

# Format SRA summary
python format_sra_summary.py --input metadata.tsv --output summary.tsv

# BioProject metadata
bash fetch_bioproject_metadata.sh -i bioproject_list.txt -o metadata/

# PRIDE metadata
python fetch_pride_metadata.py -i pxd_list.txt -o pride_metadata.tsv

# Parse PRIDE JSON
python parse_pride_json.py ./pride_jsons/ -o results.csv

# Ribo-seq SRA search
python search_sra_riboseq.py --species "Oryza sativa" --output-dir results/

# Ribo-seq BioProject retrieval
python search_riboseq_bioproject.py -s "Arabidopsis thaliana" -o results.csv

# Rfam species RNA extraction
bash extract_species_rna_rfam.sh -i rfam_dir/ -a family.txt -o output/ "Oryza sativa"
```
</details>

---

### Enrichment

> Functional enrichment analysis.

| Script | Language | Description |
|--------|----------|-------------|
| `kegg_enrichment.R` | R | KEGG pathway enrichment using clusterProfiler |

<details>
<summary>Usage examples</summary>

```bash
# KEGG enrichment
Rscript kegg_enrichment.R -i gene_list.txt -o kegg_results/ -g dosa -p 0.05 -q 0.05
```
</details>

---

### Proteomics

> MS/MS analysis, peptide properties, protein structure visualization.

| Script | Language | Description |
|--------|----------|-------------|
| `plot_ms2_spectrum.py` | Python | Annotated MS2 spectra with b/y ion labels |
| `batch_dssp.sh` | Bash | Batch DSSP secondary structure prediction |
| `dssp_summary.py` | Python | DSSP output statistics summary |
| `batch_pymol_render.py` | Python | Batch PyMOL rendering with pLDDT coloring |
| `batch_chimerax_render.py` | Python | Batch ChimeraX structure rendering |
| `plot_structure_3d.py` | Python | 3D structure with hydropathy/secondary structure coloring |
| `peptide_properties.py` | Python | Peptide physicochemical properties |
| `hydropathy_distribution.py` | Python | Sliding window hydropathy analysis |
| `protein_molecular_weight.py` | Python | Protein molecular weight calculation |

<details>
<summary>Usage examples</summary>

```bash
# MS2 spectrum plot
python plot_ms2_spectrum.py --mzid result.mzid --spectra data.mzML --protein P12345 -o spec

# Batch DSSP
bash batch_dssp.sh -i pdb_dir/ -o dssp_output/ -t 4

# DSSP summary
python dssp_summary.py -i dssp_output/ -o summary.csv

# PyMOL batch render
python batch_pymol_render.py -i pdb_files/ -o images/ --width 1200 --dpi 300

# ChimeraX batch render
python batch_chimerax_render.py -i pdb_files/ -o images/ --color bfactor

# 3D structure plot
python plot_structure_3d.py -i pdb_files/ -o images/ --color hydro

# Peptide properties
python peptide_properties.py -i peptides.fasta -o properties.tsv --visualize

# Hydropathy distribution
python hydropathy_distribution.py -i proteins.fasta -o hydropathy.csv --window-size 15

# Protein molecular weight
python protein_molecular_weight.py -i proteins.fasta -o weights.tsv
```
</details>

---

### Visualization

> Publication-quality plots for genomics and proteomics data.

| Script | Language | Description |
|--------|----------|-------------|
| `plot_heatmap_multi.R` | R | Multi-panel ComplexHeatmap with z-score scaling |
| `plot_heatmap_single.R` | R | Single expression matrix heatmap |
| `plot_base_content.R` | R | Sliding window base composition analysis |
| `plot_flower.R` | R | Flower/petal plot for set intersection sizes |
| `plot_gene_structure.R` | R | Gene structure visualization from GTF |
| `plot_upset.R` | R | UpSet intersection plots |
| `plot_gwas_manhattan.R` | R | GWAS Manhattan + QQ + regional Manhattan plots |
| `plot_qpcr_analysis.R` | R | qPCR delta-delta CT analysis and visualization |
| `plot_go_tree_heatmap.R` | R | GO semantic similarity tree + enrichment barplot |
| `plot_go_heatmap_batch.sh` | Bash | Batch GO term heatmaps from expression data |

<details>
<summary>Usage examples</summary>

```bash
# Multi-panel heatmap
Rscript plot_heatmap_multi.R -f file_list.txt -g groups.txt -o Heatmap -c blue_white_red

# Single heatmap
Rscript plot_heatmap_single.R -d expr_matrix.tsv -g groups.txt -o heatmap

# Base content
Rscript plot_base_content.R -f genome.fa -b GC -w 100 --cutoff 50

# Flower plot
Rscript plot_flower.R -i samples.tsv -o flower.pdf

# Gene structure
Rscript plot_gene_structure.R -g annotation.gtf -o gene_structure.pdf --genes GeneA,GeneB

# UpSet plot
Rscript plot_upset.R -i gene_lists.txt -o upset

# GWAS Manhattan
Rscript plot_gwas_manhattan.R -e emmax.ps -m map.txt -b genes.bed -o gwas_result

# qPCR analysis
Rscript plot_qpcr_analysis.R -d data.csv -r ACTIN -c WT -o qpcr_result

# GO tree + heatmap
Rscript plot_go_tree_heatmap.R -i go_results.tsv -o go_tree.pdf -g org.At.tair.db

# Batch GO heatmaps
bash plot_go_heatmap_batch.sh go_ids.txt expr_matrix.tsv sample_info.csv -g go_terms.gmt
```
</details>

---

### ML & Statistics

> Machine learning, clustering, and statistical analysis.

| Script | Language | Description |
|--------|----------|-------------|
| `mfuzz_soft_cluster.R` | R | Mfuzz time-series soft clustering |
| `batch_parse_sanger.R` | R | Batch parse Sanger .ab1 sequencing traces |
| `batch_sanger_blast.sh` | Bash | Batch BLAST for Sanger sequencing results |

<details>
<summary>Usage examples</summary>

```bash
# Mfuzz soft clustering
Rscript mfuzz_soft_cluster.R -i expr.tsv -k 8 -o mfuzz_result

# Batch Sanger parsing
Rscript batch_parse_sanger.R -r reference.txt -i traces/ -o results/

# Sanger BLAST
bash batch_sanger_blast.sh -d database.fa -w work_dir/ -t 8 -o sanger_BLAST
```
</details>

---

## Dependencies

### Bioinformatics Tools
| Tool | Purpose |
|------|---------|
| STAR / HISAT2 / bowtie2 | Read alignment |
| featureCounts / StringTie | Quantification |
| MACS2 | Peak calling |
| FastQC / MultiQC / Trim Galore | Quality control |
| BLAST+ | Sequence similarity |
| SAMtools / deepTools | BAM processing |
| DSSP | Secondary structure |
| PyMOL / ChimeraX | Structure visualization |
| PLINK / VCFTools | Variant analysis |
| Aspera ascp | High-speed data transfer |

### Python Packages
Biopython, pandas, numpy, matplotlib, networkx, tqdm, pyteomics

### R / Bioconductor Packages
DESeq2, clusterProfiler, ComplexHeatmap, ggplot2, rtracklayer, Mfuzz, sangerseqR, UpSetR, qqman, GOSemSim, ggtree

---

## License

This project is licensed under the [GNU General Public License v3.0](LICENSE).

---

## Author

ChengYu
