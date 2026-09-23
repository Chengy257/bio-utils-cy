# myscripts

Personal bioinformatics toolkit for next-generation sequencing (NGS) data analysis.

This repository contains standalone utility scripts covering RNA-seq, ChIP-seq, CUT&Tag, ATAC-seq, sRNA-seq, ribosome profiling, lncRNA analysis, proteomics, and general data processing.

- **`bin/`** — All executable scripts in one flat directory (Python / R / Bash), production quality with robust error handling, full CLI help, and consistent code style. Functional categories are documented in the index below, not in the file tree.
- **`config/`** — Unified project-wide configuration shared by all scripts (see [Configuration](#configuration)).
- **`tools/`** — Maintenance utilities: `doctor.sh` (dependency checks & path locator) and `smoke_test.sh` (CLI regression).
- **`docs/`** — Reference material: per-script dependency matrix (`DEPENDENCIES.md`).
- **`archives/`** — Original pre-rewrite scripts preserved as archive.

> The scripts in `bin/` were refactored and improved with the assistance of AI (Claude), adding proper argument parsing, logging, input validation, and comprehensive documentation while preserving the original algorithmic logic.

---

## Directory Structure

```
myscripts/
├── bin/                    # All executable scripts, flat — the only place to look
├── config/
│   ├── env.sh               # Unified environment / dependency configuration (committed defaults)
│   ├── env.local.sh         # Machine-specific overrides (optional, gitignored — edit this one)
│   └── env.local.sh.example # Template to copy when setting up a new machine
├── docs/
│   └── DEPENDENCIES.md      # Per-script dependency matrix (tools / Python / R)
├── tools/
│   ├── smoke_test.sh        # CLI regression: every script must answer -h
│   └── doctor.sh            # Dependency checks, path locator, migration checklist
├── archives/                # Original scripts preserved as archive
└── README.md                # Script index — functional categories live here, not in directories
```

## Quick Start

```bash
# One-time: put the toolkit on your PATH (add to ~/.bashrc)
export PATH="$HOME/myscripts/bin:$PATH"

# One-time per machine: dependency configuration
#   cp config/env.local.sh.example config/env.local.sh   then fill in tool paths
#   (use  tools/doctor.sh --locate <tool>  to find candidates on this machine)
tools/doctor.sh                     # verify — everything OK / PATH / UNSET is listed

# Python scripts — use argparse
python bin/blast_align_analysis.py -h     # or just: blast_align_analysis.py -h (if on PATH)

# R scripts — use getopt
Rscript bin/deseq2_multigroup.R -h

# Bash scripts — use getopts
bash bin/bam_to_bigwig.sh -h
```

---

## Scripts by Category

> The functional categories below exist only in this index — every script itself lives flat in `bin/`. Usage examples assume `bin/` is on your `PATH`; otherwise prefix the script name with `bin/`.

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

## Configuration

All machine-dependent settings live in one place: `config/env.local.sh` (copy from `config/env.local.sh.example`; gitignored, never leaves the machine). `config/env.sh` holds the committed defaults and the `mys_resolve_bin` helper. Every value is a **path** — there is no environment activation in this mechanism; conda-installed binaries are self-contained and work when called by absolute path.

Tool resolution follows one fixed order, everywhere:

```
CLI option  >  MYS_*_BIN path variable  >  PATH  >  error (+ hint to run tools/doctor.sh)
```

Consumed variables:

| Variable | Used by |
|----------|---------|
| `NCBI_EMAIL`, `NCBI_API_KEY` | SRA/ENA fetch & search scripts (E-utilities contact / rate limits) |
| `MYS_THREADS` | default thread count for scripts accepting `-t/--threads` |
| `MYS_PYTHON_BIN`, `MYS_RSCRIPT_BIN`, `MYS_R_LIBS` | interpreters and R package library (exported as `R_LIBS`) |
| `MYS_PREFETCH_BIN`, `MYS_FASTERQ_DUMP_BIN`, `MYS_FASTQ_DUMP_BIN`, … | one slot per external tool — see `config/env.local.sh.example` for the full list |
| `MYS_EXTRA_PATH` | escape hatch: directories prepended to PATH for scripts that look up companion binaries by name |
| `MYS_HOME` | set automatically (repository root) |

Maintenance commands:

```bash
tools/doctor.sh                     # validate everything: tool paths + Python + R packages
tools/doctor.sh --locate ascp       # scan candidate paths, print a ready-to-paste config line
tools/doctor.sh --migrate-check     # checklist of UNSET/BROKEN items (new-machine setup)
tools/smoke_test.sh                 # CLI regression: every script must answer -h
```

Per-script dependency details: [docs/DEPENDENCIES.md](docs/DEPENDENCIES.md).

---

## Script Conventions

Rules for new or modified scripts — keep the toolbox uniform:

- **Location & permissions**: every executable lives flat in `bin/` with the executable bit set. Non-CLI shared code (if it ever appears) goes in `lib/`.
- **Naming**: lowercase `snake_case`, verb-first (`plot_*`, `fetch_*`, `batch_*`, `calculate_*`); language is visible from the extension `.py` / `.R` / `.sh`.
- **CLI**: Python → argparse; R → getopt; Bash → getopts. Every script must answer `-h` with a usage message and exit cleanly.
- **Configuration**: bash scripts source `config/env.sh` (one line, see [Configuration](#configuration)). External tools are resolved through `mys_resolve_bin` (Python: an equivalent `resolve_tool` helper) with the fixed order: CLI option > `MYS_*_BIN` path variable > `PATH`. Never hardcode absolute tool paths in code.
- **Logging**: Python `logging` to stderr; Bash `[INFO]/[WARN]/[ERROR]` prefixes to stderr; R `message()`.
- **Exit codes**: `0` on success, non-zero on failure. Bash: `set -euo pipefail`. Python: `sys.exit(1)` on handled errors. R: `quit(status=1)`.
- **Headers**: description + changelog comment block; no machine-specific absolute paths in code or comments.
- **Regression**: run `tools/smoke_test.sh` after any change — every script must still pass its `-h`.

---

## Dependencies

Per-script requirements are mapped in [docs/DEPENDENCIES.md](docs/DEPENDENCIES.md); tool locations are pinned once in `config/env.local.sh` — see [Configuration](#configuration). Verify the current machine with `tools/doctor.sh`.

### Bioinformatics Tools
sra-tools (prefetch / fasterq-dump / fastq-dump) · GNU parallel · Aspera ascp · BLAST+ (blastn / blastp / makeblastdb) · SAMtools · subread featureCounts · RSeQC (infer_experiment.py) · deepTools (bamCoverage) · UCSC utils (gtf2bed / gff2bed) · bedtools · PLINK / VCFTools · DSSP (mkdssp) · PyMOL · ChimeraX · wget

### Python Packages
biopython, pandas, numpy, matplotlib, networkx, scipy, scikit-learn, python-louvain, tqdm, lxml, pyteomics, pybedtools, markov-clustering

### R / Bioconductor Packages
DESeq2, edgeR, clusterProfiler, GOSemSim, aPEAR, ComplexHeatmap, circlize, ggplot2, ggsci, cowplot, patchwork, aplot, ggtree, Gviz, plotrix, GenomicFeatures, GenomicRanges, rtracklayer, Biostrings, sangerseqR, Mfuzz, Biobase, BiocParallel, gplots, RColorBrewer, amap, magrittr, UpSetR, qqman, data.table, dplyr, ape, getopt — full per-script matrix in [docs/DEPENDENCIES.md](docs/DEPENDENCIES.md)

---

## License

This project is licensed under the [GNU General Public License v3.0](LICENSE).

---

## Author

ChengYu
