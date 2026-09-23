# Dependency Matrix

What each script in `bin/` needs at runtime. Verify against the current machine with `tools/doctor.sh`; pin tool locations in `config/env.local.sh` (see [Configuration](../README.md#configuration)). Scripts marked *(parses output)* consume files produced elsewhere and do not invoke the binary themselves.

Legend: **Tools** = external executables · **Py** = Python packages (import names) · **R** = R packages.

## Sequence Analysis

| Script | Tools | Py |
|--------|-------|----|
| `blast_align_analysis.py` | blastn, blastp | Bio, pandas |
| `aligned_fasta_similarity.py` | — | Bio, pandas |
| `sequence_similarity_network.py` | — | Bio, matplotlib, networkx, numpy, scipy |
| `blast_sequence_network.py` | blastp, makeblastdb | Bio, matplotlib, networkx, numpy |
| `reverse_complement.py` | — | — |
| `sequence_complexity.py` | — | Bio |

## Codon Analysis

| Script | Tools | Py |
|--------|-------|----|
| `codon_frequency.py` | — | — |
| `amino_acid_frequency.py` | — | — |
| `calculate_cai.py` | — | Bio |
| `kozak_similarity_score.py` | — | numpy |
| `calculate_dnds.py` | — | Bio |
| `calculate_nucleotide_diversity.py` | plink, vcftools | pandas |

## Format Conversion

| Script | Tools | Py |
|--------|-------|----|
| `genome_format_converter.sh` | gtf2bed, gff2bed; UCSC kent utils gff2gtf, gtfToGenePred, gff3ToGenePred, genePredToBed (see note) | — |
| `gtf_standardize.sh` | — (self-contained awk) | — |
| `bam_to_bigwig.sh` | samtools, bamCoverage (deepTools) | — |
| `bed12_effective_length.py` | bedtools | — |
| `fasta_to_alphafold_json.py` | — | — |

> **Note (this machine):** `gtf2bed`/`gff2bed` are pinned in `config/env.local.sh`. The kent utilities `gff2gtf`, `gtfToGenePred`, `gff3ToGenePred`, `genePredToBed` are **not installed** — the `gff2gtf`, `gtf2gp`, `gp2bed` modes of `genome_format_converter.sh` are unavailable until they are provided (install kent-utils, then pass paths via the script's own options).

## Data Processing

| Script | Tools | Py |
|--------|-------|----|
| `multi_file_join.py` | — | pandas |
| `transpose_table.sh` | — | — |
| `number_duplicates.R` | — | — |
| `tissue_specificity_tau.py` | — | numpy, pandas |

## Gene Annotation

| Script | Tools | Py |
|--------|-------|----|
| `extract_utr.py` | — | Bio |
| `genepred_utr_to_bed12.py` | — | — |
| `filter_long_introns.py` | — | tqdm |
| `maf_extract_regions.py` | — | Bio |

## RNA-seq

| Script | Tools | Py | R |
|--------|-------|----|---|
| `deseq2_multigroup.R` | — | — | DESeq2, ggplot2, BiocParallel, gplots, RColorBrewer, amap, getopt |
| `deg_wilcox_test.R` | — | — | edgeR, getopt |
| `featurecounts_pipeline.sh` | featureCounts, samtools | — | — |
| `merge_featurecounts.py` | — (parses output) | pandas | — |
| `calculate_fpkm_tpm.R` | — | — | GenomicFeatures, getopt |
| `calculate_effective_length.R` | — | — | GenomicFeatures, getopt |
| `calculate_tpm.R` | — (function library — sourced, not run) | — | — |
| `filter_expression.py` | — | numpy, pandas | — |

## Data Retrieval

| Script | Tools | Py | Credentials |
|--------|-------|----|-------------|
| `sra_prefetch_batch.sh` | prefetch, fasterq-dump / fastq-dump, parallel, fastqc (optional) | — | — |
| `ena_ascp_download.py` | ascp | — | — |
| `ena_ascp_download_batch.sh` | ascp | — | — |
| `fetch_sra_metadata_ncbi.py` | — (HTTPS) | — | NCBI_EMAIL, NCBI_API_KEY |
| `fetch_sra_metadata_ena.py` | — (HTTPS) | — | — |
| `fetch_sra_metadata_comprehensive.py` | — (HTTPS) | — | NCBI_EMAIL, NCBI_API_KEY |
| `fetch_sra_metadata_xml.py` | — (HTTPS) | — | NCBI_EMAIL, NCBI_API_KEY |
| `format_sra_summary.py` | — | — | — |
| `fetch_bioproject_metadata.sh` | wget | — | — |
| `fetch_pride_metadata.py` | — (HTTPS) | — | — |
| `parse_pride_json.py` | — | — | — |
| `search_sra_riboseq.py` | — (HTTPS) | Bio, tqdm | NCBI_EMAIL, NCBI_API_KEY |
| `search_riboseq_bioproject.py` | — (HTTPS) | Bio | NCBI_EMAIL, NCBI_API_KEY |
| `extract_species_rna_rfam.sh` | — | — | — |

## Enrichment

| Script | Tools | R |
|--------|-------|---|
| `kegg_enrichment.R` | — | clusterProfiler, aPEAR, ggplot2, magrittr, getopt |

## Proteomics

| Script | Tools | Py |
|--------|-------|----|
| `plot_ms2_spectrum.py` | — | lxml, pyteomics, matplotlib, numpy |
| `batch_dssp.sh` | dssp or mkdssp, parallel | — |
| `dssp_summary.py` | — (parses output) | pandas |
| `batch_pymol_render.py` | pymol | — |
| `batch_chimerax_render.py` | ChimeraX | — |
| `plot_structure_3d.py` | — (parses PDB/DSSP output) | matplotlib, numpy |
| `peptide_properties.py` | — | Bio, matplotlib |
| `hydropathy_distribution.py` | — | Bio, matplotlib, numpy |
| `protein_molecular_weight.py` | — | Bio |

## Visualization

| Script | Tools | R |
|--------|-------|---|
| `plot_heatmap_multi.R` | — | circlize, ComplexHeatmap, getopt |
| `plot_heatmap_single.R` | — | circlize, ComplexHeatmap, getopt, RColorBrewer |
| `plot_base_content.R` | — | Biostrings, getopt, ggplot2 |
| `plot_flower.R` | — | getopt, plotrix |
| `plot_gene_structure.R` | — | dplyr, getopt, ggplot2, rtracklayer |
| `plot_upset.R` | — | getopt, UpSetR |
| `plot_gwas_manhattan.R` | — | cowplot, data.table, GenomicRanges, getopt, ggplot2, Gviz, qqman |
| `plot_qpcr_analysis.R` | — | getopt, ggplot2, ggsci, patchwork |
| `plot_go_tree_heatmap.R` | — | ape, aplot, getopt, ggplot2, ggtree, GOSemSim |
| `plot_go_heatmap_batch.sh` | Rscript + `plot_heatmap_multi.R` (same directory) | — |

## ML & Statistics

| Script | Tools | R |
|--------|-------|---|
| `mfuzz_soft_cluster.R` | — | Biobase, getopt, Mfuzz |
| `batch_parse_sanger.R` | — | Biostrings, getopt, sangerseqR |
| `batch_sanger_blast.sh` | blastn, makeblastdb | — |

---

## Summary: machine-level requirements

- **External tools** (18): prefetch, fasterq-dump, fastq-dump, fastqc, parallel, ascp, blastn, blastp, makeblastdb, mkdssp, pymol, ChimeraX, plink, vcftools, featureCounts, samtools, bedtools, bamCoverage, gtf2bed, gff2bed — pin each in `config/env.local.sh`.
- **Python** (9): biopython, pandas, numpy, matplotlib, networkx, scipy, tqdm, lxml, pyteomics.
- **R** (33): see tables above — the Rscript interpreter and library path come from `MYS_RSCRIPT_BIN` / `MYS_R_LIBS`.
- **Credentials**: `NCBI_EMAIL`, `NCBI_API_KEY` for NCBI E-utilities scripts.
