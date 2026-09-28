# Dependency Matrix

What each script in `bin/` needs at runtime. Verify against the current machine with `tools/doctor.sh`; pin tool locations in `config/env.local.sh` (see [Configuration](../README.md#configuration)). Scripts marked *(parses output)* consume files produced elsewhere and do not invoke the binary themselves.

Legend: **Tools** = external executables · **Py** = Python packages (import names) · **R** = R packages.

## Sequence Analysis

| Script | Tools | Py |
|--------|-------|----|
| `blast_align_analysis.py` | blastn, blastp | Bio, pandas |
| `aligned_fasta_similarity.py` | — | Bio, pandas |
| `sequence_similarity_network.py` | — | Bio, matplotlib, networkx, scipy, sklearn, community (python-louvain) |
| `blast_sequence_network.py` | blastp, makeblastdb | Bio, matplotlib, markov_clustering, networkx, numpy |
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
| `bed12_effective_length.py` | — | — |
| `fasta_to_alphafold_json.py` | — | — |

> **Note (this machine):** `gtf2bed`/`gff2bed` (bedops) are pinned in `config/env.local.sh`. The UCSC kent utilities (`gtfToGenePred`, `gff3ToGenePred`, `genePredToGtf`, `genePredToBed`) live in `~/soft/ucsc-tools` and are on PATH via `~/soft/bin` symlinks, so all five conversion modes are usable; each tool can also be pinned per-path with the script's own `--xxx` options or `BUC_<TOOL>_BIN` slots.

## Data Processing

| Script | Tools | Py |
|--------|-------|----|
| `multi_file_join.py` | — | pandas |
| `transpose_table.sh` | — | — |
| `number_duplicates.R` | getopt | — |
| `tissue_specificity_tau.py` | — | numpy, pandas |

## Gene Annotation

| Script | Tools | Py |
|--------|-------|----|
| `extract_utr.py` | — | Bio |
| `genepred_utr_to_bed12.py` | — | — |
| `filter_long_introns.py` | — | — |
| `maf_extract_regions.py` | — | Bio |

## RNA-seq

| Script | Tools | Py | R |
|--------|-------|----|---|
| `deseq2_multigroup.R` | — | — | DESeq2, ggplot2, BiocParallel, pheatmap, RColorBrewer, ashr, getopt |
| `deg_wilcox_test.R` | — | — | edgeR, getopt |
| `featurecounts_pipeline.sh` | featureCounts, samtools, infer_experiment.py (RSeQC) | — | — |
| `merge_featurecounts.py` | — (parses output) | pandas | — |
| `calculate_fpkm_tpm.R` | — | — | GenomicFeatures, getopt |
| `calculate_effective_length.R` | — | — | GenomicFeatures, getopt |
| `lib/calculate_tpm.R` | — (function library — sourced, not run) | — | — |
| `filter_expression.py` | — | numpy, pandas | — |

## Data Retrieval

> Eight scripts after consolidation — one per database and function. Superseded variants (the batch ENA download wrapper, per-variant NCBI fetchers, the PRIDE JSON parser, and the BioProject filereport shell wrapper) live in the local `archives/`.

| Script | Tools | Py | Credentials |
|--------|-------|----|-------------|
| `sra_prefetch_batch.sh` | prefetch, fasterq-dump / fastq-dump, parallel, fastqc (optional) | — | — |
| `ena_ascp_download.py` | ascp | — | — |
| `fetch_sra_metadata_ncbi.py` | — (HTTPS) | — | NCBI_EMAIL, NCBI_API_KEY |
| `fetch_sra_metadata_ena.py` | — (HTTPS) | — | — |
| `pride_metadata.py` | — (HTTPS for `fetch`) | openpyxl (optional — xlsx output; falls back to CSV) | — |
| `search_sra_riboseq.py` | — (HTTPS) | Bio, tqdm | NCBI_EMAIL, NCBI_API_KEY |
| `format_sra_summary.py` | — | — | — |
| `extract_species_rna_rfam.sh` | seqkit (BUC_SEQKIT_BIN) | — | — |

## Enrichment

| Script | Tools | R |
|--------|-------|---|
| `kegg_enrichment.R` | — (KEGG, needs network) | clusterProfiler, aPEAR, ggplot2, dplyr, R.utils, magrittr, getopt |

## Proteomics

| Script | Tools | Py |
|--------|-------|----|
| `plot_ms2_spectrum.py` | — | matplotlib, numpy, pymzml (optional — only for mzML input; MGF/mzid are stdlib) |
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

- **External tools** (20): prefetch, fasterq-dump, fastq-dump, fastqc, parallel, ascp, blastn, blastp, makeblastdb, mkdssp, pymol, ChimeraX, plink, vcftools, featureCounts, samtools, infer_experiment.py, bamCoverage, gtf2bed, gff2bed — pin each in `config/env.local.sh`. (`bedtools` is no longer needed by any script; the former `BUC_BEDTOOLS_BIN` slot was removed from `config/env.sh`.)
- **Python** (10 hard requirements): biopython, pandas, numpy, matplotlib, networkx, scipy, scikit-learn, community (python-louvain), tqdm, markov-clustering. Optional: pymzml (mzML input for plot_ms2_spectrum), openpyxl (xlsx output for pride_metadata).
- **R** (34): see tables above — the Rscript interpreter and library path come from `BUC_RSCRIPT_BIN` / `BUC_R_LIBS`. If a package fails with a `dyn.load` error (typically stringi resolving `libicui18n.so.70`), set `BUC_LD_LIBRARY_PATH` in `config/env.local.sh` to the directory holding those conda libraries (e.g. `$HOME/soft/miniconda3/lib`).
- **Credentials**: `NCBI_EMAIL`, `NCBI_API_KEY` for NCBI E-utilities scripts.
