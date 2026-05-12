# Script Rewrite Plan — myscripts_raw → myscripts

## Overview

Rewrite 83 scripts from `myscripts_raw/` into `myscripts/`, organized by category.
Scripts are processed in batches of 5-8, grouped by functional similarity so shared
patterns can be established early and reused.

---

## Batch 1 — sequence_analysis (6 scripts)

Establish Python template patterns (argparse, logging, error handling) that all subsequent batches will follow.

| # | Original | New name | Notes |
|---|----------|----------|-------|
| 1 | `align_analysis.py` | `blast_align_analysis.py` | BLAST-based similarity, concurrent processing |
| 2 | `align_analysis2.py` | `aligned_fasta_similarity.py` | Direct alignment comparison, no BLAST |
| 3 | `sequence_network.py` | `sequence_similarity_network.py` | Jaccard + Louvain clustering |
| 4 | `sequence_network_blast.py` | `blast_sequence_network.py` | BLAST + MCL clustering |
| 5 | `seq_revcom.py` | `reverse_complement.py` | Simple FASTA reverse complement |
| 6 | `calc_seq_complexity.py` | `sequence_complexity.py` | LZ complexity + entropy (moved from data_processing) |

## Batch 2 — codon_analysis (6 scripts)

| # | Original | New name | Notes |
|---|----------|----------|-------|
| 7 | `codon_freq.py` | `codon_frequency.py` | Codon usage table |
| 8 | `aa_freq.py` | `amino_acid_frequency.py` | AA composition |
| 9 | `calc_cai_batch.py` | `calculate_cai.py` | Currently uses sys.argv — needs argparse |
| 10 | `Kozak_similarty_score.py` | `kozak_similarity_score.py` | Fix typo in original name |
| 11 | `batch_dNdS.py` | `calculate_dnds.py` | Pairwise dN/dS |
| 12 | `calculate_average_pi.py` | `calculate_nucleotide_diversity.py` | VCFTools wrapper |

## Batch 3 — format_conversion (5 scripts)

- #12 gff2gtf + #13 gtf2bed 合并为一个 UCSC 工具 wrapper 脚本
- #16 bam_to_bigwig 增加更多标准化方法支持（RPKM/TPM/CPM/None）

| # | Original | New name | Notes |
|---|----------|----------|-------|
| 13 | `gff2gtf.sh` + `gtf2bed.sh` | `genome_format_converter.sh` | 合并：UCSC 工具 wrapper，支持 gff2gtf / gtf2bed 等多种转换 |
| 14 | `gtf_simplified.sh` | `gtf_standardize.sh` | GTF attribute standardization |
| 15 | `calculate_bed12_length.py` | `bed12_effective_length.py` | Non-redundant length from BED12 |
| 16 | `bam2bw.sh` | `bam_to_bigwig.sh` | BAM → bigWig，支持 RPKM/TPM/CPM/None 多种标准化 |
| 17 | `fasta2AlphaFoldserver_JSON.py` | `fasta_to_alphafold_json.py` | FASTA → AlphaFold JSON |
| ~~18~~ | ~~`gtf2gff3.sh`~~ | **SKIP** | Empty file, no implementation |

## Batch 4 — data_processing (4 scripts)

- Skip #20 `merge.R` 和 #21 `njoin.sh`
- #24 `calc_seq_complexity.py` 已移至 Batch 1 sequence_analysis

| # | Original | New name | Notes |
|---|----------|----------|-------|
| 19 | `batch_join.py` | `multi_file_join.py` | Join by common column |
| 22 | `transposition.sh` | `transpose_table.sh` | Matrix transpose |
| 23 | `dup_number.R` | `number_duplicates.R` | Duplicate enumeration |
| 25 | `calc_tau.py` | `tissue_specificity_tau.py` | Currently uses sys.argv — needs argparse |

## Batch 5 — gene_annotation (4 scripts)

Establish patterns for genomics coordinate handling.

- Skip #26 `ORFannotate.pl`（Perl 版本）
- Skip #27 `ORFannotate.py`（Python 版本也跳过）
- Skip #31 #32 #33 三个 `maf_ORFsubstr` 版本（均跳过，不合并）
- 仅保留 UTR 提取、intron 过滤和 MAF region 提取

| # | Original | New name | Notes |
|---|----------|----------|-------|
| 28 | `extract_UTR.py` | `extract_utr.py` | 5'/3' UTR from GTF |
| 29 | `get5UTR_genePred2BED.py` | `genepred_utr_to_bed12.py` | genePred → BED12 UTR |
| 30 | `filter_long_intron.py` | `filter_long_introns.py` | Filter by intron length |
| 34 | `ExtractBedRegionsFromMAF.py` | `maf_extract_regions.py` | BED12-based MAF extraction |

## Batch 6 — rnaseq (8 scripts)

Establish R script template patterns.

| # | Original | New name | Notes |
|---|----------|----------|-------|
| 35 | `runDESeq2_multiGroup_230705.R` | `deseq2_multigroup.R` | Full DESeq2 pipeline |
| 36 | `DEG_UseWilcoxTest.R` | `deg_wilcox_test.R` | Wilcoxon-based DEG |
| 37 | `run_featurecount_pipeline_251114.sh` | `featurecounts_pipeline.sh` | QC + featureCounts |
| 38 | `featureCount.R_result_merge.py` | `merge_featurecounts.py` | Count matrix merging |
| 39 | `cal_expr_FPKM_TPM.R` | `calculate_fpkm_tpm.R` | Expression from STAR counts |
| 40 | `calculate_EffectiveLength.R` | `calculate_effective_length.R` | Non-redundant exon length |
| 41 | `TPM_count.R` | `calculate_tpm.R` | Count → TPM/FPKM |
| 42 | `filter_expression_matrix.py` | `filter_expression.py` | Low expression filter |

## Batch 7 — data_retrieval (12 scripts)

| # | Original | New name | Notes |
|---|----------|----------|-------|
| 43 | `batch_prefetch_251118.sh` | `sra_prefetch_batch.sh` | SRA download pipeline |
| ~~44~~ | ~~`batch_prefetch_deprecated.sh`~~ | **SKIP** | Superseded |
| 45 | `batch_ascp_ena.py` | `ena_ascp_download.py` | ENA Aspera download |
| 46 | `download_ascp.sh` | `ena_ascp_download_batch.sh` | ENA stable downloader |
| 47 | `fetch_sra_info.py` | `fetch_sra_metadata_ncbi.py` | NCBI SRA metadata |
| 48 | `fetch_sra_info_ena.py` | `fetch_sra_metadata_ena.py` | ENA metadata |
| 49 | `batch_fetch_sra_metadata.py` | `fetch_sra_metadata_comprehensive.py` | Comprehensive SRA metadata |
| 50 | `get_sra_full_metadata.py` | `fetch_sra_metadata_xml.py` | XML-based SRA metadata |
| 51 | `format_sra_metadata.py` | `format_sra_summary.py` | Metadata → standardized summary |
| 52 | `fetch_bioProject_Info.sh` | `fetch_bioproject_metadata.sh` | BioProject from ENA |
| 53 | `batch_pride_fetch.py` | `fetch_pride_metadata.py` | PRIDE proteomics metadata |
| 54 | `parse_json_pride.py` | `parse_pride_json.py` | PRIDE JSON parser |

## Batch 8 — data_retrieval (continued) + enrichment (4 scripts)

| # | Original | New name | Notes |
|---|----------|----------|-------|
| 55 | `search_sra_riboseq.py` | `search_sra_riboseq.py` | Ribo-seq SRA search |
| 56 | `search_riboSeq_bioproject_improved251117.py` | `search_riboseq_bioproject.py` | Multi-layer BioProject retrieval |
| 57 | `enrichKEGG.R` | `kegg_enrichment.R` | KEGG pathway enrichment |
| 58 | `getSpeciesRNA_fromRfam_251113.sh` | `extract_species_rna_rfam.sh` | Species RNA from Rfam |

## Batch 9 — proteomics (9 scripts)

Skip #59 `MSGFplus_results_stat.py` 和 #60 `filter_first_search_results.py`。

| # | Original | New name | Notes |
|---|----------|----------|-------|
| ~~59~~ | ~~`MSGFplus_results_stat.py`~~ | **SKIP** | |
| ~~60~~ | ~~`filter_first_search_results.py`~~ | **SKIP** | |
| 61 | `plot_MS2.py` | `plot_ms2_spectrum.py` | Annotated MS2 spectra |
| 62 | `batch_runDSSP.sh` | `batch_dssp.sh` | DSSP secondary structure |
| 63 | `dssp_output_summary.py` | `dssp_summary.py` | DSSP statistics |
| 64 | `batch_pymol_images.py` | `batch_pymol_render.py` | PyMOL batch rendering |
| 65 | `chimerax_batch_render.py` | `batch_chimerax_render.py` | ChimeraX batch rendering |
| 66 | `plot_3Dstructure.py` | `plot_structure_3d.py` | Hydropathy/SS structure plots |
| 67 | `calculate_peptide_properties.py` | `peptide_properties.py` | Physicochemical properties |
| 68 | `calculate_hydropathy_distribution.py` | `hydropathy_distribution.py` | Sliding window hydropathy |
| 69 | `protein_weight.py` | `protein_molecular_weight.py` | Protein mass calculation |

## Batch 10 — visualization (10 scripts)

| # | Original | New name | Notes |
|---|----------|----------|-------|
| 70 | `plotHeatmap.R` | `plot_heatmap_multi.R` | ComplexHeatmap multi-cluster |
| 71 | `plotHeatmap_singlePlot.R` | `plot_heatmap_single.R` | Individual heatmaps |
| 72 | `plot_BaseContent.R` | `plot_base_content.R` | Sliding window base composition |
| 73 | `plot_flower.R` | `plot_flower.R` | Circular/flower plots |
| 74 | `plot_geneStructure.R` | `plot_gene_structure.R` | Gene structure from GTF |
| 75 | `plot_Upset.R` | `plot_upset.R` | UpSet intersection plots |
| 76 | `gwas_plot.R` | `plot_gwas_manhattan.R` | Manhattan + QQ plots |
| 77 | `qPCR_analysis_plot.R` | `plot_qpcr_analysis.R` | qPCR delta-delta CT |
| 78 | `240313_GO_tree_anno.note.R` | `plot_go_tree_heatmap.R` | GO tree + heatmap |
| 79 | `run_plot_GO_Heatmap.sh` | `plot_go_heatmap_batch.sh` | Batch GO heatmaps |

## Batch 11 — ml_stats + workflow utilities (3 scripts)

Skip #80 `train_xgboost_rice_lncORF_250815.py`。

| # | Original | New name | Notes |
|---|----------|----------|-------|
| ~~80~~ | ~~`train_xgboost_rice_lncORF_250815.py`~~ | **SKIP** | |
| 81 | `run_MfuzzCluster.R` | `mfuzz_soft_cluster.R` | Time-series soft clustering |
| 82 | `batchParsingSangerSeq_231206.R` | `batch_parse_sanger.R` | AB1 Sanger sequencing |
| 83 | `batchSangerBlast.sh` | `batch_sanger_blast.sh` | Sanger BLAST pipeline |

---

## Summary

| Category | Rewritten | Skipped |
|----------|-----------|---------|
| sequence_analysis | 6 | — |
| codon_analysis | 6 | — |
| format_conversion | 5 | 1 (empty file) + 1 merge (gff2gtf+gtf2bed→1) |
| data_processing | 4 | 2 (merge.R, njoin.sh) + 1 moved (→ seq_analysis) |
| gene_annotation | 4 | 3 (ORFannotate×2, Perl) + 3 (maf_ORFsubstr×3) |
| rnaseq | 8 | — |
| data_retrieval | 12 | 1 (deprecated) |
| enrichment | 2 | — |
| proteomics | 9 | 2 (MSGFplus, filter_ms) |
| visualization | 10 | — |
| ml_stats + misc | 3 | 1 (xgboost) |
| **Total** | **69** | **~14** |

## Verification

After each batch:
1. Run `script.py -h` / `Rscript script.R --help` / `bash script.sh -h` to verify CLI
2. Confirm logging output works with `--log-level DEBUG`
3. Test with minimal input data if available
4. Check that output format matches original script behavior

Final verification:
- Compare output directory structure against plan
- Ensure no original scripts were missed
- Run all scripts with `--help` to confirm consistent CLI patterns
