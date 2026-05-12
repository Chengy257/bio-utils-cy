# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

Personal bioinformatics toolkit for next-generation sequencing (NGS) data analysis. Contains standalone utility scripts and Snakemake workflow templates covering RNA-seq, ChIP-seq, CUT&Tag, ATAC-seq, sRNA-seq, ribosome profiling, lncRNA analysis, and proteomics.

## Repository Layout

- **`myscripts_raw/`** — Original utility scripts (Python, R, Bash) — source for rewrite
- **`myscripts/`** — Rewritten, optimized scripts organized by category (output directory)
- **`workflows/snakemake/`** — Snakemake workflow templates (not part of rewrite)
- **`myscripts_new/`** — Legacy placeholder, unused

## Script Rewrite Project

### Goal

Rewrite all 83 scripts from `myscripts_raw/` into `myscripts/`, organized by functional category. Each rewritten script must be production-quality with robust error handling, comprehensive CLI, and consistent code style.

### Output directory structure

```
myscripts/
├── sequence_analysis/      # Sequence alignment, similarity networks, reverse complement
├── rnaseq/                 # DEG analysis, expression quantification, count merging
├── data_retrieval/         # SRA/ENA/PRIDE downloading and metadata fetching
├── proteomics/             # MS/MS analysis, peptide properties, structure visualization
├── gene_annotation/        # ORF annotation, UTR extraction, MAF region extraction
├── codon_analysis/         # Codon/AA frequency, CAI, dN/dS, Kozak score
├── format_conversion/      # GFF/GTF/BED/BAM format converters, FASTA utilities
├── visualization/          # Heatmaps, volcano, Manhattan, UpSet, gene structure plots
├── ml_stats/               # Machine learning, statistics, clustering
├── data_processing/        # File joining, merging, transposition, general utilities
└── enrichment/             # GO/KEGG enrichment, functional annotation
```

### Rewrite standards — Python

Every Python script must include:

**CLI (argparse)**
- `description` summarizing purpose
- `epilog` with usage examples
- `formatter_class=argparse.RawDescriptionHelpFormatter`
- Short and long flags (`-i`/`--input`)
- Proper types, defaults, and `required=True` for essential args
- `--version` flag
- Input file existence validation after parse

**Logging**
- `logging` module (never `print()` for diagnostic output)
- Configurable via `--log-level` (default INFO)
- Console + optional file handler via `--log-file`
- Format: `"%(asctime)s - %(levelname)s - %(message)s"`

**Error handling**
- Specific exceptions (never bare `except`)
- `sys.exit(1)` with descriptive message on fatal errors
- File existence and format validation before processing
- Context managers (`with`) for all file I/O

**Code structure**
- Module docstring with purpose, inputs, outputs
- `main()` function behind `if __name__ == "__main__"`
- Single-responsibility functions with docstrings
- Type hints on function signatures
- Named constants (no magic numbers)
- `concurrent.futures` for parallel processing where applicable

**Robustness**
- `--threads`/`-t` for parallelizable scripts (default: 4)
- `--force` to allow output overwrite
- Progress indication (`tqdm` or logging) for long operations
- Graceful handling of empty input / edge cases
- Exit code 0 on success, 1 on error

### Rewrite standards — R

- Use `argparse` or `getopt` with full help text
- `suppressMessages(library())` for all package loads
- File existence checks with `stop()` on missing files
- PDF/PNG output with configurable dimensions
- `tryCatch` for error handling in critical sections

### Rewrite standards — Bash

- `set -euo pipefail`
- `getopts` with `usage()` function and examples
- `--help`/`-h` support
- Trap for cleanup on exit
- All paths and tool names parameterized (no hardcoded paths)

### Process

1. Read the original script fully to understand logic and edge cases
2. Rewrite following the standards above
3. Preserve the original algorithmic logic — optimize robustness, not science
4. Test with `--help` to verify CLI
5. Mark as complete in plan document

### Naming convention

Keep original descriptive names in snake_case. Add version suffix only if multiple versions exist (e.g., `maf_extract_cds.py` replaces `maf_ORFsubstr.py`, `maf_ORFsubstr.2411.py`, `maf_ORFsubstr.2412.py`).

### Deprecated scripts to skip

- `batch_prefetch_deprecated.sh` — superseded by `batch_prefetch_251118.sh`
- `gtf2gff3.sh` — empty file, no implementation
- Duplicate versions merged into one (see plan document)

## Running Scripts

```bash
# Python — use argparse
python script.py -h
python script.py -i input.txt -o output.txt -t 8

# R — use getopt or argparse
Rscript script.R --help

# Shell — use getopts
bash script.sh -h
bash script.sh -i input -o output
```

## Snakemake Workflows

Standard workflow directory structure:

```
workflow-name/
├── main_run.sh       # Entry point — submits Snakemake to cluster
├── config.yaml       # All parameters, reference paths, thread counts
├── *.smk             # Snakemake rule files
├── rules/            # Modular rule files included by main .smk
├── scripts/          # Helper scripts called by rules
└── logs/             # Execution logs
```

Run via `main_run.sh` — submits jobs to PBS/Torque cluster via `qsub`. Configuration in `config.yaml`.

## Coding Conventions

### File headers (all languages)

```
#########################################################################
# File Name: /path/to/script.ext
# Author: ChengYu
# Description: ...
# Created Time: ...
#########################################################################
```

### Snakemake rules

- Reference config via `config["key"]`
- Extract sample list from CSV with a `get_samples()` function
- Every rule has a `log` directive
- Use `conda` directive for per-rule environments
- Use `threads: int(config["threads"])` for thread allocation

## Key Dependencies

No centralized requirements file. Dependencies managed per-workflow via conda environments.

- **Aligners**: STAR, HISAT2, bowtie2
- **Peak callers**: MACS2
- **QC**: FastQC, MultiQC, Trim Galore
- **Quantification**: featureCounts, StringTie
- **Python**: Biopython, pandas, numpy, matplotlib, networkx
- **R/Bioconductor**: DESeq2, clusterProfiler, ChIPseeker, ggplot2, rtracklayer
