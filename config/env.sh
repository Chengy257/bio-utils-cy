#!/usr/bin/env bash
# =============================================================================
# bio-utils-cy — unified environment / dependency configuration
# =============================================================================
# Single place for machine-dependent settings shared by all scripts in bin/.
#
# Bash scripts load it with one line (after their `set -` line):
#     _my_conf="$(dirname "$(readlink -f "${BASH_SOURCE[0]}")")/../config/env.sh"
#     if [ -f "${_my_conf}" ]; then . "${_my_conf}"; fi
#
# Python / R scripts do NOT source this file. They read the same variables
# from the environment (os.environ / Sys.getenv), so one configuration
# serves every language.
#
# This file holds committed, machine-independent defaults. Machine-specific
# values (tool paths, credentials) live in config/env.local.sh — gitignored,
# sourced at the end of this file (it always has the last word).
# =============================================================================

# Guard against double-sourcing.
if [ -n "${BUC_CONFIG_LOADED:-}" ]; then
    return 0 2>/dev/null || true
fi
BUC_CONFIG_LOADED=1

# --- Repository layout -------------------------------------------------------
# Absolute path of this repository, derived from this file's location.
BUC_HOME="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"

# --- Shared defaults (overridable in env.local.sh or the environment) ---------
# Default thread count for scripts that accept -t/--threads.
: "${BUC_THREADS:=8}"

# NCBI E-utilities credentials (consumed by the fetch/search data-retrieval
# scripts via os.environ).
: "${NCBI_EMAIL:=}"
: "${NCBI_API_KEY:=}"

# --- Tool path slots ----------------------------------------------------------
# Empty = fall back to PATH at resolution time (buc_resolve_bin / doctor).
# Set an absolute path per machine in config/env.local.sh to pin a binary.
: "${BUC_PYTHON_BIN:=}"
: "${BUC_RSCRIPT_BIN:=}"
: "${BUC_R_LIBS:=}"
: "${BUC_PREFETCH_BIN:=}"
: "${BUC_FASTERQ_DUMP_BIN:=}"
: "${BUC_FASTQ_DUMP_BIN:=}"
: "${BUC_FASTQC_BIN:=}"
: "${BUC_PARALLEL_BIN:=}"
: "${BUC_ASCP_BIN:=}"
: "${BUC_BLASTN_BIN:=}"
: "${BUC_BLASTP_BIN:=}"
: "${BUC_MAKEBLASTDB_BIN:=}"
: "${BUC_MKDSSP_BIN:=}"
: "${BUC_PYMOL_BIN:=}"
: "${BUC_CHIMERAX_BIN:=}"
: "${BUC_PLINK_BIN:=}"
: "${BUC_VCFTOOLS_BIN:=}"
: "${BUC_FEATURECOUNTS_BIN:=}"
: "${BUC_SAMTOOLS_BIN:=}"
: "${BUC_INFER_EXP_BIN:=}"
: "${BUC_BEDTOOLS_BIN:=}"
: "${BUC_BAMCOVERAGE_BIN:=}"
: "${BUC_GTF2BED_BIN:=}"
: "${BUC_GFF2BED_BIN:=}"

# Escape hatch for tools that look up companion binaries by name at runtime:
# directories prepended to PATH for script subprocesses. Normally empty.
: "${BUC_EXTRA_PATH:=}"

# Where tools/doctor.sh --locate searches for candidate binaries.
: "${BUC_LOCATE_DIRS:=$HOME/soft/miniconda3/envs:/opt/anaconda3/envs:/usr/local/bin:$HOME/soft}"

# --- Tool resolution ----------------------------------------------------------
# Fixed order: explicit BUC_<NAME>_BIN path > PATH.
# Usage:
#     samtools="$(buc_resolve_bin BUC_SAMTOOLS_BIN samtools)" || exit 1
buc_resolve_bin() {
    local _var="$1" _tool="$2" _cand
    _cand="${!_var:-}"
    if [ -n "${_cand}" ]; then
        if [ -x "${_cand}" ]; then
            printf '%s\n' "${_cand}"
            return 0
        fi
        printf 'buc_resolve_bin: %s="%s" is not executable (check config/env.local.sh)\n' "${_var}" "${_cand}" >&2
        return 1
    fi
    if command -v -- "${_tool}" >/dev/null 2>&1; then
        command -v -- "${_tool}"
        return 0
    fi
    printf 'buc_resolve_bin: tool "%s" not found — set %s in config/env.local.sh or run tools/doctor.sh --locate %s\n' "${_tool}" "${_var}" "${_tool}" >&2
    return 1
}
export -f buc_resolve_bin 2>/dev/null || true

# --- Machine-specific overrides (gitignored; sourced last, has the last word) --
_local="$(dirname "${BASH_SOURCE[0]}")/env.local.sh"
if [ -f "$_local" ]; then
    . "$_local"
fi
unset -v _local

# R package library: expose BUC_R_LIBS to R itself (prepend, keep any existing).
if [ -n "${BUC_R_LIBS}" ]; then
    export R_LIBS="${BUC_R_LIBS}${R_LIBS:+:${R_LIBS}}"
fi

# PATH escape hatch (see BUC_EXTRA_PATH above).
if [ -n "${BUC_EXTRA_PATH}" ]; then
    export PATH="${BUC_EXTRA_PATH%/}:${PATH}"
fi

# --- Export everything scripts should see -------------------------------------
export BUC_HOME BUC_CONFIG_LOADED BUC_THREADS NCBI_EMAIL NCBI_API_KEY \
    BUC_PYTHON_BIN BUC_RSCRIPT_BIN BUC_R_LIBS \
    BUC_PREFETCH_BIN BUC_FASTERQ_DUMP_BIN BUC_FASTQ_DUMP_BIN BUC_FASTQC_BIN \
    BUC_PARALLEL_BIN BUC_ASCP_BIN BUC_BLASTN_BIN BUC_BLASTP_BIN \
    BUC_MAKEBLASTDB_BIN BUC_MKDSSP_BIN BUC_PYMOL_BIN BUC_CHIMERAX_BIN \
    BUC_PLINK_BIN BUC_VCFTOOLS_BIN BUC_FEATURECOUNTS_BIN BUC_SAMTOOLS_BIN \
    BUC_INFER_EXP_BIN \
    BUC_BEDTOOLS_BIN BUC_BAMCOVERAGE_BIN BUC_GTF2BED_BIN BUC_GFF2BED_BIN \
    BUC_EXTRA_PATH BUC_LOCATE_DIRS
