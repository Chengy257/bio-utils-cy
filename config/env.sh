#!/usr/bin/env bash
# =============================================================================
# myscripts — unified environment / dependency configuration
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
if [ -n "${MYS_CONFIG_LOADED:-}" ]; then
    return 0 2>/dev/null || true
fi
MYS_CONFIG_LOADED=1

# --- Repository layout -------------------------------------------------------
# Absolute path of this repository, derived from this file's location.
MYS_HOME="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"

# --- Shared defaults (overridable in env.local.sh or the environment) ---------
# Default thread count for scripts that accept -t/--threads.
: "${MYS_THREADS:=8}"

# NCBI E-utilities credentials (consumed by the fetch/search data-retrieval
# scripts via os.environ).
: "${NCBI_EMAIL:=}"
: "${NCBI_API_KEY:=}"

# --- Tool path slots ----------------------------------------------------------
# Empty = fall back to PATH at resolution time (mys_resolve_bin / doctor).
# Set an absolute path per machine in config/env.local.sh to pin a binary.
: "${MYS_PYTHON_BIN:=}"
: "${MYS_RSCRIPT_BIN:=}"
: "${MYS_R_LIBS:=}"
: "${MYS_PREFETCH_BIN:=}"
: "${MYS_FASTERQ_DUMP_BIN:=}"
: "${MYS_FASTQ_DUMP_BIN:=}"
: "${MYS_FASTQC_BIN:=}"
: "${MYS_PARALLEL_BIN:=}"
: "${MYS_ASCP_BIN:=}"
: "${MYS_BLASTN_BIN:=}"
: "${MYS_BLASTP_BIN:=}"
: "${MYS_MAKEBLASTDB_BIN:=}"
: "${MYS_MKDSSP_BIN:=}"
: "${MYS_PYMOL_BIN:=}"
: "${MYS_CHIMERAX_BIN:=}"
: "${MYS_PLINK_BIN:=}"
: "${MYS_VCFTOOLS_BIN:=}"
: "${MYS_FEATURECOUNTS_BIN:=}"
: "${MYS_SAMTOOLS_BIN:=}"
: "${MYS_BEDTOOLS_BIN:=}"
: "${MYS_BAMCOVERAGE_BIN:=}"
: "${MYS_GTF2BED_BIN:=}"
: "${MYS_GFF2BED_BIN:=}"

# Escape hatch for tools that look up companion binaries by name at runtime:
# directories prepended to PATH for script subprocesses. Normally empty.
: "${MYS_EXTRA_PATH:=}"

# Where tools/doctor.sh --locate searches for candidate binaries.
: "${MYS_LOCATE_DIRS:=$HOME/soft/miniconda3/envs:/opt/anaconda3/envs:/usr/local/bin:$HOME/soft}"

# --- Tool resolution ----------------------------------------------------------
# Fixed order: explicit MYS_<NAME>_BIN path > PATH.
# Usage:
#     samtools="$(mys_resolve_bin MYS_SAMTOOLS_BIN samtools)" || exit 1
mys_resolve_bin() {
    local _var="$1" _tool="$2" _cand
    _cand="${!_var:-}"
    if [ -n "${_cand}" ]; then
        if [ -x "${_cand}" ]; then
            printf '%s\n' "${_cand}"
            return 0
        fi
        printf 'mys_resolve_bin: %s="%s" is not executable (check config/env.local.sh)\n' "${_var}" "${_cand}" >&2
        return 1
    fi
    if command -v -- "${_tool}" >/dev/null 2>&1; then
        command -v -- "${_tool}"
        return 0
    fi
    printf 'mys_resolve_bin: tool "%s" not found — set %s in config/env.local.sh or run tools/doctor.sh --locate %s\n' "${_tool}" "${_var}" "${_tool}" >&2
    return 1
}
export -f mys_resolve_bin 2>/dev/null || true

# --- Machine-specific overrides (gitignored; sourced last, has the last word) --
_local="$(dirname "${BASH_SOURCE[0]}")/env.local.sh"
if [ -f "$_local" ]; then
    . "$_local"
fi
unset -v _local

# R package library: expose MYS_R_LIBS to R itself (prepend, keep any existing).
if [ -n "${MYS_R_LIBS}" ]; then
    export R_LIBS="${MYS_R_LIBS}${R_LIBS:+:${R_LIBS}}"
fi

# PATH escape hatch (see MYS_EXTRA_PATH above).
if [ -n "${MYS_EXTRA_PATH}" ]; then
    export PATH="${MYS_EXTRA_PATH%/}:${PATH}"
fi

# --- Export everything scripts should see -------------------------------------
export MYS_HOME MYS_CONFIG_LOADED MYS_THREADS NCBI_EMAIL NCBI_API_KEY \
    MYS_PYTHON_BIN MYS_RSCRIPT_BIN MYS_R_LIBS \
    MYS_PREFETCH_BIN MYS_FASTERQ_DUMP_BIN MYS_FASTQ_DUMP_BIN MYS_FASTQC_BIN \
    MYS_PARALLEL_BIN MYS_ASCP_BIN MYS_BLASTN_BIN MYS_BLASTP_BIN \
    MYS_MAKEBLASTDB_BIN MYS_MKDSSP_BIN MYS_PYMOL_BIN MYS_CHIMERAX_BIN \
    MYS_PLINK_BIN MYS_VCFTOOLS_BIN MYS_FEATURECOUNTS_BIN MYS_SAMTOOLS_BIN \
    MYS_BEDTOOLS_BIN MYS_BAMCOVERAGE_BIN MYS_GTF2BED_BIN MYS_GFF2BED_BIN \
    MYS_EXTRA_PATH MYS_LOCATE_DIRS
