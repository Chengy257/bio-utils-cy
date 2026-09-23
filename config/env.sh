#!/usr/bin/env bash
# =============================================================================
# myscripts — unified environment / dependency configuration
# =============================================================================
# Single place for machine-dependent settings shared by all scripts in bin/.
#
# Usage in bash scripts (convention, adopted script by script):
#     source "$(dirname "${BASH_SOURCE[0]}")/../config/env.sh"
#
# Python / R scripts do NOT source this file. They read the same variables
# from the environment (os.environ / Sys.getenv), so one configuration
# serves every language.
#
# This file holds committed, machine-independent defaults. Private or
# machine-specific values (API keys, email, tool paths) go into
# config/env.local.sh, which is gitignored and sourced at the end of this
# file (it always has the last word).
# =============================================================================

# Guard against double-sourcing.
if [ -n "${MYS_CONFIG_LOADED:-}" ]; then
    return 0 2>/dev/null || true
fi
MYS_CONFIG_LOADED=1

# --- Repository layout -------------------------------------------------------
# Absolute path of this repository, derived from this file's location.
MYS_HOME="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"

# --- Shared defaults (overridable) -------------------------------------------
# Default thread count for scripts that accept -t/--threads.
: "${MYS_THREADS:=8}"

# NCBI E-utilities (already consumed by the fetch/search data-retrieval
# scripts via os.environ / Sys.getenv).
: "${NCBI_EMAIL:=}"
: "${NCBI_API_KEY:=}"

# --- Optional tool slots ------------------------------------------------------
# Most scripts locate external tools on PATH and already provide CLI options
# to override individual binaries (e.g. --ascp-bin, --parallel-bin, -r).
# Add per-tool defaults here as they get adopted, e.g.:
# : "${MYS_ASCP_BIN:=/opt/aspera/connect/bin/ascp}"
# : "${MYS_PREFETCH_BIN:=prefetch}"

# --- PATH additions -----------------------------------------------------------
# Add directories of tools that are not on the system PATH, e.g.:
# export PATH="/opt/aspera/connect/bin:$PATH"

# --- Machine-specific overrides (gitignored; sourced last, has the last word) --
_local="$(dirname "${BASH_SOURCE[0]}")/env.local.sh"
if [ -f "$_local" ]; then
    . "$_local"
fi
unset -v _local

# R package library: expose MYS_R_LIBS to R itself (prepend, keep any existing).
if [ -n "${MYS_R_LIBS:-}" ]; then
    export R_LIBS="${MYS_R_LIBS}${R_LIBS:+:${R_LIBS}}"
fi

# --- Export everything scripts should see -------------------------------------
export MYS_HOME MYS_CONFIG_LOADED MYS_THREADS NCBI_EMAIL NCBI_API_KEY
