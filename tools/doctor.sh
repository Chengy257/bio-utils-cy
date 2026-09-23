#!/usr/bin/env bash
# ==============================================================================
# doctor.sh — dependency configuration checker for myscripts
#
# Validates the tool paths declared in config/env.local.sh against what the
# scripts in bin/ actually need, and helps fill the gaps:
#
#   tools/doctor.sh                    full check (tools + Python + R packages)
#   tools/doctor.sh --locate ascp      every candidate path + ready config line
#   tools/doctor.sh --migrate-check    only items needing attention (migration)
#
# Statuses:
#   OK      declared path exists and is executable
#   PATH    not declared but found on PATH — works; consider pinning a path
#   BROKEN  declared but not a runnable file  -> fix the path in env.local.sh
#   UNSET   neither declared nor on PATH      -> fill it in config/env.local.sh
#
# Exit codes: 0 = everything OK; 1 = at least one item needs action.
# ==============================================================================
set -uo pipefail

TOOLS_DIR="$(dirname "$(readlink -f "${BASH_SOURCE[0]}")")"
REPO_ROOT="$(dirname "${TOOLS_DIR}")"

if [ -f "${REPO_ROOT}/config/env.sh" ]; then
    . "${REPO_ROOT}/config/env.sh"
else
    echo "ERROR: ${REPO_ROOT}/config/env.sh not found" >&2
    exit 1
fi

# --- Registry: VAR|command used by scripts ------------------------------------
REGISTRY=(
    "MYS_PYTHON_BIN|python3"
    "MYS_RSCRIPT_BIN|Rscript"
    "MYS_PREFETCH_BIN|prefetch"
    "MYS_FASTERQ_DUMP_BIN|fasterq-dump"
    "MYS_FASTQ_DUMP_BIN|fastq-dump"
    "MYS_FASTQC_BIN|fastqc"
    "MYS_PARALLEL_BIN|parallel"
    "MYS_ASCP_BIN|ascp"
    "MYS_BLASTN_BIN|blastn"
    "MYS_BLASTP_BIN|blastp"
    "MYS_MAKEBLASTDB_BIN|makeblastdb"
    "MYS_MKDSSP_BIN|mkdssp"
    "MYS_PYMOL_BIN|pymol"
    "MYS_CHIMERAX_BIN|chimerax"
    "MYS_PLINK_BIN|plink"
    "MYS_VCFTOOLS_BIN|vcftools"
    "MYS_FEATURECOUNTS_BIN|featureCounts"
    "MYS_SAMTOOLS_BIN|samtools"
    "MYS_INFER_EXP_BIN|infer_experiment.py"
    "MYS_BEDTOOLS_BIN|bedtools"
    "MYS_BAMCOVERAGE_BIN|bamCoverage"
    "MYS_GTF2BED_BIN|gtf2bed"
    "MYS_GFF2BED_BIN|gff2bed"
)

# Python packages imported by bin/*.py (import name)
PY_PKGS=(Bio pandas numpy matplotlib networkx scipy sklearn community tqdm lxml pyteomics pybedtools markov_clustering)

# R packages loaded by bin/*.R
R_PKGS=(getopt ggplot2 Biostrings sangerseqR GenomicFeatures edgeR DESeq2
        BiocParallel gplots RColorBrewer amap aPEAR clusterProfiler magrittr
        plotrix dplyr rtracklayer ape aplot ggtree GOSemSim cowplot data.table
        GenomicRanges Gviz qqman circlize ComplexHeatmap UpSetR ggsci patchwork
        Mfuzz Biobase)

tool_status() {  # $1=var $2=tool -> echoes "STATUS|resolved_path"
    local var="$1" tool="$2" cand
    cand="${!var:-}"
    if [ -n "${cand}" ]; then
        if [ -x "${cand}" ]; then echo "OK|${cand}"; else echo "BROKEN|${cand}"; fi
        return
    fi
    if command -v "${tool}" >/dev/null 2>&1; then
        echo "PATH|$(command -v "${tool}")"
    else
        echo "UNSET|"
    fi
}

locate_tool() {  # $1=tool name
    local tool="$1" d hit found=0
    echo "Candidates for '${tool}' (from MYS_LOCATE_DIRS):"
    for d in ${MYS_LOCATE_DIRS//:/ }; do
        d="${d/#\~/$HOME}"
        for hit in "${d}/bin/${tool}" "${d}/${tool}" "${d}"/*/bin/"${tool}"; do
            if [ -x "${hit}" ]; then
                echo "  ${hit}"
                found=1
            fi
        done
    done
    [ "${found}" -eq 1 ] || echo "  (none found)"
    local var
    for entry in "${REGISTRY[@]}"; do
        if [ "${entry#*|}" = "${tool}" ]; then
            var="${entry%%|*}"
            echo "Config line:  ${var}=/path/to/${tool}"
            return
        fi
    done
}

MIGRATE=0
case "${1:-}" in
    --locate)
        [ -n "${2:-}" ] || { echo "Usage: $0 --locate <tool>" >&2; exit 1; }
        locate_tool "${2}"
        exit 0
        ;;
    --migrate-check) MIGRATE=1 ;;
    -h|--help) sed -n '2,26p' "${BASH_SOURCE[0]}" | sed 's/^# \{0,1\}//'; exit 0 ;;
esac

# --- Tool path checks ---------------------------------------------------------
need_action=0
printf '%-8s %-16s %-24s %s\n' "STATUS" "TOOL" "VARIABLE" "PATH / HINT"
printf '%.0s-' {1..90}; echo
for entry in "${REGISTRY[@]}"; do
    var="${entry%%|*}"; tool="${entry#*|}"
    IFS='|' read -r status path <<< "$(tool_status "${var}" "${tool}")"
    case "${status}" in
        OK)     hint="${path}" ;;
        PATH)   hint="${path}  (unpinned — consider setting ${var})" ;;
        BROKEN) hint="${path}  (not executable — fix ${var})"; need_action=1 ;;
        UNSET)  hint="(not declared, not on PATH — run: $0 --locate ${tool})"; need_action=1 ;;
    esac
    if [ "${MIGRATE}" -eq 1 ] && [ "${status}" = "OK" ]; then continue; fi
    printf '%-8s %-16s %-24s %s\n' "${status}" "${tool}" "${var}" "${hint}"
done

# --- Python package checks -----------------------------------------------------
py_missing=()
if [ "${MIGRATE}" -eq 0 ] || [ "${need_action}" -eq 0 ]; then
    echo
    echo "Python packages (interpreter: ${MYS_PYTHON_BIN:-$(command -v python3 2>/dev/null) or PATH}):"
    if command -v "${MYS_PYTHON_BIN:-python3}" >/dev/null 2>&1; then
        for pkg in "${PY_PKGS[@]}"; do
            if ! "${MYS_PYTHON_BIN:-python3}" -c "import ${pkg}" >/dev/null 2>&1; then
                py_missing+=("${pkg}")
            fi
        done
        if [ ${#py_missing[@]} -eq 0 ]; then
            echo "  all ${#PY_PKGS[@]} packages importable"
        else
            echo "  MISSING: ${py_missing[*]}"
            echo "  -> install them into that interpreter, or point MYS_PYTHON_BIN at one that has them"
            need_action=1
        fi
    else
        echo "  interpreter not found — set MYS_PYTHON_BIN"
        need_action=1
    fi
fi

# --- R package checks ----------------------------------------------------------
r_missing=()
if [ -n "${MYS_RSCRIPT_BIN}" ] && [ -x "${MYS_RSCRIPT_BIN}" ]; then
    echo
    echo "R packages (Rscript: ${MYS_RSCRIPT_BIN}; R_LIBS: ${R_LIBS:-unset}):"
    pkg_list="$(IFS=,; echo "${R_PKGS[*]}")"
    while read -r pkg ok; do
        [ "${ok}" = "TRUE" ] || r_missing+=("${pkg}")
    done < <("${MYS_RSCRIPT_BIN}" --vanilla -e \
        "for (p in strsplit('${pkg_list}',',')[[1]]) cat(p, requireNamespace(p, quietly=TRUE), '\n')" 2>/dev/null)
    if [ ${#r_missing[@]} -eq 0 ]; then
        echo "  all ${#R_PKGS[@]} packages resolvable"
    else
        echo "  MISSING: ${r_missing[*]}"
        echo "  -> install into ${MYS_R_LIBS:-the R library} or adjust MYS_R_LIBS / MYS_RSCRIPT_BIN"
        need_action=1
    fi
else
    echo
    echo "R packages: SKIPPED (MYS_RSCRIPT_BIN not set or not executable)"
    need_action=1
fi

echo
if [ "${need_action}" -eq 0 ]; then
    echo "DOCTOR: all dependencies resolved."
    exit 0
else
    echo "DOCTOR: action needed — see items above (or: $0 --migrate-check)."
    exit 1
fi
