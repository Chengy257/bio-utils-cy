#!/usr/bin/env bash
# ==============================================================================
# smoke_test.sh — CLI health check for every script in bin/
#
# Runs each script with -h and verifies it answers with a usage message
# (exit 0, or output containing usage/options/help — R getopt scripts print
# usage and exit non-zero). A script that cannot print its own usage is
# broken; this is the minimum regression bar for any change in bin/.
#
# Usage:
#   tools/smoke_test.sh              # full run (Python + Bash + R)
#   tools/smoke_test.sh --quick      # skip R scripts (slow library loads)
#   tools/smoke_test.sh bin/foo.R    # check a single script
#
# Exit code: 0 if every checked script passes, 1 otherwise.
# ==============================================================================
set -uo pipefail

TOOLS_DIR="$(dirname "$(readlink -f "${BASH_SOURCE[0]}")")"
REPO_ROOT="$(dirname "${TOOLS_DIR}")"

# Pick up MYS_PYTHON_BIN / MYS_RSCRIPT_BIN when configured.
if [ -f "${REPO_ROOT}/config/env.sh" ]; then
    . "${REPO_ROOT}/config/env.sh"
fi
PYTHON_BIN="${MYS_PYTHON_BIN:-$(command -v python3 || true)}"
RSCRIPT_BIN="${MYS_RSCRIPT_BIN:-$(command -v Rscript || true)}"

QUICK=0
TARGETS=()
for arg in "$@"; do
    case "${arg}" in
        --quick) QUICK=1 ;;
        -h|--help) sed -n '2,15p' "${BASH_SOURCE[0]}" | sed 's/^# \{0,1\}//'; exit 0 ;;
        *) TARGETS+=("${arg}") ;;
    esac
done
if [ ${#TARGETS[@]} -eq 0 ]; then
    while IFS= read -r f; do TARGETS+=("$f"); done < <(find "${REPO_ROOT}/bin" -maxdepth 1 -type f | sort)
fi

pass=0; fail=0; skip=0
printf '%-38s %-7s %s\n' "SCRIPT" "STATUS" "NOTE"
printf '%.0s-' {1..80}; echo

for f in "${TARGETS[@]}"; do
    name="$(basename "$f")"
    ext="${name##*.}"
    if [ ! -f "$f" ]; then
        printf '%-38s %-7s %s\n' "$name" "SKIP" "not found"; skip=$((skip+1)); continue
    fi
    case "${ext}" in
        py) [ -n "${PYTHON_BIN}" ] || { printf '%-38s %-7s %s\n' "$name" "SKIP" "no python"; skip=$((skip+1)); continue; }
            out="$("${PYTHON_BIN}" "$f" -h 2>&1)"; rc=$? ;;
        R)  if [ "${QUICK}" -eq 1 ]; then printf '%-38s %-7s %s\n' "$name" "SKIP" "--quick"; skip=$((skip+1)); continue; fi
            if [ -z "${RSCRIPT_BIN}" ]; then printf '%-38s %-7s %s\n' "$name" "SKIP" "no Rscript"; skip=$((skip+1)); continue; fi
            out="$("${RSCRIPT_BIN}" "$f" -h 2>&1)"; rc=$? ;;
        sh) out="$(bash "$f" -h 2>&1)"; rc=$? ;;
        *)  printf '%-38s %-7s %s\n' "$name" "SKIP" "unknown type"; skip=$((skip+1)); continue ;;
    esac
    if [ "${rc}" -eq 0 ] || grep -qiE 'usage|options|help' <<<"${out}"; then
        printf '%-38s %-7s %s\n' "$name" "PASS" ""
        pass=$((pass+1))
    else
        printf '%-38s %-7s %s\n' "$name" "FAIL" "rc=${rc}; first line: $(head -n1 <<<"${out}")"
        fail=$((fail+1))
    fi
done

printf '%.0s-' {1..80}; echo
printf 'RESULT: %d pass, %d fail, %d skip\n' "${pass}" "${fail}" "${skip}"
[ "${fail}" -eq 0 ]
