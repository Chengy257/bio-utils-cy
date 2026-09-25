#!/bin/bash
#########################################################################
# File Name: transpose_table.sh
# Author: ChengYu
# Description: Transpose a tab-delimited table (rows <-> columns).
# Created Time: 2026
#
# Changelog:
#   v1.1.0  2026-09-25
#   - FIX: rows with inconsistent field counts now fail with an error
#     naming the offending line; previously they silently shifted every
#     following column.
#   - FIX: output is streamed from awk directly (stdout or via a temp file
#     + atomic rename); the result no longer round-trips through a shell
#     variable, which doubled memory, dropped trailing empty fields/lines,
#     and broke on content starting with "-e"/"-n".
#   - FIX: missing option values no longer trip `set -u` with an obscure
#     error.
#########################################################################
set -euo pipefail

# --- unified project configuration (paths/defaults; no-op if missing) ---
_my_conf="$(dirname "$(readlink -f "${BASH_SOURCE[0]}")")/../config/env.sh"
if [ -f "${_my_conf}" ]; then . "${_my_conf}"; fi
unset -v _my_conf

VERSION="1.1.0"

usage() {
    cat <<EOF
Usage: $(basename "$0") [OPTIONS] -i <INPUT>

Transpose a tab-delimited table (swap rows and columns).

Options:
  -i, --input   <FILE>   Input TSV file (required)
  -o, --output  <FILE>   Output file (default: stdout)
  -s, --sep     <CHAR>   Input field separator (default: tab; must be a
                         single character)
  -h, --help             Show this help
  -v, --version          Show version

Examples:
  $(basename "$0") -i matrix.tsv
  $(basename "$0") -i matrix.tsv -o transposed.tsv
  $(basename "$0") -i data.csv -s ","

All rows must have the same number of fields; mismatched rows are an error.
EOF
}

log_error() { echo "[ERROR] $1" >&2; exit 1; }

INPUT=""
OUTPUT=""
SEP="\t"

while [[ $# -gt 0 ]]; do
    case $1 in
        -i|--input)   [[ $# -ge 2 ]] || log_error "Missing value for $1"; INPUT="$2"; shift 2 ;;
        -o|--output)  [[ $# -ge 2 ]] || log_error "Missing value for $1"; OUTPUT="$2"; shift 2 ;;
        -s|--sep)     [[ $# -ge 2 ]] || log_error "Missing value for $1"; SEP="$2"; shift 2 ;;
        -h|--help)    usage; exit 0 ;;
        -v|--version) echo "$(basename "$0") ${VERSION}"; exit 0 ;;
        *)            log_error "Unknown option: $1" ;;
    esac
done

[[ -z "${INPUT}" ]] && log_error "Missing required option: -i/--input"
[[ ! -f "${INPUT}" ]] && log_error "Input file not found: ${INPUT}"

AWK_PROG='
NR == 1 { expected = NF }
{
    if (NF != expected) {
        printf "[ERROR] transpose_table: line %d has %d field(s), expected %d (all rows must match)\n", NR, NF, expected > "/dev/stderr"
        exit 1
    }
    for (i = 1; i <= NF; i++) {
        col[i] = col[i] (NR > 1 ? OFS : "") $i
    }
}
END {
    for (i = 1; i in col; i++) {
        print col[i]
    }
}'

if [[ -n "${OUTPUT}" ]]; then
    TMP_OUT="$(mktemp)"
    trap 'rm -f "${TMP_OUT}"' EXIT
    awk -F"${SEP}" -v OFS="${SEP}" "${AWK_PROG}" "${INPUT}" > "${TMP_OUT}"
    mv "${TMP_OUT}" "${OUTPUT}"
    trap - EXIT
else
    awk -F"${SEP}" -v OFS="${SEP}" "${AWK_PROG}" "${INPUT}"
fi
