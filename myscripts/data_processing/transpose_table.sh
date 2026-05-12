#!/bin/bash
#########################################################################
# File Name: transpose_table.sh
# Author: ChengYu
# Description: Transpose a tab-delimited table (rows <-> columns).
# Created Time: 2026
#########################################################################
set -euo pipefail

VERSION="1.0.0"

usage() {
    cat <<EOF
Usage: $(basename "$0") [OPTIONS] -i <INPUT>

Transpose a tab-delimited table (swap rows and columns).

Options:
  -i, --input   <FILE>   Input TSV file (required)
  -o, --output  <FILE>   Output file (default: stdout)
  -s, --sep     <CHAR>   Input field separator (default: tab)
  -h, --help             Show this help
  -v, --version          Show version

Examples:
  $(basename "$0") -i matrix.tsv
  $(basename "$0") -i matrix.tsv -o transposed.tsv
  $(basename "$0") -i data.csv -s ","
EOF
}

log_error() { echo "[ERROR] $1" >&2; exit 1; }

INPUT=""
OUTPUT=""
SEP="\t"

while [[ $# -gt 0 ]]; do
    case $1 in
        -i|--input)   INPUT="$2"; shift 2 ;;
        -o|--output)  OUTPUT="$2"; shift 2 ;;
        -s|--sep)     SEP="$2"; shift 2 ;;
        -h|--help)    usage; exit 0 ;;
        -v|--version) echo "$(basename "$0") ${VERSION}"; exit 0 ;;
        *)            log_error "Unknown option: $1" ;;
    esac
done

[[ -z "${INPUT}" ]] && log_error "Missing required option: -i/--input"
[[ ! -f "${INPUT}" ]] && log_error "Input file not found: ${INPUT}"

# Transpose using awk
RESULT=$(awk -F"${SEP}" -v OFS="${SEP}" '
{
    for (i = 1; i <= NF; i++) {
        col[i] = col[i] (NR > 1 ? OFS : "") $i
    }
}
END {
    for (i = 1; i in col; i++) {
        print col[i]
    }
}' "${INPUT}")

if [[ -n "${OUTPUT}" ]]; then
    echo "${RESULT}" > "${OUTPUT}"
else
    echo "${RESULT}"
fi
