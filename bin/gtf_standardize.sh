#!/bin/bash
#########################################################################
# File Name: gtf_standardize.sh
# Author: ChengYu
# Description: Standardize GTF files: filter feature types, streamline
#              attributes, backfill gene_biotype, synthesize gene lines
#              when absent, and sort the output.
# Created Time: 2026-03-20
#
# Changelog:
#   v1.1.0  2026-09-26
#   - FIX: three-argument match() is a gawk extension; the awk programs
#     are now POSIX awk (works on gawk, mawk and busybox awk).
#   - FIX: out-of-order GTFs (transcript lines after their exon/CDS
#     lines, or after the gene line) silently lost gene_id/gene_biotype/
#     gene_name backfill; both processing modes now read the file twice
#     (collect pass, then emit pass), making backfill order-independent.
#   - FIX: help promised natural chromosome sorting (-k1,1V) but the
#     code sorted lexicographically (chr10 before chr2); the code now
#     uses LC_ALL=C sort -k1,1V -k4,4n -k5,5n and the docs agree.
#   - FIX: attribute strings no longer start with a stray "; " when
#     gene_id is absent.
#   - log_info/log_verbose go to stderr; temp files are removed via an
#     EXIT trap; validate_output takes its arguments properly instead
#     of relying on dynamic scoping; duplicate rm removed; unknown
#     --attrs keys are rejected with a clear error.
#########################################################################
set -euo pipefail

# --- unified project configuration (paths/defaults; no-op if missing) ---
_my_conf="$(dirname "$(readlink -f "${BASH_SOURCE[0]}")")/../config/env.sh"
if [ -f "${_my_conf}" ]; then . "${_my_conf}"; fi
unset -v _my_conf

# Default configuration
DEFAULT_ATTRS="gene_id,transcript_id,gene_name,exon_number,gene_biotype,protein_id"
KNOWN_ATTRS="gene_id transcript_id gene_name exon_number gene_biotype protein_id"
OUTPUT_SUFFIX=".standardized.gtf"

TMP_FILES=()
cleanup() { rm -f ${TMP_FILES[@]+${TMP_FILES[@]}}; }
trap cleanup EXIT
new_tmp() {
    local t
    t=$(mktemp) || { log_error "mktemp failed"; exit 1; }
    TMP_FILES+=("$t")
    printf '%s\n' "$t"
}

show_help() {
    cat << EOF
================================================================================
                    GTF File Standardization Tool
================================================================================

Usage:
    $0 [OPTIONS] -i <INPUT_FILE> [-o <OUTPUT_FILE>]

Options:
    -h, --help              Show this help information
    -i, --input <FILE>      Input GTF file (required)
    -o, --output <FILE>     Output GTF file (optional, default adds ${OUTPUT_SUFFIX} suffix)
    -a, --attrs <ATTR_LIST> List of attributes to keep, comma-separated
                           Default: ${DEFAULT_ATTRS}
                           Available: gene_id, transcript_id, gene_name, exon_number,
                                      gene_biotype, protein_id
    -v, --verbose           Show detailed processing information
    -f, --force             Force overwrite existing output file

Description:
    This tool standardizes GTF files with the following features:
    1. Keep only gene, transcript, exon, CDS, five_prime_utr, three_prime_utr,
       start_codon, stop_codon feature types
    2. Streamline attributes, keep only specified key attributes
    3. Ensure every line contains gene_biotype information
    4. Preserve phase column (column 8) for CDS features (0, 1, 2)
    5. Automatically handle two formats:
       - Files with gene lines (e.g., RefSeq): preserve original gene lines
       - Files without gene lines (e.g., RNAcentral): automatically generate gene lines
    6. For files without gene lines, calculate coordinate range for each gene and generate gene lines
    7. Output is sorted by genomic coordinates (chromosome, start, end) using
       LC_ALL=C sort -k1,1V -k4,4n -k5,5n (natural chromosome order)

Examples:
    # Basic usage (use default configuration)
    $0 -i input.gtf

    # Specify output file
    $0 -i input.gtf -o output.gtf

    # Customize attributes to keep
    $0 -i input.gtf -a "gene_id,transcript_id,gene_biotype"

    # Show detailed processing information
    $0 -i input.gtf -v

Output Format:
    Attribute order: gene_id -> gene_name -> transcript_id -> exon_number -> protein_id -> gene_biotype

    Sorting: Output is sorted by genomic coordinates (chromosome, start, end)
             Using: LC_ALL=C sort -k1,1V -k4,4n -k5,5n

    Supported feature types and attributes:
        gene:           gene_id, gene_name, gene_biotype
        transcript:     gene_id, transcript_id, gene_name, gene_biotype
        exon:           gene_id, transcript_id, exon_number, gene_name, gene_biotype
        CDS:            gene_id, transcript_id, exon_number, gene_name, gene_biotype, protein_id
        five_prime_utr: gene_id, transcript_id, exon_number, gene_name, gene_biotype
        three_prime_utr: gene_id, transcript_id, exon_number, gene_name, gene_biotype
        start_codon:    gene_id, transcript_id, exon_number, gene_name, gene_biotype, protein_id
        stop_codon:     gene_id, transcript_id, exon_number, gene_name, gene_biotype, protein_id

    Phase column (column 8):
        CDS features:    phase preserved from input (0, 1, or 2)
        Codon features:  phase set to 0
        Other features:  phase set to "."

    RefSeq format example (with gene lines):
        gene:       gene_id "XXX"; gene_biotype "lncRNA"
        transcript: gene_id "XXX"; transcript_id "XXX"; gene_biotype "lncRNA"
        exon:       gene_id "XXX"; transcript_id "XXX"; exon_number "N"; gene_biotype "lncRNA"
        CDS:        gene_id "XXX"; transcript_id "XXX"; exon_number "N"; gene_biotype "protein_coding"; protein_id "XXX"

    RNAcentral format example (auto-generate gene lines):
        gene:       gene_id "XXX"; gene_biotype "lncRNA"      (auto-generated)
        transcript: gene_id "XXX"; transcript_id "XXX"; gene_biotype "lncRNA"
        exon:       gene_id "XXX"; transcript_id "XXX"; exon_number "N"; gene_biotype "lncRNA"

    Note: Auto-generated gene lines are calculated from coordinate ranges of all transcripts

================================================================================
EOF
}

# Logging functions (all to stderr; stdout is reserved for output files)
log_info() {
    echo "[INFO] $1" >&2
}

log_warn() {
    echo "[WARN] $1" >&2
}

log_error() {
    echo "[ERROR] $1" >&2
}

log_verbose() {
    if [[ "${VERBOSE}" == "true" ]]; then
        echo "[VERBOSE] $1" >&2
    fi
}

# Check if file exists
check_file() {
    local file=$1
    if [[ ! -f "${file}" ]]; then
        log_error "File does not exist: ${file}"
        exit 1
    fi
}

# Check if overwrite is allowed
check_overwrite() {
    local file=$1
    if [[ -f "${file}" && "${FORCE}" != "true" ]]; then
        log_error "Output file already exists: ${file}"
        log_error "Use -f/--force option to overwrite"
        exit 1
    fi
}

# Check dependencies
check_dependencies() {
    local deps=("awk" "grep" "sort")
    for dep in "${deps[@]}"; do
        if ! command -v "${dep}" &> /dev/null; then
            log_error "Missing dependency: ${dep}"
            exit 1
        fi
    done
}

# Detect if file has gene lines
detect_gene_lines() {
    local input_file=$1
    local gene_count
    gene_count=$(grep -c $'\tgene\t' "${input_file}" 2>/dev/null || true)
    if [[ "${gene_count:-0}" -gt 0 ]]; then
        echo "with_gene"
    else
        echo "without_gene"
    fi
}

# ---------------------------------------------------------------------------
# Shared awk library (POSIX awk) prepended to every program below.
#
#   attr_value   : first quoted value of a raw attribute piece, or the
#                  bare token after the key
#   parse_attrs  : split column 9 into key->value for the six known keys
#                  (exact key match; quoted or bare values)
#   add_attr     : append 'key "value"' with proper "; " separation
#   emit_attrs   : build the standardized attribute string for a feature
#   init_attrs   : parse the comma-separated --attrs list
# ---------------------------------------------------------------------------
read -r -d '' AWK_LIB <<'AWKLIB' || true
function attr_value(piece) {
    if (match(piece, /"[^"]*"/)) {
        return substr(piece, RSTART + 1, RLENGTH - 2)
    }
    sub(/^[^ \t]+[ \t]+/, "", piece)   # drop the key for bare values
    sub(/[ \t]+$/, "", piece)
    return piece
}

function parse_attrs(    i, n, parts, attr, key) {
    gene_id = ""; transcript_id = ""; gene_name = ""
    gene_biotype = ""; exon_number = ""; protein_id = ""
    n = split($9, parts, ";")
    for (i = 1; i <= n; i++) {
        attr = parts[i]
        gsub(/^[ \t]+|[ \t]+$/, "", attr)
        if (attr == "") continue
        key = attr
        sub(/[ \t].*$/, "", key)
        if (key == "gene_id")                         gene_id = attr_value(attr)
        else if (key == "transcript_id")              transcript_id = attr_value(attr)
        else if (key == "gene_name" || key == "gene") gene_name = attr_value(attr)
        else if (key == "gene_biotype")               gene_biotype = attr_value(attr)
        else if (key == "exon_number")                exon_number = attr_value(attr)
        else if (key == "protein_id")                 protein_id = attr_value(attr)
    }
}

function add_attr(result, key, value) {
    if (value == "") return result
    if (result == "") return key " \"" value "\""
    return result "; " key " \"" value "\""
}

function emit_attrs(ftype,    res) {
    res = ""
    if (wanted_attrs["gene_id"])        res = add_attr(res, "gene_id", gene_id)
    if (wanted_attrs["gene_name"])      res = add_attr(res, "gene_name", gene_name)
    if (wanted_attrs["transcript_id"] && ftype != "gene") {
        res = add_attr(res, "transcript_id", transcript_id)
    }
    if (wanted_attrs["exon_number"] && ftype != "gene" && ftype != "transcript") {
        res = add_attr(res, "exon_number", exon_number)
    }
    if (wanted_attrs["protein_id"] && (ftype == "CDS" || ftype == "start_codon" || ftype == "stop_codon")) {
        res = add_attr(res, "protein_id", protein_id)
    }
    if (wanted_attrs["gene_biotype"])   res = add_attr(res, "gene_biotype", gene_biotype)
    return res
}

function init_attrs(attrs,    i, n, list) {
    split("", wanted_attrs)   # POSIX idiom: clear the array
    n = split(attrs, list, ",")
    for (i = 1; i <= n; i++) wanted_attrs[list[i]] = 1
}
AWKLIB

# Sort: natural chromosome order, then start/end (locale-independent).
sort_output() {
    local src=$1 dst=$2
    LC_ALL=C sort -k1,1V -k4,4n -k5,5n "${src}" > "${dst}"
}

# Process GTF files with gene lines (e.g., RefSeq).
# Pass 1 collects gene-level biotype/name; pass 2 emits. Backfill is
# therefore independent of line order.
process_with_gene() {
    local input_file=$1
    local output_file=$2
    local attrs=$3

    local temp_unsorted prog
    temp_unsorted=$(new_tmp)
    prog=$(printf '%s\n%s\n' "${AWK_LIB}" 'BEGIN { init_attrs(attrs) }
FNR == NR {
    if ($3 != "gene") next
    parse_attrs()
    if (gene_id == "") next
    if (gene_biotype != "") gene_biotypes[gene_id] = gene_biotype
    if (gene_name != "" && gene_name != gene_id) gene_names[gene_id] = gene_name
    next
}
{
    if ($3 !~ /^(gene|transcript|exon|CDS|five_prime_utr|three_prime_utr|start_codon|stop_codon)$/) next
    parse_attrs()

    if (gene_biotype == "" && gene_id in gene_biotypes) gene_biotype = gene_biotypes[gene_id]
    if (gene_name == "" && gene_id in gene_names) gene_name = gene_names[gene_id]

    print $1, $2, $3, $4, $5, $6, $7, $8, emit_attrs($3)
}')

    awk -F'\t' -v attrs="${attrs}" -v OFS='\t' "${prog}" \
        "${input_file}" "${input_file}" > "${temp_unsorted}"

    sort_output "${temp_unsorted}" "${output_file}"
}

# Process GTF files without gene lines (e.g., RNAcentral): synthesize one
# gene line per gene from the extent of its transcript lines.
process_without_gene() {
    local input_file=$1
    local output_file=$2
    local attrs=$3

    local temp_unsorted prog
    temp_unsorted=$(new_tmp)
    prog=$(printf '%s\n%s\n' "${AWK_LIB}" 'BEGIN { init_attrs(attrs) }
FNR == NR {
    if ($3 != "transcript") next
    parse_attrs()
    if (gene_id == "") next
    if (!(gene_id in gene_start)) {
        gene_chrom[gene_id] = $1
        gene_source[gene_id] = $2
        gene_start[gene_id] = $4 + 0
        gene_end[gene_id] = $5 + 0
        gene_strand[gene_id] = $7
    } else {
        if ($4 + 0 < gene_start[gene_id]) gene_start[gene_id] = $4 + 0
        if ($5 + 0 > gene_end[gene_id]) gene_end[gene_id] = $5 + 0
    }
    if (gene_biotype != "" && !(gene_id in gene_biotypes)) gene_biotypes[gene_id] = gene_biotype
    if (gene_name != "" && gene_name != gene_id && !(gene_id in gene_names)) gene_names[gene_id] = gene_name
    if (transcript_id != "") transcript_to_gene[transcript_id] = gene_id
    next
}
{
    if ($3 !~ /^(transcript|exon|CDS|five_prime_utr|three_prime_utr|start_codon|stop_codon)$/) next
    parse_attrs()

    # gene line for this feature gene; only synthesizable when pass 1
    # collected transcript-derived coordinates for it
    g_id = (gene_id != "") ? gene_id : transcript_to_gene[transcript_id]
    if (g_id != "" && (g_id in gene_chrom) && !(g_id in gene_emitted)) {
        gene_emitted[g_id] = 1
        g_biotype = (g_id in gene_biotypes) ? gene_biotypes[g_id] : "unknown"
        g_name = (g_id in gene_names) ? gene_names[g_id] : ""
        gene_id = g_id; gene_biotype = g_biotype; gene_name = g_name
        print gene_chrom[g_id], gene_source[g_id], "gene", gene_start[g_id], gene_end[g_id], ".", gene_strand[g_id], ".", emit_attrs("gene")
        parse_attrs()   # restore the per-line values clobbered above
    }
    if (gene_id == "" && transcript_id != "" && transcript_id in transcript_to_gene) {
        gene_id = transcript_to_gene[transcript_id]
    }
    if (gene_biotype == "" && gene_id in gene_biotypes) gene_biotype = gene_biotypes[gene_id]
    if (gene_name == "" && gene_id in gene_names) gene_name = gene_names[gene_id]

    print $1, $2, $3, $4, $5, $6, $7, $8, emit_attrs($3)
}')

    awk -F'\t' -v attrs="${attrs}" -v OFS='\t' "${prog}" \
        "${input_file}" "${input_file}" > "${temp_unsorted}"

    sort_output "${temp_unsorted}" "${output_file}"
}

# Validate output file
validate_output() {
    local input_file=$1
    local output_file=$2
    local input_lines output_lines with_biotype missing_biotype

    input_lines=$(wc -l < "${input_file}")
    output_lines=$(wc -l < "${output_file}")

    log_verbose "Input file lines: ${input_lines}"
    log_verbose "Output file lines: ${output_lines}"

    with_biotype=$(grep -c "gene_biotype" "${output_file}" || true)
    missing_biotype=$((output_lines - with_biotype))

    if [[ "${missing_biotype}" -gt 0 ]]; then
        log_warn "${missing_biotype} lines missing gene_biotype"
    else
        log_info "All lines contain gene_biotype"
    fi

    log_info "Processing complete!"
}

# Main function
main() {
    local input_file=""
    local output_file=""
    local attrs="${DEFAULT_ATTRS}"
    local VERBOSE="false"
    local FORCE="false"

    # Parse command line arguments
    while [[ $# -gt 0 ]]; do
        case $1 in
            -h|--help)
                show_help
                exit 0
                ;;
            -i|--input)
                input_file="$2"
                shift 2
                ;;
            -o|--output)
                output_file="$2"
                shift 2
                ;;
            -a|--attrs)
                attrs="$2"
                shift 2
                ;;
            -v|--verbose)
                VERBOSE="true"
                shift
                ;;
            -f|--force)
                FORCE="true"
                shift
                ;;
            *)
                log_error "Unknown option: $1"
                echo "Use -h/--help to see help information" >&2
                exit 1
                ;;
        esac
    done

    # Check required parameters
    if [[ -z "${input_file}" ]]; then
        log_error "Missing required parameter: -i/--input"
        echo "Use -h/--help to see help information" >&2
        exit 1
    fi

    # Reject unknown attribute keys early (they would be silently dropped)
    IFS=',' read -r -a attr_items <<< "${attrs}"
    for item in "${attr_items[@]}"; do
        valid=false
        for known in ${KNOWN_ATTRS}; do
            if [[ "${item}" == "${known}" ]]; then valid=true; break; fi
        done
        if [[ "${valid}" != "true" ]]; then
            log_error "Unknown attribute in --attrs: ${item}. Valid: ${KNOWN_ATTRS// /, }"
            exit 1
        fi
    done

    # Check dependencies
    check_dependencies

    # Check input file
    check_file "${input_file}"

    # Set output file
    if [[ -z "${output_file}" ]]; then
        output_file="${input_file}${OUTPUT_SUFFIX}"
    fi

    # Check output file
    check_overwrite "${output_file}"

    # Display processing information
    log_info "Starting GTF file processing..."
    log_verbose "Input file: ${input_file}"
    log_verbose "Output file: ${output_file}"
    log_verbose "Attributes to keep: ${attrs}"

    # Detect file type and process
    local file_type
    file_type=$(detect_gene_lines "${input_file}")
    log_verbose "Detected file type: ${file_type}"

    if [[ "${file_type}" == "with_gene" ]]; then
        log_verbose "Using processing mode for files with gene lines"
        process_with_gene "${input_file}" "${output_file}" "${attrs}"
    else
        log_verbose "Using processing mode for files without gene lines"
        process_without_gene "${input_file}" "${output_file}" "${attrs}"
    fi

    # Validate output
    validate_output "${input_file}" "${output_file}"
}

# Execute main function
main "$@"
