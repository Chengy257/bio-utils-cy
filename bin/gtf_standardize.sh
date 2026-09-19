#!/bin/bash

##############################################################################
# GTF File Standardization Tool
#
# Description: Streamline GTF file attributes, keep specified key attributes,
#              ensure every line contains gene_biotype
#
# Author: Auto-generated
# Date: 2026-03-20
##############################################################################

set -euo pipefail

# Default configuration
DEFAULT_ATTRS="gene_id,transcript_id,gene_name,exon_number,gene_biotype,protein_id"
OUTPUT_SUFFIX=".standardized.gtf"

# Display help information
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
    7. Output is sorted by genomic coordinates (chromosome, start, end) using sort -k1,1V -k4,4n -k5,5n

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
             Using: sort -k1,1 -k4,4n -k5,5n

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

# Logging functions
log_info() {
    echo "[INFO] $1"
}

log_warn() {
    echo "[WARN] $1" >&2
}

log_error() {
    echo "[ERROR] $1" >&2
}

log_verbose() {
    if [[ "${VERBOSE}" == "true" ]]; then
        echo "[VERBOSE] $1"
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
    local deps=("awk" "grep" "wc")
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
    local gene_count=$(grep -c $'\tgene\t' "${input_file}" 2>/dev/null || echo "0")
    gene_count=$(echo "${gene_count}" | tr -d '\n\r' | tr -d ' ')
    if [[ "${gene_count}" -gt 0 ]]; then
        echo "with_gene"
    else
        echo "without_gene"
    fi
}

# Process GTF files with gene lines (e.g., RefSeq)
process_with_gene() {
    local input_file=$1
    local output_file=$2
    local attrs=$3

    local temp_unsorted=$(mktemp)
    awk -F'\t' -v attrs="${attrs}" '
    BEGIN {
        OFS="\t"
        # Parse attribute list
        n = split(attrs, attr_list, ",")
        for (i = 1; i <= n; i++) {
            wanted_attrs[attr_list[i]] = 1
        }
    }
    {
        if ($3 !~ /^(gene|transcript|exon|CDS|five_prime_utr|three_prime_utr|start_codon|stop_codon)$/) next

        attrs = $9
        gene_id = ""
        transcript_id = ""
        gene_name = ""
        gene_biotype = ""
        exon_number = ""
        protein_id = ""

        # Parse attributes
        n = split(attrs, parts, ";")
        for (i = 1; i <= n; i++) {
            attr = parts[i]
            gsub(/^[ \t]+|[ \t]+$/, "", attr)

            if (attr ~ /^gene_id[ \t]+"/) {
                match(attr, /gene_id[ \t]+"([^"]+)"/, m)
                if (m[1] != "") gene_id = m[1]
            } else if (attr ~ /^transcript_id[ \t]+"/ && attr !~ /^transcript_id[ \t]+""/) {
                match(attr, /transcript_id[ \t]+"([^"]+)"/, m)
                if (m[1] != "") transcript_id = m[1]
            } else if (attr ~ /^gene[ \t]+"/) {
                match(attr, /gene[ \t]+"([^"]+)"/, m)
                if (m[1] != "") gene_name = m[1]
            } else if (attr ~ /^gene_biotype[ \t]+"/) {
                match(attr, /gene_biotype[ \t]+"([^"]+)"/, m)
                if (m[1] != "") gene_biotype = m[1]
            } else if (attr ~ /^exon_number[ \t]+/) {
                if (attr ~ /^exon_number[ \t]+"/) {
                    match(attr, /exon_number[ \t]+"([^"]+)"/, m)
                    if (m[1] != "") exon_number = m[1]
                } else {
                    match(attr, /exon_number[ \t]+([^; \t]+)/, m)
                    if (m[1] != "") exon_number = m[1]
                }
            } else if (attr ~ /^protein_id[ \t]+"/) {
                match(attr, /protein_id[ \t]+"([^"]+)"/, m)
                if (m[1] != "") protein_id = m[1]
            }
        }

        # Store gene biotype and name
        if ($3 == "gene" && gene_id != "") {
            if (gene_biotype != "") gene_biotypes[gene_id] = gene_biotype
            if (gene_name != "" && gene_name != gene_id) gene_names[gene_id] = gene_name
        }

        # Get gene biotype
        if (gene_id != "" && gene_id in gene_biotypes && gene_biotype == "") {
            gene_biotype = gene_biotypes[gene_id]
        }

        # Get gene name
        if (gene_name == "" && gene_id in gene_names) {
            gene_name = gene_names[gene_id]
        }

        # Build new attribute string
        new_attrs = ""
        sep = ""

        if (wanted_attrs["gene_id"] && gene_id != "") {
            new_attrs = sep "gene_id \"" gene_id "\""
            sep = "; "
        }
        if (wanted_attrs["gene_name"] && gene_name != "" && gene_name != gene_id) {
            new_attrs = new_attrs "; gene_name \"" gene_name "\""
        }
        if (wanted_attrs["transcript_id"] && ($3 == "transcript" || $3 == "exon" || $3 == "CDS" || $3 == "five_prime_utr" || $3 == "three_prime_utr" || $3 == "start_codon" || $3 == "stop_codon") && transcript_id != "") {
            new_attrs = new_attrs "; transcript_id \"" transcript_id "\""
        }
        if (wanted_attrs["exon_number"] && ($3 == "exon" || $3 == "CDS" || $3 == "five_prime_utr" || $3 == "three_prime_utr" || $3 == "start_codon" || $3 == "stop_codon") && exon_number != "") {
            new_attrs = new_attrs "; exon_number \"" exon_number "\""
        }
        if (wanted_attrs["protein_id"] && ($3 == "CDS" || $3 == "start_codon" || $3 == "stop_codon") && protein_id != "") {
            new_attrs = new_attrs "; protein_id \"" protein_id "\""
        }
        if (wanted_attrs["gene_biotype"] && gene_biotype != "") {
            new_attrs = new_attrs "; gene_biotype \"" gene_biotype "\""
        }

        print $1,$2,$3,$4,$5,$6,$7,$8,new_attrs
    }' "${input_file}" > "${temp_unsorted}"

    # Sort by genomic coordinates
    sort -k1,1 -k4,4n -k5,5n "${temp_unsorted}" > "${output_file}"

    # Clean up temporary file
    rm -f "${temp_unsorted}"
}

# Process GTF files without gene lines (e.g., RNAcentral)
process_without_gene() {
    local input_file=$1
    local output_file=$2
    local attrs=$3

    # Step 1: Collect gene information and save to temporary file
    local temp_gene_file=$(mktemp)
    awk -F'\t' '
    BEGIN { OFS="\t" }
    {
        if ($3 != "transcript") next

        # Parse attributes
        attrs_str = $9
        gene_id = ""
        gene_biotype = ""
        gene_name = ""

        n = split(attrs_str, parts, ";")
        for (i = 1; i <= n; i++) {
            attr = parts[i]
            gsub(/^[ \t]+|[ \t]+$/, "", attr)

            if (attr ~ /^gene_id[ \t]+"/) {
                match(attr, /gene_id[ \t]+"([^"]+)"/, m)
                if (m[1] != "") gene_id = m[1]
            } else if (attr ~ /^gene_biotype[ \t]+"/) {
                match(attr, /gene_biotype[ \t]+"([^"]+)"/, m)
                if (m[1] != "") gene_biotype = m[1]
            } else if (attr ~ /^gene[ \t]+"/) {
                match(attr, /gene[ \t]+"([^"]+)"/, m)
                if (m[1] != "") gene_name = m[1]
            }
        }

        # Collect gene information
        if (gene_id != "") {
            # Check if gene already exists
            if (!(gene_id in gene_start)) {
                gene_chrom[gene_id] = $1
                gene_source[gene_id] = $2
                gene_start[gene_id] = $4 + 0
                gene_end[gene_id] = $5 + 0
                gene_strand[gene_id] = $7
                gene_biotype_val[gene_id] = gene_biotype
                if (gene_name != "" && gene_name != gene_id) {
                    gene_name_val[gene_id] = gene_name
                }
            } else {
                # Update coordinate range
                if ($4 + 0 < gene_start[gene_id]) gene_start[gene_id] = $4 + 0
                if ($5 + 0 > gene_end[gene_id]) gene_end[gene_id] = $5 + 0
                # Add gene_name if available and not present
                if (gene_name != "" && gene_name != gene_id && !(gene_id in gene_name_val)) {
                    gene_name_val[gene_id] = gene_name
                }
            }
        }
    }
    END {
        for (gene_id in gene_chrom) {
            biotype = (gene_id in gene_biotype_val && gene_biotype_val[gene_id] != "") ? gene_biotype_val[gene_id] : "unknown"
            gene_name = (gene_id in gene_name_val) ? gene_name_val[gene_id] : ""
            print gene_chrom[gene_id], gene_source[gene_id], gene_start[gene_id], gene_end[gene_id], gene_strand[gene_id], gene_id, biotype, gene_name
        }
    }' "${input_file}" > "${temp_gene_file}"

    # Step 2: Output gene lines, then process transcript/exon lines, then sort
    local temp_unsorted=$(mktemp)
    {
        # First output all gene lines (read from temp file and format)
        awk -F'\t' -v attrs="${attrs}" '
        BEGIN {
            OFS="\t"
            # Parse attribute list
            n = split(attrs, attr_list, ",")
            for (i = 1; i <= n; i++) {
                wanted_attrs[attr_list[i]] = 1
            }
        }
        {
            chrom = $1
            source = $2
            start = $3
            end = $4
            strand = $5
            gene_id = $6
            biotype = $7
            gene_name = $8

            # Build gene line attributes
            new_attrs = ""
            sep = ""

            if (wanted_attrs["gene_id"]) {
                new_attrs = sep "gene_id \"" gene_id "\""
                sep = "; "
            }
            if (wanted_attrs["gene_name"] && gene_name != "") {
                new_attrs = new_attrs "; gene_name \"" gene_name "\""
            }
            if (wanted_attrs["gene_biotype"] && biotype != "") {
                new_attrs = new_attrs "; gene_biotype \"" biotype "\""
            }

            print chrom, source, "gene", start, end, ".", strand, ".", new_attrs
        }' "${temp_gene_file}"

        # Then process transcript and exon lines
        awk -F'\t' -v attrs="${attrs}" '
        BEGIN {
            OFS="\t"
            # Parse attribute list
            n = split(attrs, attr_list, ",")
            for (i = 1; i <= n; i++) {
                wanted_attrs[attr_list[i]] = 1
            }
        }
        {
            if ($3 !~ /^(transcript|exon|CDS|five_prime_utr|three_prime_utr|start_codon|stop_codon)$/) next

            attrs_str = $9
            gene_id = ""
            transcript_id = ""
            gene_biotype = ""
            exon_number = ""
            gene_name = ""
            protein_id = ""

            # Parse attributes
            n = split(attrs_str, parts, ";")
            for (i = 1; i <= n; i++) {
                attr = parts[i]
                gsub(/^[ \t]+|[ \t]+$/, "", attr)

                if (attr ~ /^gene_id[ \t]+"/) {
                    match(attr, /gene_id[ \t]+"([^"]+)"/, m)
                    if (m[1] != "") gene_id = m[1]
                } else if (attr ~ /^transcript_id[ \t]+"/) {
                    match(attr, /transcript_id[ \t]+"([^"]+)"/, m)
                    if (m[1] != "") transcript_id = m[1]
                } else if (attr ~ /^gene_biotype[ \t]+"/) {
                    match(attr, /gene_biotype[ \t]+"([^"]+)"/, m)
                    if (m[1] != "") gene_biotype = m[1]
                } else if (attr ~ /^gene[ \t]+"/) {
                    match(attr, /gene[ \t]+"([^"]+)"/, m)
                    if (m[1] != "") gene_name = m[1]
                } else if (attr ~ /^exon_number[ \t]+/) {
                    if (attr ~ /^exon_number[ \t]+"/) {
                        match(attr, /exon_number[ \t]+"([^"]+)"/, m)
                        if (m[1] != "") exon_number = m[1]
                    } else {
                        match(attr, /exon_number[ \t]+([^; \t]+)/, m)
                        if (m[1] != "") exon_number = m[1]
                    }
                } else if (attr ~ /^protein_id[ \t]+"/) {
                    match(attr, /protein_id[ \t]+"([^"]+)"/, m)
                    if (m[1] != "") protein_id = m[1]
                }
            }

            # Store transcript to gene mapping and biotype
            if ($3 == "transcript") {
                if (gene_id != "" && gene_biotype != "") {
                    gene_biotypes[gene_id] = gene_biotype
                }
                if (transcript_id != "" && gene_id != "") {
                    transcript_to_gene[transcript_id] = gene_id
                }
            }

            # Exon lines need to find gene_id through transcript_id
            if ($3 == "exon" && transcript_id != "" && transcript_id in transcript_to_gene) {
                gene_id = transcript_to_gene[transcript_id]
            }

            # Get gene biotype
            if (gene_id != "" && gene_id in gene_biotypes) {
                gene_biotype = gene_biotypes[gene_id]
            }

            # Get gene name
            if (gene_name == "" && gene_id in gene_names) {
                gene_name = gene_names[gene_id]
            }

            # Build new attribute string
            new_attrs = ""
            sep = ""

            if (wanted_attrs["gene_id"] && gene_id != "") {
                new_attrs = sep "gene_id \"" gene_id "\""
                sep = "; "
            }
            if (wanted_attrs["gene_name"] && gene_name != "" && gene_name != gene_id) {
                new_attrs = new_attrs "; gene_name \"" gene_name "\""
            }
            if (wanted_attrs["transcript_id"] && ($3 == "transcript" || $3 == "exon" || $3 == "CDS" || $3 == "five_prime_utr" || $3 == "three_prime_utr" || $3 == "start_codon" || $3 == "stop_codon") && transcript_id != "") {
                new_attrs = new_attrs "; transcript_id \"" transcript_id "\""
            }
            if (wanted_attrs["exon_number"] && ($3 == "exon" || $3 == "CDS" || $3 == "five_prime_utr" || $3 == "three_prime_utr" || $3 == "start_codon" || $3 == "stop_codon") && exon_number != "") {
                new_attrs = new_attrs "; exon_number \"" exon_number "\""
            }
            if (wanted_attrs["protein_id"] && ($3 == "CDS" || $3 == "start_codon" || $3 == "stop_codon") && protein_id != "") {
                new_attrs = new_attrs "; protein_id \"" protein_id "\""
            }
            if (wanted_attrs["gene_biotype"] && gene_biotype != "") {
                new_attrs = new_attrs "; gene_biotype \"" gene_biotype "\""
            }

            print $1,$2,$3,$4,$5,$6,$7,$8,new_attrs
        }' "${input_file}"
    } > "${temp_unsorted}"

    # Step 3: Sort by genomic coordinates
    sort -k1,1 -k4,4n -k5,5n "${temp_unsorted}" > "${output_file}"

    # Clean up temporary files
    rm -f "${temp_gene_file}" "${temp_unsorted}"

    # Clean up temporary file
    rm -f "${temp_gene_file}"
}

# Validate output file
validate_output() {
    local output_file=$2
    local input_lines=$(wc -l < "${input_file}")
    local output_lines=$(wc -l < "${output_file}")

    log_verbose "Input file lines: ${input_lines}"
    log_verbose "Output file lines: ${output_lines}"

    # Check gene_biotype - use grep -c to count
    local total_lines=$(wc -l < "${output_file}")
    local with_biotype=$(grep -c "gene_biotype" "${output_file}" || echo "0")
    local missing_biotype=$((total_lines - with_biotype))

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
                echo "Use -h/--help to see help information"
                exit 1
                ;;
        esac
    done

    # Check required parameters
    if [[ -z "${input_file}" ]]; then
        log_error "Missing required parameter: -i/--input"
        echo "Use -h/--help to see help information"
        exit 1
    fi

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
    local file_type=$(detect_gene_lines "${input_file}")
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
