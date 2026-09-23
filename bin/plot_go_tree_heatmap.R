#!/usr/bin/env Rscript
#########################################################################
# File Name: bin/plot_go_tree_heatmap.R
# Author: ChengYu
# Description: Combine a GO semantic similarity tree (NJ) with an
#              enrichment barplot (NES coloured by p-value) and GO
#              description labels into a single publication-quality
#              figure.  Uses GOSemSim for pairwise similarity, ape for
#              neighbour-joining tree construction, ggtree for tree
#              rendering, and aplot for layout assembly.
# Created Time: 2026
#########################################################################

suppressMessages(library(getopt))
suppressMessages(library(GOSemSim))
suppressMessages(library(ape))
suppressMessages(library(ggtree))
suppressMessages(library(ggplot2))
suppressMessages(library(aplot))

# ── CLI specification ──────────────────────────────────────────────────
spec <- matrix(c(
    "input",    "i", 1, "character", "Path to input TSV. Required columns: GO_ID, Description, Cluster, NES, pvalue. Required.",
    "output",   "o", 1, "character", "Output PDF file path. Required.",
    "ontology", "O", 2, "character", "GO ontology: BP, MF, or CC. Default: 'BP'.",
    "organism", "g", 2, "character", "OrgDb name for GOSemSim (e.g. 'org.At.tair.db', 'org.Hs.eg.db'). Default: 'org.At.tair.db'.",
    "measure",  "m", 2, "character", "Semantic similarity measure ('Wang', 'Rel', 'Jiang', etc.). Default: 'Wang'.",
    "keytype",  "k", 2, "character", "Gene ID keytype passed to GOSemSim godata(). Default: 'GID'.",
    "height",   "H", 2, "numeric",   "Plot height in inches. Default: 8.",
    "width",    "W", 2, "numeric",   "Plot width in inches. Default: 14.",
    "help",     "h", 0, "logical",   "Print this help message and exit."
), byrow = TRUE, ncol = 5)
colnames(spec) <- c("long", "short", "argflag", "type", "help")

opt <- getopt(spec)

# ── Help ───────────────────────────────────────────────────────────────
if (!is.null(opt$help)) {
    cat("Usage: Rscript plot_go_tree_heatmap.R -i <go_results.tsv> -o <output.pdf> [options]\n\n")
    cat("Combine a GO semantic similarity tree with an enrichment barplot.\n\n")
    cat("Options:\n")
    for (i in seq_len(nrow(spec))) {
        flag <- paste0("-", spec[i, "short"], ", --", spec[i, "long"])
        cat(sprintf("  %-30s %s\n", flag, spec[i, "help"]))
    }
    cat("\nInput TSV must have columns (tab-delimited, header row):\n")
    cat("  GO_ID       GO term accession (e.g. GO:0008150)\n")
    cat("  Description Human-readable GO term name\n")
    cat("  Cluster     Group/cluster label for facet wrapping\n")
    cat("  NES         Normalised enrichment score (numeric)\n")
    cat("  pvalue      Adjusted p-value for colour mapping (numeric)\n")
    cat("\nExamples:\n")
    cat("  Rscript plot_go_tree_heatmap.R -i go_results.tsv -o go_tree.pdf\n")
    cat("  Rscript plot_go_tree_heatmap.R -i go.tsv -o out.pdf -g org.Hs.eg.db -O BP\n")
    cat("  Rscript plot_go_tree_heatmap.R -i go.tsv -o out.pdf --measure Rel --keytype ENTREZID\n")
    quit(status = 0)
}

# ── Validate required arguments ────────────────────────────────────────
if (is.null(opt$input)) {
    stop("Error: -i/--input is required. Use -h for help.")
}
if (is.null(opt$output)) {
    stop("Error: -o/--output is required. Use -h for help.")
}
if (!file.exists(opt$input)) {
    stop(paste0("Error: input file not found: ", opt$input))
}

# ── Set defaults ───────────────────────────────────────────────────────
ontology   <- if (is.null(opt$ontology)) "BP"            else toupper(opt$ontology)
organism   <- if (is.null(opt$organism)) "org.At.tair.db" else opt$organism
measure    <- if (is.null(opt$measure))  "Wang"           else opt$measure
keytype    <- if (is.null(opt$keytype))  "GID"            else opt$keytype
plot_height <- if (is.null(opt$height))  8                else opt$height
plot_width  <- if (is.null(opt$width))   14               else opt$width

# Validate ontology
valid_ont <- c("BP", "MF", "CC")
if (!(ontology %in% valid_ont)) {
    stop(paste0("Error: --ontology must be one of ", paste(valid_ont, collapse = "/"),
                ". Got: ", ontology))
}

# ── Read input data ────────────────────────────────────────────────────
df <- tryCatch(
    read.table(opt$input, header = TRUE, sep = "\t", quote = "",
               comment.char = "", stringsAsFactors = FALSE, check.names = FALSE),
    error = function(e) {
        stop(paste0("Error reading input file: ", e$message))
    }
)

required_cols <- c("GO_ID", "Description", "Cluster", "NES", "pvalue")
missing_cols <- setdiff(required_cols, colnames(df))
if (length(missing_cols) > 0) {
    stop(paste0("Error: input file is missing required columns: ",
                paste(missing_cols, collapse = ", ")))
}
if (nrow(df) == 0) {
    stop("Error: input file contains no data rows.")
}

# Ensure numeric columns
df$NES    <- suppressWarnings(as.numeric(df$NES))
df$pvalue <- suppressWarnings(as.numeric(df$pvalue))
if (anyNA(df$NES)) {
    stop("Error: column 'NES' contains non-numeric values.")
}
if (anyNA(df$pvalue)) {
    stop("Error: column 'pvalue' contains non-numeric values.")
}

unique_go    <- unique(df$GO_ID)
n_go         <- length(unique_go)
n_clusters   <- length(unique(df$Cluster))
n_total_rows <- nrow(df)

cat(sprintf("Input summary:\n"))
cat(sprintf("  Total rows:      %d\n", n_total_rows))
cat(sprintf("  Unique GO terms: %d\n", n_go))
cat(sprintf("  Clusters:        %d (%s)\n", n_clusters, paste(unique(df$Cluster), collapse = ", ")))
cat(sprintf("  Ontology:        %s\n", ontology))
cat(sprintf("  Organism:        %s\n", organism))
cat(sprintf("  Similarity:      %s\n", measure))

# ── Build GO semantic similarity and NJ tree ──────────────────────────
# Load the OrgDb package and construct GOSemSim data
org_pkg <- organism
if (!requireNamespace(org_pkg, quietly = TRUE)) {
    stop(paste0("Error: organism package '", org_pkg, "' is not installed. ",
                "Install it with BiocManager::install('", org_pkg, "')."))
}

go_data <- tryCatch(
    godata(OrgDb = org_pkg, ont = ontology, keytype = keytype),
    error = function(e) {
        stop(paste0("Error creating GOSemSim godata object: ", e$message))
    }
)

# Compute pairwise GO term similarity
go_sim <- tryCatch(
    mgoSim(unique_go, unique_go, semData = go_data, measure = measure, combine = NULL),
    error = function(e) {
        stop(paste0("Error computing pairwise GO similarity: ", e$message))
    }
)

# Replace any NA in similarity matrix with 0 (can happen for orphan terms)
go_sim[is.na(go_sim)] <- 0

# ── Decide: tree + barplot or fallback barplot only ────────────────────
MIN_TERMS_FOR_TREE <- 3
use_tree <- n_go >= MIN_TERMS_FOR_TREE

if (use_tree) {
    # Distance matrix: 1 - similarity
    go_dist <- as.dist(1 - go_sim)

    # Build neighbour-joining tree
    tree <- tryCatch(
        nj(go_dist),
        error = function(e) {
            warning(paste0("NJ tree construction failed (", e$message,
                           "). Falling back to barplot only."))
            return(NULL)
        }
    )

    if (is.null(tree)) {
        use_tree <- FALSE
    }
}

# ── Prepare barplot data ──────────────────────────────────────────────
# Reorder factor levels to match the tree tip order when tree is used
if (use_tree) {
    tip_order <- tree$tip.label
    df$GO_ID <- factor(df$GO_ID, levels = rev(tip_order))
} else {
    df$GO_ID <- factor(df$GO_ID, levels = rev(sort(unique(df$GO_ID))))
}

p_bar <- ggplot(df, aes(x = GO_ID, y = NES, fill = pvalue)) +
    geom_bar(stat = "identity", width = 0.7) +
    coord_flip() +
    facet_wrap(~ Cluster, nrow = 1, scales = "free_x") +
    scale_fill_gradient(low = "red", high = "blue", name = "p-value") +
    labs(y = "NES", x = NULL) +
    theme_minimal(base_size = 10) +
    theme(
        axis.text.y   = element_blank(),
        axis.ticks.y  = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank(),
        strip.text    = element_text(face = "bold")
    )

# ── Build tree plot ───────────────────────────────────────────────────
if (use_tree) {
    p_tree <- ggtree(tree, branch.length = "none") +
        geom_tiplab(size = 3) +
        coord_cartesian(xlim = c(-0.1, 1.3), clip = "off") +
        theme_tree2()
} else {
    cat(sprintf("Note: fewer than %d GO terms — skipping tree, producing barplot only.\n",
                MIN_TERMS_FOR_TREE))
}

# ── Build description text panel ──────────────────────────────────────
# Match GO_ID order to the barplot y-axis
desc_df <- unique(df[, c("GO_ID", "Description")])
if (use_tree) {
    desc_df$GO_ID <- factor(desc_df$GO_ID, levels = rev(tip_order))
} else {
    desc_df$GO_ID <- factor(desc_df$GO_ID, levels = rev(sort(unique(as.character(desc_df$GO_ID)))))
}
desc_df <- desc_df[order(desc_df$GO_ID), ]

p_text <- ggplot(desc_df, aes(x = 1, y = GO_ID, label = Description)) +
    geom_text(aes(hjust = 0), size = 3, family = "sans", fontface = "italic") +
    ylim(0.5, n_go + 0.5) +
    theme_void() +
    coord_flip() +
    scale_x_continuous(limits = c(0.5, 1.5))

# ── Assemble combined plot ────────────────────────────────────────────
if (use_tree) {
    p_combined <- p_bar %>%
        insert_left(p_tree, width = 0.2) %>%
        insert_right(p_text, width = 0.8)
} else {
    # Barplot only + description on the right
    p_combined <- p_bar %>%
        insert_right(p_text, width = 0.8)
}

# ── Write PDF ─────────────────────────────────────────────────────────
pdf(file = opt$output, height = plot_height, width = plot_width)
tryCatch(
    {
        print(p_combined)
        cat(sprintf("\nPlot saved to %s\n", opt$output))
        if (use_tree) {
            cat(sprintf("  Mode:            tree + barplot + description\n"))
        } else {
            cat(sprintf("  Mode:            barplot + description (no tree)\n"))
        }
        cat(sprintf("  GO terms plotted: %d\n", n_go))
    },
    error = function(e) {
        stop(paste0("Error rendering plot: ", e$message))
    },
    finally = {
        dev.off()
    }
)
