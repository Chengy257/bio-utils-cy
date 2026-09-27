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
# Changelog:
#   v1.1.0  2026-09-28
#   - FIX: the description panel mapped a discrete y (GO_ID factor)
#     onto ylim(0.5, n+0.5), a continuous scale — every render failed
#     with "Discrete value supplied to a continuous scale" and no
#     output of any kind could be produced. The bogus ylim is removed.
#   - FIX: the OrgDb package (GOSemSim godata) was loaded even for
#     inputs with fewer than 3 GO terms, where no tree is built and no
#     similarity is needed. The OrgDb is now required only in tree
#     mode; without it the script degrades to barplot-only with a
#     warning (same fallback policy as a failed NJ tree).
#   - CHANGE: bar fill now encodes -log10(p) (field standard) instead
#     of the raw p-value, which visually compressed significance
#     differences. Colour semantics are preserved: red = significant,
#     blue = not. p <= 0 is clipped to the smallest positive double
#     with a warning.
#   - FIX: when one GO_ID carries several descriptions, the label panel
#     stacked one text label per description on the same y position
#     (overlapping text). The first description is used; conflicts are
#     reported with a warning.
#   - FIX: aplot/patchwork assembly shared the y axis through the
#     NES mapping (coord_flip'ed plots), which aplot cannot align; the
#     GO_ID is now the y aesthetic of both panels (horizontal bars,
#     no coord_flip), the only layout aplot can align.
#   - FIX: aplot <= 0.2.3 cannot render legends inside assembled plots
#     with ggplot2 4.x ("object is not a unit"). The figure is first
#     drawn with the legend; on failure it is re-drawn legend-less with
#     a warning instead of producing nothing.
#   - CHANGE: INFO/summary logs go to stderr via message(); stdout
#     stays reserved for data. dev.off() no longer prints "null
#     device" to stdout.
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
    "organism", "g", 2, "character", "OrgDb name for GOSemSim (e.g. 'org.At.tair.db', 'org.Hs.eg.db'). Required only when >= 3 GO terms (tree mode). Default: 'org.At.tair.db'.",
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
    cat("  pvalue      Adjusted p-value; bars are coloured by -log10(p)\n")
    cat("              (red = significant, blue = not)\n")
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

# -log10(p) colour variable; clip non-positive p so the log stays finite.
if (any(df$pvalue <= 0)) {
    warning(sprintf("%d p-value(s) <= 0 clipped to %.3e for the -log10 colour scale.",
                    sum(df$pvalue <= 0), .Machine$double.xmin))
}
df$neg_log10_p <- -log10(pmax(df$pvalue, .Machine$double.xmin))

unique_go    <- unique(df$GO_ID)
n_go         <- length(unique_go)
n_clusters   <- length(unique(df$Cluster))
n_total_rows <- nrow(df)

message(sprintf("Input summary:"))
message(sprintf("  Total rows:      %d", n_total_rows))
message(sprintf("  Unique GO terms: %d", n_go))
message(sprintf("  Clusters:        %d (%s)", n_clusters, paste(unique(df$Cluster), collapse = ", ")))
message(sprintf("  Ontology:        %s", ontology))
message(sprintf("  Organism:        %s", organism))
message(sprintf("  Similarity:      %s", measure))

# ── Build GO semantic similarity and NJ tree (tree mode only) ─────────
MIN_TERMS_FOR_TREE <- 3
use_tree <- n_go >= MIN_TERMS_FOR_TREE

if (use_tree) {
    # The OrgDb is only needed to build the semantic-similarity tree;
    # barplot-only output does not depend on it.
    if (!requireNamespace(organism, quietly = TRUE)) {
        warning(paste0("organism package '", organism, "' is not installed — ",
                       "falling back to barplot-only output. Install it with ",
                       "BiocManager::install('", organism, "') to get the tree."))
        use_tree <- FALSE
    }
}

if (use_tree) {
    go_data <- tryCatch(
        godata(OrgDb = organism, ont = ontology, keytype = keytype),
        error = function(e) {
            warning(paste0("GOSemSim godata failed (", e$message,
                           "). Falling back to barplot only."))
            NULL
        }
    )
    if (is.null(go_data)) {
        use_tree <- FALSE
    } else {
        go_sim <- tryCatch(
            mgoSim(unique_go, unique_go, semData = go_data, measure = measure, combine = NULL),
            error = function(e) {
                warning(paste0("Error computing pairwise GO similarity (", e$message,
                               "). Falling back to barplot only."))
                NULL
            }
        )
        if (is.null(go_sim)) {
            use_tree <- FALSE
        } else {
            # NA can occur for orphan terms
            go_sim[is.na(go_sim)] <- 0
            go_dist <- as.dist(1 - go_sim)
            tree <- tryCatch(
                nj(go_dist),
                error = function(e) {
                    warning(paste0("NJ tree construction failed (", e$message,
                                   "). Falling back to barplot only."))
                    NULL
                }
            )
            if (is.null(tree)) {
                use_tree <- FALSE
            }
        }
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

p_bar <- ggplot(df, aes(x = NES, y = GO_ID, fill = neg_log10_p)) +
    geom_col(width = 0.7) +
    facet_wrap(~ Cluster, nrow = 1, scales = "free_x") +
    scale_fill_gradient(low = "blue", high = "red", name = "-log10(p)") +
    labs(x = "NES", y = NULL) +
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
    message(sprintf("Note: fewer than %d GO terms (or tree unavailable) — producing barplot only.",
                    MIN_TERMS_FOR_TREE))
}

# ── Build description text panel ──────────────────────────────────────
# One label per GO_ID: keep the first description, warn on conflicts.
dup_ids <- unique(df$GO_ID[duplicated(df$GO_ID)])
if (length(dup_ids) > 0) {
    conflicting <- dup_ids[vapply(dup_ids, function(g) {
        length(unique(df$Description[df$GO_ID == g])) > 1
    }, logical(1))]
    if (length(conflicting) > 0) {
        warning(paste0("GO terms with multiple descriptions (first one used): ",
                       paste(conflicting, collapse = ", ")))
    }
}
desc_df <- df[!duplicated(df$GO_ID), c("GO_ID", "Description")]
# Match GO_ID order to the barplot y-axis
if (use_tree) {
    desc_df$GO_ID <- factor(as.character(desc_df$GO_ID), levels = rev(tip_order))
} else {
    desc_df$GO_ID <- factor(as.character(desc_df$GO_ID),
                            levels = rev(sort(unique(as.character(desc_df$GO_ID)))))
}
desc_df <- desc_df[order(desc_df$GO_ID), ]

p_text <- ggplot(desc_df, aes(x = 1, y = GO_ID, label = Description)) +
    geom_text(aes(hjust = 0), size = 3, family = "sans", fontface = "italic") +
    theme_void() +
    scale_x_continuous(limits = c(0.5, 1.5))

# ── Write PDF ─────────────────────────────────────────────────────────
assemble <- function(pb) {
    if (use_tree) {
        pb %>% insert_left(p_tree, width = 0.2) %>% insert_right(p_text, width = 0.8)
    } else {
        pb %>% insert_right(p_text, width = 0.8)
    }
}

draw_combined <- function(pb, err_out) {
    pdf(file = opt$output, height = plot_height, width = plot_width)
    ok <- tryCatch(
        { print(assemble(pb)); TRUE },
        error = function(e) { assign("msg", conditionMessage(e), envir = err_out); FALSE }
    )
    invisible(dev.off())
    ok
}

err_env <- new.env(parent = emptyenv())
err_env$msg <- ""
if (!draw_combined(p_bar, err_env)) {
    # aplot <= 0.2.3 cannot render legends inside assembled plots with
    # ggplot2 4.x ("object is not a unit"). Degrade: redraw legend-less.
    warning(paste0("Legend rendering failed (", err_env$msg, "; aplot ",
                   packageVersion("aplot"), " + ggplot2 ",
                   packageVersion("ggplot2"), " legend incompatibility). ",
                   "Re-drawing without the -log10(p) colour legend — ",
                   "upgrade aplot to restore it."))
    if (!draw_combined(p_bar + theme(legend.position = "none"), err_env)) {
        stop(paste0("Error rendering plot: ", err_env$msg))
    }
}

message(sprintf("Plot saved to %s", opt$output))
if (use_tree) {
    message("  Mode:            tree + barplot + description")
} else {
    message("  Mode:            barplot + description (no tree)")
}
message(sprintf("  GO terms plotted: %d", n_go))
