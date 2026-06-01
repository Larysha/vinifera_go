#!/usr/bin/env Rscript
################################################################################
#                  BUILD ANNOTATION DATA FOR vinifera_go                       #
#                                                                              #
# One-time (re-runnable) setup script. Converts the Vitis vinifera T2T v5.1    #
# blast2go GMT file into three compact .rds objects that the Shiny app loads   #
# at start-up:                                                                 #
#                                                                              #
#   data/vitis_v5_term2gene.rds      data.frame(term, gene)                    #
#   data/vitis_v5_term2name.rds      data.frame(term, name, ontology)          #
#   data/vitis_v5_gene_universe.rds  character vector of annotated genes       #
#                                                                              #
# The ontology (BP / MF / CC) is resolved here, once, via GO.db so that the    #
# deployed app never needs GO.db at runtime (smaller bundle, faster start-up). #
#                                                                              #
# Re-run this when Grapedia releases a new annotation: download the new        #
# blast2go file, convert it to GMT, drop it in data-raw/, and run this script. #
#                                                                              #
# GO annotations sourced from:                                                 #
#   https://grapedia.org/wp-content/uploads/2024/07/T2T_5.1_blast2go.zip       #
# converted to GMT format with a custom script.                                #
#                                                                              #
# Usage:                                                                       #
#   Rscript data-raw/build_annotation.R [path/to/annotations.gmt]              #
################################################################################

suppressPackageStartupMessages({
  library(dplyr)
  library(GO.db)
})

# ---- Resolve input / output paths -------------------------------------------

args <- commandArgs(trailingOnly = TRUE)
gmt_file <- if (length(args) >= 1 && nzchar(args[1])) {
  args[1]
} else {
  "data-raw/blast2go_t2t_5.1.gmt"
}

out_dir <- "data"
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

if (!file.exists(gmt_file)) {
  stop("GMT file not found: ", gmt_file)
}

cat("Building annotation data for vinifera_go\n")
cat(strrep("=", 70), "\n", sep = "")
cat("  Input GMT: ", gmt_file, "\n", sep = "")
cat("  Output:    ", out_dir, "/\n\n", sep = "")

# ---- Parse the GMT file -----------------------------------------------------
# GMT format: GO_ID <tab> GO_NAME <tab> gene1 <tab> gene2 <tab> ...

cat("Parsing GMT file...\n")
gmt_lines <- readLines(gmt_file)

parsed <- lapply(gmt_lines, function(line) {
  parts <- strsplit(line, "\t", fixed = TRUE)[[1]]
  if (length(parts) < 3) return(NULL)
  genes <- parts[3:length(parts)]
  genes <- genes[nzchar(genes)]
  if (length(genes) == 0) return(NULL)
  list(term = parts[1], name = parts[2], genes = genes)
})
parsed <- parsed[!vapply(parsed, is.null, logical(1))]

# One row per gene-term pair
term2gene <- do.call(rbind, lapply(parsed, function(x) {
  data.frame(term = x$term, gene = x$genes, stringsAsFactors = FALSE)
}))

# One row per term
term2name <- data.frame(
  term = vapply(parsed, function(x) x$term, character(1)),
  name = vapply(parsed, function(x) x$name, character(1)),
  stringsAsFactors = FALSE
)
term2name <- term2name[!duplicated(term2name$term), ]

cat("  Parsed", nrow(term2name), "GO terms,",
    nrow(term2gene), "gene-term associations\n\n")

# ---- Annotate each term with its ontology (BP / MF / CC) via GO.db ----------

cat("Resolving ontologies from GO.db (this is the only GO.db step)...\n")
term2name$ontology <- vapply(term2name$term, function(go_id) {
  tryCatch(Ontology(GOTERM[[go_id]]), error = function(e) NA_character_)
}, character(1))

# Drop terms whose ontology could not be resolved (obsolete / unknown IDs)
valid <- term2name$ontology %in% c("BP", "MF", "CC")
n_dropped <- sum(!valid)
if (n_dropped > 0) {
  cat("  Dropped", n_dropped, "terms with unknown / obsolete ontology\n")
}
term2name <- term2name[valid, ]
term2gene <- term2gene[term2gene$term %in% term2name$term, ]

# ---- Gene universe = all genes that carry at least one GO annotation --------

gene_universe <- sort(unique(term2gene$gene))

# ---- Save -------------------------------------------------------------------

saveRDS(term2gene, file.path(out_dir, "vitis_v5_term2gene.rds"), compress = "xz")
saveRDS(term2name, file.path(out_dir, "vitis_v5_term2name.rds"), compress = "xz")
saveRDS(gene_universe, file.path(out_dir, "vitis_v5_gene_universe.rds"),
        compress = "xz")

# ---- Summary ----------------------------------------------------------------

cat("\nDone. Summary:\n")
cat("  GO terms:                ", nrow(term2name), "\n", sep = "")
cat("    BP (Biological Process): ",
    sum(term2name$ontology == "BP"), "\n", sep = "")
cat("    MF (Molecular Function): ",
    sum(term2name$ontology == "MF"), "\n", sep = "")
cat("    CC (Cellular Component): ",
    sum(term2name$ontology == "CC"), "\n", sep = "")
cat("  Gene-term associations:  ", nrow(term2gene), "\n", sep = "")
cat("  Annotated genes (universe): ", length(gene_universe), "\n", sep = "")
cat("\nFiles written to ", out_dir, "/:\n", sep = "")
for (f in c("vitis_v5_term2gene.rds", "vitis_v5_term2name.rds",
            "vitis_v5_gene_universe.rds")) {
  size_kb <- round(file.info(file.path(out_dir, f))$size / 1024, 1)
  cat("  ", f, " (", size_kb, " KB)\n", sep = "")
}
