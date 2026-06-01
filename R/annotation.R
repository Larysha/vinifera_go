################################################################################
# Annotation loading
#
# Reads the pre-built .rds annotation objects (produced by
# data-raw/build_annotation.R) once, at app start-up. No GO.db needed here:
# the ontology of every term is already baked into term2name.
################################################################################

#' Load the Vitis T2T v5.1 GO annotation
#'
#' @param data_dir Directory holding the three vitis_v5_*.rds files
#' @return List with:
#'   - term2gene: data.frame(term, gene)
#'   - term2name: data.frame(term, name, ontology)
#'   - universe:  character vector of all annotated genes
#'   - n_terms, n_genes, ontology_counts: summary values for the UI
load_annotation <- function(data_dir = "data") {
  files <- c(
    term2gene = "vitis_v5_term2gene.rds",
    term2name = "vitis_v5_term2name.rds",
    universe  = "vitis_v5_gene_universe.rds"
  )
  paths <- stats::setNames(file.path(data_dir, files), names(files))
  missing <- !file.exists(paths)
  if (any(missing)) {
    stop("Missing annotation file(s): ",
         paste(files[missing], collapse = ", "),
         ". Run data-raw/build_annotation.R first.")
  }

  term2gene <- readRDS(paths[["term2gene"]])
  term2name <- readRDS(paths[["term2name"]])
  universe  <- readRDS(paths[["universe"]])

  list(
    term2gene = term2gene,
    term2name = term2name,
    universe  = universe,
    n_terms   = nrow(term2name),
    n_genes   = length(universe),
    ontology_counts = table(term2name$ontology)
  )
}
