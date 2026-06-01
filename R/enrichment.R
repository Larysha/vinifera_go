################################################################################
# Enrichment core
#
# GO over-representation (ORA) and GSEA on Vitis vinifera T2T v5.1 gene IDs.
# Adapted from ~/phd/go_analysis_tut/go_enrichment.R and
# ~/phd/mocha/.../09_motif_go_analysis/04_gsea_importance.R, with the
# config-file machinery removed so the functions are callable directly.
################################################################################

#' Subset a TERM2GENE table to a single ontology
#'
#' @param annotation List from load_annotation() (term2gene, term2name, universe)
#' @param ont Ontology code: "BP", "MF" or "CC"
#' @return TERM2GENE data.frame (term, gene) restricted to that ontology
.term2gene_for_ontology <- function(annotation, ont) {
  terms <- annotation$term2name$term[annotation$term2name$ontology == ont]
  annotation$term2gene[annotation$term2gene$term %in% terms, , drop = FALSE]
}

#' Numeric GeneRatio / BgRatio from clusterProfiler "a/b" strings
.ratio_to_numeric <- function(x) {
  vapply(x, function(r) {
    p <- as.numeric(strsplit(r, "/", fixed = TRUE)[[1]])
    p[1] / p[2]
  }, numeric(1))
}

#' Run ORA for one ontology
#'
#' Thin wrapper around clusterProfiler::enricher(). Returns a tidy data.frame
#' (or NULL when nothing can be tested / nothing is enriched).
#'
#' @param genes Character vector of gene IDs of interest
#' @param annotation Annotation list from load_annotation()
#' @param ont Ontology code ("BP", "MF", "CC")
#' @param universe Character vector to use as background
#' @param p_adjust Adjustment method: "BH", "bonferroni", "holm" or "none"
#' @param p_cutoff Adjusted p-value cutoff applied inside enricher()
#' @param min_size,max_size Min / max gene-set size to test
#' @return tidy data.frame or NULL
run_ora_ontology <- function(genes, annotation, ont, universe,
                             p_adjust = "BH", p_cutoff = 0.05,
                             min_size = 10, max_size = 500) {
  t2g <- .term2gene_for_ontology(annotation, ont)
  if (nrow(t2g) == 0) return(NULL)

  t2n <- annotation$term2name[annotation$term2name$ontology == ont,
                              c("term", "name")]

  # Only genes present in the background can contribute
  genes_in <- intersect(genes, universe)
  if (length(genes_in) < 1) return(NULL)

  res <- tryCatch(
    clusterProfiler::enricher(
      gene          = genes_in,
      universe      = universe,
      TERM2GENE     = t2g,
      TERM2NAME     = t2n,
      pAdjustMethod = p_adjust,
      pvalueCutoff  = 1,        # keep everything; threshold applied downstream
      qvalueCutoff  = 1,
      minGSSize     = min_size,
      maxGSSize     = max_size
    ),
    error = function(e) NULL
  )

  if (is.null(res) || nrow(as.data.frame(res)) == 0) return(NULL)

  df <- as.data.frame(res, stringsAsFactors = FALSE)
  rownames(df) <- NULL
  df$Ontology <- ont
  df$GeneRatioNum <- .ratio_to_numeric(df$GeneRatio)
  df$BgRatioNum <- .ratio_to_numeric(df$BgRatio)
  df
}

#' Run ORA across the selected ontologies
#'
#' @param genes Character vector of gene IDs of interest
#' @param annotation Annotation list from load_annotation()
#' @param ontologies Character vector subset of c("BP","MF","CC")
#' @param background Optional character vector; defaults to the full annotated
#'   universe. Intersected with the annotated universe before use.
#' @inheritParams run_ora_ontology
#' @return Combined tidy data.frame (possibly 0 rows), never NULL
run_ora <- function(genes, annotation, ontologies = c("BP", "MF", "CC"),
                    background = NULL, p_adjust = "BH", p_cutoff = 0.05,
                    min_size = 10, max_size = 500) {
  universe <- if (is.null(background)) {
    annotation$universe
  } else {
    intersect(background, annotation$universe)
  }

  parts <- lapply(ontologies, function(ont) {
    run_ora_ontology(genes, annotation, ont, universe,
                     p_adjust = p_adjust, p_cutoff = p_cutoff,
                     min_size = min_size, max_size = max_size)
  })
  parts <- parts[!vapply(parts, is.null, logical(1))]

  if (length(parts) == 0) return(empty_ora_result())
  do.call(rbind, parts)
}

#' Run ORA separately for up- and down-regulated genes
#'
#' @param genes Character vector of gene IDs
#' @param expression Numeric vector (same length / order as genes), e.g. log2FC
#' @param ... Passed to run_ora()
#' @return Combined data.frame with an extra "Direction" column ("Up"/"Down")
run_ora_split <- function(genes, expression, annotation, ...) {
  up <- genes[expression > 0]
  down <- genes[expression < 0]

  res <- list()
  if (length(up) > 0) {
    r <- run_ora(up, annotation, ...)
    if (nrow(r) > 0) { r$Direction <- "Up"; res$up <- r }
  }
  if (length(down) > 0) {
    r <- run_ora(down, annotation, ...)
    if (nrow(r) > 0) { r$Direction <- "Down"; res$down <- r }
  }

  if (length(res) == 0) {
    out <- empty_ora_result()
    out$Direction <- character(0)
    return(out)
  }
  do.call(rbind, res)
}

#' Run GSEA across the selected ontologies on a ranked gene list
#'
#' @param ranked Named numeric vector (names = gene IDs), need not be sorted
#' @param annotation Annotation list from load_annotation()
#' @param ontologies Character vector subset of c("BP","MF","CC")
#' @param p_adjust,p_cutoff,min_size,max_size As for run_ora
#' @return Combined tidy GSEA data.frame (possibly 0 rows). Carries the per-
#'   ontology GSEA objects in attr(x, "gsea_objects") for enrichment plots.
run_gsea <- function(ranked, annotation, ontologies = c("BP", "MF", "CC"),
                     p_adjust = "BH", p_cutoff = 0.05,
                     min_size = 10, max_size = 500) {
  ranked <- sort(ranked, decreasing = TRUE)

  objs <- list()
  parts <- list()
  for (ont in ontologies) {
    t2g <- .term2gene_for_ontology(annotation, ont)
    if (nrow(t2g) == 0) next
    t2n <- annotation$term2name[annotation$term2name$ontology == ont,
                                c("term", "name")]

    obj <- tryCatch(
      clusterProfiler::GSEA(
        geneList      = ranked,
        TERM2GENE     = t2g,
        TERM2NAME     = t2n,
        pvalueCutoff  = 1,
        pAdjustMethod = p_adjust,
        minGSSize     = min_size,
        maxGSSize     = max_size,
        by            = "fgsea",
        seed          = TRUE,
        verbose       = FALSE
      ),
      error = function(e) NULL
    )

    if (is.null(obj) || nrow(as.data.frame(obj)) == 0) next

    objs[[ont]] <- obj
    df <- as.data.frame(obj, stringsAsFactors = FALSE)
    rownames(df) <- NULL
    df$Ontology <- ont
    parts[[ont]] <- df
  }

  if (length(parts) == 0) {
    out <- empty_gsea_result()
    attr(out, "gsea_objects") <- objs
    return(out)
  }
  out <- do.call(rbind, parts)
  attr(out, "gsea_objects") <- objs
  out
}

#' Remove redundant GO terms (ShinyGO-style), using gene-set overlap only
#'
#' Walks terms from most to least significant; drops a term if its gene set
#' overlaps an already-kept (more significant) term beyond `overlap_cutoff`.
#' Overlap is containment: |A ∩ B| / min(|A|, |B|), which catches parent/child
#' GO terms that share most of their genes. Uses only the gene IDs in the
#' result, so no GO.db / semantic-similarity data is needed at runtime.
#'
#' @param df ORA result data.frame (needs geneID, p.adjust)
#' @param overlap_cutoff Numeric in (0, 1]; higher = less aggressive
#' @return Filtered data.frame (same columns)
remove_redundancy <- function(df, overlap_cutoff = 0.5) {
  if (nrow(df) <= 1) return(df)

  ord <- order(df$p.adjust)
  df_ord <- df[ord, , drop = FALSE]
  gene_sets <- strsplit(df_ord$geneID, "/", fixed = TRUE)

  keep <- logical(nrow(df_ord))
  kept_sets <- list()
  for (i in seq_len(nrow(df_ord))) {
    gi <- gene_sets[[i]]
    redundant <- FALSE
    for (gk in kept_sets) {
      inter <- length(intersect(gi, gk))
      if (inter / min(length(gi), length(gk)) > overlap_cutoff) {
        redundant <- TRUE
        break
      }
    }
    if (!redundant) {
      keep[i] <- TRUE
      kept_sets[[length(kept_sets) + 1]] <- gi
    }
  }

  out <- df_ord[keep, , drop = FALSE]
  rownames(out) <- NULL
  out
}

#' Empty ORA result with the canonical column set (for graceful "no results")
empty_ora_result <- function() {
  data.frame(
    ID = character(0), Description = character(0), GeneRatio = character(0),
    BgRatio = character(0), pvalue = numeric(0), p.adjust = numeric(0),
    qvalue = numeric(0), geneID = character(0), Count = integer(0),
    Ontology = character(0), GeneRatioNum = numeric(0), BgRatioNum = numeric(0),
    stringsAsFactors = FALSE
  )
}

#' Empty GSEA result with the canonical column set
empty_gsea_result <- function() {
  data.frame(
    ID = character(0), Description = character(0), setSize = integer(0),
    enrichmentScore = numeric(0), NES = numeric(0), pvalue = numeric(0),
    p.adjust = numeric(0), qvalue = numeric(0), core_enrichment = character(0),
    Ontology = character(0), stringsAsFactors = FALSE
  )
}
