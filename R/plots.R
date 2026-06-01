################################################################################
# Plots
#
# ggplot2 builders that preserve the style of go_analysis_tut/go_enrichment.R:
# faceted dotplot (size = Count, colour = -log10(FDR)) and a horizontal barplot
# of -log10(FDR), in the burgundy / blue / teal ontology palette. Plus a GSEA
# enrichment-score plot via enrichplot::gseaplot2().
#
# Every builder returns a ggplot (or recordedplot/grob), so the same object can
# be shown in the app and saved to PNG / PDF.
################################################################################

ONTOLOGY_COLOURS <- c(
  BP = "#671436",  # burgundy  - Biological Process
  MF = "#052b67",  # dark blue - Molecular Function
  CC = "#176272"   # teal      - Cellular Component
)

ONTOLOGY_LABELS <- c(
  BP = "BP (Biological Process)",
  MF = "MF (Molecular Function)",
  CC = "CC (Cellular Component)"
)

#' Wrap a long string onto (roughly) two balanced lines
.wrap_label <- function(x, width = 40) {
  vapply(x, function(s) {
    if (nchar(s) <= width) return(s)
    words <- strsplit(s, " ", fixed = TRUE)[[1]]
    if (length(words) < 2) return(s)
    mid <- ceiling(length(words) / 2)
    paste0(paste(words[1:mid], collapse = " "), "\n",
           paste(words[(mid + 1):length(words)], collapse = " "))
  }, character(1))
}

#' Keep the top-N most significant terms per ontology (and per direction, when
#' a split-by-direction result is supplied)
.top_per_ontology <- function(df, top_n) {
  df <- df[order(df$p.adjust), , drop = FALSE]
  keys <- if ("Direction" %in% names(df)) {
    paste(df$Ontology, df$Direction, sep = "|")
  } else {
    df$Ontology
  }
  do.call(rbind, lapply(split(df, keys), utils::head, top_n))
}

#' ORA dot plot, faceted by ontology
#'
#' @param df ORA result data.frame (run_ora output)
#' @param top_n Number of terms to show per ontology
#' @return ggplot object (or NULL if df is empty)
go_dotplot <- function(df, top_n = 20) {
  if (nrow(df) == 0) return(NULL)
  d <- .top_per_ontology(df, top_n)
  d$Description_wrapped <- .wrap_label(d$Description)
  d$OntologyLabel <- ONTOLOGY_LABELS[d$Ontology]
  split_mode <- "Direction" %in% names(d)

  facet <- if (split_mode) {
    ggplot2::facet_grid(OntologyLabel ~ Direction, scales = "free_y",
                        space = "free_y")
  } else {
    ggplot2::facet_grid(OntologyLabel ~ ., scales = "free_y",
                        space = "free_y")
  }
  subtitle <- if (split_mode) {
    "Up- and down-regulated genes analysed separately"
  } else {
    NULL
  }

  ggplot2::ggplot(
    d,
    ggplot2::aes(x = .data$GeneRatioNum,
                 y = stats::reorder(.data$Description_wrapped,
                                    .data$GeneRatioNum))
  ) +
    ggplot2::geom_point(ggplot2::aes(size = .data$Count,
                                     colour = -log10(.data$p.adjust))) +
    ggplot2::scale_colour_gradient(low = "#3b5bbf", high = "#d62828",
                                   name = "-log10(FDR)") +
    ggplot2::scale_size_continuous(name = "Gene count", range = c(2.5, 9)) +
    facet +
    ggplot2::labs(x = "Gene ratio", y = NULL,
                  title = "GO enrichment (over-representation)",
                  subtitle = subtitle) +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(
      strip.text = ggplot2::element_text(face = "bold"),
      axis.text.y = ggplot2::element_text(size = 9),
      plot.title = ggplot2::element_text(face = "bold")
    )
}

#' ORA bar plot of -log10(FDR), coloured by ontology
#'
#' @param df ORA result data.frame
#' @param top_n Terms per ontology
#' @param fdr_cutoff Draws a reference line at -log10(cutoff)
#' @return ggplot object (or NULL if df is empty)
go_barplot <- function(df, top_n = 20, fdr_cutoff = 0.05) {
  if (nrow(df) == 0) return(NULL)
  d <- .top_per_ontology(df, top_n)
  d$neglogFDR <- -log10(d$p.adjust)
  d$label <- paste0("[", d$Ontology, "] ", d$Description)
  d$label <- .wrap_label(d$label, width = 55)
  d <- d[order(d$neglogFDR), , drop = FALSE]
  # Labels can repeat across the Up / Down panels, so make them unique per row
  # while keeping the displayed text via a named factor.
  d$row_id <- factor(seq_len(nrow(d)), levels = seq_len(nrow(d)))
  split_mode <- "Direction" %in% names(d)

  p <- ggplot2::ggplot(d, ggplot2::aes(x = .data$neglogFDR, y = .data$row_id,
                                       fill = .data$Ontology)) +
    ggplot2::geom_col(colour = "grey20", width = 0.75) +
    ggplot2::geom_vline(xintercept = -log10(fdr_cutoff),
                        linetype = "dashed", colour = "#dc3030") +
    ggplot2::scale_y_discrete(labels = stats::setNames(d$label, d$row_id)) +
    ggplot2::scale_fill_manual(values = ONTOLOGY_COLOURS,
                               labels = ONTOLOGY_LABELS, name = "Ontology") +
    ggplot2::labs(x = "-log10(FDR)", y = NULL,
                  title = "GO enrichment (top terms by FDR)") +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(
      axis.text.y = ggplot2::element_text(size = 9),
      plot.title = ggplot2::element_text(face = "bold"),
      legend.position = "bottom"
    )

  if (split_mode) {
    p <- p +
      ggplot2::facet_grid(Direction ~ ., scales = "free_y", space = "free_y") +
      ggplot2::labs(subtitle = "Up- and down-regulated genes analysed separately")
  }
  p
}

#' GSEA dot plot (NES vs term, sized by set size)
#'
#' @param df GSEA result data.frame (run_gsea output)
#' @param top_n Terms per ontology
#' @return ggplot object (or NULL if df is empty)
gsea_dotplot <- function(df, top_n = 20) {
  if (nrow(df) == 0) return(NULL)
  d <- df[order(df$Ontology, df$p.adjust), , drop = FALSE]
  d <- do.call(rbind, lapply(split(d, d$Ontology), utils::head, top_n))
  d$Description_wrapped <- .wrap_label(d$Description)
  d$OntologyLabel <- ONTOLOGY_LABELS[d$Ontology]

  ggplot2::ggplot(
    d,
    ggplot2::aes(x = .data$NES,
                 y = stats::reorder(.data$Description_wrapped, .data$NES))
  ) +
    ggplot2::geom_point(ggplot2::aes(size = .data$setSize,
                                     colour = -log10(.data$p.adjust))) +
    ggplot2::geom_vline(xintercept = 0, linetype = "dashed",
                        colour = "grey50") +
    ggplot2::scale_colour_gradient(low = "#3b5bbf", high = "#d62828",
                                   name = "-log10(FDR)") +
    ggplot2::scale_size_continuous(name = "Gene-set size", range = c(2.5, 9)) +
    ggplot2::facet_grid(OntologyLabel ~ ., scales = "free_y",
                        space = "free_y") +
    ggplot2::labs(x = "Normalised enrichment score (NES)", y = NULL,
                  title = "GSEA enrichment by ontology") +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(
      strip.text = ggplot2::element_text(face = "bold"),
      axis.text.y = ggplot2::element_text(size = 9),
      plot.title = ggplot2::element_text(face = "bold")
    )
}

#' GSEA enrichment-score (running-score) plot for a single term
#'
#' @param gsea_result run_gsea() output (carries gsea_objects attribute)
#' @param term_id GO ID to plot
#' @return ggplot/patchwork object from enrichplot::gseaplot2(), or NULL
gsea_enrichment_plot <- function(gsea_result, term_id) {
  objs <- attr(gsea_result, "gsea_objects")
  if (is.null(objs) || nrow(gsea_result) == 0) return(NULL)

  ont <- gsea_result$Ontology[match(term_id, gsea_result$ID)]
  if (is.na(ont) || is.null(objs[[ont]])) return(NULL)

  desc <- gsea_result$Description[match(term_id, gsea_result$ID)]
  nes <- gsea_result$NES[match(term_id, gsea_result$ID)]
  fdr <- gsea_result$p.adjust[match(term_id, gsea_result$ID)]

  tryCatch(
    enrichplot::gseaplot2(
      objs[[ont]], geneSetID = term_id, ES_geom = "line",
      title = sprintf("%s\nNES = %.2f, FDR = %.3g", desc, nes, fdr)
    ),
    error = function(e) NULL
  )
}

#' Save a plot to PNG or PDF
#'
#' @param plot A ggplot (or recordedplot / grob)
#' @param file Output path
#' @param fmt "png" or "pdf"
#' @param width,height Inches
save_plot <- function(plot, file, fmt = c("png", "pdf"),
                      width = 10, height = 8) {
  fmt <- match.arg(fmt)
  args <- list(filename = file, plot = plot, device = fmt,
               width = width, height = height, limitsize = FALSE)
  if (fmt == "png") args$dpi <- 200
  do.call(ggplot2::ggsave, args)
}
