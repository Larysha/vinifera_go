################################################################################
# Input parsing and validation
#
# Reads uploaded gene lists (CSV / TSV / whitespace / plain one-per-line),
# tolerating comment lines, and validates IDs against the T2T v5.1 format and
# the annotated universe.
################################################################################

# Gene IDs look like Vitvi05_01chr01g00350. Validation is anchored on the
# assembly prefix so scaffold / unplaced genes still pass.
VITIS_V5_ID_PATTERN <- "^Vitvi05_01"

#' Read an uploaded gene table
#'
#' Auto-detects the delimiter, skips leading comment lines (starting with #),
#' and returns a data.frame. A single-column / one-gene-per-line file yields a
#' one-column frame. Does not assume a header.
#'
#' @param path File path
#' @return data.frame with at least one column
read_gene_table <- function(path) {
  all_lines <- readLines(path, warn = FALSE)
  # Drop comment and blank lines
  keep <- !grepl("^\\s*#", all_lines) & nzchar(trimws(all_lines))
  lines <- all_lines[keep]
  if (length(lines) == 0) stop("The uploaded file contains no data.")

  # Detect delimiter from the first data line
  first <- lines[1]
  delim <- if (grepl("\t", first)) {
    "\t"
  } else if (grepl(",", first)) {
    ","
  } else {
    ""  # whitespace / single column
  }

  txt <- paste(lines, collapse = "\n")
  if (identical(delim, "")) {
    # One token per line: treat as a single unnamed column
    df <- utils::read.table(text = txt, header = FALSE,
                            stringsAsFactors = FALSE, fill = TRUE)
  } else {
    df <- utils::read.table(text = txt, header = TRUE, sep = delim,
                            stringsAsFactors = FALSE, fill = TRUE,
                            quote = "\"", check.names = TRUE)
  }
  df
}

#' Guess which column holds the gene IDs
#'
#' Prefers a column literally named "gene_id" / "gene"; otherwise the column
#' with the most values matching the V5 ID pattern; falls back to the first.
#'
#' @param df data.frame from read_gene_table()
#' @return Column name (character)
guess_gene_column <- function(df) {
  nms <- names(df)
  named <- nms[tolower(nms) %in% c("gene_id", "geneid", "gene", "genes")]
  if (length(named) > 0) return(named[1])

  match_rate <- vapply(df, function(col) {
    col <- as.character(col)
    if (length(col) == 0) return(0)
    mean(grepl(VITIS_V5_ID_PATTERN, col))
  }, numeric(1))
  if (any(match_rate > 0)) return(nms[which.max(match_rate)])
  nms[1]
}

#' Numeric columns that could carry expression / log2FC values
#'
#' @param df data.frame
#' @param exclude Column name(s) to exclude (e.g. the gene column)
#' @return Character vector of candidate column names (may be empty)
numeric_columns <- function(df, exclude = character(0)) {
  cand <- setdiff(names(df), exclude)
  cand[vapply(cand, function(col) {
    suppressWarnings(is.numeric(df[[col]]) ||
      !all(is.na(as.numeric(as.character(df[[col]])))))
  }, logical(1))]
}

#' Validate a set of gene IDs against format and the annotated universe
#'
#' @param genes Character vector of gene IDs
#' @param universe Character vector of all annotated genes
#' @return List of counts / percentages and the unmatched IDs
validate_gene_ids <- function(genes, universe) {
  genes <- unique(genes[nzchar(genes)])
  n <- length(genes)
  fmt_ok <- grepl(VITIS_V5_ID_PATTERN, genes)
  in_universe <- genes %in% universe

  list(
    n              = n,
    n_format_ok    = sum(fmt_ok),
    pct_format_ok  = if (n) 100 * sum(fmt_ok) / n else 0,
    n_in_universe  = sum(in_universe),
    pct_in_universe = if (n) 100 * sum(in_universe) / n else 0,
    bad_format     = genes[!fmt_ok],
    not_annotated  = genes[fmt_ok & !in_universe]
  )
}
