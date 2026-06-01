library(testthat)

# Locate and source the app's R/ helpers regardless of the working directory
# testthat uses when running this file.
.find_R_dir <- function() {
  for (p in c("R", "../../R", "../R")) {
    if (dir.exists(p) && file.exists(file.path(p, "enrichment.R"))) {
      return(normalizePath(p))
    }
  }
  stop("Could not locate the R/ directory")
}
.r_dir <- .find_R_dir()
source(file.path(.r_dir, "enrichment.R"))
source(file.path(.r_dir, "io.R"))

# ---- A small synthetic annotation where one term is guaranteed enriched -----
make_fixture <- function() {
  universe <- sprintf("Vitvi05_01chr01g%05d", 1:1000)

  # Target term: 30 genes, of which 25 are in our gene-of-interest list
  target_genes <- universe[1:30]
  # A few background terms with no special overlap
  other_terms <- lapply(1:20, function(i) {
    sample(universe, 25)
  })

  term2gene <- rbind(
    data.frame(term = "GO:TEST01", gene = target_genes,
               stringsAsFactors = FALSE),
    do.call(rbind, Map(function(g, i) {
      data.frame(term = sprintf("GO:BG%03d", i), gene = g,
                 stringsAsFactors = FALSE)
    }, other_terms, seq_along(other_terms)))
  )

  terms <- unique(term2gene$term)
  term2name <- data.frame(
    term = terms,
    name = paste("synthetic term", terms),
    ontology = "BP",
    stringsAsFactors = FALSE
  )

  list(term2gene = term2gene, term2name = term2name, universe = universe)
}

test_that("run_ora recovers a planted over-represented term", {
  set.seed(42)
  annotation <- make_fixture()

  # 25 of the 30 target genes, plus 25 unrelated genes
  genes <- c(annotation$universe[1:25], annotation$universe[500:524])

  res <- run_ora(genes, annotation, ontologies = "BP", min_size = 5,
                 max_size = 500, p_adjust = "BH")

  expect_s3_class(res, "data.frame")
  expect_true(nrow(res) >= 1)
  expect_true("GO:TEST01" %in% res$ID)

  hit <- res[res$ID == "GO:TEST01", ]
  expect_lt(hit$p.adjust, 0.05)
  expect_equal(hit$Count, 25)            # 25 of the input genes are in the term
  expect_true(all(c("Ontology", "GeneRatioNum", "BgRatioNum") %in% names(res)))
})

test_that("run_ora returns an empty frame when nothing is enriched", {
  annotation <- make_fixture()
  # Genes drawn from a region with no concentrated annotation
  genes <- annotation$universe[800:810]
  res <- run_ora(genes, annotation, ontologies = "BP", min_size = 5)
  expect_s3_class(res, "data.frame")
  # Either zero rows, or no term passes a strict threshold
  if (nrow(res) > 0) expect_true(all(res$ID != "GO:TEST01"))
})

test_that("remove_redundancy drops a near-duplicate term", {
  df <- data.frame(
    ID = c("GO:A", "GO:B"),
    Description = c("a", "b"),
    p.adjust = c(0.001, 0.01),
    geneID = c("g1/g2/g3/g4/g5", "g1/g2/g3/g4"),  # B ⊂ A
    stringsAsFactors = FALSE
  )
  out <- remove_redundancy(df, overlap_cutoff = 0.5)
  expect_equal(nrow(out), 1)
  expect_equal(out$ID, "GO:A")  # the more significant term is kept
})

test_that("validate_gene_ids reports format and universe match correctly", {
  universe <- c("Vitvi05_01chr01g00001", "Vitvi05_01chr01g00002")
  genes <- c("Vitvi05_01chr01g00001",   # valid + annotated
             "Vitvi05_01chr09g09999",   # valid format, not annotated
             "AT1G01010")               # wrong format
  v <- validate_gene_ids(genes, universe)

  expect_equal(v$n, 3)
  expect_equal(v$n_format_ok, 2)
  expect_equal(v$n_in_universe, 1)
  expect_true("AT1G01010" %in% v$bad_format)
  expect_true("Vitvi05_01chr09g09999" %in% v$not_annotated)
})
