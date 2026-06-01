#!/usr/bin/env Rscript
################################################################################
# Deploy vinifera_go to shinyapps.io
#
# Prerequisites (run once, interactively):
#   install.packages("rsconnect")
#   rsconnect::setAccountInfo(name = "<account>", token = "<token>",
#                             secret = "<secret>")
#
# Then deploy with:
#   source("deploy.R")
################################################################################

library(rsconnect)

# clusterProfiler / enrichplot / DOSE / are on Bioconductor; tell rsconnect where
# to find them so the server can install the dependency tree.
options(repos = BiocManager::repositories())

# Upload only what the running app needs. The 2.4 MB GMT (data-raw/), the tests
# and any scratch files are deliberately excluded to keep the bundle small.
app_files <- c(
  "app.R",
  "guide.md",
  "README.md",
  list.files("R", full.names = TRUE),
  file.path("www", "favicon.png"),
  file.path("data", c("vitis_v5_term2gene.rds",
                      "vitis_v5_term2name.rds",
                      "vitis_v5_gene_universe.rds",
                      "example_genes.csv"))
)

rsconnect::deployApp(
  appDir   = ".",
  appName  = "vinifera_go",
  appTitle = "vinifera_go - GO enrichment for Vitis vinifera T2T v5.1",
  appFiles = app_files,
  forceUpdate = TRUE
)
