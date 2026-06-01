# vinifera_go

A lightweight Shiny app for gene ontology (GO) enrichment analysis on the latest
telomere-to-telomere *Vitis vinifera* assembly, PN40024 T2T v5.1 (Grapedia
annotation, gene IDs of the form `Vitvi05_01chr01g00350`).

It exists because ShinyGO only carries the older 12X / V1 annotation and cannot map
current T2T v5.1 gene IDs. This app does the equivalent analysis against the current
assembly, in the browser.


## What it does

Over-representation analysis (ORA) via `clusterProfiler::enricher()` is the default.
When the user supplies an expression or log2 fold-change column, two further modes
become available: ORA split into up- and down-regulated genes, and GSEA on the ranked
list via `clusterProfiler::GSEA()`. Results come as a filterable, downloadable table, a
clusterProfiler-style dot plot, a bar plot, and GSEA running-score plots, with the
plots exportable as PNG or PDF. The adjustable parameters (FDR cutoff and method,
number of pathways shown, minimum and maximum pathway size, and a redundancy filter)
sit beneath the results, in the style of ShinyGO.

## Data provenance

The GO annotations were downloaded from Grapedia as the
[T2T_5.1_blast2go.zip](https://grapedia.org/wp-content/uploads/2024/07/T2T_5.1_blast2go.zip)
release and converted to a gene-set (GMT) file. These are orthology-based (blast2go)
transfers covering roughly 65% of the genome, so not every gene is annotated. The
default background is every annotated gene; users can upload their own background
instead.

## Project layout

```
vinifera_go/
├── app.R                     the Shiny app (UI and server)
├── guide.md                  user-facing landing page (rendered on the Home tab)
├── R/
│   ├── annotation.R          loads the pre-built .rds at start-up
│   ├── io.R                  reads and validates uploaded gene lists
│   ├── enrichment.R          ORA, split, GSEA, redundancy filter
│   └── plots.R               dot, bar, and GSEA plot builders
├── data/                     pre-built annotation (committed and deployed)
│   ├── vitis_v5_term2gene.rds
│   ├── vitis_v5_term2name.rds
│   └── vitis_v5_gene_universe.rds
├── data-raw/
│   ├── build_annotation.R    one-time setup: GMT to .rds (re-runnable)
│   └── blast2go_t2t_5.1.gmt  the GO source (not deployed)
├── tests/testthat/           automated tests
├── deploy.R                  shinyapps.io deployment
└── README.md
```

## Local installation

Requires R (4.2 or newer). The clusterProfiler dependencies live on Bioconductor:

```r
install.packages(c("shiny", "bslib", "DT", "ggplot2", "dplyr",
                   "BiocManager", "testthat", "rsconnect"))
BiocManager::install(c("clusterProfiler", "enrichplot", "DOSE", "GO.db"))
```

`GO.db` is only needed to build the annotation, not to run the app. The ontology of
each term is baked into the `.rds` files at build time, so the running app never loads
`GO.db`. That keeps start-up fast and the deployment small.

## Building the annotation data

The `.rds` files in `data/` are already built and committed. Rebuild them when Grapedia
releases a new annotation:

```bash
# Download the new blast2go annotation, convert it to GMT format
# (GO_ID <tab> name <tab> gene1 <tab> gene2 ...), drop it in data-raw/, then:
Rscript data-raw/build_annotation.R data-raw/your_new_annotation.gmt
```

With no argument it defaults to `data-raw/blast2go_t2t_5.1.gmt`.

## Running the app locally

```r
shiny::runApp(".")
```

A sample gene list is provided at `tests/testthat/sample_gene_list.csv` (50 genes),
with a matching background at `tests/testthat/sample_background_genes.csv`. Input can be
CSV, TSV, or plain text: one gene per line, or a table with a gene column (auto-detected,
a column named `gene_id` or `gene` is preferred) and an optional numeric expression /
log2 fold-change column. Lines starting with `#` are ignored.

## Running the tests

```bash
Rscript -e 'testthat::test_dir("tests/testthat")'
```

The tests check that `run_ora()` recovers a planted over-represented term, that the
redundancy filter drops near-duplicate terms, and that gene-ID validation reports the
format and universe matches correctly. They use a synthetic fixture, so they do not
require the full annotation build.

## Deploying to shinyapps.io

Create a free account at https://www.shinyapps.io and copy the token from your account
settings, then configure rsconnect once:

```r
rsconnect::setAccountInfo(name = "<account>", token = "<token>", secret = "<secret>")
```

Deploy with:

```r
source("deploy.R")
```

`deploy.R` uploads only what the running app needs (the app, `guide.md`, the `R/`
helpers, and the small `.rds` data), excludes the GMT and the tests, and points
rsconnect at the Bioconductor repositories so clusterProfiler resolves on the server.
The bundle is a few MB, well inside the shinyapps.io free-tier limits.

## Notes
The annotation is loaded
once at start-up and shared across sessions, so analyses stay fast and memory use
modest.
