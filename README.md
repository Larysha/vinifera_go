# <img src="www/favicon.png" alt="vinifera_go logo" height="40" valign="middle"> vinifera_go

A lightweight Shiny app for gene ontology (GO) enrichment analysis on the latest
telomere-to-telomere *Vitis vinifera* assembly, PN40024 T2T v5.1 (Grapedia
annotation, gene IDs of the form `Vitvi05_01chr01g00350`).

It exists because ShinyGO only carries the older 12X / V1 annotation and cannot map
current T2T v5.1 gene IDs. This app does the equivalent analysis against the current
assembly, in the browser.

Live app: https://biopod.shinyapps.io/vinifera_go/


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


## Notes
The annotation is loaded
once at start-up and shared across sessions, so analyses stay fast and memory use
modest.
