################################################################################
# vinifera_go  -  GO enrichment for the Vitis vinifera T2T v5.1 assembly
#
# A lightweight, single-purpose Shiny app that runs GO over-representation
# analysis (ORA) and optional GSEA on T2T v5.1 gene IDs (Vitvi05_01...), filling
# the gap left by ShinyGO, which only supports the older 12X / V1 annotation.
#
# Annotation data is pre-built by data-raw/build_annotation.R and loaded once at
# start-up. Run locally with shiny::runApp(); deploy with deploy.R.
################################################################################

suppressPackageStartupMessages({
  library(shiny)
  library(bslib)
  library(DT)
  library(ggplot2)
})

# Null-coalescing helper (used throughout the server)
`%||%` <- function(a, b) if (is.null(a)) b else a

# Select + rename result columns for display / download. Shared by the on-screen
# table and both CSV downloads so the filtered and full exports look identical.
format_result_table <- function(df, type) {
  if (is.null(df) || nrow(df) == 0) return(df)

  if (type == "gsea") {
    cols <- c("ID", "Description", "Ontology", "setSize", "NES",
              "enrichmentScore", "pvalue", "p.adjust", "core_enrichment")
    df <- df[, intersect(cols, names(df)), drop = FALSE]
    names(df) <- c("GO ID", "Term", "Ontology", "Set size", "NES",
                   "Enrichment score", "p-value", "Adjusted p",
                   "Core genes")[seq_along(names(df))]
  } else {
    keep <- c("ID", "Description", "Ontology",
              if ("Direction" %in% names(df)) "Direction",
              "Count", "GeneRatio", "BgRatio", "pvalue", "p.adjust", "geneID")
    df <- df[, intersect(keep, names(df)), drop = FALSE]
    nice <- c(ID = "GO ID", Description = "Term", Ontology = "Ontology",
              Direction = "Direction", Count = "Gene count",
              GeneRatio = "Gene ratio", BgRatio = "Bg ratio",
              pvalue = "p-value", p.adjust = "Adjusted p", geneID = "Genes")
    names(df) <- nice[names(df)]
  }
  df
}

# Source helpers explicitly (also auto-sourced by Shiny, but this keeps the app
# runnable from any working directory / via testthat).
for (f in c("annotation.R", "io.R", "enrichment.R", "plots.R")) {
  source(file.path("R", f), local = FALSE)
}

# ---- Load annotation ONCE, shared across all sessions -----------------------
ANNOTATION <- load_annotation("data")

# Example dataset for the "Load example" button: a deduped gene list with
# synthetic log2FC values, so all three modes (ORA, split, GSEA) can be tried.
EXAMPLE_GENES <- tryCatch(
  utils::read.csv("data/example_genes.csv", stringsAsFactors = FALSE),
  error = function(e) NULL
)
EXAMPLE_N <- if (is.null(EXAMPLE_GENES)) 0L else nrow(EXAMPLE_GENES)

# GSEA needs a ranked list spanning many genes; warn below this size.
GSEA_MIN_RECOMMENDED <- 200

GENOME_LABEL <- "T2T v5.1"
GMT_SOURCE_URL <-
  "https://grapedia.org/wp-content/uploads/2024/07/T2T_5.1_blast2go.zip"
SOURCE_CODE_URL <- "https://github.com/Larysha/vinifera_go"

ONT_CHOICES <- c("Biological Process (BP)" = "BP",
                 "Molecular Function (MF)" = "MF",
                 "Cellular Component (CC)" = "CC")

FDR_METHODS <- c("Benjamini-Hochberg (BH)" = "BH",
                 "Bonferroni" = "bonferroni",
                 "Holm" = "holm",
                 "None" = "none")

# ---- Landing-page content: the user guide, with live counts injected --------
# guide.md is the public-facing how-to; {{...}} placeholders are filled in here
# so the term / gene numbers always match the loaded annotation.
guide_md <- local({
  txt <- paste(readLines("guide.md", warn = FALSE), collapse = "\n")
  repl <- c(
    "\\{\\{N_TERMS\\}\\}" = format(ANNOTATION$n_terms, big.mark = ","),
    "\\{\\{N_BP\\}\\}"    = format(ANNOTATION$ontology_counts[["BP"]],
                                   big.mark = ","),
    "\\{\\{N_MF\\}\\}"    = format(ANNOTATION$ontology_counts[["MF"]],
                                   big.mark = ","),
    "\\{\\{N_CC\\}\\}"    = format(ANNOTATION$ontology_counts[["CC"]],
                                   big.mark = ","),
    "\\{\\{N_GENES\\}\\}" = format(ANNOTATION$n_genes, big.mark = ",")
  )
  for (p in names(repl)) txt <- gsub(p, repl[[p]], txt)
  txt
})

about_md <- sprintf('
## About and citation

vinifera_go was built to run GO enrichment on the current *Vitis vinifera*
PN40024 %s assembly, which [ShinyGO](https://bioinformatics.sdstate.edu/go/)
does not support.

Genome assembly: *Vitis vinifera* PN40024 %s.

GO annotations: Grapedia blast2go release, [T2T_5.1_blast2go.zip](%s), converted
to a gene-set (GMT) file. Roughly 65%% genome coverage.

Enrichment is performed with the clusterProfiler R package (Wu et al., 2021,
*The Innovation*), using enricher() for over-representation and GSEA() for gene
set enrichment. Plots use ggplot2 and enrichplot.
',
  GENOME_LABEL, GENOME_LABEL, GMT_SOURCE_URL)

# =============================================================================
# UI
# =============================================================================

# Grapevine palette: a burgundy / purple banner with a light-pink highlight.
GRAPE_BURGUNDY <- "#6d1f44"
GRAPE_PINK     <- "#f6a8cb"

app_theme <- bs_theme(
  version = 5,
  bootswatch = "flatly",
  primary = GRAPE_BURGUNDY,
  secondary = "#f6c6da"
)

# Navbar colours and the pink active / hover highlight (set explicitly so they
# do not depend on how bootswatch computes navbar contrast).
app_css <- sprintf("
  .navbar { background-color: %1$s !important; }
  .navbar .navbar-brand { color: #ffffff !important; font-weight: 600; }
  .navbar .nav-link { color: rgba(255, 255, 255, 0.8) !important; }
  .navbar .nav-link:hover,
  .navbar .nav-link:focus { color: #ffffff !important; }
  .navbar .nav-link.active {
    color: #ffffff !important;
    border-bottom: 3px solid %2$s;
  }
  .nav-tabs .nav-link { color: %1$s !important; }
  .nav-tabs .nav-link.active {
    color: %1$s !important;
    border-bottom: 2px solid %2$s !important;
    font-weight: 600;
  }
  .nav-tabs .nav-link:hover { color: %1$s !important; }
  a { color: #8a2a57; }
  /* Greyed-out (unavailable) analysis-mode options */
  .mode-disabled, .mode-disabled span { color: #adb5bd !important; cursor: not-allowed; }
  .mode-disabled input { cursor: not-allowed; }
  .form-control:focus, .form-select:focus {
    border-color: %2$s;
    box-shadow: 0 0 0 0.2rem rgba(246, 168, 203, 0.5);
  }
  /* DataTables pagination: replace the bootswatch teal with the pink / burgundy
     palette. Bootstrap 5 pagination is driven by CSS variables, so overriding
     them on .pagination is the reliable way to recolour every state. */
  .pagination {
    --bs-pagination-color: %1$s;
    --bs-pagination-bg: #ffffff;
    --bs-pagination-border-color: %2$s;
    --bs-pagination-hover-color: %1$s;
    --bs-pagination-hover-bg: #f6c6da;
    --bs-pagination-hover-border-color: %2$s;
    --bs-pagination-focus-color: %1$s;
    --bs-pagination-focus-bg: #f6c6da;
    --bs-pagination-focus-box-shadow: 0 0 0 0.25rem rgba(246, 168, 203, 0.5);
    --bs-pagination-active-color: %1$s;
    --bs-pagination-active-bg: %2$s;
    --bs-pagination-active-border-color: %2$s;
  }
  /* Fallback for the non-Bootstrap DataTables pagination markup */
  .dataTables_wrapper .dataTables_paginate .paginate_button { color: %1$s !important; }
  .dataTables_wrapper .dataTables_paginate .paginate_button:hover {
    background: #f6c6da !important; color: %1$s !important; border-color: %2$s !important;
  }
  .dataTables_wrapper .dataTables_paginate .paginate_button.current,
  .dataTables_wrapper .dataTables_paginate .paginate_button.current:hover {
    background: %2$s !important; color: %1$s !important; border-color: %2$s !important;
  }
", GRAPE_BURGUNDY, GRAPE_PINK)

ui <- page_navbar(
  title = tags$span(
    tags$img(src = "favicon.png", height = "30px",
             style = "margin-right: 8px; vertical-align: middle;"),
    "vinifera_go"
  ),
  theme = app_theme,
  header = tags$head(
    tags$link(rel = "icon", type = "image/png", href = "favicon.png"),
    tags$link(rel = "apple-touch-icon", href = "favicon.png"),
    tags$style(HTML(app_css))
  ),
  fillable = FALSE,

  # ---- Analysis (landing page) ----
  nav_panel(
    "Analysis",
    layout_sidebar(
      sidebar = sidebar(
        width = 340,
        title = "Inputs",

        uiOutput("gene_file_ui"),
        helpText("One gene per line, or a table with a gene column",
                 "and an optional expression / log2FC column."),

        textAreaInput("gene_paste", "Or paste a gene list",
                      placeholder = "Vitvi05_01chr01g00350\nVitvi05_01chr01g01160\n...",
                      rows = 4, resize = "vertical"),
        helpText("One gene per line. If anything is pasted here it is used",
                 "in place of an uploaded file."),

        actionButton("clear_all", "Clear all inputs",
                     class = "btn-outline-danger btn-sm", width = "100%"),
        helpText("Remove any uploaded file, pasted genes, or example data and",
                 "start from a fresh slate."),

        actionButton("load_example", "Load example data",
                     class = "btn-outline-secondary btn-sm",
                     width = "100%"),
        helpText(sprintf(paste("A %d-gene list with synthetic log2FC values,",
                               "for trying out all three analysis modes."),
                         EXAMPLE_N)),

        uiOutput("gene_col_ui"),
        uiOutput("expr_col_ui"),

        radioButtons("bg_mode", "Background (universe)",
                     choices = c(
                       "Full annotated genome (default)" = "default",
                       "Upload custom background" = "custom"),
                     selected = "default"),
        conditionalPanel(
          "input.bg_mode == 'custom'",
          uiOutput("bg_file_ui"),
          helpText("Recommended: upload only the genes that were testable /",
                   "expressed in your experiment.")
        ),

        checkboxGroupInput("ontologies", "GO categories",
                           choices = ONT_CHOICES,
                           selected = c("BP", "MF", "CC")),

        uiOutput("mode_ui"),
        uiOutput("mode_legend"),

        actionButton("run", "Run enrichment", class = "btn-primary",
                     width = "100%"),
        br(), br(),
        uiOutput("validation_box")
      ),

      # ---- Outputs ----
      navset_card_tab(
        id = "results_tabs",
        nav_panel("Results table", DT::DTOutput("results_table")),
        nav_panel("Dot plot", plotOutput("dotplot", height = "650px")),
        nav_panel("Bar plot", plotOutput("barplot", height = "650px")),
        nav_panel(
          "GSEA",
          conditionalPanel(
            "output.has_gsea == true",
            plotOutput("gsea_dotplot", height = "650px"),
            hr(),
            uiOutput("gsea_term_ui"),
            plotOutput("gsea_running", height = "450px")
          ),
          conditionalPanel(
            "output.has_gsea != true",
            div(class = "text-muted", style = "padding: 1rem;",
                "Provide an expression / log2FC column and select the GSEA",
                "mode to see gene-set enrichment results here.")
          )
        )
      ),

      # ---- ShinyGO-style adjustable parameters, beneath the outputs ----
      card(
        card_header("Parameters (adjust and the results update)"),
        layout_columns(
          col_widths = c(4, 4, 4),
          numericInput("p_cutoff", "FDR (adjusted-p) cutoff",
                       value = 0.05, min = 0, max = 1, step = 0.01),
          selectInput("p_adjust", "FDR method", choices = FDR_METHODS,
                      selected = "BH"),
          numericInput("top_n", "Number of pathways to show",
                       value = 20, min = 1, max = 100, step = 1)
        ),
        layout_columns(
          col_widths = c(4, 4, 4),
          numericInput("min_size", "Min pathway size",
                       value = 10, min = 1, max = 1000, step = 1),
          numericInput("max_size", "Max pathway size",
                       value = 500, min = 1, max = 5000, step = 1),
          div(
            checkboxInput("remove_redundant", "Remove redundant terms",
                          value = FALSE),
            conditionalPanel(
              "input.remove_redundant == true",
              sliderInput("overlap_cutoff", "Redundancy overlap cutoff",
                          min = 0.1, max = 0.95, value = 0.5, step = 0.05)
            )
          )
        ),
        helpText("Note: changing min / max pathway size or the FDR method",
                 "re-runs the enrichment test. The other controls re-filter the",
                 "existing results instantly.")
      ),

      # ---- Downloads ----
      card(
        card_header("Downloads"),
        layout_columns(
          col_widths = c(6, 6),
          downloadButton("dl_csv", "Significant results (CSV)"),
          downloadButton("dl_csv_all", "All results, unfiltered (CSV)")
        ),
        helpText("\"Significant results\" matches the table above (FDR cutoff",
                 "and redundancy filter applied). \"All results\" includes every",
                 "tested GO term with its p-value, even non-significant ones."),
        layout_columns(
          col_widths = c(4, 4, 4),
          downloadButton("dl_dot_png", "Dot plot (PNG)"),
          downloadButton("dl_dot_pdf", "Dot plot (PDF)"),
          downloadButton("dl_bar_pdf", "Bar plot (PDF)")
        )
      )
    )
  ),

  # ---- Docs (user guide + about, appended) ----
  nav_panel(
    "Docs",
    tags$style(HTML("
      .docs-content { font-size: 0.95rem; }
      .docs-content h1 { font-size: 1.5rem; margin-top: 1.5rem; }
      .docs-content h2 { font-size: 1.25rem; margin-top: 1.5rem; }
      .docs-content h3 { font-size: 1.05rem; margin-top: 1.2rem;
                         font-weight: 600; }
      .docs-content h4 { font-size: 0.95rem; font-weight: 600; }
    ")),
    div(class = "container docs-content",
        style = "max-width: 820px; padding-top: 1rem;",
        markdown(paste(guide_md, about_md, sep = "\n\n")))
  )
)

# =============================================================================
# Server
# =============================================================================
server <- function(input, output, session) {

  store <- reactiveValues(example = NULL, reset = 0L)

  # Enrichment output, held in a reactiveVal so the Clear button can blank it.
  results <- reactiveVal(NULL)

  # File inputs are rendered server-side so the Clear button can reset them
  # (bumping store$reset rebuilds them empty).
  output$gene_file_ui <- renderUI({
    store$reset
    fileInput("gene_file", "Gene list (CSV / TSV / text)",
              accept = c(".csv", ".tsv", ".txt", ".tab"))
  })
  output$bg_file_ui <- renderUI({
    store$reset
    fileInput("bg_file", "Background gene list",
              accept = c(".csv", ".tsv", ".txt", ".tab"))
  })

  observeEvent(input$load_example, {
    if (is.null(EXAMPLE_GENES)) {
      showNotification("Example dataset is not available.", type = "error")
      return()
    }
    store$example <- EXAMPLE_GENES
    showNotification(
      sprintf("Loaded example data: %d genes with synthetic log2FC values.",
              nrow(EXAMPLE_GENES)),
      type = "message")
  })

  # Clear everything: uploaded files, pasted genes, example data, and results.
  observeEvent(input$clear_all, {
    store$example <- NULL
    store$reset <- store$reset + 1L
    results(NULL)
    updateTextAreaInput(session, "gene_paste", value = "")
    updateRadioButtons(session, "bg_mode", selected = "default")
    showNotification("Cleared all inputs. Ready for a fresh start.",
                     type = "message")
  })

  # ---- Resolve the gene table -----------------------------------------------
  # Precedence: pasted text, then an uploaded file, then the example dataset.
  gene_table <- reactive({
    pasted <- input$gene_paste
    if (!is.null(pasted) && nzchar(trimws(pasted))) {
      return(tryCatch(read_gene_text(pasted),
                      error = function(e) {
                        showNotification(paste("Could not parse pasted genes:",
                                               e$message), type = "error")
                        NULL
                      }))
    }
    if (!is.null(input$gene_file)) {
      return(tryCatch(read_gene_table(input$gene_file$datapath),
                      error = function(e) {
                        showNotification(paste("Could not read file:",
                                               e$message), type = "error")
                        NULL
                      }))
    }
    store$example  # NULL until the example button is used
  })

  # Column pickers (populated from the uploaded file)
  output$gene_col_ui <- renderUI({
    df <- gene_table(); req(df)
    selectInput("gene_col", "Gene ID column",
                choices = names(df), selected = guess_gene_column(df))
  })

  output$expr_col_ui <- renderUI({
    df <- gene_table(); req(df, input$gene_col)
    num <- numeric_columns(df, exclude = input$gene_col)
    # Auto-select a log2FC / expression-looking column so the split and GSEA
    # modes are offered without the user having to pick it manually.
    guess <- num[grepl("log.?2?.?fc|logfc|fold|expr", tolower(num))]
    sel <- if (length(guess)) guess[1] else ""
    selectInput("expr_col", "Expression / log2FC column (optional)",
                choices = c("None" = "", num), selected = sel)
  })

  # Analysis mode. ORA is always available. Split-by-direction and GSEA need an
  # expression / log2FC column, so when none is selected they stay visible but
  # greyed out (rendered as disabled placeholders) rather than disappearing.
  output$mode_ui <- renderUI({
    has_expr <- !is.null(input$expr_col) && nzchar(input$expr_col)

    if (has_expr) {
      choices <- c("Over-representation (ORA)" = "ora",
                   "ORA split by direction (up / down)" = "split",
                   "GSEA (ranked list)" = "gsea")
      sel <- if (!is.null(input$mode) && input$mode %in% choices) {
        input$mode
      } else {
        "ora"
      }
      radioButtons("mode", "Analysis mode", choices = choices, selected = sel)
    } else {
      # Only ORA is selectable; the other two shown disabled for context.
      disabled_opt <- function(label) {
        div(class = "radio mode-disabled",
            tags$label(
              tags$input(type = "radio", disabled = NA), " ", tags$span(label)))
      }
      tagList(
        radioButtons("mode", "Analysis mode",
                     choices = c("Over-representation (ORA)" = "ora"),
                     selected = "ora"),
        disabled_opt("ORA split by direction (up / down)"),
        disabled_opt("GSEA (ranked list)")
      )
    }
  })

  output$mode_legend <- renderUI({
    helpText(
      "Over-representation (ORA) works on any gene list. ",
      "Split by direction and GSEA require an expression / log2FC column; ",
      "GSEA also needs a long, ranked list (ideally most of the genome), ",
      "so it is best suited to genome-wide results rather than a short hit list."
    )
  })

  # ---- Resolve genes, expression and background -----------------------------
  genes_input <- reactive({
    df <- gene_table(); req(df, input$gene_col)
    g <- as.character(df[[input$gene_col]])
    g[nzchar(g) & !is.na(g)]
  })

  expr_input <- reactive({
    df <- gene_table()
    if (is.null(df) || is.null(input$expr_col) || !nzchar(input$expr_col)) {
      return(NULL)
    }
    suppressWarnings(as.numeric(as.character(df[[input$expr_col]])))
  })

  background_input <- reactive({
    if (identical(input$bg_mode, "default")) return(NULL)
    req(input$bg_file)
    bg <- tryCatch(read_gene_table(input$bg_file$datapath),
                   error = function(e) NULL)
    req(bg)
    col <- guess_gene_column(bg)
    as.character(bg[[col]])
  })

  # ---- Live validation summary ----------------------------------------------
  output$validation_box <- renderUI({
    g <- tryCatch(genes_input(), error = function(e) NULL)
    if (is.null(g) || length(g) == 0) return(NULL)
    v <- validate_gene_ids(g, ANNOTATION$universe)
    warn <- if (v$pct_format_ok < 90) {
      tags$p(class = "text-danger",
             sprintf("Warning: only %.0f%% of IDs match the expected %s format (Vitvi05_01...).",
                     v$pct_format_ok, GENOME_LABEL))
    }
    div(class = "alert alert-light border small",
        tags$strong("Input check"), tags$br(),
        sprintf("%s unique genes provided.", v$n), tags$br(),
        sprintf("%.0f%% match the %s ID format.", v$pct_format_ok,
                GENOME_LABEL), tags$br(),
        sprintf("%.0f%% (%s) are present in the annotated universe.",
                v$pct_in_universe, v$n_in_universe),
        warn)
  })

  # ---- Run enrichment on demand, storing into results() ---------------------
  # The stored result depends on inputs + min/max size + FDR method. Display
  # filters (cutoff, top_n, redundancy) are applied separately and cheaply.
  observeEvent(input$run, {
    g <- tryCatch(genes_input(), error = function(e) NULL)
    if (is.null(g) || length(g) == 0) {
      showNotification("Please provide a gene list first.", type = "error")
      return()
    }
    if (length(input$ontologies) == 0) {
      showNotification("Please select at least one GO category.", type = "error")
      return()
    }

    mode <- input$mode %||% "ora"
    bg <- background_input()

    res <- withProgress(message = "Running enrichment...", value = 0.5, {
      if (identical(mode, "gsea")) {
        expr <- expr_input()
        if (is.null(expr)) {
          showNotification("GSEA needs an expression / log2FC column.",
                           type = "error")
          return(NULL)
        }
        ranked <- stats::setNames(expr, g)
        ranked <- ranked[!is.na(ranked)]
        if (length(ranked) < GSEA_MIN_RECOMMENDED) {
          showNotification(
            sprintf(paste("GSEA works best on a ranked list spanning many",
                          "genes (ideally most of the genome). You provided",
                          "%d, so the result may be unreliable. For a short",
                          "gene list, over-representation analysis (ORA) is the",
                          "better choice."), length(ranked)),
            type = "warning", duration = 12)
        }
        out <- run_gsea(ranked, ANNOTATION, ontologies = input$ontologies,
                        p_adjust = input$p_adjust, min_size = input$min_size,
                        max_size = input$max_size)
        list(type = "gsea", data = out)
      } else if (identical(mode, "split")) {
        expr <- expr_input()
        if (is.null(expr)) {
          showNotification("Split mode needs an expression / log2FC column.",
                           type = "error")
          return(NULL)
        }
        out <- run_ora_split(g, expr, ANNOTATION,
                             ontologies = input$ontologies, background = bg,
                             p_adjust = input$p_adjust,
                             min_size = input$min_size,
                             max_size = input$max_size)
        list(type = "split", data = out)
      } else {
        out <- run_ora(g, ANNOTATION, ontologies = input$ontologies,
                       background = bg, p_adjust = input$p_adjust,
                       min_size = input$min_size, max_size = input$max_size)
        list(type = "ora", data = out)
      }
    })
    results(res)
  })

  # ---- Apply cheap display filters: FDR cutoff + redundancy -----------------
  filtered_result <- reactive({
    br <- results(); req(br)
    df <- br$data
    if (nrow(df) == 0) return(br)

    df <- df[df$p.adjust <= input$p_cutoff, , drop = FALSE]

    if (isTRUE(input$remove_redundant) && br$type != "gsea" && nrow(df) > 1) {
      df <- remove_redundancy(df, overlap_cutoff = input$overlap_cutoff)
    }
    out <- br
    out$data <- df
    out
  })

  output$has_gsea <- reactive({
    br <- results()
    !is.null(br) && identical(br$type, "gsea") && nrow(br$data) > 0
  })
  outputOptions(output, "has_gsea", suspendWhenHidden = FALSE)

  # ---- Results table --------------------------------------------------------
  # Displayed table: respects the FDR cutoff and redundancy filter.
  table_data <- reactive({
    res <- filtered_result()
    format_result_table(res$data, res$type)
  })

  # Every tested term, ignoring the display filters (FDR cutoff, redundancy).
  # Used by the "All results" download so non-significant terms are included.
  full_table_data <- reactive({
    br <- results(); req(br)
    format_result_table(br$data, br$type)
  })

  output$results_table <- DT::renderDT({
    df <- table_data()
    validate(need(nrow(df) > 0, "No enriched terms at the current settings."))
    num_cols <- intersect(c("p-value", "Adjusted p", "NES", "Enrichment score",
                            "Gene ratio", "Bg ratio"), names(df))
    dt <- DT::datatable(df, rownames = FALSE, filter = "top",
                        options = list(pageLength = 25, scrollX = TRUE),
                        class = "compact stripe")
    for (cl in c("p-value", "Adjusted p")) {
      if (cl %in% names(df)) dt <- DT::formatSignif(dt, cl, digits = 3)
    }
    for (cl in intersect(c("NES", "Enrichment score"), names(df))) {
      dt <- DT::formatRound(dt, cl, digits = 3)
    }
    dt
  })

  # ---- Plots ----------------------------------------------------------------
  dot_plot <- reactive({
    res <- filtered_result()
    if (res$type == "gsea") return(NULL)
    go_dotplot(res$data, top_n = input$top_n)
  })
  bar_plot <- reactive({
    res <- filtered_result()
    if (res$type == "gsea") return(NULL)
    go_barplot(res$data, top_n = input$top_n, fdr_cutoff = input$p_cutoff)
  })

  output$dotplot <- renderPlot({
    p <- dot_plot()
    validate(need(!is.null(p), "No results to plot (or GSEA mode is active)."))
    p
  })
  output$barplot <- renderPlot({
    p <- bar_plot()
    validate(need(!is.null(p), "No results to plot (or GSEA mode is active)."))
    p
  })

  # ---- GSEA-specific outputs ------------------------------------------------
  output$gsea_dotplot <- renderPlot({
    res <- filtered_result()
    validate(need(res$type == "gsea" && nrow(res$data) > 0,
                  "No GSEA results."))
    gsea_dotplot(res$data, top_n = input$top_n)
  })

  output$gsea_term_ui <- renderUI({
    res <- filtered_result()
    if (res$type != "gsea" || nrow(res$data) == 0) return(NULL)
    choices <- stats::setNames(res$data$ID,
                               paste0(res$data$Description,
                                      " (", res$data$ID, ")"))
    selectInput("gsea_term", "Running-score plot for term:",
                choices = choices, width = "100%")
  })

  output$gsea_running <- renderPlot({
    res <- filtered_result()
    validate(need(res$type == "gsea" && !is.null(input$gsea_term),
                  "Select a term."))
    p <- gsea_enrichment_plot(res$data, input$gsea_term)
    validate(need(!is.null(p), "Could not draw the enrichment plot."))
    p
  })

  # ---- Downloads ------------------------------------------------------------
  # Significant terms only (what the table currently shows).
  output$dl_csv <- downloadHandler(
    filename = function() {
      paste0("vinifera_go_results_", Sys.Date(), ".csv")
    },
    content = function(file) {
      utils::write.csv(table_data(), file, row.names = FALSE)
    }
  )

  # Every tested term, including non-significant ones (no FDR / redundancy
  # filtering). The full set so users can see weakly-enriched terms too.
  output$dl_csv_all <- downloadHandler(
    filename = function() {
      paste0("vinifera_go_results_all_", Sys.Date(), ".csv")
    },
    content = function(file) {
      utils::write.csv(full_table_data(), file, row.names = FALSE)
    }
  )

  plot_for_download <- function(which) {
    res <- filtered_result()
    if (res$type == "gsea") {
      gsea_dotplot(res$data, top_n = input$top_n)
    } else if (which == "dot") {
      dot_plot()
    } else {
      bar_plot()
    }
  }

  output$dl_dot_png <- downloadHandler(
    filename = function() paste0("vinifera_go_dotplot_", Sys.Date(), ".png"),
    content = function(file) {
      p <- plot_for_download("dot"); req(p)
      save_plot(p, file, "png", width = 11, height = 9)
    }
  )
  output$dl_dot_pdf <- downloadHandler(
    filename = function() paste0("vinifera_go_dotplot_", Sys.Date(), ".pdf"),
    content = function(file) {
      p <- plot_for_download("dot"); req(p)
      save_plot(p, file, "pdf", width = 11, height = 9)
    }
  )
  output$dl_bar_pdf <- downloadHandler(
    filename = function() paste0("vinifera_go_barplot_", Sys.Date(), ".pdf"),
    content = function(file) {
      p <- plot_for_download("bar"); req(p)
      save_plot(p, file, "pdf", width = 11, height = 9)
    }
  )
}

shinyApp(ui, server)
