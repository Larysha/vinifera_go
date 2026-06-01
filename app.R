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

# Source helpers explicitly (also auto-sourced by Shiny, but this keeps the app
# runnable from any working directory / via testthat).
for (f in c("annotation.R", "io.R", "enrichment.R", "plots.R")) {
  source(file.path("R", f), local = FALSE)
}

# ---- Load annotation ONCE, shared across all sessions -----------------------
ANNOTATION <- load_annotation("data")

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

Source code and issues: [%s](%s).
',
  GENOME_LABEL, GENOME_LABEL, GMT_SOURCE_URL, SOURCE_CODE_URL, SOURCE_CODE_URL)

# =============================================================================
# UI
# =============================================================================
ui <- page_navbar(
  title = "vinifera_go",
  theme = bs_theme(version = 5, bootswatch = "flatly"),
  fillable = FALSE,

  # ---- Home / user guide ----
  nav_panel(
    "Home",
    div(class = "container", style = "max-width: 860px; padding-top: 1rem;",
        markdown(guide_md))
  ),

  # ---- Analysis ----
  nav_panel(
    "Analysis",
    layout_sidebar(
      sidebar = sidebar(
        width = 340,
        title = "Inputs",

        fileInput("gene_file", "Gene list (CSV / TSV / text)",
                  accept = c(".csv", ".tsv", ".txt", ".tab")),
        helpText("One gene per line, or a table with a gene column",
                 "and an optional expression / log2FC column."),

        uiOutput("gene_col_ui"),
        uiOutput("expr_col_ui"),

        radioButtons("bg_mode", "Background (universe)",
                     choices = c(
                       "Full annotated genome (default)" = "default",
                       "Upload custom background" = "custom"),
                     selected = "default"),
        conditionalPanel(
          "input.bg_mode == 'custom'",
          fileInput("bg_file", "Background gene list",
                    accept = c(".csv", ".tsv", ".txt", ".tab")),
          helpText("Recommended: upload only the genes that were testable /",
                   "expressed in your experiment.")
        ),

        checkboxGroupInput("ontologies", "GO categories",
                           choices = ONT_CHOICES,
                           selected = c("BP", "MF", "CC")),

        uiOutput("mode_ui"),

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
          col_widths = c(3, 3, 3, 3),
          downloadButton("dl_csv", "Results (CSV)"),
          downloadButton("dl_dot_png", "Dot plot (PNG)"),
          downloadButton("dl_dot_pdf", "Dot plot (PDF)"),
          downloadButton("dl_bar_pdf", "Bar plot (PDF)")
        )
      )
    )
  ),

  # ---- About ----
  nav_panel(
    "About",
    div(class = "container", style = "max-width: 860px; padding-top: 1rem;",
        markdown(about_md))
  )
)

# =============================================================================
# Server
# =============================================================================
server <- function(input, output, session) {

  # ---- Read the uploaded gene table -----------------------------------------
  gene_table <- reactive({
    req(input$gene_file)
    tryCatch(read_gene_table(input$gene_file$datapath),
             error = function(e) {
               showNotification(paste("Could not read file:", e$message),
                                type = "error")
               NULL
             })
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
    selectInput("expr_col", "Expression / log2FC column (optional)",
                choices = c("None" = "", num), selected = "")
  })

  # Analysis mode: ORA only, unless an expression column is chosen, which also
  # unlocks split-by-direction and GSEA.
  output$mode_ui <- renderUI({
    has_expr <- !is.null(input$expr_col) && nzchar(input$expr_col)
    choices <- if (has_expr) {
      c("Over-representation (ORA)" = "ora",
        "ORA split by direction (up / down)" = "split",
        "GSEA (ranked list)" = "gsea")
    } else {
      c("Over-representation (ORA)" = "ora")
    }
    sel <- if (!is.null(input$mode) && input$mode %in% choices) {
      input$mode
    } else {
      "ora"
    }
    radioButtons("mode", "Analysis mode", choices = choices, selected = sel)
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

  # ---- Run enrichment (only the steps that need the test re-run) ------------
  # The "base" result depends on inputs + min/max size + FDR method. Display
  # filters (cutoff, top_n, redundancy) are applied separately and cheaply.
  base_result <- eventReactive(input$run, {
    g <- genes_input()
    validate(need(length(g) > 0, "Please upload a gene list."))
    validate(need(length(input$ontologies) > 0,
                  "Please select at least one GO category."))

    mode <- input$mode %||% "ora"
    bg <- background_input()

    withProgress(message = "Running enrichment...", value = 0.5, {
      if (identical(mode, "gsea")) {
        expr <- expr_input()
        validate(need(!is.null(expr),
                      "GSEA needs an expression / log2FC column."))
        ranked <- stats::setNames(expr, g)
        ranked <- ranked[!is.na(ranked)]
        res <- run_gsea(ranked, ANNOTATION, ontologies = input$ontologies,
                        p_adjust = input$p_adjust, min_size = input$min_size,
                        max_size = input$max_size)
        list(type = "gsea", data = res)
      } else if (identical(mode, "split")) {
        expr <- expr_input()
        validate(need(!is.null(expr),
                      "Split mode needs an expression / log2FC column."))
        res <- run_ora_split(g, expr, ANNOTATION,
                             ontologies = input$ontologies, background = bg,
                             p_adjust = input$p_adjust,
                             min_size = input$min_size,
                             max_size = input$max_size)
        list(type = "split", data = res)
      } else {
        res <- run_ora(g, ANNOTATION, ontologies = input$ontologies,
                       background = bg, p_adjust = input$p_adjust,
                       min_size = input$min_size, max_size = input$max_size)
        list(type = "ora", data = res)
      }
    })
  })

  # ---- Apply cheap display filters: FDR cutoff + redundancy -----------------
  filtered_result <- reactive({
    br <- base_result()
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
    br <- tryCatch(base_result(), error = function(e) NULL)
    !is.null(br) && identical(br$type, "gsea") && nrow(br$data) > 0
  })
  outputOptions(output, "has_gsea", suspendWhenHidden = FALSE)

  # ---- Results table --------------------------------------------------------
  table_data <- reactive({
    res <- filtered_result()
    df <- res$data
    if (nrow(df) == 0) return(df)

    if (res$type == "gsea") {
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
  output$dl_csv <- downloadHandler(
    filename = function() {
      paste0("vinifera_go_results_", Sys.Date(), ".csv")
    },
    content = function(file) {
      utils::write.csv(table_data(), file, row.names = FALSE)
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
