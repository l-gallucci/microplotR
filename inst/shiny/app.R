# microplotr Shiny app: upload -> validate -> plot.
#
# Dev use: devtools::load_all() from the package root, then
#   shiny::runApp("inst/shiny")
# Installed use: library(microplotr); shiny::runApp(system.file("shiny", package = "microplotr"))

library(shiny)
library(bslib)

# Shiny defaults to a 5 MB upload cap -- easily blown past by a real 16S
# feature/taxonomy table (thousands of ASVs, especially with full-sequence
# Feature_IDs from dada2 rather than short OTU/ASV IDs). 500 MB comfortably
# covers that without leaving the limit unbounded.
options(shiny.maxRequestSize = 500 * 1024^2)

# ---- shared helpers -------------------------------------------------------

.UPLOAD_ACCEPT <- c(".tsv", ".csv", ".rds")

.tmp_path_for <- function(fileinfo) {
  # Shiny's uploaded file lands at a random-named temp path (fileinfo$datapath)
  # that does NOT preserve the original extension, but mp_read_table()/
  # mp_read_data() pick their parser (tsv/csv/rds) from the path's extension
  # -- so copy to a temp file named with the *original* (fileinfo$name)
  # extension first.
  if (is.null(fileinfo)) return(NULL)
  ext <- tolower(tools::file_ext(fileinfo$name))
  tmp <- tempfile(fileext = paste0(".", ext))
  file.copy(fileinfo$datapath, tmp, overwrite = TRUE)
  tmp
}

.read_table_safe <- function(fileinfo, id_col = NULL, rename = NULL) {
  # id_col promotes an .rds matrix's rownames into the right ID column
  # (e.g. a phyloseq otu_table()) instead of erroring; rename lets a
  # column under a custom name (e.g. "SampleName") be renamed to the
  # canonical one (e.g. "Sample_ID") the validators/plots expect. Only
  # single-file reads -- the taxonomy module uses mp_read_data() directly
  # instead, since that's also where transposed-feature-table detection
  # lives (it needs all three files' ID sets to check orientation against).
  tmp <- .tmp_path_for(fileinfo)
  if (is.null(tmp)) return(NULL)
  mp_read_table(tmp, id_col = id_col, rename = rename)
}

.effective_id_column <- function(actual, default) {
  if (is.null(actual) || !nzchar(trimws(actual))) default else trimws(actual)
}

.validation_ui <- function(report) {
  if (is.null(report)) {
    return(div(class = "alert alert-info", "Upload the required file(s) to validate."))
  }
  errs <- mp_errors(report)
  warns <- mp_warnings(report)
  # Findings can run long (dada2-style full-sequence IDs, thousands of
  # features) even after mp_validate()'s own per-message truncation -- cap
  # this box's height and let it scroll internally so it can never push the
  # plot-type/rank/etc. controls further down out of the sidebar.
  div(
    style = "max-height: 300px; overflow-y: auto;",
    tagList(
      if (nrow(errs) > 0) {
        div(class = "alert alert-danger",
            tags$strong("Errors (fix before plotting):"),
            tags$ul(lapply(seq_len(nrow(errs)), function(i) {
              tags$li(sprintf("[%s] %s: %s", errs$file[i], errs$field[i], errs$message[i]))
            })))
      },
      if (nrow(warns) > 0) {
        div(class = "alert alert-warning",
            tags$strong("Warnings (plotting still allowed):"),
            tags$ul(lapply(seq_len(nrow(warns)), function(i) {
              tags$li(sprintf("[%s] %s: %s", warns$file[i], warns$field[i], warns$message[i]))
            })))
      },
      if (nrow(errs) == 0 && nrow(warns) == 0) div(class = "alert alert-success", "Input looks good.")
    )
  )
}

.numeric_cols <- function(df, exclude) {
  cols <- setdiff(names(df), exclude)
  Filter(function(col) !any(is.na(suppressWarnings(as.numeric(df[[col]])))), cols)
}

`%||%` <- function(a, b) if (is.null(a) || identical(a, "(none)")) b else a

.download_controls_ui <- function(ns) {
  tagList(
    fluidRow(
      column(6, numericInput(ns("dl_width"), "Width (in)", value = 8, min = 2, step = 0.5)),
      column(6, numericInput(ns("dl_height"), "Height (in)", value = 6, min = 2, step = 0.5))
    ),
    radioButtons(ns("dl_format"), "Format", choices = c("PNG" = "png", "SVG" = "svg"), inline = TRUE),
    downloadButton(ns("download"), "Download plot")
  )
}

.make_download_handler <- function(input, plot_obj_fn, filename_prefix_fn) {
  downloadHandler(
    filename = function() paste0(filename_prefix_fn(), ".", input$dl_format),
    content = function(file) {
      ggplot2::ggsave(file, plot = plot_obj_fn(), width = input$dl_width, height = input$dl_height,
                       dpi = 300, device = input$dl_format)
    }
  )
}

# ---- Taxonomy module --------------------------------------------------------

.parse_exclude_list <- function(txt) {
  if (is.null(txt) || trimws(txt) == "") return(NULL)
  terms <- trimws(strsplit(txt, ",")[[1]])
  terms[nzchar(terms)]
}

.filter_controls_ui <- function(ns) {
  tagList(
    textInput(ns("exclude_taxa"), "Exclude taxa (comma-separated, any rank)",
              placeholder = "e.g. Chloroplast, Mitochondria"),
    checkboxInput(ns("remove_singletons"), "Remove singleton features (total count = 1)", value = FALSE),
    numericInput(ns("min_total_count"), "Min total read count per feature", value = 0, min = 0),
    uiOutput(ns("filter_summary"))
  )
}

.fmt_pct <- function(x) {
  # A fixed 1-decimal format rounds anything under 0.05% down to "0.0%",
  # which reads as "nothing happened" even when features really were
  # removed (contaminants/singletons are often a tiny fraction of total
  # reads by nature, not a sign the filter silently no-op'd).
  if (x == 0) return("0%")
  if (x < 1) return(sprintf("%.3f%%", x))
  sprintf("%.1f%%", x)
}

.filter_summary_ui <- function(summary) {
  total <- summary[summary$reason == "total", ]
  if (total$n_features_removed == 0) {
    return(div(class = "alert alert-secondary", "No features removed."))
  }
  div(
    class = "alert alert-info",
    sprintf(
      "Removed %d feature(s) (%s of features, %s of reads).",
      total$n_features_removed, .fmt_pct(total$pct_features_removed), .fmt_pct(total$pct_reads_removed)
    )
  )
}

taxonomy_ui <- function(id) {
  ns <- NS(id)
  sidebarLayout(
    sidebarPanel(
      fileInput(ns("feature_table"), "Feature table (.tsv/.csv/.rds)", accept = .UPLOAD_ACCEPT),
      fileInput(ns("taxonomy"), "Taxonomy (.tsv/.csv/.rds)", accept = .UPLOAD_ACCEPT),
      fileInput(ns("metadata"), "Metadata (.tsv/.csv/.rds)", accept = .UPLOAD_ACCEPT),
      textInput(ns("feature_id_column"), "Feature ID column name (if not \"Feature_ID\")",
                placeholder = "Feature_ID"),
      textInput(ns("sample_id_column"), "Sample ID column name (if not \"Sample_ID\")",
                placeholder = "Sample_ID"),
      selectInput(ns("feature_table_orientation"), "Feature table orientation",
                  choices = c("Auto-detect" = "auto",
                              "Features are rows, samples are columns" = "asis",
                              "Samples are rows, features are columns (transposed)" = "transpose"),
                  selected = "auto"),
      hr(),
      .filter_controls_ui(ns),
      hr(),
      uiOutput(ns("validation")),
      hr(),
      uiOutput(ns("controls"))
    ),
    mainPanel(
      plotOutput(ns("plot"), height = "600px"),
      .download_controls_ui(ns)
    )
  )
}

taxonomy_server <- function(id) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns

    raw_data <- reactive({
      req(input$feature_table, input$taxonomy, input$metadata)
      transpose_override <- switch(input$feature_table_orientation %||% "auto",
                                    asis = FALSE, transpose = TRUE, NULL)
      mp_read_data(
        .tmp_path_for(input$feature_table),
        .tmp_path_for(input$taxonomy),
        .tmp_path_for(input$metadata),
        feature_id_column = .effective_id_column(input$feature_id_column, "Feature_ID"),
        sample_id_column = .effective_id_column(input$sample_id_column, "Sample_ID"),
        feature_table_transpose = transpose_override
      )
    })

    filtered <- reactive({
      mp_filter_features(
        raw_data(),
        exclude_taxa = .parse_exclude_list(input$exclude_taxa),
        remove_singletons = isTRUE(input$remove_singletons),
        min_total_count = input$min_total_count %||% 0
      )
    })

    data <- reactive(filtered()$data)

    output$filter_summary <- renderUI({
      req(raw_data())
      .filter_summary_ui(filtered()$summary)
    })

    report <- reactive(mp_validate(data()))
    is_valid <- reactive(mp_is_valid(report()))

    # Melt+join data() into the long (Feature_ID, Sample_ID, Count, ranks,
    # metadata) shape every taxa plot needs -- once per validated upload,
    # not once per rank/top_n/plot-type tweak. Each mp_taxa_*() call below
    # accepts this precomputed table via `long =` instead of re-deriving it
    # internally, which is what actually makes the caching take effect --
    # the biggest win once a dataset runs into the thousands of features.
    # Internal (non-exported) helper, hence `:::`.
    long_data <- reactive({
      req(is_valid())
      microplotr:::.mp_long_abundance(data())
    })

    output$validation <- renderUI({
      if (is.null(input$feature_table) || is.null(input$taxonomy) || is.null(input$metadata)) {
        return(.validation_ui(NULL))
      }
      .validation_ui(report())
    })

    tax_cols <- reactive(setdiff(names(data()$taxonomy), "Feature_ID"))
    meta_cols <- reactive(setdiff(names(data()$metadata), "Sample_ID"))
    numeric_meta_cols <- reactive(.numeric_cols(data()$metadata, "Sample_ID"))

    output$controls <- renderUI({
      req(is_valid())
      validate(need(is_valid(), "Fix validation errors above first."))

      tc <- tax_cols()
      mc <- meta_cols()
      nmc <- numeric_meta_cols()

      tagList(
        selectInput(ns("plot_type"), "Plot type", choices = c(
          "Barplot" = "barplot", "Heatmap" = "heatmap", "Bubbleplot" = "bubbleplot",
          "Gradient" = "gradient", "Alpha diversity" = "alpha", "Beta diversity" = "beta",
          "Treemap" = "treemap"
        )),
        conditionalPanel(
          condition = sprintf("input['%s'] != 'alpha' && input['%s'] != 'beta'", ns("plot_type"), ns("plot_type")),
          selectInput(ns("rank"), "Rank", choices = tc, selected = if ("Genus" %in% tc) "Genus" else tc[1]),
          selectInput(ns("group_rank"), "Group rank (upper level)", choices = c("(none)", tc),
                      selected = if ("Phylum" %in% tc) "Phylum" else "(none)"),
          numericInput(ns("top_n"), "Top N taxa (blank = all)", value = 10, min = 1),
          numericInput(ns("min_rel_abund"), "Min mean relative abundance (%)", value = 0, min = 0, step = 0.1),
          numericInput(ns("min_prevalence"), "Min prevalence (fraction of samples, 0-1)",
                       value = 0, min = 0, max = 1, step = 0.05),
          numericInput(ns("detection"), "Detection threshold for prevalence (rel. abundance %)",
                       value = 0, min = 0, step = 0.1)
        ),
        conditionalPanel(
          condition = sprintf("input['%s'] == 'gradient'", ns("plot_type")),
          selectInput(ns("gradient_var"), "Gradient variable (numeric)", choices = nmc),
          selectizeInput(ns("gradient_taxa"), "Specific taxa to plot (optional, overrides Top N)",
                          choices = {
                            current_rank <- input$rank %||% (if ("Genus" %in% tc) "Genus" else tc[1])
                            tax <- data()$taxonomy
                            if (current_rank %in% names(tax)) {
                              sort(unique(stats::na.omit(tax[[current_rank]])))
                            } else {
                              character(0)
                            }
                          },
                          multiple = TRUE,
                          options = list(placeholder = "Leave empty to use Top N ranking"))
        ),
        conditionalPanel(
          condition = sprintf("input['%s'] == 'alpha' || input['%s'] == 'beta'", ns("plot_type"), ns("plot_type")),
          selectInput(ns("group_var"), "Group variable", choices = mc)
        ),
        conditionalPanel(
          condition = sprintf("input['%s'] == 'beta'", ns("plot_type")),
          selectInput(ns("method"), "Dissimilarity", choices = c("bray", "jaccard")),
          selectInput(ns("ordination"), "Ordination", choices = c("pcoa", "nmds"))
        )
      )
    })

    plot_obj <- reactive({
      req(is_valid(), input$plot_type)
      d <- data()
      gr <- input$group_rank %||% NULL
      mra <- input$min_rel_abund %||% 0
      mp <- input$min_prevalence %||% 0
      det <- input$detection %||% 0

      switch(input$plot_type,
        barplot = mp_taxa_barplot(d, rank = input$rank, group_rank = gr, top_n = input$top_n,
                                   min_rel_abund = mra, min_prevalence = mp, detection = det,
                                   long = long_data()),
        heatmap = mp_taxa_heatmap(d, rank = input$rank, top_n = input$top_n,
                                   min_rel_abund = mra, min_prevalence = mp, detection = det,
                                   long = long_data()),
        bubbleplot = mp_taxa_bubbleplot(d, rank = input$rank, group_rank = gr, top_n = input$top_n,
                                         min_rel_abund = mra, min_prevalence = mp, detection = det,
                                         long = long_data()),
        gradient = {
          req(input$gradient_var)
          gradient_taxa <- if (length(input$gradient_taxa) > 0) input$gradient_taxa else NULL
          mp_asv_gradient_plot(d, gradient_var = input$gradient_var, rank = input$rank,
                                group_rank = gr, top_n = input$top_n, taxa = gradient_taxa,
                                min_rel_abund = mra, min_prevalence = mp, detection = det,
                                long = long_data())
        },
        alpha = {
          req(input$group_var)
          mp_alpha_diversity_plot(d, group_var = input$group_var)
        },
        beta = {
          req(input$group_var)
          mp_beta_diversity_plot(d, group_var = input$group_var, method = input$method,
                                  ordination = input$ordination)
        },
        treemap = mp_taxa_treemap(d, rank = input$rank, group_rank = gr, top_n = input$top_n,
                                   min_rel_abund = mra, min_prevalence = mp, detection = det,
                                   long = long_data())
      )
    })

    output$plot <- renderPlot(plot_obj())

    output$download <- .make_download_handler(input, plot_obj, function() paste0("microplotr_", input$plot_type))
  })
}

# ---- Functional profile module ---------------------------------------------

function_ui <- function(id) {
  ns <- NS(id)
  sidebarLayout(
    sidebarPanel(
      fileInput(ns("gene_counts"), "Gene/KO count table (.tsv/.csv/.rds)", accept = .UPLOAD_ACCEPT),
      fileInput(ns("annotation"), "Function annotation (.tsv/.csv/.rds)", accept = .UPLOAD_ACCEPT),
      uiOutput(ns("validation")),
      hr(),
      uiOutput(ns("controls"))
    ),
    mainPanel(
      plotOutput(ns("plot"), height = "600px"),
      .download_controls_ui(ns)
    )
  )
}

function_server <- function(id) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns

    tables <- reactive({
      req(input$gene_counts, input$annotation)
      list(
        feature_table = .read_table_safe(input$gene_counts, id_col = "Feature_ID"),
        taxonomy = .read_table_safe(input$annotation, id_col = "Feature_ID")
      )
    })

    report <- reactive({
      t <- tables()
      stub_metadata <- data.frame(Sample_ID = setdiff(names(t$feature_table), "Feature_ID"))
      mp_validate(list(feature_table = t$feature_table, taxonomy = t$taxonomy, metadata = stub_metadata),
                  required_ranks = character(0))
    })
    is_valid <- reactive(mp_is_valid(report()))

    output$validation <- renderUI({
      if (is.null(input$gene_counts) || is.null(input$annotation)) return(.validation_ui(NULL))
      .validation_ui(report())
    })

    annot_cols <- reactive(setdiff(names(tables()$taxonomy), "Feature_ID"))

    output$controls <- renderUI({
      req(is_valid())
      ac <- annot_cols()
      tagList(
        selectInput(ns("plot_type"), "Plot type", choices = c(
          "Barplot" = "barplot", "Heatmap" = "heatmap", "Treemap" = "treemap"
        )),
        selectInput(ns("rank"), "Rank", choices = ac, selected = if ("KEGG_ko" %in% ac) "KEGG_ko" else ac[1]),
        selectInput(ns("group_rank"), "Group rank (upper level)", choices = c("(none)", ac),
                    selected = if ("COG_category" %in% ac) "COG_category" else "(none)"),
        numericInput(ns("top_n"), "Top N (blank = all)", value = 10, min = 1)
      )
    })

    plot_obj <- reactive({
      req(is_valid(), input$plot_type)
      t <- tables()
      gr <- input$group_rank %||% NULL
      nested_legend <- !is.null(gr)

      switch(input$plot_type,
        barplot = mp_function_barplot(t$feature_table, t$taxonomy, rank = input$rank, group_rank = gr,
                                       top_n = input$top_n, nested_legend = nested_legend),
        heatmap = mp_function_heatmap(t$feature_table, t$taxonomy, rank = input$rank, top_n = input$top_n),
        treemap = mp_function_treemap(t$feature_table, t$taxonomy, rank = input$rank, group_rank = gr,
                                       top_n = input$top_n)
      )
    })

    output$plot <- renderPlot(plot_obj())
    output$download <- .make_download_handler(input, plot_obj, function() paste0("microplotr_function_", input$plot_type))
  })
}

# ---- MAG quality module -----------------------------------------------------

mag_ui <- function(id) {
  ns <- NS(id)
  sidebarLayout(
    sidebarPanel(
      fileInput(ns("mag_table"), "MAG quality table (.tsv/.csv/.rds, CheckM/CheckM2)", accept = .UPLOAD_ACCEPT),
      uiOutput(ns("validation")),
      hr(),
      uiOutput(ns("controls"))
    ),
    mainPanel(
      plotOutput(ns("plot"), height = "600px"),
      .download_controls_ui(ns)
    )
  )
}

mag_server <- function(id) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns

    mag_table <- reactive({
      req(input$mag_table)
      .read_table_safe(input$mag_table, id_col = "Name")
    })
    report <- reactive(mp_validate_mag(mag_table()))
    is_valid <- reactive(mp_is_valid(report()))

    output$validation <- renderUI({
      if (is.null(input$mag_table)) return(.validation_ui(NULL))
      .validation_ui(report())
    })

    mag_cols <- reactive(setdiff(names(mag_table()), c("Name", "Completeness", "Contamination")))

    output$controls <- renderUI({
      req(is_valid())
      mc <- mag_cols()
      tagList(
        selectInput(ns("plot_type"), "Plot type", choices = c(
          "Completeness vs contamination" = "scatter", "Quality distribution" = "distribution"
        )),
        conditionalPanel(
          condition = sprintf("input['%s'] == 'scatter'", ns("plot_type")),
          selectInput(ns("size_col"), "Size by", choices = c("(none)", mc),
                      selected = if ("Genome_Size" %in% mc) "Genome_Size" else "(none)"),
          selectInput(ns("color_col"), "Color by", choices = c("(none)", mc))
        )
      )
    })

    plot_obj <- reactive({
      req(is_valid(), input$plot_type)
      mt <- mag_table()
      if (input$plot_type == "scatter") {
        mp_mag_quality_plot(mt, size_col = input$size_col %||% NULL, color_col = input$color_col %||% NULL)
      } else {
        mp_mag_quality_distribution(mt)
      }
    })

    output$plot <- renderPlot(plot_obj())
    output$download <- .make_download_handler(input, plot_obj, function() paste0("microplotr_mag_", input$plot_type))
  })
}

# ---- Assembly QC module ------------------------------------------------------

assembly_ui <- function(id) {
  ns <- NS(id)
  sidebarLayout(
    sidebarPanel(
      selectInput(ns("plot_type"), "Plot type", choices = c(
        "Nx curve" = "nx", "Summary statistic" = "summary"
      )),
      conditionalPanel(
        condition = sprintf("input['%s'] == 'nx'", ns("plot_type")),
        fileInput(ns("contig_lengths"), "Contig lengths (.tsv/.csv/.rds)", accept = .UPLOAD_ACCEPT)
      ),
      conditionalPanel(
        condition = sprintf("input['%s'] == 'summary'", ns("plot_type")),
        fileInput(ns("assembly_summary"), "Assembly summary (.tsv/.csv/.rds, QUAST report)", accept = .UPLOAD_ACCEPT)
      ),
      uiOutput(ns("validation")),
      hr(),
      uiOutput(ns("controls"))
    ),
    mainPanel(
      plotOutput(ns("plot"), height = "600px"),
      .download_controls_ui(ns)
    )
  )
}

assembly_server <- function(id) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns

    contigs <- reactive({ req(input$contig_lengths); .read_table_safe(input$contig_lengths) })
    summary_tbl <- reactive({ req(input$assembly_summary); .read_table_safe(input$assembly_summary, id_col = "Assembly_ID") })

    report <- reactive({
      if (input$plot_type == "nx") mp_validate_contig_lengths(contigs()) else mp_validate_assembly_summary(summary_tbl())
    })
    is_valid <- reactive(mp_is_valid(report()))

    output$validation <- renderUI({
      if (input$plot_type == "nx" && is.null(input$contig_lengths)) return(.validation_ui(NULL))
      if (input$plot_type == "summary" && is.null(input$assembly_summary)) return(.validation_ui(NULL))
      .validation_ui(report())
    })

    output$controls <- renderUI({
      req(is_valid())
      if (input$plot_type == "summary") {
        cols <- setdiff(names(summary_tbl()), "Assembly_ID")
        selectInput(ns("stat_col"), "Statistic", choices = cols, selected = if ("N50" %in% cols) "N50" else cols[1])
      } else {
        NULL
      }
    })

    plot_obj <- reactive({
      req(is_valid())
      if (input$plot_type == "nx") {
        mp_assembly_nx_plot(contigs())
      } else {
        req(input$stat_col)
        mp_assembly_summary_barplot(summary_tbl(), stat_col = input$stat_col)
      }
    })

    output$plot <- renderPlot(plot_obj())
    output$download <- .make_download_handler(input, plot_obj, function() paste0("microplotr_assembly_", input$plot_type))
  })
}

# ---- app --------------------------------------------------------------------

ui <- page_navbar(
  title = "microplotr",
  theme = bs_theme(bootswatch = "flatly"),
  nav_panel("Taxonomy", taxonomy_ui("tax")),
  nav_panel("Functional profile", function_ui("func")),
  nav_panel("MAG quality", mag_ui("mag")),
  nav_panel("Assembly QC", assembly_ui("asm"))
)

server <- function(input, output, session) {
  taxonomy_server("tax")
  function_server("func")
  mag_server("mag")
  assembly_server("asm")
}

shinyApp(ui, server)
