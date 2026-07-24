# microplotr Shiny app: upload -> validate -> plot.
#
# Dev use: devtools::load_all() from the package root, then
#   shiny::runApp("inst/shiny")
# Installed use: library(microplotr); shiny::runApp(system.file("shiny", package = "microplotr"))

library(shiny)
library(bslib)

# ---- shared helpers -------------------------------------------------------

.read_tsv_safe <- function(fileinfo) {
  if (is.null(fileinfo)) return(NULL)
  readr::read_delim(fileinfo$datapath, delim = "\t", show_col_types = FALSE)
}

.validation_ui <- function(report) {
  if (is.null(report)) {
    return(div(class = "alert alert-info", "Upload the required file(s) to validate."))
  }
  errs <- mp_errors(report)
  warns <- mp_warnings(report)
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

taxonomy_ui <- function(id) {
  ns <- NS(id)
  sidebarLayout(
    sidebarPanel(
      fileInput(ns("feature_table"), "Feature table (.tsv)", accept = ".tsv"),
      fileInput(ns("taxonomy"), "Taxonomy (.tsv)", accept = ".tsv"),
      fileInput(ns("metadata"), "Metadata (.tsv)", accept = ".tsv"),
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

    data <- reactive({
      req(input$feature_table, input$taxonomy, input$metadata)
      list(
        feature_table = .read_tsv_safe(input$feature_table),
        taxonomy = .read_tsv_safe(input$taxonomy),
        metadata = .read_tsv_safe(input$metadata)
      )
    })

    report <- reactive(mp_validate(data()))
    is_valid <- reactive(mp_is_valid(report()))

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
          numericInput(ns("top_n"), "Top N taxa (blank = all)", value = 10, min = 1)
        ),
        conditionalPanel(
          condition = sprintf("input['%s'] == 'gradient'", ns("plot_type")),
          selectInput(ns("gradient_var"), "Gradient variable (numeric)", choices = nmc)
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

      switch(input$plot_type,
        barplot = mp_taxa_barplot(d, rank = input$rank, group_rank = gr, top_n = input$top_n),
        heatmap = mp_taxa_heatmap(d, rank = input$rank, top_n = input$top_n),
        bubbleplot = mp_taxa_bubbleplot(d, rank = input$rank, group_rank = gr, top_n = input$top_n),
        gradient = {
          req(input$gradient_var)
          mp_asv_gradient_plot(d, gradient_var = input$gradient_var, rank = input$rank,
                                group_rank = gr, top_n = input$top_n)
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
        treemap = mp_taxa_treemap(d, rank = input$rank, group_rank = gr, top_n = input$top_n)
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
      fileInput(ns("gene_counts"), "Gene/KO count table (.tsv)", accept = ".tsv"),
      fileInput(ns("annotation"), "Function annotation (.tsv)", accept = ".tsv"),
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
      list(feature_table = .read_tsv_safe(input$gene_counts), taxonomy = .read_tsv_safe(input$annotation))
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
      fileInput(ns("mag_table"), "MAG quality table (.tsv, CheckM/CheckM2)", accept = ".tsv"),
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
      .read_tsv_safe(input$mag_table)
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
        fileInput(ns("contig_lengths"), "Contig lengths (.tsv)", accept = ".tsv")
      ),
      conditionalPanel(
        condition = sprintf("input['%s'] == 'summary'", ns("plot_type")),
        fileInput(ns("assembly_summary"), "Assembly summary (.tsv, QUAST report)", accept = ".tsv")
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

    contigs <- reactive({ req(input$contig_lengths); .read_tsv_safe(input$contig_lengths) })
    summary_tbl <- reactive({ req(input$assembly_summary); .read_tsv_safe(input$assembly_summary) })

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
