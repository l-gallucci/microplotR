example_dir <- function(name) {
  d <- system.file("extdata", name, package = "microplotr")
  if (d == "") d <- testthat::test_path("..", "..", "inst", "extdata", name)
  d
}

load_example <- function(name) {
  d <- example_dir(name)
  mp_read_data(
    file.path(d, "feature_table.tsv"),
    file.path(d, "taxonomy.tsv"),
    file.path(d, "metadata.tsv")
  )
}

test_that("default heatmap (with dendrograms) builds as a patchwork", {
  data <- load_example("example_valid")
  p <- mp_taxa_heatmap(data, rank = "Genus", top_n = 10)
  expect_s3_class(p, "patchwork")
})

test_that("show_dendrogram = FALSE returns a plain ggplot", {
  data <- load_example("example_valid")
  p <- mp_taxa_heatmap(data, rank = "Genus", top_n = 10, show_dendrogram = FALSE)
  expect_s3_class(p, "ggplot")
  expect_false(inherits(p, "patchwork"))
})

test_that("top_n limits the number of rows shown", {
  data <- load_example("example_valid")
  p <- mp_taxa_heatmap(data, rank = "Genus", top_n = 5, show_dendrogram = FALSE)
  built <- ggplot2::ggplot_build(p)
  expect_equal(length(unique(built$data[[1]]$y)), 5)
})

test_that("log10 transform is accepted and changes the legend title", {
  data <- load_example("example_valid")
  p_clr <- mp_taxa_heatmap(data, rank = "Genus", top_n = 8, transform = "clr", show_dendrogram = FALSE)
  p_log <- mp_taxa_heatmap(data, rank = "Genus", top_n = 8, transform = "log10", show_dendrogram = FALSE)
  expect_false(identical(p_clr$scales$get_scales("fill")$name, p_log$scales$get_scales("fill")$name))
})

test_that("cluster_rows/cluster_cols = FALSE skips clustering without error", {
  data <- load_example("example_valid")
  p <- mp_taxa_heatmap(data, rank = "Genus", top_n = 8,
                        cluster_rows = FALSE, cluster_cols = FALSE, show_dendrogram = FALSE)
  expect_s3_class(p, "ggplot")
})

test_that("fix_taxonomy resolves missing genus before plotting", {
  data <- load_example("example_valid")
  p <- mp_taxa_heatmap(data, rank = "Genus", top_n = 20, show_dendrogram = FALSE, fix_taxonomy = TRUE)
  built <- ggplot2::ggplot_build(p)
  labels <- ggplot2::layer_scales(p)$y$get_labels()
  expect_true(any(grepl("^Unclassified_", labels)))
})

test_that("min_prevalence excludes taxa detected in too few samples", {
  data <- load_example("example_valid")
  p_all <- mp_taxa_heatmap(data, rank = "Genus", top_n = NULL, show_dendrogram = FALSE)
  p_filtered <- mp_taxa_heatmap(data, rank = "Genus", top_n = NULL, min_prevalence = 0.99, show_dendrogram = FALSE)
  n_all <- length(ggplot2::layer_scales(p_all)$y$get_labels())
  n_filtered <- length(ggplot2::layer_scales(p_filtered)$y$get_labels())
  expect_true(n_filtered <= n_all)
})
