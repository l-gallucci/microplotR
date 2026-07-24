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

load_tsv <- function(name, file) {
  d <- example_dir(name)
  readr::read_delim(file.path(d, file), delim = "\t", show_col_types = FALSE)
}

test_that("mp_taxa_treemap builds with grouping", {
  data <- load_example("example_valid")
  p <- mp_taxa_treemap(data, rank = "Genus", group_rank = "Phylum")
  expect_s3_class(p, "ggplot")
  built <- ggplot2::ggplot_build(p)
  expect_true(length(built$data) > 0)
})

test_that("mp_taxa_treemap areas sum to ~100", {
  data <- load_example("example_valid")
  p <- mp_taxa_treemap(data, rank = "Genus", group_rank = "Phylum")
  expect_equal(sum(p$data$rel_abund), 100, tolerance = 1e-6)
})

test_that("mp_taxa_treemap works without group_rank", {
  data <- load_example("example_valid")
  p <- mp_taxa_treemap(data, rank = "Genus", group_rank = NULL, top_n = 8)
  expect_s3_class(p, "ggplot")
})

test_that("top_n pools remainder into other_label", {
  data <- load_example("example_valid")
  p <- mp_taxa_treemap(data, rank = "Genus", group_rank = "Phylum", top_n = 3)
  expect_true("Other" %in% p$data$.taxon)
})

test_that("fix_taxonomy resolves missing genus", {
  data <- load_example("example_valid")
  p <- mp_taxa_treemap(data, rank = "Genus", group_rank = "Phylum", top_n = NULL, fix_taxonomy = TRUE)
  expect_true(any(grepl("^Unclassified_", p$data$.taxon)))
})

test_that("mp_function_treemap builds on eggNOG-mapper-shaped data", {
  ft <- load_tsv("example_function_profile", "gene_count_table.tsv")
  fa <- load_tsv("example_function_profile", "function_annotation.tsv")
  p <- mp_function_treemap(ft, fa, rank = "KEGG_ko", group_rank = "COG_category")
  expect_s3_class(p, "ggplot")
})
