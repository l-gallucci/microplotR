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

# Every taxa plot function accepts a precomputed `long` (the melt+join
# .mp_long_abundance() would derive internally) so a caller -- e.g. the
# Shiny app -- can compute it once per dataset and reuse it across
# rank/top_n/plot-type tweaks instead of re-melting on every call. These
# lock in that passing `long` explicitly is equivalent to letting each
# function derive it internally.

test_that("mp_taxa_barplot: long= matches internal derivation", {
  data <- load_example("example_valid")
  long <- microplotr:::.mp_long_abundance(data)
  p1 <- mp_taxa_barplot(data, rank = "Genus", group_rank = "Phylum", top_n = 5)
  p2 <- mp_taxa_barplot(data, rank = "Genus", group_rank = "Phylum", top_n = 5, long = long)
  expect_identical(p1$data, p2$data)
})

test_that("mp_taxa_heatmap: long= matches internal derivation", {
  data <- load_example("example_valid")
  long <- microplotr:::.mp_long_abundance(data)
  p1 <- mp_taxa_heatmap(data, rank = "Genus", top_n = 5)
  p2 <- mp_taxa_heatmap(data, rank = "Genus", top_n = 5, long = long)
  expect_identical(p1$data, p2$data)
})

test_that("mp_taxa_bubbleplot: long= matches internal derivation", {
  data <- load_example("example_valid")
  long <- microplotr:::.mp_long_abundance(data)
  p1 <- mp_taxa_bubbleplot(data, rank = "Genus", group_rank = "Phylum", top_n = 5)
  p2 <- mp_taxa_bubbleplot(data, rank = "Genus", group_rank = "Phylum", top_n = 5, long = long)
  expect_identical(p1$data, p2$data)
})

test_that("mp_asv_gradient_plot: long= matches internal derivation", {
  data <- load_example("example_valid")
  long <- microplotr:::.mp_long_abundance(data)
  gradient_var <- setdiff(names(data$metadata), "Sample_ID")[1]
  p1 <- mp_asv_gradient_plot(data, gradient_var = gradient_var, rank = "Genus", top_n = 5)
  p2 <- mp_asv_gradient_plot(data, gradient_var = gradient_var, rank = "Genus", top_n = 5, long = long)
  expect_identical(p1$data, p2$data)
})

test_that("mp_taxa_treemap: long= matches internal derivation", {
  data <- load_example("example_valid")
  long <- microplotr:::.mp_long_abundance(data)
  p1 <- mp_taxa_treemap(data, rank = "Genus", group_rank = "Phylum", top_n = 5)
  p2 <- mp_taxa_treemap(data, rank = "Genus", group_rank = "Phylum", top_n = 5, long = long)
  expect_identical(p1$data, p2$data)
})

test_that("long= is reused across different rank/top_n calls without re-deriving", {
  data <- load_example("example_valid")
  long <- microplotr:::.mp_long_abundance(data)
  # Same `long` object works for multiple different plot configurations --
  # exactly the reuse pattern the Shiny app relies on.
  p_genus <- mp_taxa_barplot(data, rank = "Genus", group_rank = NULL, nested_legend = FALSE, top_n = 3, long = long)
  p_phylum <- mp_taxa_barplot(data, rank = "Phylum", group_rank = NULL, nested_legend = FALSE, top_n = 3, long = long)
  expect_s3_class(p_genus, "ggplot")
  expect_s3_class(p_phylum, "ggplot")
})
