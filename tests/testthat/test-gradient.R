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

test_that("builds with nested grouping and smoothing", {
  data <- load_example("example_valid")
  p <- mp_asv_gradient_plot(data, gradient_var = "Depth_m", rank = "Genus", group_rank = "Phylum", top_n = 8)
  expect_s3_class(p, "ggplot")
  built <- ggplot2::ggplot_build(p)
  expect_true(length(built$data) >= 2) # points + smooth layers
})

test_that("gradient column is used as numeric x", {
  data <- load_example("example_valid")
  p <- mp_asv_gradient_plot(data, gradient_var = "Depth_m", rank = "Genus", group_rank = "Phylum", top_n = 8)
  expect_true(is.numeric(p$data$Depth_m))
})

test_that("top_n caps taxa shown, rest pooled into Other", {
  data <- load_example("example_valid")
  p <- mp_asv_gradient_plot(data, gradient_var = "Depth_m", rank = "Genus", group_rank = "Phylum", top_n = 3)
  expect_equal(length(levels(p$data$.taxon)), 4) # 3 kept + Other
  expect_true("Other" %in% levels(p$data$.taxon))
})

test_that("group_rank = NULL gives a flat (non-nested) legend", {
  data <- load_example("example_valid")
  p <- mp_asv_gradient_plot(data, gradient_var = "Depth_m", rank = "Genus", group_rank = NULL, top_n = 5)
  scale <- p$scales$get_scales("colour")
  expect_false(any(grepl("\\*\\*", scale$labels)))
})

test_that("smooth = FALSE omits the smoothing layer", {
  data <- load_example("example_valid")
  p_smooth <- mp_asv_gradient_plot(data, gradient_var = "Depth_m", top_n = 5, smooth = TRUE)
  p_flat <- mp_asv_gradient_plot(data, gradient_var = "Depth_m", top_n = 5, smooth = FALSE)
  expect_true(length(p_smooth$layers) > length(p_flat$layers))
})

test_that("facet = TRUE facets by taxon and drops the legend", {
  data <- load_example("example_valid")
  p <- mp_asv_gradient_plot(data, gradient_var = "Depth_m", rank = "Genus", group_rank = "Phylum",
                             top_n = 5, facet = TRUE)
  expect_true(inherits(p$facet, "FacetWrap"))
})

test_that("fix_taxonomy resolves missing genus", {
  data <- load_example("example_valid")
  p <- mp_asv_gradient_plot(data, gradient_var = "Depth_m", rank = "Genus", group_rank = "Phylum", top_n = 20)
  expect_true(any(grepl("^Unclassified_", levels(p$data$.taxon))))
})

test_that("min_rel_abund filters out low-abundance taxa into Other", {
  data <- load_example("example_valid")
  p <- mp_asv_gradient_plot(data, gradient_var = "Depth_m", rank = "Genus", group_rank = "Phylum",
                             top_n = NULL, min_rel_abund = 8)
  expect_true("Other" %in% levels(p$data$.taxon))
})
