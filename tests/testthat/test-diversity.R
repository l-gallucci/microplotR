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

# --- alpha diversity ---

test_that("alpha diversity plot builds and facets by metric", {
  data <- load_example("example_valid")
  p <- mp_alpha_diversity_plot(data, group_var = "Group")
  expect_s3_class(p, "ggplot")
  expect_true(inherits(p$facet, "FacetWrap"))
})

test_that("Observed richness does not get fake y-jitter (constant metric stays constant)", {
  data <- load_example("example_valid")
  # Observed richness is identical across every sample in this example dataset,
  # so ggpubr's stat_compare_means() has nothing to compare -- expect its warning.
  p <- mp_alpha_diversity_plot(data, group_var = "Group", metrics = "Observed")
  built <- suppressWarnings(ggplot2::ggplot_build(p))
  jitter_layer <- built$data[[2]] # boxplot is layer 1, jitter is layer 2
  expect_equal(length(unique(jitter_layer$y)), 1)
})

test_that("metrics subset is respected", {
  data <- load_example("example_valid")
  p <- mp_alpha_diversity_plot(data, group_var = "Group", metrics = c("Shannon", "Simpson"))
  expect_equal(levels(p$data$metric), c("Shannon", "Simpson"))
})

test_that("test = NULL disables the stat annotation layer", {
  data <- load_example("example_valid")
  p_test <- mp_alpha_diversity_plot(data, group_var = "Group", test = "wilcoxon")
  p_none <- mp_alpha_diversity_plot(data, group_var = "Group", test = NULL)
  expect_true(length(p_test$layers) > length(p_none$layers))
})

test_that("a sample with NA group_var is dropped (with a warning), not left as a phantom group", {
  data <- load_example("example_valid")
  data$metadata$Group[1] <- NA
  expect_warning(
    p <- mp_alpha_diversity_plot(data, group_var = "Group"),
    "Dropped"
  )
  expect_false(anyNA(p$data$Group))
  dropped_sample <- data$metadata$Sample_ID[1]
  expect_false(dropped_sample %in% p$data$Sample_ID)
})

# --- beta diversity ---

test_that("beta diversity PCoA builds with PERMANOVA subtitle", {
  data <- load_example("example_valid")
  p <- mp_beta_diversity_plot(data, group_var = "Group", method = "bray", ordination = "pcoa")
  expect_s3_class(p, "ggplot")
  expect_true(grepl("PERMANOVA", p$labels$subtitle))
})

test_that("beta diversity NMDS builds with stress annotated", {
  data <- load_example("example_valid")
  p <- mp_beta_diversity_plot(data, group_var = "Group", method = "bray", ordination = "nmds")
  expect_s3_class(p, "ggplot")
  expect_true(grepl("stress", p$labels$subtitle))
})

test_that("jaccard method runs without error", {
  data <- load_example("example_valid")
  p <- mp_beta_diversity_plot(data, group_var = "Group", method = "jaccard", ordination = "pcoa")
  expect_s3_class(p, "ggplot")
})

test_that("show_ellipse = FALSE omits the ellipse layer", {
  data <- load_example("example_valid")
  p_ellipse <- mp_beta_diversity_plot(data, group_var = "Group", show_ellipse = TRUE, permanova = FALSE)
  p_flat <- mp_beta_diversity_plot(data, group_var = "Group", show_ellipse = FALSE, permanova = FALSE)
  expect_true(length(p_ellipse$layers) > length(p_flat$layers))
})

test_that("permanova = FALSE and show_ellipse = FALSE gives no subtitle", {
  data <- load_example("example_valid")
  p <- mp_beta_diversity_plot(data, group_var = "Group", permanova = FALSE, show_ellipse = FALSE)
  expect_null(p$labels$subtitle)
})

test_that("a sample with NA group_var is dropped (with a warning) instead of erroring adonis2", {
  data <- load_example("example_valid")
  na_sample <- data$metadata$Sample_ID[1]
  data$metadata$Group[1] <- NA
  expect_warning(
    p <- mp_beta_diversity_plot(data, group_var = "Group", method = "bray", ordination = "pcoa"),
    "Dropped"
  )
  expect_s3_class(p, "ggplot")
  expect_false(anyNA(p$data$Group))
  expect_false(na_sample %in% p$data$Sample_ID)
})
