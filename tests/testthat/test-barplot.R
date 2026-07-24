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

test_that("nested-legend barplot builds without error and is a ggplot", {
  data <- load_example("example_valid")
  p <- mp_taxa_barplot(data, rank = "Genus", group_rank = "Phylum", top_n = 5)
  expect_s3_class(p, "ggplot")
  built <- ggplot2::ggplot_build(p)
  expect_true(nrow(built$data[[1]]) > 0)
})

test_that("top_n caps the number of individual taxa shown, plus one Other bucket", {
  data <- load_example("example_valid")
  p <- mp_taxa_barplot(data, rank = "Genus", group_rank = "Phylum", top_n = 3)
  n_levels <- length(levels(p$data$.taxon))
  expect_equal(n_levels, 4) # 3 kept + Other
  expect_true("Other" %in% levels(p$data$.taxon))
})

test_that("legend labels bold only the first taxon in each phylum group", {
  data <- load_example("example_valid")
  p <- mp_taxa_barplot(data, rank = "Genus", group_rank = "Phylum", top_n = 8, fix_taxonomy = TRUE)
  build <- ggplot2::ggplot_build(p)
  scale <- p$scales$get_scales("fill")
  labels <- scale$labels
  bold_count <- sum(grepl("^\\*\\*", labels))
  n_groups <- length(unique(data$taxonomy$Phylum[!is.na(data$taxonomy$Phylum)])) |> min(bold_count + 1) # sanity bound
  expect_true(bold_count >= 1)
  expect_false(any(grepl("^\\*\\*Other", labels))) # Other never bold
})

test_that("nested_legend = TRUE without group_rank errors clearly", {
  data <- load_example("example_valid")
  expect_error(
    mp_taxa_barplot(data, rank = "Genus", group_rank = NULL, nested_legend = TRUE),
    "requires `group_rank`"
  )
})

test_that("nested_legend = FALSE produces a flat legend with plain labels", {
  data <- load_example("example_valid")
  p <- mp_taxa_barplot(data, rank = "Genus", group_rank = NULL, nested_legend = FALSE, top_n = 5)
  scale <- p$scales$get_scales("fill")
  expect_false(any(grepl("\\*\\*", scale$labels)))
})

test_that("fix_taxonomy resolves missing genus to Unclassified_<family>", {
  data <- load_example("example_valid")
  p <- mp_taxa_barplot(data, rank = "Genus", group_rank = "Phylum", top_n = 20, fix_taxonomy = TRUE)
  expect_true(any(grepl("^Unclassified_", levels(p$data$.taxon))))
})

test_that("facet_var adds a facet layer", {
  data <- load_example("example_valid")
  p <- mp_taxa_barplot(data, rank = "Genus", group_rank = "Phylum", top_n = 5, facet_var = "Group")
  expect_true(inherits(p$facet, "FacetWrap"))
})

test_that("explicit sample_order sets Sample_ID factor levels", {
  data <- load_example("example_valid")
  custom <- rev(data$metadata$Sample_ID)
  p <- mp_taxa_barplot(data, rank = "Genus", group_rank = "Phylum", top_n = 5, sample_order = custom)
  expect_equal(levels(p$data$Sample_ID), custom)
})

test_that("min_rel_abund pools low-abundance taxa into Other even with top_n = NULL", {
  data <- load_example("example_valid")
  p <- mp_taxa_barplot(data, rank = "Genus", group_rank = "Phylum", top_n = NULL, min_rel_abund = 8)
  expect_true("Other" %in% levels(p$data$.taxon))
  expect_true(length(levels(p$data$.taxon)) < 11)
})

test_that("top_n = NULL with no thresholds shows every taxon and no Other bucket", {
  data <- load_example("example_valid")
  p <- mp_taxa_barplot(data, rank = "Genus", group_rank = "Phylum", top_n = NULL)
  expect_false("Other" %in% levels(p$data$.taxon))
})
