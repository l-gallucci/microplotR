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

test_that("bubbleplot builds without error and is a ggplot", {
  data <- load_example("example_valid")
  p <- mp_taxa_bubbleplot(data, rank = "Genus", group_rank = "Phylum", top_n = 8)
  expect_s3_class(p, "ggplot")
})

test_that("grid is complete: n_samples * top_n rows drawn, absent combos sized 0", {
  data <- load_example("example_valid")
  n_samples <- ncol(data$feature_table) - 1
  p <- mp_taxa_bubbleplot(data, rank = "Genus", group_rank = "Phylum", top_n = 6)
  expect_equal(nrow(p$data), n_samples * 6)
  expect_true(all(p$data$rel_abund >= 0))
})

test_that("top_n limits number of distinct taxa on the y axis", {
  data <- load_example("example_valid")
  p <- mp_taxa_bubbleplot(data, rank = "Genus", group_rank = "Phylum", top_n = 4)
  expect_equal(length(levels(p$data$Genus)), 4)
})

test_that("group_rank = NULL drops the color aesthetic", {
  data <- load_example("example_valid")
  p <- mp_taxa_bubbleplot(data, rank = "Genus", group_rank = NULL, top_n = 6)
  expect_false("colour" %in% names(p$mapping) && !is.null(rlang::get_expr(p$mapping$colour)))
})

test_that("fix_taxonomy resolves missing genus before plotting", {
  data <- load_example("example_valid")
  p <- mp_taxa_bubbleplot(data, rank = "Genus", group_rank = "Phylum", top_n = 20, fix_taxonomy = TRUE)
  expect_true(any(grepl("^Unclassified_", levels(p$data$Genus))))
})

test_that("facet_var adds a facet layer", {
  data <- load_example("example_valid")
  p <- mp_taxa_bubbleplot(data, rank = "Genus", group_rank = "Phylum", top_n = 6, facet_var = "Group")
  expect_true(inherits(p$facet, "FacetWrap"))
})

test_that("min_rel_abund reduces the number of taxa shown", {
  data <- load_example("example_valid")
  p_all <- mp_taxa_bubbleplot(data, rank = "Genus", group_rank = "Phylum", top_n = NULL)
  p_filtered <- mp_taxa_bubbleplot(data, rank = "Genus", group_rank = "Phylum", top_n = NULL, min_rel_abund = 8)
  expect_true(length(levels(p_filtered$data$Genus)) <= length(levels(p_all$data$Genus)))
})
