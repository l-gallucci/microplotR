example_dir <- function(name) {
  system.file("extdata", name, package = "microplotr")
}

load_example <- function(name) {
  d <- example_dir(name)
  if (d == "") {
    # not yet installed as a package -- fall back to relative path for local dev/testthat runs
    d <- file.path(testthat::test_path("..", "..", "inst", "extdata", name))
  }
  mp_read_data(
    file.path(d, "feature_table.tsv"),
    file.path(d, "taxonomy.tsv"),
    file.path(d, "metadata.tsv")
  )
}

test_that("valid dataset has no errors", {
  data <- load_example("example_valid")
  report <- mp_validate(data)
  expect_true(mp_is_valid(report))
  expect_equal(nrow(mp_errors(report)), 0)
})

test_that("valid dataset passes gradient and group checks", {
  data <- load_example("example_valid")
  report <- mp_validate(data, gradient_column = "Depth_m", group_column = "Group")
  expect_true(mp_is_valid(report))
})

test_that("id mismatch detected in both directions", {
  data <- load_example("example_broken_id_mismatch")
  report <- mp_validate(data)
  expect_false(mp_is_valid(report))
  msgs <- paste(mp_errors(report)$message, collapse = " ")
  expect_true(grepl("S9", msgs))
  expect_true(grepl("S8", msgs))
})

test_that("negative counts detected", {
  data <- load_example("example_broken_negative_counts")
  report <- mp_validate(data)
  expect_false(mp_is_valid(report))
  expect_true(any(grepl("Negative abundance", mp_errors(report)$message)))
})

test_that("missing taxonomy column detected", {
  data <- load_example("example_broken_missing_columns")
  report <- mp_validate(data)
  expect_false(mp_is_valid(report))
  expect_true("Phylum" %in% mp_errors(report)$field)
})

test_that("duplicate feature ids detected", {
  data <- load_example("example_broken_duplicate_ids")
  report <- mp_validate(data)
  expect_false(mp_is_valid(report))
  expect_true(any(grepl("Duplicate Feature_ID", mp_errors(report)$message)))
})

test_that("missing group column is an error", {
  data <- load_example("example_valid")
  report <- mp_validate(data, group_column = "NotAColumn")
  expect_false(mp_is_valid(report))
})

test_that("single-level group column is an error", {
  data <- load_example("example_valid")
  data$metadata$Constant <- "only_one_value"
  report <- mp_validate(data, group_column = "Constant")
  expect_false(mp_is_valid(report))
})

test_that("required_ranks can be relaxed for non-taxonomic (e.g. functional) annotation tables", {
  ft <- data.frame(Feature_ID = c("g1", "g2"), S1 = c(10, 20), S2 = c(5, 15))
  fa <- data.frame(Feature_ID = c("g1", "g2"), COG_category = c("C", "J"), KEGG_ko = c("K00001", "K00002"))
  meta <- data.frame(Sample_ID = c("S1", "S2"))
  data <- list(feature_table = ft, taxonomy = fa, metadata = meta)

  report_default <- mp_validate(data)
  expect_false(mp_is_valid(report_default)) # Phylum/Genus missing by default

  report_relaxed <- mp_validate(data, required_ranks = character(0))
  expect_true(mp_is_valid(report_relaxed))
})
