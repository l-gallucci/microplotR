example_dir <- function(name) {
  d <- system.file("extdata", name, package = "microplotr")
  if (d == "") d <- testthat::test_path("..", "..", "inst", "extdata", name)
  d
}

load_mag <- function(name) {
  d <- example_dir(name)
  readr::read_delim(file.path(d, "mag_quality.tsv"), delim = "\t", show_col_types = FALSE)
}

test_that("valid MAG table passes validation", {
  mag <- load_mag("example_mag_quality")
  report <- mp_validate_mag(mag)
  expect_true(mp_is_valid(report))
})

test_that("completeness out of [0,100] range is an error", {
  mag <- load_mag("example_mag_broken_out_of_range")
  report <- mp_validate_mag(mag)
  expect_false(mp_is_valid(report))
  expect_true(any(grepl("outside the valid", mp_errors(report)$message)))
})

test_that("negative contamination is an error", {
  mag <- load_mag("example_mag_broken_negative")
  report <- mp_validate_mag(mag)
  expect_false(mp_is_valid(report))
  expect_true(any(grepl("cannot be negative", mp_errors(report)$message)))
})

test_that("duplicate MAG names detected", {
  mag <- load_mag("example_mag_broken_duplicate_ids")
  report <- mp_validate_mag(mag)
  expect_false(mp_is_valid(report))
  expect_true(any(grepl("Duplicate MAG Name", mp_errors(report)$message)))
})

test_that("missing required column detected", {
  mag <- load_mag("example_mag_broken_missing_columns")
  report <- mp_validate_mag(mag)
  expect_false(mp_is_valid(report))
  expect_true("Contamination" %in% mp_errors(report)$field)
})

test_that("CheckM1-style column names ('Bin Id') are normalized and accepted", {
  mag <- load_mag("example_mag_quality")
  checkm1_style <- mag
  names(checkm1_style)[names(checkm1_style) == "Name"] <- "Bin Id"
  names(checkm1_style)[names(checkm1_style) == "Genome_Size"] <- "Genome size (bp)"
  report <- mp_validate_mag(checkm1_style)
  expect_true(mp_is_valid(report))
})

test_that("mp_mag_quality_plot builds without error and applies MIMAG shading", {
  mag <- load_mag("example_mag_quality")
  p <- mp_mag_quality_plot(mag, size_col = "Genome_Size", color_col = "Phylum")
  expect_s3_class(p, "ggplot")
  built <- ggplot2::ggplot_build(p)
  expect_true(length(built$data) >= 2) # rect/hline/vline layers + points
})

test_that("mp_mag_quality_plot works with no size/color columns", {
  mag <- load_mag("example_mag_quality")
  p <- mp_mag_quality_plot(mag, size_col = NULL, color_col = NULL, show_mimag_thresholds = FALSE)
  expect_s3_class(p, "ggplot")
})

test_that("mp_mag_quality_plot errors on unknown color_col", {
  mag <- load_mag("example_mag_quality")
  expect_error(mp_mag_quality_plot(mag, color_col = "NotAColumn"), "not found")
})

test_that("mp_mag_quality_distribution builds a two-facet histogram", {
  mag <- load_mag("example_mag_quality")
  p <- mp_mag_quality_distribution(mag)
  expect_s3_class(p, "ggplot")
  expect_true(inherits(p$facet, "FacetWrap"))
})
