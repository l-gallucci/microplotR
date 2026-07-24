example_dir <- function(name) {
  d <- system.file("extdata", name, package = "microplotr")
  if (d == "") d <- testthat::test_path("..", "..", "inst", "extdata", name)
  d
}

load_tsv <- function(name, file) {
  d <- example_dir(name)
  readr::read_delim(file.path(d, file), delim = "\t", show_col_types = FALSE)
}

test_that("valid contig_lengths table passes validation", {
  contigs <- load_tsv("example_assembly", "contig_lengths.tsv")
  expect_true(mp_is_valid(mp_validate_contig_lengths(contigs)))
})

test_that("valid assembly_summary table passes validation", {
  summary_tbl <- load_tsv("example_assembly", "assembly_summary.tsv")
  expect_true(mp_is_valid(mp_validate_assembly_summary(summary_tbl)))
})

test_that("negative contig length detected", {
  bad <- load_tsv("example_assembly_broken_negative_length", "contig_lengths.tsv")
  report <- mp_validate_contig_lengths(bad)
  expect_false(mp_is_valid(report))
  expect_true(any(grepl("must be positive", mp_errors(report)$message)))
})

test_that("duplicate contig id within assembly detected", {
  bad <- load_tsv("example_assembly_broken_duplicate_contig_ids", "contig_lengths.tsv")
  report <- mp_validate_contig_lengths(bad)
  expect_false(mp_is_valid(report))
  expect_true(any(grepl("Duplicate Contig_ID", mp_errors(report)$message)))
})

test_that("missing required summary column detected", {
  bad <- load_tsv("example_assembly_broken_missing_columns", "assembly_summary.tsv")
  report <- mp_validate_assembly_summary(bad)
  expect_false(mp_is_valid(report))
  expect_true("N50" %in% mp_errors(report)$field)
})

test_that("duplicate assembly id detected", {
  bad <- load_tsv("example_assembly_broken_duplicate_assembly_ids", "assembly_summary.tsv")
  report <- mp_validate_assembly_summary(bad)
  expect_false(mp_is_valid(report))
  expect_true(any(grepl("Duplicate Assembly_ID", mp_errors(report)$message)))
})

test_that("QUAST-style column names are normalized and accepted", {
  summary_tbl <- load_tsv("example_assembly", "assembly_summary.tsv")
  quast_style <- summary_tbl
  names(quast_style)[names(quast_style) == "Assembly_ID"] <- "Assembly"
  names(quast_style)[names(quast_style) == "Total_length"] <- "Total length"
  names(quast_style)[names(quast_style) == "GC_percent"] <- "GC (%)"
  report <- mp_validate_assembly_summary(quast_style)
  expect_true(mp_is_valid(report))
})

test_that("mp_assembly_nx_plot builds one step-line per assembly", {
  contigs <- load_tsv("example_assembly", "contig_lengths.tsv")
  p <- mp_assembly_nx_plot(contigs)
  expect_s3_class(p, "ggplot")
  built <- ggplot2::ggplot_build(p)
  expect_equal(length(unique(built$data[[1]]$group)), length(unique(contigs$Assembly_ID)))
})

test_that("mp_assembly_summary_barplot uses N50 by default", {
  summary_tbl <- load_tsv("example_assembly", "assembly_summary.tsv")
  p <- mp_assembly_summary_barplot(summary_tbl)
  expect_s3_class(p, "ggplot")
  expect_equal(rlang::as_label(p$mapping$y), "N50")
})

test_that("mp_assembly_summary_barplot errors on unknown stat_col", {
  summary_tbl <- load_tsv("example_assembly", "assembly_summary.tsv")
  expect_error(mp_assembly_summary_barplot(summary_tbl, stat_col = "NotAColumn"), "not found")
})
