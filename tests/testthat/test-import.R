example_dir <- function(name) {
  d <- system.file("extdata", name, package = "microplotr")
  if (d == "") d <- testthat::test_path("..", "..", "inst", "extdata", name)
  d
}

test_that("mp_read_table reads tsv", {
  path <- file.path(example_dir("example_valid"), "metadata.tsv")
  df <- mp_read_table(path)
  expect_true(is.data.frame(df))
  expect_true("Sample_ID" %in% names(df))
  expect_true(is.character(df$Sample_ID))
})

test_that("mp_read_table reads csv", {
  tsv_path <- file.path(example_dir("example_valid"), "metadata.tsv")
  ref <- readr::read_delim(tsv_path, delim = "\t", show_col_types = FALSE)

  csv_path <- tempfile(fileext = ".csv")
  write.csv(ref, csv_path, row.names = FALSE)

  df <- mp_read_table(csv_path)
  expect_true(is.data.frame(df))
  expect_setequal(names(df), names(ref))
  expect_equal(nrow(df), nrow(ref))
})

test_that("mp_read_table reads rds and coerces to character", {
  df_in <- data.frame(Sample_ID = c("S1", "S2"), Depth_m = c(5L, 10L), stringsAsFactors = FALSE)
  rds_path <- tempfile(fileext = ".rds")
  saveRDS(df_in, rds_path)

  df_out <- mp_read_table(rds_path)
  expect_true(is.data.frame(df_out))
  expect_true(is.character(df_out$Depth_m))
  expect_equal(df_out$Depth_m, c("5", "10"))
})

test_that("mp_read_table errors informatively on a non-data-frame rds", {
  rds_path <- tempfile(fileext = ".rds")
  saveRDS(list(a = 1, b = 2), rds_path)
  expect_error(mp_read_table(rds_path), "does not contain a data frame")
})

test_that("mp_read_table promotes matrix rownames to id_col when given", {
  mat <- matrix(1:6, nrow = 2, dimnames = list(c("F1", "F2"), c("S1", "S2", "S3")))
  rds_path <- tempfile(fileext = ".rds")
  saveRDS(mat, rds_path)

  df <- mp_read_table(rds_path, id_col = "Feature_ID")
  expect_true(is.data.frame(df))
  expect_true("Feature_ID" %in% names(df))
  expect_equal(df$Feature_ID, c("F1", "F2"))
  expect_equal(df$S1, c("1", "2"))
})

test_that("mp_read_table leaves a matrix's rownames alone when id_col is NULL", {
  mat <- matrix(1:6, nrow = 2, dimnames = list(c("F1", "F2"), c("S1", "S2", "S3")))
  rds_path <- tempfile(fileext = ".rds")
  saveRDS(mat, rds_path)

  df <- mp_read_table(rds_path)
  expect_true(is.data.frame(df))
  expect_false("Feature_ID" %in% names(df))
  expect_equal(nrow(df), 2)
})

test_that("mp_read_table does not duplicate id_col if the matrix already has it as a column", {
  mat <- matrix(c("F1", "F2", "1", "2"), nrow = 2, dimnames = list(NULL, c("Feature_ID", "S1")))
  rds_path <- tempfile(fileext = ".rds")
  saveRDS(mat, rds_path)

  df <- mp_read_table(rds_path, id_col = "Feature_ID")
  expect_equal(sum(names(df) == "Feature_ID"), 1)
})

test_that("mp_read_data reads a phyloseq-style otu_table matrix end to end", {
  d <- example_dir("example_valid")
  ft_ref <- readr::read_delim(file.path(d, "feature_table.tsv"), delim = "\t", show_col_types = FALSE)
  tax_ref <- readr::read_delim(file.path(d, "taxonomy.tsv"), delim = "\t", show_col_types = FALSE)

  ft_mat <- as.matrix(ft_ref[setdiff(names(ft_ref), "Feature_ID")])
  rownames(ft_mat) <- ft_ref$Feature_ID
  ft_rds <- tempfile(fileext = ".rds")
  saveRDS(ft_mat, ft_rds)

  data <- mp_read_data(ft_rds, file.path(d, "taxonomy.tsv"), file.path(d, "metadata.tsv"))
  expect_true("Feature_ID" %in% names(data$feature_table))
  report <- mp_validate(data)
  expect_true(mp_is_valid(report))
})

test_that("mp_read_data works end to end with a mix of tsv/csv/rds inputs", {
  d <- example_dir("example_valid")
  ft_ref <- readr::read_delim(file.path(d, "feature_table.tsv"), delim = "\t", show_col_types = FALSE)
  tax_ref <- readr::read_delim(file.path(d, "taxonomy.tsv"), delim = "\t", show_col_types = FALSE)
  meta_ref <- readr::read_delim(file.path(d, "metadata.tsv"), delim = "\t", show_col_types = FALSE)

  csv_path <- tempfile(fileext = ".csv")
  write.csv(tax_ref, csv_path, row.names = FALSE)
  rds_path <- tempfile(fileext = ".rds")
  saveRDS(meta_ref, rds_path)

  data <- mp_read_data(file.path(d, "feature_table.tsv"), csv_path, rds_path)
  report <- mp_validate(data)
  expect_true(mp_is_valid(report))
})
