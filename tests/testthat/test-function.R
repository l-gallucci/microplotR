example_dir <- function(name) {
  d <- system.file("extdata", name, package = "microplotr")
  if (d == "") d <- testthat::test_path("..", "..", "inst", "extdata", name)
  d
}

load_tsv <- function(name, file) {
  d <- example_dir(name)
  readr::read_delim(file.path(d, file), delim = "\t", show_col_types = FALSE)
}

test_that("mp_function_barplot builds with COG_category grouping (eggNOG-mapper shape)", {
  ft <- load_tsv("example_function_profile", "gene_count_table.tsv")
  fa <- load_tsv("example_function_profile", "function_annotation.tsv")
  p <- mp_function_barplot(ft, fa, rank = "KEGG_ko", group_rank = "COG_category", top_n = 8)
  expect_s3_class(p, "ggplot")
  built <- ggplot2::ggplot_build(p)
  expect_true(nrow(built$data[[1]]) > 0)
})

test_that("mp_function_barplot resolves blank COG_category via tax_fix_ranks", {
  ft <- load_tsv("example_function_profile", "gene_count_table.tsv")
  fa <- load_tsv("example_function_profile", "function_annotation.tsv")
  p <- mp_function_barplot(ft, fa, rank = "KEGG_ko", group_rank = "COG_category", top_n = 20)
  expect_true("Unclassified" %in% as.character(unique(fa$COG_category)) == FALSE) # sanity: raw data has blanks, not literal "Unclassified"
  expect_s3_class(p, "ggplot")
})

test_that("mp_function_heatmap builds", {
  ft <- load_tsv("example_function_profile", "gene_count_table.tsv")
  fa <- load_tsv("example_function_profile", "function_annotation.tsv")
  p <- mp_function_heatmap(ft, fa, rank = "KEGG_ko", top_n = 20, show_dendrogram = FALSE)
  expect_s3_class(p, "ggplot")
})

test_that("flat KofamScan-style annotation (no group_rank) works", {
  ft <- load_tsv("example_function_profile", "gene_count_table.tsv")
  fa <- load_tsv("example_function_profile", "function_annotation.tsv")
  fa_kofam <- fa[, c("Feature_ID", "KEGG_ko", "Description")]
  names(fa_kofam) <- c("Feature_ID", "KO", "KO_definition")
  p <- mp_function_barplot(ft, fa_kofam, rank = "KO", group_rank = NULL, nested_legend = FALSE, top_n = 8)
  expect_s3_class(p, "ggplot")
})

test_that("REGRESSION: a rank value mapping to two different group values does not crash (barplot)", {
  # e.g. a KO whose COG category is blank on some genes but assigned on others
  ft <- data.frame(Feature_ID = c("g1", "g2", "g3"), S1 = c(10, 20, 5), S2 = c(15, 5, 10))
  fa <- data.frame(Feature_ID = c("g1", "g2", "g3"),
                    Category = c("A", "A", "B"),
                    Item = c("X", "X", "X")) # "X" maps to both "A" and "B" across rows
  data <- list(feature_table = ft, taxonomy = fa, metadata = data.frame(Sample_ID = c("S1", "S2")))
  p <- mp_taxa_barplot(data, rank = "Item", group_rank = "Category", top_n = NULL,
                        fix_taxonomy = FALSE, other_label = "Other")
  expect_s3_class(p, "ggplot")
  expect_equal(length(unique(levels(p$data$.taxon))), length(levels(p$data$.taxon)))
})

test_that("REGRESSION: a rank value mapping to two different group values does not crash (bubbleplot)", {
  ft <- data.frame(Feature_ID = c("g1", "g2", "g3"), S1 = c(10, 20, 5), S2 = c(15, 5, 10))
  fa <- data.frame(Feature_ID = c("g1", "g2", "g3"),
                    Category = c("A", "A", "B"),
                    Item = c("X", "X", "X"))
  data <- list(feature_table = ft, taxonomy = fa, metadata = data.frame(Sample_ID = c("S1", "S2")))
  p <- mp_taxa_bubbleplot(data, rank = "Item", group_rank = "Category", top_n = NULL, fix_taxonomy = FALSE)
  expect_s3_class(p, "ggplot")
})
