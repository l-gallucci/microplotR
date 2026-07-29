make_data <- function() {
  # F1: abundant, kept. F2: Mitochondria contaminant. F3: Chloroplast
  # contaminant (lowercase in the file, to check case-insensitivity).
  # F4: singleton (total count 1). F5: explicitly excluded by ID.
  feature_table <- data.frame(
    Feature_ID = c("F1", "F2", "F3", "F4", "F5"),
    S1 = c(50, 30, 20, 1, 5),
    S2 = c(60, 10, 15, 0, 5),
    stringsAsFactors = FALSE
  )
  taxonomy <- data.frame(
    Feature_ID = c("F1", "F2", "F3", "F4", "F5"),
    Phylum = c("Proteobacteria", "Proteobacteria", "Cyanobacteria", "Proteobacteria", "Proteobacteria"),
    Family = c("Enterobacteriaceae", "Mitochondria", "Streptophyta", "Enterobacteriaceae", "Enterobacteriaceae"),
    Order = c("Enterobacterales", "Rickettsiales", "chloroplast", "Enterobacterales", "Enterobacterales"),
    Genus = c("Escherichia", "Escherichia", "Chloroplast", "Escherichia", "Escherichia"),
    stringsAsFactors = FALSE
  )
  metadata <- data.frame(Sample_ID = c("S1", "S2"), stringsAsFactors = FALSE)
  list(feature_table = feature_table, taxonomy = taxonomy, metadata = metadata)
}

test_that("no filters (defaults) keep everything and report zero removed", {
  out <- mp_filter_features(make_data())
  expect_equal(nrow(out$data$feature_table), 5)
  expect_equal(nrow(out$data$taxonomy), 5)
  total <- out$summary[out$summary$reason == "total", ]
  expect_equal(total$n_features_removed, 0)
  expect_equal(total$pct_features_removed, 0)
})

test_that("exclude_taxa drops features by taxonomy match across ranks, case-insensitively", {
  out <- mp_filter_features(make_data(), exclude_taxa = c("Mitochondria", "Chloroplast"))
  expect_setequal(out$data$feature_table$Feature_ID, c("F1", "F4", "F5"))
  expect_setequal(out$data$taxonomy$Feature_ID, c("F1", "F4", "F5"))
  row <- out$summary[out$summary$reason == "exclude_taxa", ]
  expect_equal(row$n_features_removed, 2)
})

test_that("exclude_ranks restricts which taxonomy columns are checked", {
  # "Rickettsiales" only appears in Order; restricting to Family should miss it
  out <- mp_filter_features(make_data(), exclude_taxa = "Rickettsiales", exclude_ranks = "Family")
  expect_setequal(out$data$feature_table$Feature_ID, c("F1", "F2", "F3", "F4", "F5"))
})

test_that("exclude_features drops by exact Feature_ID regardless of taxonomy", {
  out <- mp_filter_features(make_data(), exclude_features = "F5")
  expect_false("F5" %in% out$data$feature_table$Feature_ID)
  expect_false("F5" %in% out$data$taxonomy$Feature_ID)
  row <- out$summary[out$summary$reason == "exclude_features", ]
  expect_equal(row$n_features_removed, 1)
})

test_that("remove_singletons drops only the total-count-1 feature", {
  out <- mp_filter_features(make_data(), remove_singletons = TRUE)
  expect_false("F4" %in% out$data$feature_table$Feature_ID)
  expect_setequal(out$data$feature_table$Feature_ID, c("F1", "F2", "F3", "F5"))
})

test_that("min_total_count generalizes beyond singletons", {
  out <- mp_filter_features(make_data(), min_total_count = 10)
  # F4 (total 1) and F5 (total 10, not < 10) -- only F4 removed
  expect_false("F4" %in% out$data$feature_table$Feature_ID)
  expect_true("F5" %in% out$data$feature_table$Feature_ID)
})

test_that("reasons combine and total isn't double-counted for features hit by more than one", {
  out <- mp_filter_features(
    make_data(),
    exclude_taxa = c("Mitochondria", "Chloroplast"),
    exclude_features = "F2",  # already caught by exclude_taxa too
    remove_singletons = TRUE
  )
  expect_setequal(out$data$feature_table$Feature_ID, c("F1", "F5"))
  total <- out$summary[out$summary$reason == "total", ]
  expect_equal(total$n_features_removed, 3)  # F2, F3, F4 -- not 4
})

test_that("summary reports percentage of reads removed against the original total", {
  out <- mp_filter_features(make_data(), exclude_features = "F1")
  # F1 total = 110 out of grand total 196 (50+30+20+1+5+60+10+15+0+5)
  total <- out$summary[out$summary$reason == "total", ]
  expect_equal(total$n_reads_removed, 110)
  expect_equal(round(total$pct_reads_removed, 1), round(100 * 110 / 196, 1))
})

test_that("filtered output is still valid input to mp_validate", {
  out <- mp_filter_features(make_data(), exclude_taxa = c("Mitochondria", "Chloroplast"))
  report <- mp_validate(out$data)
  expect_true(mp_is_valid(report))
})

test_that("require_ranks drops features missing a value at any given rank", {
  data <- make_data()
  data$taxonomy$Phylum[data$taxonomy$Feature_ID == "F1"] <- NA
  data$taxonomy$Phylum[data$taxonomy$Feature_ID == "F2"] <- "  "  # blank after trim

  out <- mp_filter_features(data, require_ranks = "Phylum")
  expect_false("F1" %in% out$data$feature_table$Feature_ID)
  expect_false("F2" %in% out$data$feature_table$Feature_ID)
  expect_true(all(c("F3", "F4", "F5") %in% out$data$feature_table$Feature_ID))
  row <- out$summary[out$summary$reason == "missing_rank", ]
  expect_equal(row$n_features_removed, 2)
})

test_that("require_ranks with multiple ranks drops a feature missing any one of them", {
  data <- make_data()
  data$taxonomy$Genus[data$taxonomy$Feature_ID == "F5"] <- NA  # Phylum fine, Genus missing

  out <- mp_filter_features(data, require_ranks = c("Phylum", "Genus"))
  expect_false("F5" %in% out$data$feature_table$Feature_ID)
  expect_true(all(c("F1", "F2", "F3", "F4") %in% out$data$feature_table$Feature_ID))
})
