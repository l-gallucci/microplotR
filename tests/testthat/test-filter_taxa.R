make_agg <- function() {
  # 4 samples, 3 taxa: A abundant+prevalent, B abundant but only in 1 sample,
  # C low abundance but prevalent everywhere.
  data.frame(
    Sample_ID = rep(c("S1", "S2", "S3", "S4"), times = 3),
    taxon = rep(c("A", "B", "C"), each = 4),
    rel_abund = c(20, 25, 22, 18, 40, 0, 0, 0, 1, 1, 1, 1)
  )
}

test_that("no filters (defaults) return every taxon", {
  agg <- make_agg()
  out <- mp_filter_taxa(agg, taxon_col = "taxon")
  expect_setequal(out, c("A", "B", "C"))
})

test_that("min_rel_abund excludes low-mean-abundance taxa", {
  agg <- make_agg()
  out <- mp_filter_taxa(agg, taxon_col = "taxon", min_rel_abund = 5)
  expect_setequal(out, c("A", "B"))
})

test_that("min_prevalence as a fraction excludes taxa detected in too few samples", {
  agg <- make_agg()
  out <- mp_filter_taxa(agg, taxon_col = "taxon", min_prevalence = 0.5)
  expect_setequal(out, c("A", "C"))
})

test_that("min_prevalence > 1 is treated as an absolute sample count", {
  agg <- make_agg()
  out <- mp_filter_taxa(agg, taxon_col = "taxon", min_prevalence = 3)
  expect_setequal(out, c("A", "C"))
})

test_that("both filters combine with AND logic", {
  agg <- make_agg()
  out <- mp_filter_taxa(agg, taxon_col = "taxon", min_rel_abund = 5, min_prevalence = 0.5)
  expect_setequal(out, "A")
})

test_that("detection threshold changes what counts as 'present' for prevalence", {
  agg <- make_agg()
  # C is always 1%; raising detection above 1 makes it never "detected"
  out <- mp_filter_taxa(agg, taxon_col = "taxon", min_prevalence = 0.5, detection = 1)
  expect_false("C" %in% out)
})
