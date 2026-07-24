test_that("known values pass through unchanged", {
  tax <- data.frame(
    Feature_ID = "F1", Phylum = "Firmicutes", Family = "Lactobacillaceae",
    Genus = "Lactobacillus", stringsAsFactors = FALSE
  )
  fixed <- mp_tax_fix(tax)
  expect_equal(fixed$Genus, "Lactobacillus")
})

test_that("single unknown falls back to last known ancestor", {
  tax <- data.frame(
    Feature_ID = "F1", Phylum = "Firmicutes", Family = "Lactobacillaceae",
    Genus = NA_character_, stringsAsFactors = FALSE
  )
  fixed <- mp_tax_fix(tax)
  expect_equal(fixed$Genus, "Unclassified_Lactobacillaceae")
})

test_that("cascading unknowns anchor to the same true ancestor, not to each other", {
  tax <- data.frame(
    Feature_ID = "F1", Phylum = "Firmicutes", Family = "Lactobacillaceae",
    Genus = "unclassified", Species = "uncultured", stringsAsFactors = FALSE
  )
  fixed <- mp_tax_fix(tax)
  expect_equal(fixed$Genus, "Unclassified_Lactobacillaceae")
  expect_equal(fixed$Species, "Unclassified_Lactobacillaceae")
})

test_that("unknown with no known ancestor becomes plain 'Unclassified'", {
  tax <- data.frame(
    Feature_ID = "F1", Domain = NA_character_, Phylum = NA_character_,
    stringsAsFactors = FALSE
  )
  fixed <- mp_tax_fix(tax)
  expect_equal(fixed$Domain, "Unclassified")
  expect_equal(fixed$Phylum, "Unclassified")
})

test_that("unknown-string matching is case-insensitive and whitespace-trimmed", {
  tax <- data.frame(
    Feature_ID = c("F1", "F2", "F3"), Phylum = "Bacteroidota",
    Genus = c("UNCULTURED", "  ", "Unknown"), stringsAsFactors = FALSE
  )
  fixed <- mp_tax_fix(tax)
  expect_true(all(fixed$Genus == "Unclassified_Bacteroidota"))
})

test_that("errors informatively when no recognized rank columns are present", {
  tax <- data.frame(Feature_ID = "F1", Notes = "x", stringsAsFactors = FALSE)
  expect_error(mp_tax_fix(tax), "No recognized taxonomy rank columns")
})
