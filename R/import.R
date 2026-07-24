#' Read the raw object behind a path, with no shaping applied yet
#' (matrix/data.frame/tibble as-is for `.rds`, parsed table for tsv/csv).
#' @keywords internal
#' @noRd
.mp_read_raw <- function(path) {
  if (grepl("\\.rds$", path, ignore.case = TRUE)) {
    return(readRDS(path))
  }
  delim <- if (grepl("\\.(tsv|txt)$", path, ignore.case = TRUE)) "\t" else ","
  readr::read_delim(
    path,
    delim = delim,
    col_types = readr::cols(.default = readr::col_character()),
    show_col_types = FALSE
  )
}

#' Rename whichever of `rename`'s keys are present as columns to their
#' mapped canonical name (skips a key if the canonical name is already a
#' column, so it never clobbers or double-applies).
#' @keywords internal
#' @noRd
.mp_apply_rename <- function(df, rename) {
  if (is.null(rename)) return(df)
  for (actual in names(rename)) {
    canonical <- rename[[actual]]
    if (actual %in% names(df) && !(canonical %in% names(df))) {
      names(df)[names(df) == actual] <- canonical
    }
  }
  df
}

#' Build a `rename` vector (see [mp_read_table()]) for a single custom ID
#' column name, or `NULL` if it already matches the canonical name.
#' @keywords internal
#' @noRd
.mp_id_rename <- function(actual, canonical) {
  if (identical(actual, canonical)) return(NULL)
  stats::setNames(canonical, actual)
}

#' Turn a raw `.mp_read_raw()` result into the finalized character data
#' frame `mp_read_table()`/`mp_read_data()` return: matrix -> data frame,
#' column rename, rowname -> `id_col` promotion, character coercion.
#' @keywords internal
#' @noRd
.mp_finalize_table <- function(obj, path, id_col = NULL, rename = NULL) {
  if (is.matrix(obj)) {
    obj <- as.data.frame(obj, stringsAsFactors = FALSE, check.names = FALSE)
  }
  if (!is.data.frame(obj)) {
    stop(sprintf("'%s' is an .rds file but does not contain a data frame or matrix (got class: %s).",
                  path, paste(class(obj), collapse = "/")))
  }

  obj <- .mp_apply_rename(obj, rename)

  if (!is.null(id_col) && !(id_col %in% names(obj))) {
    rn <- rownames(obj)
    has_named_rows <- !is.null(rn) && !identical(rn, as.character(seq_len(nrow(obj))))
    if (has_named_rows) {
      obj <- cbind(setNames(data.frame(rn, stringsAsFactors = FALSE), id_col), obj)
      rownames(obj) <- NULL
    }
  }

  as.data.frame(lapply(obj, as.character), stringsAsFactors = FALSE, check.names = FALSE)
}

#' If a matrix's rows/columns are swapped relative to known feature/sample
#' IDs (e.g. dada2's `seqtab.rds`, which comes out samples x ASVs -- the
#' opposite of the features x samples this package expects), flip it back.
#' No-op for anything but a matrix with both dimnames set, or when neither
#' orientation is a clearly better fit (e.g. IDs not yet known, or a
#' fresh/tiny table with no overlap either way).
#' @keywords internal
#' @noRd
.mp_maybe_transpose <- function(obj, feature_ids, sample_ids) {
  if (!is.matrix(obj) || is.null(rownames(obj)) || is.null(colnames(obj))) return(obj)
  if (length(feature_ids) == 0 || length(sample_ids) == 0) return(obj)

  rn <- rownames(obj)
  cn <- colnames(obj)
  as_is_score <- length(intersect(rn, feature_ids)) + length(intersect(cn, sample_ids))
  flipped_score <- length(intersect(rn, sample_ids)) + length(intersect(cn, feature_ids))

  if (flipped_score > as_is_score) t(obj) else obj
}

#' Read tidy microbiome input files
#'
#' Reads the three tidy input files as character-typed tibbles (no
#' validation, no type coercion) so [mp_validate()] can distinguish
#' "not numeric" from "numeric" without prior silent coercion.
#'
#' Two common real-world mismatches are handled automatically:
#' - **Custom ID column names.** If your files use something other than
#'   `Feature_ID`/`Sample_ID` (e.g. a sample sheet with a `SampleName`
#'   column), pass its actual name via `feature_id_column`/
#'   `sample_id_column` and it's renamed to the canonical name on read --
#'   no need to edit the file itself.
#' - **Transposed feature tables.** Tools like dada2 (`seqtab.rds`) or a
#'   `phyloseq::otu_table()` sometimes save samples x features instead of
#'   features x samples. When `feature_table_path` is a `.rds` matrix, its
#'   orientation is checked against the (already-read) taxonomy/metadata
#'   ID sets and transposed back if that's a clearly better fit.
#'
#' @param feature_table_path Path to feature_table.tsv/.csv/.rds.
#' @param taxonomy_path Path to taxonomy.tsv/.csv/.rds.
#' @param metadata_path Path to metadata.tsv/.csv/.rds.
#' @param feature_id_column Actual name of the feature-ID column in
#'   `feature_table_path`/`taxonomy_path`, if not already `"Feature_ID"`.
#' @param sample_id_column Actual name of the sample-ID column in
#'   `metadata_path` (and `feature_table_path`'s header row, for tsv/csv),
#'   if not already `"Sample_ID"`.
#' @return A list with `feature_table`, `taxonomy`, `metadata` tibbles.
#' @export
mp_read_data <- function(feature_table_path, taxonomy_path, metadata_path,
                          feature_id_column = "Feature_ID", sample_id_column = "Sample_ID") {
  feature_rename <- .mp_id_rename(feature_id_column, "Feature_ID")
  sample_rename <- .mp_id_rename(sample_id_column, "Sample_ID")

  taxonomy <- mp_read_table(taxonomy_path, id_col = "Feature_ID", rename = feature_rename)
  metadata <- mp_read_table(metadata_path, id_col = "Sample_ID", rename = sample_rename)

  ft_raw <- .mp_maybe_transpose(
    .mp_read_raw(feature_table_path),
    feature_ids = taxonomy$Feature_ID,
    sample_ids = metadata$Sample_ID
  )
  feature_table <- .mp_finalize_table(ft_raw, feature_table_path, id_col = "Feature_ID", rename = feature_rename)

  list(feature_table = feature_table, taxonomy = taxonomy, metadata = metadata)
}

#' Read a single tidy table (.tsv, .csv, or .rds)
#'
#' Delimiter/format is picked from the file extension: `.tsv`/`.txt` reads
#' tab-delimited, `.csv` reads comma-delimited, `.rds` reads a serialized R
#' object via [readRDS()] (a data frame/tibble, or a matrix -- e.g. a
#' `phyloseq::otu_table()` -- which is coerced to one, see `id_col` below).
#' Values always come back character-typed (no numeric coercion) so
#' [mp_validate()] can distinguish "not numeric" from "NA" -- an `.rds`
#' input with real numeric/integer columns is converted to character on
#' read for the same reason. Shared by [mp_read_data()] and the Shiny
#' apps, so upload handling and library loading never drift apart.
#'
#' @param path Path to the file. For a Shiny `fileInput`, pass the
#'   `datapath` but check/derive the extension from the original `name`
#'   (Shiny renames the temp file), e.g.
#'   `mp_read_table(input$file$datapath)` only works if `datapath` itself
#'   ends in the right extension -- see `inst/shiny/app.R` for the pattern
#'   used when that's not the case.
#' @param id_col Name of the ID column this table is expected to have
#'   (e.g. `"Feature_ID"`, `"Sample_ID"`). Only used for `.rds` input that
#'   turns out to be a `matrix` (common when exporting from `phyloseq` or
#'   similar, where feature/sample IDs live in `rownames()` rather than a
#'   column): the matrix is converted to a data frame and, if it doesn't
#'   already have a column named `id_col`, its rownames are promoted into
#'   one. Leave `NULL` to skip that promotion (rownames are dropped as
#'   usual, matching plain `as.data.frame()` behavior).
#' @param rename Optional named character vector, `c(actual_name =
#'   "canonical_name")`, for renaming one or more columns on read -- e.g.
#'   `c(SampleName = "Sample_ID")` if your metadata's ID column isn't
#'   already called `Sample_ID`. Applied before `id_col` promotion, and
#'   never overwrites a column that's already named `canonical_name`.
#' @return A tibble (or data frame, for `.rds`) with character columns.
#' @export
mp_read_table <- function(path, id_col = NULL, rename = NULL) {
  .mp_finalize_table(.mp_read_raw(path), path, id_col = id_col, rename = rename)
}
