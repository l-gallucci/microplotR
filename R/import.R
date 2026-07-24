#' Read tidy microbiome input files
#'
#' Reads the three tidy input files as character-typed tibbles (no
#' validation, no type coercion) so [mp_validate()] can distinguish
#' "not numeric" from "numeric" without prior silent coercion.
#'
#' @param feature_table_path Path to feature_table.tsv/.csv.
#' @param taxonomy_path Path to taxonomy.tsv/.csv.
#' @param metadata_path Path to metadata.tsv/.csv.
#' @return A list with `feature_table`, `taxonomy`, `metadata` tibbles.
#' @export
mp_read_data <- function(feature_table_path, taxonomy_path, metadata_path) {
  list(
    feature_table = .mp_read_table(feature_table_path),
    taxonomy = .mp_read_table(taxonomy_path),
    metadata = .mp_read_table(metadata_path)
  )
}

.mp_read_table <- function(path) {
  delim <- if (grepl("\\.(tsv|txt)$", path, ignore.case = TRUE)) "\t" else ","
  readr::read_delim(
    path,
    delim = delim,
    col_types = readr::cols(.default = readr::col_character()),
    show_col_types = FALSE
  )
}
