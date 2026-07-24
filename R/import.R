#' Read tidy microbiome input files
#'
#' Reads the three tidy input files as character-typed tibbles (no
#' validation, no type coercion) so [mp_validate()] can distinguish
#' "not numeric" from "numeric" without prior silent coercion.
#'
#' @param feature_table_path Path to feature_table.tsv/.csv/.rds.
#' @param taxonomy_path Path to taxonomy.tsv/.csv/.rds.
#' @param metadata_path Path to metadata.tsv/.csv/.rds.
#' @return A list with `feature_table`, `taxonomy`, `metadata` tibbles.
#' @export
mp_read_data <- function(feature_table_path, taxonomy_path, metadata_path) {
  list(
    feature_table = mp_read_table(feature_table_path),
    taxonomy = mp_read_table(taxonomy_path),
    metadata = mp_read_table(metadata_path)
  )
}

#' Read a single tidy table (.tsv, .csv, or .rds)
#'
#' Delimiter/format is picked from the file extension: `.tsv`/`.txt` reads
#' tab-delimited, `.csv` reads comma-delimited, `.rds` reads a serialized R
#' object via [readRDS()] (must be a data frame/tibble). Values always come
#' back character-typed (no numeric coercion) so [mp_validate()] can
#' distinguish "not numeric" from "NA" -- an `.rds` data frame with real
#' numeric/integer columns is converted to character on read for the same
#' reason. Shared by [mp_read_data()] and the Shiny apps, so upload
#' handling and library loading never drift apart.
#'
#' @param path Path to the file. For a Shiny `fileInput`, pass the
#'   `datapath` but check/derive the extension from the original `name`
#'   (Shiny renames the temp file), e.g.
#'   `mp_read_table(input$file$datapath)` only works if `datapath` itself
#'   ends in the right extension -- see `inst/shiny/app.R` for the pattern
#'   used when that's not the case.
#' @return A tibble (or data frame, for `.rds`) with character columns.
#' @export
mp_read_table <- function(path) {
  if (grepl("\\.rds$", path, ignore.case = TRUE)) {
    obj <- readRDS(path)
    if (!is.data.frame(obj)) {
      stop(sprintf("'%s' is an .rds file but does not contain a data frame (got class: %s).",
                    path, paste(class(obj), collapse = "/")))
    }
    return(as.data.frame(lapply(obj, as.character), stringsAsFactors = FALSE, check.names = FALSE))
  }

  delim <- if (grepl("\\.(tsv|txt)$", path, ignore.case = TRUE)) "\t" else ","
  readr::read_delim(
    path,
    delim = delim,
    col_types = readr::cols(.default = readr::col_character()),
    show_col_types = FALSE
  )
}
