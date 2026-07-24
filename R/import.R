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
    feature_table = mp_read_table(feature_table_path, id_col = "Feature_ID"),
    taxonomy = mp_read_table(taxonomy_path, id_col = "Feature_ID"),
    metadata = mp_read_table(metadata_path, id_col = "Sample_ID")
  )
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
#' @return A tibble (or data frame, for `.rds`) with character columns.
#' @export
mp_read_table <- function(path, id_col = NULL) {
  if (grepl("\\.rds$", path, ignore.case = TRUE)) {
    obj <- readRDS(path)
    if (is.matrix(obj)) {
      obj <- as.data.frame(obj, stringsAsFactors = FALSE, check.names = FALSE)
      has_named_rows <- !is.null(rownames(obj)) && !identical(rownames(obj), as.character(seq_len(nrow(obj))))
      if (!is.null(id_col) && has_named_rows && !(id_col %in% names(obj))) {
        obj <- cbind(setNames(data.frame(rownames(obj), stringsAsFactors = FALSE), id_col), obj)
        rownames(obj) <- NULL
      }
    }
    if (!is.data.frame(obj)) {
      stop(sprintf("'%s' is an .rds file but does not contain a data frame or matrix (got class: %s).",
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
