#' Column name aliases accepted from QUAST's `report.tsv`, mapped to the
#' canonical names used internally.
#' @keywords internal
#' @noRd
.mp_assembly_summary_alias_map <- list(
  Assembly_ID = c("Assembly_ID", "Assembly"),
  N_contigs = c("N_contigs", "# contigs"),
  Largest_contig = c("Largest_contig", "Largest contig"),
  Total_length = c("Total_length", "Total length"),
  N50 = c("N50"),
  L50 = c("L50"),
  GC_percent = c("GC_percent", "GC (%)", "GC")
)

#' @keywords internal
#' @noRd
.mp_assembly_normalize <- function(df, alias_map) {
  for (canonical in names(alias_map)) {
    if (canonical %in% names(df)) next
    alias <- Find(function(a) a %in% names(df), alias_map[[canonical]])
    if (!is.null(alias)) names(df)[names(df) == alias] <- canonical
  }
  df
}

#' Validate a per-contig length table
#'
#' Checks a tidy `Assembly_ID`, `Contig_ID`, `Length` table (one row per
#' contig) -- the input [mp_assembly_nx_plot()] needs. Not the QUAST
#' `report.tsv` summary shape (see [mp_validate_assembly_summary()]).
#'
#' @param contig_lengths A data frame with `Assembly_ID`, `Contig_ID`,
#'   `Length` columns.
#' @return An `mp_validation_report` (see [mp_validate()]).
#' @export
mp_validate_contig_lengths <- function(contig_lengths) {
  findings <- list()
  add <- function(level, field, message) {
    findings[[length(findings) + 1]] <<- list(level = level, file = "contig_lengths.tsv", field = field, message = message)
  }

  required <- c("Assembly_ID", "Contig_ID", "Length")
  for (col in required) {
    if (!col %in% names(contig_lengths)) add("error", col, sprintf("Required column '%s' not found.", col))
  }
  if (any(vapply(findings, function(f) f$level == "error", logical(1)))) return(.mp_report(findings))

  dup <- contig_lengths[duplicated(contig_lengths[c("Assembly_ID", "Contig_ID")]), c("Assembly_ID", "Contig_ID")]
  if (nrow(dup) > 0) {
    add("error", "Contig_ID", sprintf("Duplicate Contig_ID within an Assembly_ID: %s.",
                                       .mp_vec(paste(dup$Assembly_ID, dup$Contig_ID, sep = "/"))))
  }

  length_num <- suppressWarnings(as.numeric(contig_lengths$Length))
  if (any(is.na(length_num) & !is.na(contig_lengths$Length))) {
    add("error", "Length", "Non-numeric value(s) in 'Length'.")
  } else if (any(length_num <= 0, na.rm = TRUE)) {
    add("error", "Length", "Value(s) must be positive (contig length in bp).")
  }

  .mp_report(findings)
}

#' Validate an assembly summary table
#'
#' Checks a QUAST `report.tsv`-shaped table (one row per assembly). Column
#' names from QUAST's own report (e.g. `"# contigs"`, `"Total length"`) are
#' accepted and normalized internally.
#'
#' @param assembly_summary A data frame with an assembly-identifier column
#'   (`Assembly_ID` or QUAST's `"Assembly"`), and `N50`, `Total_length`
#'   (or QUAST's `"Total length"`).
#' @return An `mp_validation_report` (see [mp_validate()]).
#' @export
mp_validate_assembly_summary <- function(assembly_summary) {
  assembly_summary <- .mp_assembly_normalize(assembly_summary, .mp_assembly_summary_alias_map)
  findings <- list()
  add <- function(level, field, message) {
    findings[[length(findings) + 1]] <<- list(level = level, file = "assembly_summary.tsv", field = field, message = message)
  }

  required <- c("Assembly_ID", "N50", "Total_length")
  for (col in required) {
    if (!col %in% names(assembly_summary)) add("error", col, sprintf("Required column '%s' not found.", col))
  }
  if (any(vapply(findings, function(f) f$level == "error", logical(1)))) return(.mp_report(findings))

  dup <- unique(assembly_summary$Assembly_ID[duplicated(assembly_summary$Assembly_ID)])
  if (length(dup) > 0) {
    add("error", "Assembly_ID", sprintf("Duplicate Assembly_ID value(s) found: %s.", .mp_vec(dup)))
  }

  for (col in c("N50", "Total_length")) {
    vals <- suppressWarnings(as.numeric(assembly_summary[[col]]))
    if (any(is.na(vals) & !is.na(assembly_summary[[col]]))) {
      add("error", col, sprintf("Non-numeric value(s) in '%s'.", col))
    } else if (any(vals <= 0, na.rm = TRUE)) {
      add("error", col, sprintf("Value(s) in '%s' must be positive.", col))
    }
  }

  .mp_report(findings)
}
