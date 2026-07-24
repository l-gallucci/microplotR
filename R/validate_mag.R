#' Column name aliases accepted from CheckM (v1) and CheckM2, mapped to the
#' canonical names used internally. CheckM1's `qa`/`bin_stats` table uses
#' "Bin Id"; CheckM2's `quality_report.tsv` uses "Name". Both tools agree on
#' "Completeness"/"Contamination" already.
#' @keywords internal
#' @noRd
.mp_mag_alias_map <- list(
  Name = c("Name", "Bin Id", "Bin_Id", "bin_id"),
  Genome_Size = c("Genome_Size", "Genome size (bp)", "Genome_size"),
  Contig_N50 = c("Contig_N50", "N50 (contigs)", "N50"),
  GC_Content = c("GC_Content", "GC", "GC (%)")
)

#' Rename whichever alias of each known MAG-table column is present to its
#' canonical name (CheckM1 and CheckM2 column-naming differences).
#' @keywords internal
#' @noRd
.mp_mag_normalize <- function(mag_table) {
  for (canonical in names(.mp_mag_alias_map)) {
    if (canonical %in% names(mag_table)) next
    alias <- Find(function(a) a %in% names(mag_table), .mp_mag_alias_map[[canonical]])
    if (!is.null(alias)) names(mag_table)[names(mag_table) == alias] <- canonical
  }
  mag_table
}

#' Validate a MAG (metagenome-assembled genome) quality table
#'
#' Checks a single flat table shaped like CheckM2's `quality_report.tsv` or
#' CheckM (v1)'s `qa`/`bin_stats` table (one row per bin/MAG) -- not the
#' 3-file taxonomy/metadata schema used by [mp_validate()]. Column names from
#' either tool are accepted (e.g. CheckM1's `"Bin Id"` or CheckM2's `"Name"`)
#' and normalized internally. Mirrors the same `{level, file, field,
#' message}` finding shape.
#'
#' @param mag_table A data frame with at least a bin-identifier column
#'   (`Name` or CheckM1's `"Bin Id"`), `Completeness`, `Contamination`,
#'   and a genome-size column (`Genome_Size` or CheckM1's `"Genome size (bp)"`).
#' @return An `mp_validation_report` (see [mp_validate()]).
#' @export
mp_validate_mag <- function(mag_table) {
  mag_table <- .mp_mag_normalize(mag_table)
  findings <- list()
  add <- function(level, field, message) {
    findings[[length(findings) + 1]] <<- list(level = level, file = "mag_quality.tsv", field = field, message = message)
  }

  required <- c("Name", "Completeness", "Contamination", "Genome_Size")
  for (col in required) {
    if (!col %in% names(mag_table)) {
      add("error", col, sprintf("Required column '%s' not found.", col))
    }
  }
  if (any(vapply(findings, function(f) f$level == "error", logical(1)))) {
    return(.mp_report(findings))
  }

  dup_names <- unique(mag_table$Name[duplicated(mag_table$Name)])
  if (length(dup_names) > 0) {
    add("error", "Name", sprintf("Duplicate MAG Name value(s) found: %s.", .mp_vec(dup_names)))
  }

  completeness <- suppressWarnings(as.numeric(mag_table$Completeness))
  if (any(is.na(completeness) & !is.na(mag_table$Completeness))) {
    add("error", "Completeness", "Non-numeric value(s) in 'Completeness'.")
  } else if (any(completeness < 0 | completeness > 100, na.rm = TRUE)) {
    add("error", "Completeness", "Value(s) outside the valid [0, 100] range.")
  }

  contamination <- suppressWarnings(as.numeric(mag_table$Contamination))
  if (any(is.na(contamination) & !is.na(mag_table$Contamination))) {
    add("error", "Contamination", "Non-numeric value(s) in 'Contamination'.")
  } else if (any(contamination < 0, na.rm = TRUE)) {
    add("error", "Contamination", "Negative value(s) found; contamination cannot be negative.")
  }

  genome_size <- suppressWarnings(as.numeric(mag_table$Genome_Size))
  if (any(is.na(genome_size) & !is.na(mag_table$Genome_Size))) {
    add("error", "Genome_Size", "Non-numeric value(s) in 'Genome_Size'.")
  } else if (any(genome_size <= 0, na.rm = TRUE)) {
    add("error", "Genome_Size", "Value(s) must be positive (genome size in bp).")
  }

  .mp_report(findings)
}
