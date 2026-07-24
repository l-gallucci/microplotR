#' @importFrom rlang .data
NULL

#' Melt feature table into long format, joined with taxonomy and metadata
#'
#' Internal helper shared by the taxa-level plot functions.
#'
#' @param data List with `feature_table`, `taxonomy`, `metadata` (as from
#'   [mp_read_data()]).
#' @param fix_taxonomy If `TRUE`, run [mp_tax_fix()] on `data$taxonomy` first.
#' @param tax_fix_ranks Passed through to [mp_tax_fix()]'s `ranks` argument.
#'   `NULL` (default) auto-detects standard 16S rank columns; pass explicitly
#'   when `data$taxonomy` holds a non-taxonomic hierarchy (e.g. a functional
#'   annotation table) whose column names aren't in the fixed 16S rank set.
#' @return A long tibble: `Feature_ID`, `Sample_ID`, `Count`, taxonomy rank
#'   columns, metadata columns.
#' @keywords internal
#' @noRd
.mp_long_abundance <- function(data, fix_taxonomy = TRUE, tax_fix_ranks = NULL) {
  taxonomy <- if (fix_taxonomy) mp_tax_fix(data$taxonomy, ranks = tax_fix_ranks) else data$taxonomy
  ft <- data$feature_table
  sample_cols <- setdiff(names(ft), "Feature_ID")

  long <- tidyr::pivot_longer(ft, cols = dplyr::all_of(sample_cols),
                               names_to = "Sample_ID", values_to = "Count")
  long$Count <- as.numeric(long$Count)
  long <- dplyr::left_join(long, taxonomy, by = "Feature_ID")
  long <- dplyr::left_join(long, data$metadata, by = "Sample_ID")
  long
}
