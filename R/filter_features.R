#' Filter features (ASVs/OTUs) from a MicrobiomeData object
#'
#' Removes whole features from `feature_table`/`taxonomy` before any
#' plot-level aggregation -- e.g. dropping chloroplast/mitochondrial
#' contaminants, singleton ASVs, or an explicit exclude list -- and
#' reports how much was removed (feature count and read count, in
#' absolute terms and percentage) so filtering decisions aren't made
#' blind. Complements [mp_filter_taxa()], which filters *after*
#' aggregation to a chosen rank for a specific plot; this one operates on
#' raw features and feeds a cleaned dataset into [mp_validate()] or the
#' plot functions.
#'
#' @param data List with `feature_table`, `taxonomy`, `metadata` (as from
#'   [mp_read_data()]).
#' @param exclude_taxa Character vector of taxon names to drop wherever
#'   they appear (exact match, case/whitespace-insensitive) in any of
#'   `exclude_ranks` -- e.g. `c("Chloroplast", "Mitochondria")` to remove
#'   the classic 16S plastid/organelle contaminants regardless of which
#'   rank column a given classifier put them in.
#' @param exclude_ranks Which `data$taxonomy` columns to check for
#'   `exclude_taxa` matches. `NULL` (default) checks every column except
#'   `Feature_ID`.
#' @param exclude_features Character vector of exact `Feature_ID` values
#'   to drop directly, regardless of taxonomy (e.g. a denoiser's own
#'   contaminant/chimera flag list).
#' @param min_total_count Minimum total read count (summed across all
#'   samples) a feature must have to be kept. `0` (default) disables this
#'   filter.
#' @param remove_singletons If `TRUE`, also drop features with a total
#'   read count of exactly 1 (the classic amplicon "singleton") --
#'   equivalent to raising `min_total_count` to at least `2`.
#' @return A list with:
#'   - `data`: the filtered `feature_table`/`taxonomy` (same shape/columns,
#'     `metadata` untouched), ready for [mp_validate()] or the plot
#'     functions.
#'   - `summary`: a tibble, one row per active filter reason plus a
#'     `"total"` row (removed by any reason, i.e. not simply the sum, so
#'     features hit by more than one reason aren't double-counted), with
#'     `n_features_removed`/`pct_features_removed` and
#'     `n_reads_removed`/`pct_reads_removed`. Percentages are always
#'     against the *original* (pre-filter) totals, so they reflect what
#'     filtering actually took out of the input.
#' @export
mp_filter_features <- function(data,
                                exclude_taxa = NULL,
                                exclude_ranks = NULL,
                                exclude_features = NULL,
                                min_total_count = 0,
                                remove_singletons = FALSE) {
  ft <- data$feature_table
  tax <- data$taxonomy
  sample_cols <- setdiff(names(ft), "Feature_ID")

  counts <- as.data.frame(lapply(ft[sample_cols], as.numeric), check.names = FALSE)
  total_count <- rowSums(counts, na.rm = TRUE)
  total_reads_all <- sum(total_count)
  n_features_all <- nrow(ft)

  reasons <- list()

  if (!is.null(exclude_features)) {
    reasons[["exclude_features"]] <- ft$Feature_ID %in% exclude_features
  }

  if (!is.null(exclude_taxa)) {
    ranks <- if (is.null(exclude_ranks)) setdiff(names(tax), "Feature_ID") else intersect(exclude_ranks, names(tax))
    needle <- tolower(trimws(exclude_taxa))
    hit_by_id <- stats::setNames(rep(FALSE, nrow(tax)), tax$Feature_ID)
    for (r in ranks) {
      hit_by_id <- hit_by_id | (tolower(trimws(tax[[r]])) %in% needle)
    }
    reasons[["exclude_taxa"]] <- unname(hit_by_id[ft$Feature_ID])
  }

  threshold <- if (remove_singletons) max(min_total_count, 2) else min_total_count
  if (threshold > 0) {
    reasons[["min_total_count"]] <- total_count < threshold
  }

  remove <- if (length(reasons) == 0) rep(FALSE, nrow(ft)) else Reduce(`|`, reasons)
  remove[is.na(remove)] <- FALSE

  .row <- function(reason, hit) {
    hit[is.na(hit)] <- FALSE
    tibble::tibble(
      reason = reason,
      n_features_removed = sum(hit),
      pct_features_removed = if (n_features_all > 0) 100 * sum(hit) / n_features_all else NA_real_,
      n_reads_removed = sum(total_count[hit]),
      pct_reads_removed = if (total_reads_all > 0) 100 * sum(total_count[hit]) / total_reads_all else NA_real_
    )
  }

  summary <- dplyr::bind_rows(
    lapply(names(reasons), function(nm) .row(nm, reasons[[nm]])),
    .row("total", remove)
  )

  keep_ids <- ft$Feature_ID[!remove]
  list(
    data = list(
      feature_table = ft[!remove, , drop = FALSE],
      taxonomy = tax[tax$Feature_ID %in% keep_ids, , drop = FALSE],
      metadata = data$metadata
    ),
    summary = summary
  )
}
