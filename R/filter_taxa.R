#' Filter taxa by minimum abundance and/or prevalence
#'
#' Shared by [mp_taxa_barplot()], [mp_taxa_heatmap()], [mp_taxa_bubbleplot()].
#' Narrows the candidate taxa *before* `top_n` is applied: `top_n` then picks
#' the most abundant among survivors, so a taxon is excluded from the "top"
#' set either because it failed a threshold here or because it wasn't in the
#' top N of what passed -- both cases are treated identically by the calling
#' plot function (pooled into `other_label` for the barplot, simply dropped
#' for the heatmap/bubbleplot).
#'
#' @param agg A per-(Sample_ID, taxon) aggregated relative-abundance table
#'   (as built internally by the plot functions), with columns `Sample_ID`,
#'   `taxon_col`, `rel_abund`. Every taxon must have one row per sample
#'   (0 where absent) for prevalence/mean to be computed correctly.
#' @param taxon_col Name of the taxon column in `agg`.
#' @param min_rel_abund Minimum mean relative abundance (%) across all
#'   samples. `0` (default) disables this filter.
#' @param min_prevalence Minimum fraction of samples (`0-1`) -- or, if a
#'   value `> 1` is given, an absolute sample count -- in which the taxon
#'   must be "detected" (see `detection`). `0` (default) disables this
#'   filter.
#' @param detection Relative-abundance (%) threshold above which a taxon
#'   counts as "detected" in a sample, for `min_prevalence` purposes.
#'   Default `0` (any nonzero abundance counts as present).
#' @return Character vector of taxa (from `taxon_col`) passing both filters.
#' @export
mp_filter_taxa <- function(agg, taxon_col, min_rel_abund = 0, min_prevalence = 0, detection = 0) {
  n_samples <- length(unique(agg$Sample_ID))
  prevalence_threshold <- if (min_prevalence > 1) min_prevalence else min_prevalence * n_samples

  summary <- agg |>
    dplyr::group_by(dplyr::across(dplyr::all_of(taxon_col))) |>
    dplyr::summarise(
      mean_rel_abund = mean(.data$rel_abund),
      prevalence = sum(.data$rel_abund > detection),
      .groups = "drop"
    )

  summary[[taxon_col]][summary$mean_rel_abund >= min_rel_abund & summary$prevalence >= prevalence_threshold]
}
