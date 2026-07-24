#' Taxa bubble (dot) plot
#'
#' One row per taxon at `rank` (top-N by total relative abundance), one
#' column per sample; dot size = relative abundance (%), dot color = the
#' upper taxonomy level `group_rank` -- the `phyloseq::psmelt()` + ggplot
#' bubble / MicrobiomeAnalyst dot-plot convention. Every sample x taxon
#' combination is drawn (size 0 where absent) so presence/absence patterns
#' read directly off the grid. Taxa (rows) are ordered grouped by
#' `group_rank` block (blocks and within-block taxa both ordered by
#' descending total abundance), matching the ordering convention used by
#' [mp_taxa_barplot()].
#'
#' By default (`fix_taxonomy = TRUE`), unknown/missing values at `rank` are
#' first resolved via [mp_tax_fix()].
#'
#' @param data List with `feature_table`, `taxonomy`, `metadata`.
#' @param rank Taxonomy rank to plot (rows). Default `"Genus"`.
#' @param group_rank Upper taxonomy rank used to color and order `rank`.
#'   Default `"Phylum"`. Set to `NULL` for a single-color plot.
#' @param top_n Number of `rank` taxa to show (ranked by total relative
#'   abundance among taxa passing `min_rel_abund`/`min_prevalence`). Set
#'   `NULL` to show every taxon that passes the filters, uncapped.
#' @param min_rel_abund Minimum mean relative abundance (%) across all
#'   samples for a taxon to be shown. `0` (default) disables this filter.
#' @param min_prevalence Minimum fraction of samples (`0-1`) -- or, if
#'   `> 1`, an absolute sample count -- in which a taxon must be detected
#'   (relative abundance `> detection`) to be shown. `0` (default) disables
#'   this filter.
#' @param detection Relative-abundance (%) threshold above which a taxon
#'   counts as "detected" in a sample, for `min_prevalence`. Default `0`.
#' @param facet_var Optional metadata column to facet by.
#' @param sample_order Either `NULL`, a single metadata column name to sort
#'   samples by, or an explicit character vector of `Sample_ID`.
#' @param fix_taxonomy Run [mp_tax_fix()] on the taxonomy table first.
#' @param tax_fix_ranks Passed through to [mp_tax_fix()]'s `ranks` argument.
#'   `NULL` (default) auto-detects standard 16S rank columns; pass
#'   explicitly when `data$taxonomy` is a non-taxonomic hierarchy.
#' @param max_size Maximum dot size (passed to `ggplot2::scale_size_area()`).
#' @param long Advanced: a precomputed long-format table (same shape
#'   `.mp_long_abundance()` returns) to reuse instead of re-melting/rejoining
#'   `data` from scratch. `NULL` (default) derives it from `data` as before;
#'   `fix_taxonomy`/`tax_fix_ranks` are ignored if `long` is supplied.
#' @return A ggplot2 object.
#' @export
mp_taxa_bubbleplot <- function(data,
                                rank = "Genus",
                                group_rank = "Phylum",
                                top_n = 25,
                                min_rel_abund = 0,
                                min_prevalence = 0,
                                detection = 0,
                                facet_var = NULL,
                                sample_order = NULL,
                                fix_taxonomy = TRUE,
                                tax_fix_ranks = NULL,
                                max_size = 10,
                                long = NULL) {
  if (is.null(long)) long <- .mp_long_abundance(data, fix_taxonomy = fix_taxonomy, tax_fix_ranks = tax_fix_ranks)
  long <- long |>
    dplyr::group_by(.data$Sample_ID) |>
    dplyr::mutate(rel_abund = 100 * .data$Count / sum(.data$Count)) |>
    dplyr::ungroup()

  agg_cols <- unique(c("Sample_ID", rank, group_rank, facet_var))
  agg <- long |>
    dplyr::group_by(dplyr::across(dplyr::all_of(agg_cols))) |>
    dplyr::summarise(rel_abund = sum(.data$rel_abund), .groups = "drop")

  survivors <- if (min_rel_abund > 0 || min_prevalence > 0) {
    mp_filter_taxa(agg, taxon_col = rank, min_rel_abund = min_rel_abund,
                   min_prevalence = min_prevalence, detection = detection)
  } else {
    unique(agg[[rank]])
  }

  taxon_totals <- agg[agg[[rank]] %in% survivors, , drop = FALSE] |>
    dplyr::group_by(dplyr::across(dplyr::all_of(c(rank, group_rank)))) |>
    dplyr::summarise(total = sum(.data$rel_abund), .groups = "drop") |>
    dplyr::arrange(dplyr::desc(.data$total))
  keep <- if (is.null(top_n)) taxon_totals[[rank]] else utils::head(taxon_totals[[rank]], top_n)
  taxon_totals <- taxon_totals[taxon_totals[[rank]] %in% keep, , drop = FALSE]

  if (!is.null(group_rank) && anyDuplicated(taxon_totals[[rank]])) {
    # A taxon should map to exactly one group; if the data assigns it
    # inconsistent groups across rows (e.g. messy/partial annotation),
    # collapse to the group with the highest partial abundance, but keep the
    # taxon's full cross-group total so ranking still reflects it.
    taxon_totals <- taxon_totals |>
      dplyr::group_by(dplyr::across(dplyr::all_of(rank))) |>
      dplyr::mutate(.full_total = sum(.data$total)) |>
      dplyr::slice_max(order_by = .data$total, n = 1, with_ties = FALSE) |>
      dplyr::ungroup() |>
      dplyr::mutate(total = .data$.full_total) |>
      dplyr::select(-".full_total")
  }

  if (!is.null(group_rank)) {
    spec <- mp_nested_order_palette(
      taxa = taxon_totals[[rank]], groups = taxon_totals[[group_rank]], totals = taxon_totals$total,
      other_label = "none" # sentinel that will never match, so no "Other" bucketing happens here
    )
    taxon_order <- spec$order
    group_of_taxon <- stats::setNames(taxon_totals[[group_rank]], taxon_totals[[rank]])
    groups <- unique(group_of_taxon[taxon_order])
    group_palette <- stats::setNames(colorspace::qualitative_hcl(length(groups), palette = "Dark 3"), groups)
  } else {
    taxon_order <- taxon_totals[[rank]][order(-taxon_totals$total)]
    group_palette <- NULL
  }

  # complete grid: every sample x kept-taxon combination, 0 where absent
  samples <- unique(agg$Sample_ID)
  grid <- tidyr::expand_grid(Sample_ID = samples, .taxon_key = taxon_order)
  names(grid)[names(grid) == ".taxon_key"] <- rank

  agg_small <- agg[, c("Sample_ID", rank, "rel_abund")]
  plot_df <- dplyr::left_join(grid, agg_small, by = c("Sample_ID", rank))
  plot_df$rel_abund[is.na(plot_df$rel_abund)] <- 0

  if (!is.null(group_rank)) {
    plot_df[[group_rank]] <- group_of_taxon[plot_df[[rank]]]
  }
  if (!is.null(facet_var)) {
    facet_lookup <- stats::setNames(data$metadata[[facet_var]], as.character(data$metadata$Sample_ID))
    plot_df[[facet_var]] <- facet_lookup[as.character(plot_df$Sample_ID)]
  }

  plot_df[[rank]] <- factor(plot_df[[rank]], levels = rev(taxon_order))

  if (!is.null(sample_order)) {
    if (length(sample_order) == 1 && sample_order %in% names(data$metadata)) {
      sample_levels <- as.character(data$metadata$Sample_ID[order(data$metadata[[sample_order]])])
    } else {
      sample_levels <- sample_order
    }
    plot_df$Sample_ID <- factor(plot_df$Sample_ID, levels = sample_levels)
  }

  mapping <- if (!is.null(group_rank)) {
    ggplot2::aes(x = .data$Sample_ID, y = .data[[rank]], size = .data$rel_abund, color = .data[[group_rank]])
  } else {
    ggplot2::aes(x = .data$Sample_ID, y = .data[[rank]], size = .data$rel_abund)
  }

  p <- ggplot2::ggplot(plot_df, mapping) +
    ggplot2::geom_point() +
    ggplot2::scale_size_area(name = "Relative abundance (%)", max_size = max_size) +
    ggplot2::labs(x = NULL, y = NULL) +
    mp_theme_pub() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))

  if (!is.null(group_rank)) {
    p <- p + ggplot2::scale_color_manual(values = group_palette, name = NULL)
  }

  if (!is.null(facet_var)) {
    p <- p + ggplot2::facet_wrap(ggplot2::vars(.data[[facet_var]]), scales = "free_x")
  }

  p
}
