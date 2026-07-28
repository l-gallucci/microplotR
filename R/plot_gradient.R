#' ASV/taxon abundance over a continuous gradient
#'
#' Relative abundance (%) of the top-N taxa at `rank` plotted against a
#' continuous metadata gradient (depth, pH, time...) -- the standard way to
#' show which taxa respond to an environmental gradient (e.g.
#' `phyloseq::psmelt()` + `ggplot(aes(x = gradient, y = Abundance, color =
#' Taxon))`, common in depth-profile / time-series microbiome papers).
#'
#' Same taxa filtering (`top_n`/`min_rel_abund`/`min_prevalence`/`detection`)
#' and nested-legend grouping (`group_rank`) conventions as
#' [mp_taxa_barplot()], so a barplot and gradient plot of the same data
#' select/color/order taxa consistently.
#'
#' By default (`fix_taxonomy = TRUE`), unknown/missing values at `rank` are
#' first resolved via [mp_tax_fix()].
#'
#' @param data List with `feature_table`, `taxonomy`, `metadata`.
#' @param gradient_var Metadata column holding the continuous gradient
#'   (e.g. `"Depth_m"`, `"pH"`). Must be numeric.
#' @param rank Taxonomy rank to plot. Default `"Genus"`.
#' @param group_rank Upper taxonomy rank used to group/color/order `rank`
#'   (nested-legend convention, see [mp_taxa_barplot()]). Default `"Phylum"`.
#'   Set `NULL` for a flat legend.
#' @param top_n Number of `rank` taxa to show (ranked by total relative
#'   abundance among taxa passing `min_rel_abund`/`min_prevalence`); the
#'   rest are pooled into `other_label`. `NULL` shows every taxon that
#'   passes the filters, uncapped. Ignored if `taxa` is given.
#' @param taxa Optional character vector of exact `rank` values to show
#'   (e.g. `"Escherichia"` for a single taxon's gradient, or a short
#'   explicit list) -- bypasses `top_n`/`min_rel_abund`/`min_prevalence`
#'   ranking entirely and shows exactly these taxa (any not present in the
#'   data are silently ignored). `NULL` (default) falls back to the
#'   `top_n`-ranked selection.
#' @param min_rel_abund,min_prevalence,detection See [mp_taxa_barplot()].
#' @param facet Facet into one panel per taxon (small multiples) instead of
#'   overlaying all taxa in one panel with a color legend.
#' @param smooth Add a per-taxon trend line (`ggplot2::geom_smooth()`).
#' @param smooth_method Smoothing method passed to `geom_smooth()`. Default
#'   `"loess"`.
#' @param fix_taxonomy Run [mp_tax_fix()] on the taxonomy table first.
#' @param tax_fix_ranks Passed through to [mp_tax_fix()]'s `ranks` argument;
#'   see [mp_taxa_barplot()].
#' @param other_label Label for the pooled "everything else" bucket.
#' @param long Advanced: a precomputed long-format table (same shape
#'   `.mp_long_abundance()` returns) to reuse instead of re-melting/rejoining
#'   `data` from scratch. `NULL` (default) derives it from `data` as before;
#'   `fix_taxonomy`/`tax_fix_ranks` are ignored if `long` is supplied.
#' @return A ggplot2 object.
#' @export
mp_asv_gradient_plot <- function(data,
                                  gradient_var,
                                  rank = "Genus",
                                  group_rank = "Phylum",
                                  top_n = 10,
                                  taxa = NULL,
                                  min_rel_abund = 0,
                                  min_prevalence = 0,
                                  detection = 0,
                                  facet = FALSE,
                                  smooth = TRUE,
                                  smooth_method = "loess",
                                  fix_taxonomy = TRUE,
                                  tax_fix_ranks = NULL,
                                  other_label = "Other",
                                  long = NULL) {
  if (is.null(long)) long <- .mp_long_abundance(data, fix_taxonomy = fix_taxonomy, tax_fix_ranks = tax_fix_ranks)
  long <- long |>
    dplyr::group_by(.data$Sample_ID) |>
    dplyr::mutate(rel_abund = 100 * .data$Count / sum(.data$Count)) |>
    dplyr::ungroup()

  agg_cols <- unique(c("Sample_ID", rank, group_rank))
  agg <- long |>
    dplyr::group_by(dplyr::across(dplyr::all_of(agg_cols))) |>
    dplyr::summarise(rel_abund = sum(.data$rel_abund), .groups = "drop")

  survivors <- if (min_rel_abund > 0 || min_prevalence > 0) {
    mp_filter_taxa(agg, taxon_col = rank, min_rel_abund = min_rel_abund,
                   min_prevalence = min_prevalence, detection = detection)
  } else {
    unique(agg[[rank]])
  }

  totals <- agg[agg[[rank]] %in% survivors, , drop = FALSE] |>
    dplyr::group_by(dplyr::across(dplyr::all_of(rank))) |>
    dplyr::summarise(total = sum(.data$rel_abund), .groups = "drop") |>
    dplyr::arrange(dplyr::desc(.data$total))
  keep <- if (!is.null(taxa)) {
    intersect(taxa, unique(agg[[rank]]))
  } else if (is.null(top_n)) {
    totals[[rank]]
  } else {
    utils::head(totals[[rank]], top_n)
  }

  agg$.taxon <- ifelse(agg[[rank]] %in% keep, agg[[rank]], other_label)
  if (!is.null(group_rank)) {
    agg$.group <- ifelse(agg$.taxon == other_label, other_label, agg[[group_rank]])
  }

  bucket_cols <- unique(c("Sample_ID", ".taxon", if (!is.null(group_rank)) ".group"))
  agg <- agg |>
    dplyr::group_by(dplyr::across(dplyr::all_of(bucket_cols))) |>
    dplyr::summarise(rel_abund = sum(.data$rel_abund), .groups = "drop")

  taxon_totals <- agg |>
    dplyr::group_by(dplyr::across(dplyr::all_of(c(".taxon", if (!is.null(group_rank)) ".group")))) |>
    dplyr::summarise(total = sum(.data$rel_abund), .groups = "drop")

  if (!is.null(group_rank)) {
    # A taxon should map to exactly one group; collapse to the group with
    # the highest partial abundance if the data is inconsistent -- see
    # mp_taxa_barplot() for the same defensive fix.
    taxon_totals <- taxon_totals |>
      dplyr::group_by(.data$.taxon) |>
      dplyr::mutate(.full_total = sum(.data$total)) |>
      dplyr::slice_max(order_by = .data$total, n = 1, with_ties = FALSE) |>
      dplyr::ungroup() |>
      dplyr::mutate(total = .data$.full_total) |>
      dplyr::select(-".full_total")

    spec <- mp_nested_order_palette(
      taxa = taxon_totals$.taxon, groups = taxon_totals$.group, totals = taxon_totals$total,
      other_label = other_label
    )
  } else {
    ord <- taxon_totals[order(taxon_totals$.taxon == other_label, -taxon_totals$total), ".taxon", drop = TRUE]
    n <- length(ord)
    is_other <- ord == other_label
    cols <- character(n)
    cols[!is_other] <- colorspace::qualitative_hcl(sum(!is_other), palette = "Dark 3")
    cols[is_other] <- "grey60"
    spec <- list(order = ord, palette = stats::setNames(cols, ord), labels = stats::setNames(ord, ord))
  }

  agg$.taxon <- factor(agg$.taxon, levels = spec$order)
  agg <- dplyr::left_join(agg, data$metadata[, c("Sample_ID", gradient_var)], by = "Sample_ID")
  agg[[gradient_var]] <- suppressWarnings(as.numeric(agg[[gradient_var]]))

  p <- ggplot2::ggplot(agg, ggplot2::aes(x = .data[[gradient_var]], y = .data$rel_abund, color = .data$.taxon)) +
    ggplot2::geom_point(alpha = 0.7, size = 2) +
    ggplot2::scale_color_manual(values = spec$palette, breaks = spec$order, labels = spec$labels, name = NULL) +
    ggplot2::labs(x = gradient_var, y = "Relative abundance (%)") +
    mp_theme_pub()

  if (!is.null(group_rank)) {
    # Same 2-line "**Group**<br>taxon" legend labels as mp_taxa_barplot() --
    # needs the taller key height for the same reason (see there).
    p <- p + ggplot2::theme(
      legend.text = ggtext::element_markdown(size = ggplot2::rel(0.8), lineheight = 1.1),
      legend.key.height = grid::unit(1.8, "lines"),
      legend.spacing.y = grid::unit(0.15, "lines")
    )
  }

  if (smooth) {
    p <- p + ggplot2::geom_smooth(method = smooth_method, formula = y ~ x, se = FALSE, linewidth = 0.8)
  }

  if (facet) {
    p <- p + ggplot2::facet_wrap(ggplot2::vars(.data$.taxon), scales = "free_y") +
      ggplot2::theme(legend.position = "none")
  }

  p
}
