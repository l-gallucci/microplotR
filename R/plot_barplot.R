#' Stacked relative-abundance taxa barplot, optionally with a nested legend
#'
#' Standard 16S/marker-gene stacked barplot: one bar per sample, segments are
#' relative abundance (%) of the top-N taxa at `rank` (everything else pooled
#' into `other_label`). With `nested_legend = TRUE` (default), taxa are
#' ordered and colored grouped by `group_rank` -- each group gets a base hue,
#' its member taxa get shades of that hue, and the legend shows a bold
#' `**group_rank**` header above the first taxon of each group (rendered via
#' `ggtext::element_markdown()`), matching the grouped-legend convention used
#' e.g. in the `ggnested` package for microbiome barplots. Bars/legend are
#' ordered by abundance at `rank` (default Genus), grouped by `group_rank`
#' (default Phylum) abundance.
#'
#' By default (`fix_taxonomy = TRUE`), unknown/missing values at `rank` are
#' first resolved via [mp_tax_fix()] so no bar segment is left unlabeled.
#'
#' The plot carries no title/subtitle -- only axis titles and the legend --
#' so it's ready to drop into a figure panel as-is.
#'
#' @param data List with `feature_table`, `taxonomy`, `metadata` (as from
#'   [mp_read_data()]).
#' @param rank Taxonomy rank to plot (bar segments). Default `"Genus"`.
#' @param group_rank Upper taxonomy rank used to group/color/order `rank`
#'   and (if `nested_legend = TRUE`) to build the bold legend headers.
#'   Default `"Phylum"`. Set to `NULL` to disable grouping entirely (implies
#'   `nested_legend = FALSE`).
#' @param top_n Number of `rank` taxa to show individually (ranked by total
#'   relative abundance among taxa passing `min_rel_abund`/`min_prevalence`);
#'   the rest are pooled into `other_label`. Set `NULL` to show every taxon
#'   that passes the filters, uncapped.
#' @param min_rel_abund Minimum mean relative abundance (%) across all
#'   samples for a taxon to be eligible; taxa below this are pooled into
#'   `other_label` regardless of `top_n`. `0` (default) disables this filter.
#' @param min_prevalence Minimum fraction of samples (`0-1`) -- or, if
#'   `> 1`, an absolute sample count -- in which a taxon must be detected
#'   (relative abundance `> detection`) to be eligible. `0` (default)
#'   disables this filter.
#' @param detection Relative-abundance (%) threshold above which a taxon
#'   counts as "detected" in a sample, for `min_prevalence`. Default `0`.
#' @param facet_var Optional metadata column to facet by (e.g. a treatment
#'   group), one panel per level.
#' @param sample_order Either `NULL` (keep the sample order from
#'   `feature_table`), a single metadata column name to sort samples by, or
#'   an explicit character vector of `Sample_ID` in the desired order.
#' @param fix_taxonomy Run [mp_tax_fix()] on the taxonomy table first.
#' @param tax_fix_ranks Passed through to [mp_tax_fix()]'s `ranks` argument.
#'   `NULL` (default) auto-detects standard 16S rank columns; pass
#'   explicitly (e.g. `c(group_rank, rank)`) when `data$taxonomy` is a
#'   non-taxonomic hierarchy (e.g. a functional annotation table).
#' @param nested_legend Use the grouped/bold-header legend. Requires
#'   `group_rank`.
#' @param other_label Label for the pooled "everything else" bucket (taxa
#'   failing `min_rel_abund`/`min_prevalence`, or outside `top_n`).
#' @return A ggplot2 object.
#' @export
mp_taxa_barplot <- function(data,
                             rank = "Genus",
                             group_rank = "Phylum",
                             top_n = 10,
                             min_rel_abund = 0,
                             min_prevalence = 0,
                             detection = 0,
                             facet_var = NULL,
                             sample_order = NULL,
                             fix_taxonomy = TRUE,
                             tax_fix_ranks = NULL,
                             nested_legend = TRUE,
                             other_label = "Other") {
  if (nested_legend && is.null(group_rank)) {
    stop("`nested_legend = TRUE` requires `group_rank` (e.g. \"Phylum\"); ",
         "pass group_rank = NULL together with nested_legend = FALSE for a flat legend.")
  }

  long <- .mp_long_abundance(data, fix_taxonomy = fix_taxonomy, tax_fix_ranks = tax_fix_ranks)
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

  totals <- agg[agg[[rank]] %in% survivors, , drop = FALSE] |>
    dplyr::group_by(dplyr::across(dplyr::all_of(rank))) |>
    dplyr::summarise(total = sum(.data$rel_abund), .groups = "drop") |>
    dplyr::arrange(dplyr::desc(.data$total))
  keep <- if (is.null(top_n)) totals[[rank]] else utils::head(totals[[rank]], top_n)

  agg$.taxon <- ifelse(agg[[rank]] %in% keep, agg[[rank]], other_label)
  if (!is.null(group_rank)) {
    agg$.group <- ifelse(agg$.taxon == other_label, other_label, agg[[group_rank]])
  }

  bucket_cols <- unique(c("Sample_ID", ".taxon", if (!is.null(group_rank)) ".group", facet_var))
  agg <- agg |>
    dplyr::group_by(dplyr::across(dplyr::all_of(bucket_cols))) |>
    dplyr::summarise(rel_abund = sum(.data$rel_abund), .groups = "drop")

  taxon_totals <- agg |>
    dplyr::group_by(dplyr::across(dplyr::all_of(c(".taxon", if (!is.null(group_rank)) ".group")))) |>
    dplyr::summarise(total = sum(.data$rel_abund), .groups = "drop")

  if (!is.null(group_rank)) {
    # A taxon should map to exactly one group; if the data assigns it
    # inconsistent groups across rows (e.g. messy/partial annotation),
    # collapse to the group with the highest partial abundance, but keep the
    # taxon's full cross-group total so ranking/top_n still reflects it.
    taxon_totals <- taxon_totals |>
      dplyr::group_by(.data$.taxon) |>
      dplyr::mutate(.full_total = sum(.data$total)) |>
      dplyr::slice_max(order_by = .data$total, n = 1, with_ties = FALSE) |>
      dplyr::ungroup() |>
      dplyr::mutate(total = .data$.full_total) |>
      dplyr::select(-".full_total")
  }

  if (nested_legend) {
    spec <- mp_nested_order_palette(
      taxa = taxon_totals$.taxon,
      groups = taxon_totals$.group,
      totals = taxon_totals$total,
      other_label = other_label
    )
  } else {
    ord <- taxon_totals[order(taxon_totals$.taxon == other_label, -taxon_totals$total), ".taxon", drop = TRUE]
    n <- length(ord)
    is_other <- ord == other_label
    cols <- character(n)
    cols[!is_other] <- colorspace::qualitative_hcl(sum(!is_other), palette = "Dark 3")
    cols[is_other] <- "grey80"
    spec <- list(order = ord, palette = stats::setNames(cols, ord), labels = stats::setNames(ord, ord))
  }

  agg$.taxon <- factor(agg$.taxon, levels = spec$order)

  if (!is.null(sample_order)) {
    if (length(sample_order) == 1 && sample_order %in% names(data$metadata)) {
      sample_levels <- as.character(data$metadata$Sample_ID[order(data$metadata[[sample_order]])])
    } else {
      sample_levels <- sample_order
    }
    agg$Sample_ID <- factor(agg$Sample_ID, levels = sample_levels)
  }

  p <- ggplot2::ggplot(agg, ggplot2::aes(x = .data$Sample_ID, y = .data$rel_abund, fill = .data$.taxon)) +
    ggplot2::geom_col(width = 0.9) +
    ggplot2::scale_fill_manual(values = spec$palette, breaks = spec$order, labels = spec$labels) +
    ggplot2::labs(x = NULL, y = "Relative abundance (%)", fill = NULL) +
    ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0, 0.02))) +
    mp_theme_pub() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))

  p <- p + if (nested_legend) {
    ggplot2::theme(legend.text = ggtext::element_markdown(size = ggplot2::rel(0.8)))
  } else {
    ggplot2::theme(legend.text = ggplot2::element_text(size = ggplot2::rel(0.8)))
  }

  if (!is.null(facet_var)) {
    p <- p + ggplot2::facet_wrap(ggplot2::vars(.data[[facet_var]]), scales = "free_x")
  }

  p
}
