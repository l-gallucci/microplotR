#' Taxa (or functional) treemap
#'
#' Hierarchical composition as nested rectangles -- a static, publication-
#' ready alternative to an interactive Krona chart: area = relative
#' abundance (mean across samples, so areas sum to ~100%), grouped/colored
#' by `group_rank` with a subgroup border, each rectangle labeled directly
#' (no external legend needed). Works generically on taxonomic or functional
#' data (see [mp_function_barplot()] for the same generalization on the
#' barplot/heatmap).
#'
#' By default (`fix_taxonomy = TRUE`), unknown/missing values at `rank` are
#' first resolved via [mp_tax_fix()].
#'
#' @param data List with `feature_table`, `taxonomy`, `metadata`.
#' @param rank Taxonomy/functional level to plot (rectangles). Default `"Genus"`.
#' @param group_rank Upper level used to group/color `rank` and draw subgroup
#'   borders. Default `"Phylum"`. Set `NULL` to color by `rank` directly
#'   instead (no subgroup borders).
#' @param top_n Number of `rank` taxa to show individually (ranked by mean
#'   relative abundance among taxa passing `min_rel_abund`/`min_prevalence`);
#'   the rest are pooled into `other_label`. `NULL` (default) shows every
#'   taxon that passes the filters, uncapped.
#' @param min_rel_abund,min_prevalence,detection See [mp_taxa_barplot()].
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
mp_taxa_treemap <- function(data,
                             rank = "Genus",
                             group_rank = "Phylum",
                             top_n = NULL,
                             min_rel_abund = 0,
                             min_prevalence = 0,
                             detection = 0,
                             fix_taxonomy = TRUE,
                             tax_fix_ranks = NULL,
                             other_label = "Other",
                             long = NULL) {
  if (is.null(long)) long <- .mp_long_abundance(data, fix_taxonomy = fix_taxonomy, tax_fix_ranks = tax_fix_ranks)
  long <- long |>
    dplyr::group_by(.data$Sample_ID) |>
    dplyr::mutate(rel_abund = 100 * .data$Count / sum(.data$Count)) |>
    dplyr::ungroup()

  group_cols <- unique(c("Sample_ID", rank, group_rank))
  agg <- long |>
    dplyr::group_by(dplyr::across(dplyr::all_of(group_cols))) |>
    dplyr::summarise(rel_abund = sum(.data$rel_abund), .groups = "drop")

  survivors <- if (min_rel_abund > 0 || min_prevalence > 0) {
    mp_filter_taxa(agg, taxon_col = rank, min_rel_abund = min_rel_abund,
                   min_prevalence = min_prevalence, detection = detection)
  } else {
    unique(agg[[rank]])
  }

  # Mean per-sample relative abundance per taxon -- sums to ~100% across taxa
  # (linearity of the mean), so treemap area is directly comparable to the
  # barplot's per-sample relative abundance (%).
  taxon_cols <- unique(c(rank, group_rank))
  taxon_mean <- agg[agg[[rank]] %in% survivors, , drop = FALSE] |>
    dplyr::group_by(dplyr::across(dplyr::all_of(taxon_cols))) |>
    dplyr::summarise(rel_abund = mean(.data$rel_abund), .groups = "drop") |>
    dplyr::arrange(dplyr::desc(.data$rel_abund))

  keep <- if (is.null(top_n)) taxon_mean[[rank]] else utils::head(taxon_mean[[rank]], top_n)

  taxon_mean$.taxon <- ifelse(taxon_mean[[rank]] %in% keep, taxon_mean[[rank]], other_label)
  if (!is.null(group_rank)) {
    taxon_mean$.group <- ifelse(taxon_mean$.taxon == other_label, other_label, taxon_mean[[group_rank]])
  }

  bucket_cols <- unique(c(".taxon", if (!is.null(group_rank)) ".group"))
  final_df <- taxon_mean |>
    dplyr::group_by(dplyr::across(dplyr::all_of(bucket_cols))) |>
    dplyr::summarise(rel_abund = sum(.data$rel_abund), .groups = "drop")

  if (!is.null(group_rank)) {
    groups <- unique(final_df$.group)
    n_groups <- length(groups)
    group_colors <- stats::setNames(colorspace::qualitative_hcl(n_groups, palette = "Dark 3"), groups)
    if (other_label %in% names(group_colors)) group_colors[other_label] <- "grey70"

    p <- ggplot2::ggplot(final_df, ggplot2::aes(area = .data$rel_abund, fill = .data$.group,
                                                 subgroup = .data$.group, label = .data$.taxon)) +
      treemapify::geom_treemap() +
      treemapify::geom_treemap_subgroup_border(colour = "white", size = 2) +
      treemapify::geom_treemap_subgroup_text(place = "topleft", grow = FALSE, alpha = 0.7,
                                              colour = "white", fontface = "bold") +
      treemapify::geom_treemap_text(colour = "white", place = "center", grow = FALSE, reflow = TRUE) +
      ggplot2::scale_fill_manual(values = group_colors, name = NULL)
  } else {
    n <- nrow(final_df)
    is_other <- final_df$.taxon == other_label
    cols <- character(n)
    cols[!is_other] <- colorspace::qualitative_hcl(sum(!is_other), palette = "Dark 3")
    cols[is_other] <- "grey70"
    taxon_colors <- stats::setNames(cols, final_df$.taxon)

    p <- ggplot2::ggplot(final_df, ggplot2::aes(area = .data$rel_abund, fill = .data$.taxon, label = .data$.taxon)) +
      treemapify::geom_treemap() +
      treemapify::geom_treemap_text(colour = "white", place = "center", grow = FALSE, reflow = TRUE) +
      ggplot2::scale_fill_manual(values = taxon_colors, name = NULL)
  }

  p + mp_theme_pub() +
    ggplot2::theme(legend.position = "none", axis.line = ggplot2::element_blank(),
                   axis.ticks = ggplot2::element_blank(), axis.text = ggplot2::element_blank())
}

#' Functional treemap
#'
#' Thin wrapper around [mp_taxa_treemap()] for functional (not taxonomic)
#' data -- see [mp_function_barplot()] for the input shape (eggNOG-mapper /
#' KofamScan).
#'
#' @param feature_table Per-sample gene/KO/pathway abundance table.
#' @param function_annotation `Feature_ID` + functional hierarchy columns.
#' @param rank,group_rank,top_n,min_rel_abund,min_prevalence,detection,fix_taxonomy,other_label See [mp_taxa_treemap()] -- identical parameters and defaults.
#' @return A ggplot2 object.
#' @export
mp_function_treemap <- function(feature_table,
                                 function_annotation,
                                 rank = "KEGG_ko",
                                 group_rank = "COG_category",
                                 top_n = NULL,
                                 min_rel_abund = 0,
                                 min_prevalence = 0,
                                 detection = 0,
                                 fix_taxonomy = TRUE,
                                 other_label = "Other") {
  data <- list(feature_table = feature_table, taxonomy = function_annotation,
               metadata = .mp_function_metadata_stub(feature_table, NULL))
  mp_taxa_treemap(data, rank = rank, group_rank = group_rank, top_n = top_n,
                   min_rel_abund = min_rel_abund, min_prevalence = min_prevalence, detection = detection,
                   fix_taxonomy = fix_taxonomy, tax_fix_ranks = c(group_rank, rank), other_label = other_label)
}
