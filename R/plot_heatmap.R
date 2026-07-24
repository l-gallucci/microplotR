#' Clustered taxa heatmap
#'
#' Top-N taxa at `rank` x samples, CLR- or log10-transformed relative
#' abundance, with hierarchical clustering (UPGMA) on both axes -- samples
#' clustered by Bray-Curtis dissimilarity on relative abundance (the
#' community-ecology standard, e.g. `vegan::vegdist`), taxa clustered by
#' Euclidean distance on the transformed values shown -- matching the
#' `ampvis2::amp_heatmap` / `pheatmap` / `ComplexHeatmap` convention for
#' marker-gene heatmaps. Dendrograms are drawn on both margins by default
#' (`show_dendrogram = TRUE`); set to `FALSE` for a plain ordered tile grid.
#'
#' By default (`fix_taxonomy = TRUE`), unknown/missing values at `rank` are
#' first resolved via [mp_tax_fix()].
#'
#' @param data List with `feature_table`, `taxonomy`, `metadata`.
#' @param rank Taxonomy rank to plot (rows). Default `"Genus"`.
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
#' @param transform `"clr"` (centered log-ratio, computed per sample across
#'   the shown taxa) or `"log10"` (`log10(x + pseudocount)`).
#' @param pseudocount Added before log/CLR transform to handle zeros.
#' @param cluster_rows,cluster_cols Cluster taxa / samples.
#' @param hclust_method Agglomeration method passed to [stats::hclust()].
#' @param show_dendrogram Draw dendrograms on the margins (requires at least
#'   one of `cluster_rows`/`cluster_cols`).
#' @param fix_taxonomy Run [mp_tax_fix()] on the taxonomy table first.
#' @param tax_fix_ranks Passed through to [mp_tax_fix()]'s `ranks` argument.
#'   `NULL` (default) auto-detects standard 16S rank columns; pass
#'   explicitly when `data$taxonomy` is a non-taxonomic hierarchy (e.g. a
#'   functional annotation table).
#' @return A ggplot2 object (`show_dendrogram = FALSE`) or a patchwork object
#'   (`show_dendrogram = TRUE`).
#' @export
mp_taxa_heatmap <- function(data,
                             rank = "Genus",
                             top_n = 25,
                             min_rel_abund = 0,
                             min_prevalence = 0,
                             detection = 0,
                             transform = c("clr", "log10"),
                             pseudocount = NULL,
                             cluster_rows = TRUE,
                             cluster_cols = TRUE,
                             hclust_method = "average",
                             show_dendrogram = TRUE,
                             fix_taxonomy = TRUE,
                             tax_fix_ranks = NULL) {
  transform <- match.arg(transform)

  long <- .mp_long_abundance(data, fix_taxonomy = fix_taxonomy, tax_fix_ranks = tax_fix_ranks)
  long <- long |>
    dplyr::group_by(.data$Sample_ID) |>
    dplyr::mutate(rel_abund = 100 * .data$Count / sum(.data$Count)) |>
    dplyr::ungroup()

  agg <- long |>
    dplyr::group_by(.data$Sample_ID, .data[[rank]]) |>
    dplyr::summarise(rel_abund = sum(.data$rel_abund), .groups = "drop")

  survivors <- if (min_rel_abund > 0 || min_prevalence > 0) {
    mp_filter_taxa(agg, taxon_col = rank, min_rel_abund = min_rel_abund,
                   min_prevalence = min_prevalence, detection = detection)
  } else {
    unique(agg[[rank]])
  }

  totals <- agg[agg[[rank]] %in% survivors, , drop = FALSE] |>
    dplyr::group_by(.data[[rank]]) |>
    dplyr::summarise(total = sum(.data$rel_abund), .groups = "drop") |>
    dplyr::arrange(dplyr::desc(.data$total))
  keep <- if (is.null(top_n)) totals[[rank]] else utils::head(totals[[rank]], top_n)
  agg <- agg[agg[[rank]] %in% keep, , drop = FALSE]

  samples <- sort(unique(agg$Sample_ID))
  taxa <- keep
  rel_mat <- matrix(0, nrow = length(taxa), ncol = length(samples),
                     dimnames = list(taxa, samples))
  for (i in seq_len(nrow(agg))) {
    rel_mat[agg[[rank]][i], agg$Sample_ID[i]] <- agg$rel_abund[i]
  }

  pc <- pseudocount %||% (min(rel_mat[rel_mat > 0], na.rm = TRUE) / 2)
  if (transform == "log10") {
    val_mat <- log10(rel_mat + pc)
    legend_title <- "log10(rel. abundance %)"
  } else {
    log_mat <- log(rel_mat + pc)
    val_mat <- sweep(log_mat, 2, colMeans(log_mat), FUN = "-")
    legend_title <- "CLR abundance"
  }

  row_order <- seq_len(nrow(val_mat))
  col_order <- seq_len(ncol(val_mat))
  hc_rows <- hc_cols <- NULL

  if (cluster_rows && nrow(val_mat) > 2) {
    hc_rows <- stats::hclust(stats::dist(val_mat, method = "euclidean"), method = hclust_method)
    row_order <- hc_rows$order
  }
  if (cluster_cols && ncol(val_mat) > 2) {
    bray <- vegan::vegdist(t(rel_mat), method = "bray")
    hc_cols <- stats::hclust(bray, method = hclust_method)
    col_order <- hc_cols$order
  }

  val_mat <- val_mat[row_order, col_order, drop = FALSE]
  taxa_ord <- rownames(val_mat)
  sample_ord <- colnames(val_mat)

  tile_df <- as.data.frame(as.table(val_mat), stringsAsFactors = FALSE)
  names(tile_df) <- c("taxon", "sample", "value")
  tile_df$x <- match(tile_df$sample, sample_ord)
  tile_df$y <- match(tile_df$taxon, taxa_ord)

  n_row <- length(taxa_ord)
  n_col <- length(sample_ord)

  heat <- ggplot2::ggplot(tile_df, ggplot2::aes(x = .data$x, y = .data$y, fill = .data$value)) +
    ggplot2::geom_tile() +
    ggplot2::scale_fill_viridis_c(name = legend_title) +
    ggplot2::scale_x_continuous(breaks = seq_len(n_col), labels = sample_ord,
                                 limits = c(0.5, n_col + 0.5), expand = c(0, 0)) +
    ggplot2::scale_y_continuous(breaks = seq_len(n_row), labels = taxa_ord,
                                 limits = c(0.5, n_row + 0.5), expand = c(0, 0)) +
    ggplot2::labs(x = NULL, y = NULL) +
    mp_theme_pub() +
    ggplot2::theme(axis.line = ggplot2::element_blank(),
                   axis.ticks = ggplot2::element_blank(),
                   axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))

  if (!show_dendrogram || (is.null(hc_rows) && is.null(hc_cols))) {
    return(heat)
  }

  blank_dendro_theme <- ggplot2::theme_void() +
    ggplot2::theme(plot.margin = ggplot2::margin(0, 0, 0, 0))

  top <- patchwork::plot_spacer()
  if (!is.null(hc_cols)) {
    seg <- ggdendro::dendro_data(hc_cols)$segments
    top <- ggplot2::ggplot(seg) +
      ggplot2::geom_segment(ggplot2::aes(x = .data$x, y = .data$y, xend = .data$xend, yend = .data$yend)) +
      ggplot2::scale_x_continuous(limits = c(0.5, n_col + 0.5), expand = c(0, 0)) +
      blank_dendro_theme
  }

  left <- patchwork::plot_spacer()
  if (!is.null(hc_rows)) {
    seg <- ggdendro::dendro_data(hc_rows)$segments
    left <- ggplot2::ggplot(seg) +
      ggplot2::geom_segment(ggplot2::aes(x = .data$x, y = .data$y, xend = .data$xend, yend = .data$yend)) +
      ggplot2::coord_flip() +
      ggplot2::scale_x_continuous(limits = c(0.5, n_row + 0.5), expand = c(0, 0)) +
      ggplot2::scale_y_reverse() +
      blank_dendro_theme
  }

  layout <- "
#TTT
LHHH
LHHH
LHHH
"
  patchwork::wrap_plots(T = top, L = left, H = heat, design = layout) +
    patchwork::plot_layout(heights = c(1, 3, 3, 3), widths = c(1, 3, 3, 3))
}

`%||%` <- function(a, b) if (is.null(a)) b else a
