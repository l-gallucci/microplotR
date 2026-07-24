#' Assembly Nx curve
#'
#' Cumulative contig-length distribution per assembly -- the standard
#' QUAST-style plot: x = Nx percentile (0-100, cumulative % of total
#' assembly length), y = contig length (log10 scale), one line per assembly.
#'
#' @param contig_lengths A data frame with `Assembly_ID`, `Contig_ID`,
#'   `Length` columns (as validated by [mp_validate_contig_lengths()]).
#' @return A ggplot2 object.
#' @export
mp_assembly_nx_plot <- function(contig_lengths) {
  contig_lengths$Length <- as.numeric(contig_lengths$Length)

  nx <- do.call(rbind, lapply(split(contig_lengths, contig_lengths$Assembly_ID), function(d) {
    lengths <- sort(d$Length, decreasing = TRUE)
    total <- sum(lengths)
    cum_pct <- cumsum(lengths) / total * 100
    data.frame(Assembly_ID = d$Assembly_ID[1], x = cum_pct, y = lengths)
  }))

  ggplot2::ggplot(nx, ggplot2::aes(x = .data$x, y = .data$y, color = .data$Assembly_ID)) +
    ggplot2::geom_step() +
    ggplot2::scale_y_log10() +
    ggplot2::labs(x = "Nx (% of total assembly length)", y = "Contig length (bp)", color = NULL) +
    mp_theme_pub()
}

#' Assembly summary barplot
#'
#' A single summary statistic (default `N50`) per assembly, from a
#' QUAST-`report.tsv`-shaped table.
#'
#' @param assembly_summary A data frame with an `Assembly_ID` column and the
#'   column named by `stat_col` (as validated by
#'   [mp_validate_assembly_summary()]).
#' @param stat_col Which column to plot. Default `"N50"`.
#' @return A ggplot2 object.
#' @export
mp_assembly_summary_barplot <- function(assembly_summary, stat_col = "N50") {
  assembly_summary <- .mp_assembly_normalize(assembly_summary, .mp_assembly_summary_alias_map)
  if (!stat_col %in% names(assembly_summary)) {
    stop(sprintf("stat_col '%s' not found in assembly_summary.", stat_col))
  }
  assembly_summary[[stat_col]] <- as.numeric(assembly_summary[[stat_col]])

  ggplot2::ggplot(assembly_summary, ggplot2::aes(x = .data$Assembly_ID, y = .data[[stat_col]])) +
    ggplot2::geom_col(fill = "steelblue", width = 0.7) +
    ggplot2::labs(x = NULL, y = gsub("_", " ", stat_col)) +
    ggplot2::scale_y_continuous(labels = scales::label_comma()) +
    mp_theme_pub() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))
}
