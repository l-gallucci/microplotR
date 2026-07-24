#' MAG (metagenome-assembled genome) quality plot
#'
#' One dot per MAG/bin, x = Contamination (%), y = Completeness (%) -- the
#' standard first figure in any MAG-recovery paper. Accepts CheckM2's
#' `quality_report.tsv` or CheckM (v1)'s `qa`/`bin_stats` table directly
#' (column-name differences between the two tools are normalized
#' internally). Draws dashed reference lines and a shaded region at the
#' MIMAG quality-tier thresholds by default (Bowers et al. 2017, *Nature
#' Biotechnology* -- "Minimum information about a single amplified genome
#' and a metagenome-assembled genome"): high quality >=90% completeness &
#' <5% contamination; medium quality >=50% completeness & <10% contamination.
#'
#' @param mag_table A data frame with a bin-identifier column (`Name` or
#'   CheckM1's `"Bin Id"`), `Completeness`, `Contamination`, and optionally
#'   a genome-size column (`Genome_Size` or CheckM1's `"Genome size (bp)"`).
#' @param size_col Column to map to dot size (e.g. genome size). `NULL`
#'   disables sizing (uniform dots). Default `"Genome_Size"` (silently
#'   ignored if not present in `mag_table`).
#' @param color_col Optional column to map to dot color (e.g. a taxon or
#'   group column).
#' @param show_mimag_thresholds Draw the MIMAG quality-tier reference lines
#'   and shaded high-quality region.
#' @return A ggplot2 object.
#' @export
mp_mag_quality_plot <- function(mag_table,
                                 size_col = "Genome_Size",
                                 color_col = NULL,
                                 show_mimag_thresholds = TRUE) {
  mag_table <- .mp_mag_normalize(mag_table)
  mag_table$Completeness <- suppressWarnings(as.numeric(mag_table$Completeness))
  mag_table$Contamination <- suppressWarnings(as.numeric(mag_table$Contamination))

  if (!is.null(size_col) && !size_col %in% names(mag_table)) size_col <- NULL
  if (!is.null(color_col) && !color_col %in% names(mag_table)) {
    stop(sprintf("color_col '%s' not found in mag_table.", color_col))
  }

  p <- ggplot2::ggplot(mag_table, ggplot2::aes(x = .data$Contamination, y = .data$Completeness))

  if (show_mimag_thresholds) {
    p <- p +
      ggplot2::annotate("rect", xmin = -Inf, xmax = 5, ymin = 90, ymax = Inf,
                         fill = "grey85", alpha = 0.4) +
      ggplot2::geom_hline(yintercept = c(50, 90), linetype = "dashed", color = "grey50") +
      ggplot2::geom_vline(xintercept = c(5, 10), linetype = "dashed", color = "grey50")
  }

  p <- p + if (!is.null(size_col) && !is.null(color_col)) {
    ggplot2::geom_point(ggplot2::aes(size = .data[[size_col]], color = .data[[color_col]]), alpha = 0.85)
  } else if (!is.null(size_col)) {
    ggplot2::geom_point(ggplot2::aes(size = .data[[size_col]]), alpha = 0.85, color = "steelblue")
  } else if (!is.null(color_col)) {
    ggplot2::geom_point(ggplot2::aes(color = .data[[color_col]]), alpha = 0.85, size = 2.5)
  } else {
    ggplot2::geom_point(alpha = 0.85, size = 2.5, color = "steelblue")
  }

  if (!is.null(size_col)) {
    p <- p + ggplot2::scale_size_area(name = gsub("_", " ", size_col), max_size = 10,
                                       labels = scales::label_comma())
  }

  p + ggplot2::labs(x = "Contamination (%)", y = "Completeness (%)", color = NULL) +
    mp_theme_pub()
}

#' MAG completeness / contamination distribution
#'
#' Faceted histogram of Completeness and Contamination across all MAGs in
#' `mag_table` -- the standard supplementary figure showing the overall
#' quality spread of a recovered MAG set. Same MIMAG threshold lines as
#' [mp_mag_quality_plot()], drawn per facet. Accepts CheckM2 or CheckM (v1)
#' column names (normalized internally).
#'
#' @param mag_table A data frame with a bin-identifier column, `Completeness`,
#'   `Contamination`.
#' @param bins Number of histogram bins per facet.
#' @param show_mimag_thresholds Draw the MIMAG quality-tier reference lines.
#' @return A ggplot2 object.
#' @export
mp_mag_quality_distribution <- function(mag_table, bins = 20, show_mimag_thresholds = TRUE) {
  mag_table <- .mp_mag_normalize(mag_table)
  completeness <- suppressWarnings(as.numeric(mag_table$Completeness))
  contamination <- suppressWarnings(as.numeric(mag_table$Contamination))

  long <- data.frame(
    metric = rep(c("Completeness (%)", "Contamination (%)"), each = length(completeness)),
    value = c(completeness, contamination)
  )

  p <- ggplot2::ggplot(long, ggplot2::aes(x = .data$value)) +
    ggplot2::geom_histogram(bins = bins, fill = "steelblue", color = "white") +
    ggplot2::facet_wrap(ggplot2::vars(.data$metric), scales = "free_x") +
    ggplot2::labs(x = NULL, y = "Number of MAGs") +
    mp_theme_pub()

  if (show_mimag_thresholds) {
    vlines <- data.frame(
      metric = c("Completeness (%)", "Completeness (%)", "Contamination (%)", "Contamination (%)"),
      xintercept = c(50, 90, 5, 10)
    )
    p <- p + ggplot2::geom_vline(data = vlines, ggplot2::aes(xintercept = .data$xintercept),
                                  linetype = "dashed", color = "grey40")
  }

  p
}
