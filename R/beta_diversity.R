#' Beta diversity ordination plot
#'
#' Bray-Curtis or Jaccard dissimilarity (on relative abundance), PCoA
#' (`stats::cmdscale()`) or NMDS (`vegan::metaMDS()`) ordination, colored by
#' a grouping metadata variable with 95% confidence ellipses, PERMANOVA
#' (`vegan::adonis2()`) annotated as a plot subtitle -- the standard beta-
#' diversity ordination figure in amplicon/marker-gene studies. No
#' phylogenetic (UniFrac) distances -- Bray-Curtis/Jaccard only.
#'
#' @param data List with `feature_table`, `taxonomy`, `metadata`.
#' @param group_var Metadata column (categorical) to color/test by.
#' @param method Dissimilarity method: `"bray"` (default) or `"jaccard"`.
#' @param ordination `"pcoa"` (default) or `"nmds"`.
#' @param show_ellipse Draw 95% confidence ellipses per group
#'   (`ggplot2::stat_ellipse()`).
#' @param permanova Run and annotate PERMANOVA (`vegan::adonis2()`).
#' @param permutations Number of permutations for PERMANOVA. Default 999.
#' @return A ggplot2 object.
#' @export
mp_beta_diversity_plot <- function(data,
                                    group_var,
                                    method = c("bray", "jaccard"),
                                    ordination = c("pcoa", "nmds"),
                                    show_ellipse = TRUE,
                                    permanova = TRUE,
                                    permutations = 999) {
  method <- match.arg(method)
  ordination <- match.arg(ordination)

  ft <- data$feature_table
  sample_cols <- setdiff(names(ft), "Feature_ID")
  mat <- apply(as.matrix(ft[, sample_cols]), 2, as.numeric)
  mat <- sweep(mat, 2, colSums(mat), "/") * 100
  comm <- t(mat)
  rownames(comm) <- sample_cols

  dist <- vegan::vegdist(comm, method = method)
  meta <- data$metadata[match(sample_cols, data$metadata$Sample_ID), , drop = FALSE]

  stress_label <- NULL
  if (ordination == "pcoa") {
    pco <- stats::cmdscale(dist, k = 2, eig = TRUE)
    points <- as.data.frame(pco$points)
    names(points) <- c("Axis1", "Axis2")
    var_explained <- pco$eig[1:2] / sum(pco$eig[pco$eig > 0]) * 100
    xlab <- sprintf("PCoA1 (%.1f%%)", var_explained[1])
    ylab <- sprintf("PCoA2 (%.1f%%)", var_explained[2])
  } else {
    nmds <- vegan::metaMDS(dist, k = 2, trace = FALSE)
    points <- as.data.frame(nmds$points)
    names(points) <- c("Axis1", "Axis2")
    xlab <- "NMDS1"
    ylab <- "NMDS2"
    stress_label <- sprintf("stress = %.3f", nmds$stress)
  }

  points$Sample_ID <- sample_cols
  points <- dplyr::left_join(points, meta[, c("Sample_ID", group_var)], by = "Sample_ID")

  p <- ggplot2::ggplot(points, ggplot2::aes(x = .data$Axis1, y = .data$Axis2, color = .data[[group_var]])) +
    ggplot2::geom_point(size = 2.5, alpha = 0.85) +
    ggplot2::labs(x = xlab, y = ylab, color = NULL) +
    mp_theme_pub() +
    ggplot2::theme(plot.subtitle = ggplot2::element_text(size = ggplot2::rel(0.85))) # PERMANOVA/stress is
    # essential statistical content here, not decorative chrome -- re-enable
    # the subtitle that mp_theme_pub() otherwise blanks.

  if (show_ellipse) {
    p <- p + ggplot2::stat_ellipse(level = 0.95)
  }

  caption_parts <- character(0)
  if (permanova) {
    adonis_df <- data.frame(group = points[[group_var]])
    perm <- vegan::adonis2(dist ~ group, data = adonis_df, permutations = permutations)
    r2 <- perm$R2[1]
    pval <- perm$`Pr(>F)`[1]
    pval_str <- if (pval < 0.001) "< 0.001" else sprintf("= %.3f", pval)
    caption_parts <- c(caption_parts, sprintf("PERMANOVA R\u00b2 = %.3f, p %s", r2, pval_str))
  }
  if (!is.null(stress_label)) caption_parts <- c(caption_parts, stress_label)

  if (length(caption_parts) > 0) {
    p <- p + ggplot2::labs(subtitle = paste(caption_parts, collapse = "; "))
  }

  p
}
