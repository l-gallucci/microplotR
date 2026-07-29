#' Alpha diversity boxplot
#'
#' Observed richness, Shannon, and Simpson diversity per sample (via
#' `vegan::specnumber()`/`vegan::diversity()`), boxplot + jitter by a
#' grouping metadata variable, faceted one panel per metric, with a
#' Wilcoxon (2 groups) or Kruskal-Wallis (>2 groups) test annotated via
#' `ggpubr::stat_compare_means()` -- the standard alpha-diversity comparison
#' figure in amplicon/marker-gene studies.
#'
#' @param data List with `feature_table`, `taxonomy`, `metadata`.
#' @param group_var Metadata column (categorical, >=2 levels) to compare.
#' @param metrics Diversity metrics to compute/plot, any of `"Observed"`,
#'   `"Shannon"`, `"Simpson"`. Default all three.
#' @param test `"auto"` (default: Wilcoxon for 2 groups, Kruskal-Wallis for
#'   more than 2), `"wilcoxon"`, `"kruskal"`, or `NULL` to disable the test
#'   annotation.
#' @return A ggplot2 object, faceted by metric.
#' @export
mp_alpha_diversity_plot <- function(data,
                                     group_var,
                                     metrics = c("Observed", "Shannon", "Simpson"),
                                     test = "auto") {
  ft <- data$feature_table
  sample_cols <- setdiff(names(ft), "Feature_ID")
  mat <- apply(as.matrix(ft[, sample_cols]), 2, as.numeric)
  comm <- t(mat) # samples x features

  df <- data.frame(Sample_ID = sample_cols, stringsAsFactors = FALSE)
  if ("Observed" %in% metrics) df$Observed <- vegan::specnumber(comm)
  if ("Shannon" %in% metrics) df$Shannon <- vegan::diversity(comm, index = "shannon")
  if ("Simpson" %in% metrics) df$Simpson <- vegan::diversity(comm, index = "simpson")

  metrics_present <- intersect(metrics, names(df))
  long <- tidyr::pivot_longer(df, cols = dplyr::all_of(metrics_present),
                               names_to = "metric", values_to = "value")
  long$metric <- factor(long$metric, levels = metrics_present)
  long <- dplyr::left_join(long, data$metadata[, c("Sample_ID", group_var)], by = "Sample_ID")

  # Samples with no group_var value can't be grouped/tested -- drop them
  # rather than let a phantom "NA" group show up as its own box/facet and
  # silently skew the group count the Wilcoxon/Kruskal-Wallis choice below
  # is based on.
  is_missing <- is.na(long[[group_var]]) | trimws(as.character(long[[group_var]])) == ""
  if (any(is_missing)) {
    warning(sprintf("Dropped %d sample-metric row(s) with missing '%s'.", sum(is_missing), group_var))
    long <- long[!is_missing, , drop = FALSE]
  }
  long[[group_var]] <- droplevels(as.factor(long[[group_var]]))

  n_groups <- length(unique(long[[group_var]]))
  test_method <- test
  if (identical(test, "auto")) {
    test_method <- if (n_groups == 2) "wilcox.test" else "kruskal.test"
  } else if (identical(test, "wilcoxon")) {
    test_method <- "wilcox.test"
  } else if (identical(test, "kruskal")) {
    test_method <- "kruskal.test"
  }

  p <- ggplot2::ggplot(long, ggplot2::aes(x = .data[[group_var]], y = .data$value)) +
    ggplot2::geom_boxplot(outlier.shape = NA, fill = "grey90") +
    ggplot2::geom_jitter(width = 0.15, height = 0, alpha = 0.7, size = 1.8) +
    ggplot2::facet_wrap(ggplot2::vars(.data$metric), scales = "free_y") +
    ggplot2::labs(x = NULL, y = "Diversity value") +
    mp_theme_pub()

  if (!is.null(test)) {
    p <- p + ggpubr::stat_compare_means(method = test_method, label = "p.format")
  }

  p
}
