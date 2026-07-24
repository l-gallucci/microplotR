#' Build ordering, nested color palette, and nested legend labels for a
#' taxon-colored plot
#'
#' Orders taxa grouped by an upper taxonomy level (each group's taxa kept
#' contiguous, groups and within-group taxa both ordered by descending
#' abundance), assigns each group a base hue and shades its member taxa from
#' that hue (so related taxa read as a family of colors), and builds legend
#' labels where the first taxon in each group carries a bold `**group**`
#' markdown header (rendered via `ggtext::element_markdown()`).
#' `other_label` (e.g. "Other") is always ordered last and colored grey, with
#' a plain (non-bold) label.
#'
#' @param taxa Character vector, one entry per unique taxon (e.g. genus).
#' @param groups Character vector, same length as `taxa`, the upper taxonomy
#'   level (e.g. phylum) each taxon belongs to.
#' @param totals Numeric vector, same length as `taxa`, overall abundance
#'   used to rank groups and within-group taxa.
#' @param other_label The sentinel value used for the aggregated "everything
#'   else" bucket, if present in `taxa`.
#' @param other_color Color for `other_label`.
#' @return A list with `order` (character vector, plotting/legend order),
#'   `palette` (named character vector, taxon -> hex color), `labels` (named
#'   character vector, taxon -> markdown legend label).
#' @export
mp_nested_order_palette <- function(taxa, groups, totals,
                                     other_label = "Other", other_color = "grey80") {
  stopifnot(length(taxa) == length(groups), length(taxa) == length(totals))

  is_other <- taxa == other_label
  main_taxa <- taxa[!is_other]
  main_groups <- groups[!is_other]
  main_totals <- totals[!is_other]

  group_totals <- stats::aggregate(main_totals, by = list(group = main_groups), FUN = sum)
  group_order <- group_totals$group[order(-group_totals$x)]

  ord <- order(match(main_groups, group_order), -main_totals)
  main_taxa <- main_taxa[ord]
  main_groups <- main_groups[ord]

  order <- main_taxa
  if (any(is_other)) order <- c(order, other_label)

  n_groups <- length(group_order)
  base_hues <- if (n_groups > 0) colorspace::qualitative_hcl(n_groups, palette = "Dark 3") else character(0)
  names(base_hues) <- group_order

  palette <- stats::setNames(character(length(order)), order)
  labels <- stats::setNames(character(length(order)), order)

  for (g in group_order) {
    members <- main_taxa[main_groups == g]
    n <- length(members)
    shades <- if (n == 1) base_hues[[g]] else colorspace::lighten(base_hues[[g]], amount = seq(0, 0.55, length.out = n))
    palette[members] <- shades
    labels[members] <- members
    labels[members[1]] <- sprintf("**%s**<br>%s", g, members[1])
  }

  if (any(is_other)) {
    palette[other_label] <- other_color
    labels[other_label] <- other_label
  }

  list(order = order, palette = palette, labels = labels)
}
