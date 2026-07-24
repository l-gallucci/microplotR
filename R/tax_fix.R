#' Fix unknown/missing taxonomy labels by falling back to the last known rank
#'
#' Genus-level (or any rank) plots otherwise show blank/NA/"uncultured" bars,
#' which is common in 16S taxonomy tables and misleading in a stacked
#' barplot. This fills an unknown value at a given rank with the last known
#' ancestor rank so every feature still gets a meaningful label -- the same
#' idea as `tax_fix()` in the microViz package (Barnett et al. 2021,
#' Bioinformatics). Cascading unknowns (e.g. both Genus and Species missing)
#' all fall back to the same true ancestor rather than chaining onto each
#' other's fabricated labels.
#'
#' @param taxonomy A taxonomy data frame/tibble with `Feature_ID` and rank columns.
#' @param ranks Rank columns, high-to-low order. Default: auto-detected from
#'   `c("Domain","Kingdom","Phylum","Class","Order","Family","Genus","Species")`
#'   intersected with `names(taxonomy)`.
#' @param unknown_strings Values (case-insensitive, whitespace-trimmed) treated
#'   as "unknown".
#' @param label_format Template for the fallback label; `{last}` is replaced
#'   with the last known ancestor value.
#' @return `taxonomy` with unknown values replaced.
#' @export
mp_tax_fix <- function(taxonomy,
                        ranks = NULL,
                        unknown_strings = c("", "na", "unknown", "unclassified",
                                             "uncultured", "unidentified", "metagenome"),
                        label_format = "Unclassified_{last}") {
  all_possible <- c("Domain", "Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")
  if (is.null(ranks)) ranks <- intersect(all_possible, names(taxonomy))
  if (length(ranks) == 0) {
    stop("No recognized taxonomy rank columns found in `taxonomy`. ",
         "Expected one or more of: ", paste(all_possible, collapse = ", "))
  }

  is_unknown <- function(x) is.na(x) | tolower(trimws(x)) %in% tolower(unknown_strings)

  n <- nrow(taxonomy)
  last_good <- rep(NA_character_, n)
  out <- taxonomy

  for (rank in ranks) {
    values <- as.character(taxonomy[[rank]])
    unk <- is_unknown(values)
    fixed <- values

    has_ancestor <- unk & !is.na(last_good)
    no_ancestor <- unk & is.na(last_good)

    if (any(has_ancestor)) {
      fixed[has_ancestor] <- vapply(last_good[has_ancestor], function(l) {
        gsub("{last}", l, label_format, fixed = TRUE)
      }, character(1), USE.NAMES = FALSE)
    }
    fixed[no_ancestor] <- "Unclassified"

    out[[rank]] <- fixed
    # Advance the ancestor pointer only where the real value was known --
    # a fixed/fabricated label never becomes the next rank's ancestor.
    last_good[!unk] <- values[!unk]
  }

  out
}
