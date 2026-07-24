#' Validate tidy microbiome input (feature table + taxonomy + metadata)
#'
#' Mirrors the checks documented in data-format.md and implemented on the
#' Python side by `validate()` in microplotpy -- keep both in sync.
#'
#' @param data List as returned by [mp_read_data()], with `feature_table`,
#'   `taxonomy`, `metadata` tibbles.
#' @param gradient_column Optional metadata column to check is numeric
#'   (required by the ASV-level gradient plot).
#' @param group_column Optional metadata column to check has >= 2 levels
#'   (required by alpha/beta diversity group comparisons).
#' @param required_ranks Column names that must be present in `data$taxonomy`.
#'   Default `c("Phylum", "Genus")` for the standard 16S taxonomy shape.
#'   Pass `character(0)` when `data$taxonomy` is actually a non-taxonomic
#'   hierarchy (e.g. a functional annotation table from eggNOG-mapper/
#'   KofamScan, see [mp_function_barplot()]) whose column names aren't ranks.
#' @return An `mp_validation_report`: a list with `findings`, a data frame
#'   with columns `level`, `file`, `field`, `message`.
#' @export
mp_validate <- function(data, gradient_column = NULL, group_column = NULL,
                         required_ranks = c("Phylum", "Genus")) {
  ft <- data$feature_table
  tax <- data$taxonomy
  meta <- data$metadata
  findings <- list()

  add <- function(level, file, field, message) {
    findings[[length(findings) + 1]] <<- list(level = level, file = file, field = field, message = message)
  }

  if (!"Feature_ID" %in% names(ft)) {
    add("error", "feature_table.tsv", "Feature_ID", "Required column 'Feature_ID' not found.")
  }
  if (!"Feature_ID" %in% names(tax)) {
    add("error", "taxonomy.tsv", "Feature_ID", "Required column 'Feature_ID' not found.")
  }
  if (!"Sample_ID" %in% names(meta)) {
    add("error", "metadata.tsv", "Sample_ID", "Required column 'Sample_ID' not found.")
  }
  for (rank in required_ranks) {
    if (!rank %in% names(tax)) {
      add("error", "taxonomy.tsv", rank, sprintf("Required taxonomy rank column '%s' not found.", rank))
    }
  }

  if (any(vapply(findings, function(f) f$level == "error", logical(1)))) {
    return(.mp_report(findings))
  }

  # --- uniqueness ---
  dup_features_ft <- unique(ft$Feature_ID[duplicated(ft$Feature_ID)])
  if (length(dup_features_ft) > 0) {
    add("error", "feature_table.tsv", "Feature_ID",
        sprintf("Duplicate Feature_ID values found: %s.", .mp_vec(dup_features_ft)))
  }
  dup_features_tax <- unique(tax$Feature_ID[duplicated(tax$Feature_ID)])
  if (length(dup_features_tax) > 0) {
    add("error", "taxonomy.tsv", "Feature_ID",
        sprintf("Duplicate Feature_ID values found: %s.", .mp_vec(dup_features_tax)))
  }
  dup_samples <- unique(meta$Sample_ID[duplicated(meta$Sample_ID)])
  if (length(dup_samples) > 0) {
    add("error", "metadata.tsv", "Sample_ID",
        sprintf("Duplicate Sample_ID values found: %s.", .mp_vec(dup_samples)))
  }

  # --- sample ID match: feature_table columns vs metadata Sample_ID ---
  ft_samples <- setdiff(names(ft), "Feature_ID")
  meta_samples <- as.character(meta$Sample_ID)
  missing_in_ft <- sort(setdiff(meta_samples, ft_samples))
  missing_in_meta <- sort(setdiff(ft_samples, meta_samples))
  if (length(missing_in_ft) > 0) {
    add("error", "metadata.tsv", "Sample_ID",
        sprintf("Sample(s) %s found in metadata.tsv but missing from feature_table.tsv columns.",
                .mp_vec(missing_in_ft)))
  }
  if (length(missing_in_meta) > 0) {
    add("error", "feature_table.tsv", "Sample_ID",
        sprintf("Sample column(s) %s found in feature_table.tsv but missing from metadata.tsv Sample_ID.",
                .mp_vec(missing_in_meta)))
  }

  # --- feature ID match: feature_table vs taxonomy ---
  ft_features <- as.character(ft$Feature_ID)
  tax_features <- as.character(tax$Feature_ID)
  missing_in_tax <- sort(setdiff(ft_features, tax_features))
  missing_in_ft_feat <- sort(setdiff(tax_features, ft_features))
  if (length(missing_in_tax) > 0) {
    add("error", "taxonomy.tsv", "Feature_ID",
        sprintf("Feature(s) %s found in feature_table.tsv but missing from taxonomy.tsv.",
                .mp_vec(missing_in_tax)))
  }
  if (length(missing_in_ft_feat) > 0) {
    add("error", "feature_table.tsv", "Feature_ID",
        sprintf("Feature(s) %s found in taxonomy.tsv but missing from feature_table.tsv.",
                .mp_vec(missing_in_ft_feat)))
  }

  # --- numeric, non-negative abundance values ---
  sample_cols <- setdiff(names(ft), "Feature_ID")
  numeric_ft <- as.data.frame(lapply(ft[sample_cols], function(x) suppressWarnings(as.numeric(x))))
  present_mask <- as.data.frame(lapply(ft[sample_cols], function(x) !is.na(x)))
  non_numeric_mask <- is.na(numeric_ft) & present_mask
  if (any(unlist(non_numeric_mask))) {
    bad_cells <- .mp_bad_cells(ft$Feature_ID, sample_cols, non_numeric_mask)
    add("error", "feature_table.tsv", "values",
        sprintf("Non-numeric abundance value(s) at: %s.", .mp_vec(bad_cells)))
  } else {
    negative_mask <- numeric_ft < 0
    negative_mask[is.na(negative_mask)] <- FALSE
    if (any(unlist(negative_mask))) {
      bad_cells <- .mp_bad_cells(ft$Feature_ID, sample_cols, negative_mask)
      add("error", "feature_table.tsv", "values",
          sprintf("Negative abundance value(s) at: %s.", .mp_vec(bad_cells)))
    }

    # --- all-zero rows / columns (warning only) ---
    row_sums <- rowSums(numeric_ft, na.rm = TRUE)
    zero_features <- ft$Feature_ID[row_sums == 0]
    if (length(zero_features) > 0) {
      add("warning", "feature_table.tsv", "Feature_ID",
          sprintf("Feature(s) %s are all-zero across every sample.", .mp_vec(zero_features)))
    }
    col_sums <- colSums(numeric_ft, na.rm = TRUE)
    zero_samples <- sample_cols[col_sums == 0]
    if (length(zero_samples) > 0) {
      add("warning", "feature_table.tsv", "Sample_ID",
          sprintf("Sample(s) %s are all-zero across every feature.", .mp_vec(zero_samples)))
    }
  }

  # --- missing taxonomy values (warning, treated as Unclassified) ---
  for (rank in setdiff(names(tax), "Feature_ID")) {
    values <- tax[[rank]]
    missing_mask <- is.na(values) | trimws(values) == ""
    if (any(missing_mask)) {
      feats <- tax$Feature_ID[missing_mask]
      add("warning", "taxonomy.tsv", rank,
          sprintf("Missing '%s' for feature(s) %s; will be labeled 'Unclassified'.", rank, .mp_vec(feats)))
    }
  }

  # --- gradient column check ---
  if (!is.null(gradient_column)) {
    if (!gradient_column %in% names(meta)) {
      add("error", "metadata.tsv", gradient_column,
          sprintf("Gradient column '%s' not found in metadata.tsv.", gradient_column))
    } else {
      coerced <- suppressWarnings(as.numeric(meta[[gradient_column]]))
      present <- !is.na(meta[[gradient_column]])
      if (any(is.na(coerced) & present)) {
        add("error", "metadata.tsv", gradient_column,
            sprintf("Gradient column '%s' contains non-numeric value(s).", gradient_column))
      }
    }
  }

  # --- group column check ---
  if (!is.null(group_column)) {
    if (!group_column %in% names(meta)) {
      add("error", "metadata.tsv", group_column,
          sprintf("Group column '%s' not found in metadata.tsv.", group_column))
    } else {
      levels <- unique(stats::na.omit(meta[[group_column]]))
      if (length(levels) < 2) {
        add("error", "metadata.tsv", group_column,
            sprintf("Group column '%s' has %d level(s); need at least 2 for group comparison.",
                    group_column, length(levels)))
      }
    }
  }

  .mp_report(findings)
}

.mp_vec <- function(x, limit = 10, max_chars = 40) {
  # Truncates both the number of values shown and each value's own length --
  # unbounded here makes findings unreadable (and can hang the Shiny
  # sidebar) once Feature_ID/Sample_ID hold full DNA sequences (e.g. dada2
  # ASV tables) or datasets run into the thousands of features. The
  # underlying finding still carries the full untruncated vector; only this
  # display string is shortened.
  shown <- ifelse(nchar(x) > max_chars, paste0(substr(x, 1, max_chars), "..."), x)
  if (length(shown) > limit) shown <- c(shown[seq_len(limit)], sprintf("(and %d more)", length(x) - limit))
  paste0("[", paste(shown, collapse = ", "), "]")
}

.mp_bad_cells <- function(feature_ids, sample_cols, mask, limit = 10) {
  cells <- character(0)
  for (j in seq_along(sample_cols)) {
    rows <- which(mask[[j]])
    if (length(rows) > 0) {
      cells <- c(cells, sprintf("%s/%s", feature_ids[rows], sample_cols[j]))
    }
  }
  if (length(cells) > limit) cells <- c(cells[seq_len(limit)], "(truncated)")
  cells
}

.mp_report <- function(findings) {
  df <- if (length(findings) == 0) {
    data.frame(level = character(0), file = character(0), field = character(0), message = character(0))
  } else {
    do.call(rbind, lapply(findings, as.data.frame, stringsAsFactors = FALSE))
  }
  structure(list(findings = df), class = "mp_validation_report")
}

#' Errors from a validation report
#'
#' @param report An `mp_validation_report` from [mp_validate()].
#' @return A data frame of the `error`-level findings only.
#' @export
mp_errors <- function(report) report$findings[report$findings$level == "error", , drop = FALSE]

#' Warnings from a validation report
#'
#' @param report An `mp_validation_report` from [mp_validate()].
#' @return A data frame of the `warning`-level findings only.
#' @export
mp_warnings <- function(report) report$findings[report$findings$level == "warning", , drop = FALSE]

#' Whether a validation report has no errors
#'
#' @param report An `mp_validation_report` from [mp_validate()].
#' @return `TRUE` if there are no `error`-level findings (warnings are fine).
#' @export
mp_is_valid <- function(report) nrow(mp_errors(report)) == 0
