#' Build a trivial one-column metadata stub (just Sample_ID) when the caller
#' doesn't supply one -- function-profile plots usually don't need grouping
#' metadata, but the taxa engine expects a metadata table to exist.
#' @keywords internal
#' @noRd
.mp_function_metadata_stub <- function(feature_table, metadata) {
  if (!is.null(metadata)) return(metadata)
  data.frame(Sample_ID = setdiff(names(feature_table), "Feature_ID"), stringsAsFactors = FALSE)
}

#' Functional profile barplot
#'
#' Stacked relative-abundance barplot for functional (not taxonomic) data --
#' gene family / KO / pathway abundance per sample, grouped/colored by a
#' broader functional category. Thin wrapper around [mp_taxa_barplot()]:
#' identical engine, different input names and defaults suited to
#' functional-annotation tools.
#'
#' Realistic input shape: a `feature_table` of per-sample gene/KO counts
#' (from any read quantifier, e.g. featureCounts/Salmon on predicted ORFs)
#' joined by `Feature_ID` to a `function_annotation` table from a gene
#' functional-annotation tool:
#' - **eggNOG-mapper** (`*.emapper.annotations`): use `rank = "KEGG_ko"`,
#'   `group_rank = "COG_category"` (defaults) -- eggNOG-mapper's own column
#'   names, unchanged.
#' - **KofamScan** (KO-only, no broader category): pass a `function_annotation`
#'   with just `Feature_ID`/`KO`/`KO_definition`, set `rank = "KO"`,
#'   `group_rank = NULL`, `nested_legend = FALSE` (flat legend, no grouping
#'   tier to nest under).
#'
#' Missing/blank category values (e.g. a gene with no COG hit) are resolved
#' by [mp_tax_fix()] the same way an unclassified genus is on the taxonomic
#' side -- see `fix_taxonomy`.
#'
#' @param feature_table Per-sample gene/KO/pathway abundance table
#'   (`Feature_ID` + one column per sample).
#' @param function_annotation `Feature_ID` + functional hierarchy columns
#'   (e.g. eggNOG-mapper's `COG_category`, `KEGG_ko`, `Description`).
#' @param metadata Optional sample metadata (only needed for `facet_var`/
#'   `sample_order` by metadata column); auto-stubbed from `feature_table`'s
#'   sample columns if `NULL`.
#' @param rank Functional level to plot (bar segments). Default `"KEGG_ko"`.
#' @param group_rank Broader functional category used to group/color/order
#'   `rank`. Default `"COG_category"`. Set `NULL` to disable (implies
#'   `nested_legend = FALSE`).
#' @param top_n,min_rel_abund,min_prevalence,detection,facet_var,sample_order,fix_taxonomy,nested_legend,other_label See [mp_taxa_barplot()] -- identical parameters and defaults.
#' @return A ggplot2 object.
#' @export
mp_function_barplot <- function(feature_table,
                                 function_annotation,
                                 metadata = NULL,
                                 rank = "KEGG_ko",
                                 group_rank = "COG_category",
                                 top_n = 10,
                                 min_rel_abund = 0,
                                 min_prevalence = 0,
                                 detection = 0,
                                 facet_var = NULL,
                                 sample_order = NULL,
                                 fix_taxonomy = TRUE,
                                 nested_legend = TRUE,
                                 other_label = "Other") {
  data <- list(feature_table = feature_table, taxonomy = function_annotation,
               metadata = .mp_function_metadata_stub(feature_table, metadata))
  mp_taxa_barplot(data, rank = rank, group_rank = group_rank, top_n = top_n,
                   min_rel_abund = min_rel_abund, min_prevalence = min_prevalence, detection = detection,
                   facet_var = facet_var, sample_order = sample_order, fix_taxonomy = fix_taxonomy,
                   tax_fix_ranks = c(group_rank, rank), nested_legend = nested_legend, other_label = other_label)
}

#' Functional profile heatmap
#'
#' Clustered heatmap for functional (not taxonomic) data. Thin wrapper
#' around [mp_taxa_heatmap()] -- identical engine, different input names and
#' defaults suited to functional-annotation tools (eggNOG-mapper/KofamScan,
#' see [mp_function_barplot()]).
#'
#' @param feature_table Per-sample gene/KO/pathway abundance table
#'   (`Feature_ID` + one column per sample).
#' @param function_annotation `Feature_ID` + functional hierarchy columns
#'   (e.g. eggNOG-mapper's `COG_category`, `KEGG_ko`, `Description`).
#' @param rank Functional level to plot (rows). Default `"KEGG_ko"`.
#' @param top_n,min_rel_abund,min_prevalence,detection,transform,pseudocount,cluster_rows,cluster_cols,hclust_method,show_dendrogram,fix_taxonomy See [mp_taxa_heatmap()] -- identical parameters and defaults.
#' @return A ggplot2 object (`show_dendrogram = FALSE`) or a patchwork object.
#' @export
mp_function_heatmap <- function(feature_table,
                                 function_annotation,
                                 rank = "KEGG_ko",
                                 top_n = 25,
                                 min_rel_abund = 0,
                                 min_prevalence = 0,
                                 detection = 0,
                                 transform = "clr",
                                 pseudocount = NULL,
                                 cluster_rows = TRUE,
                                 cluster_cols = TRUE,
                                 hclust_method = "average",
                                 show_dendrogram = TRUE,
                                 fix_taxonomy = TRUE) {
  data <- list(feature_table = feature_table, taxonomy = function_annotation,
               metadata = .mp_function_metadata_stub(feature_table, NULL))
  mp_taxa_heatmap(data, rank = rank, top_n = top_n, min_rel_abund = min_rel_abund,
                   min_prevalence = min_prevalence, detection = detection, transform = transform,
                   pseudocount = pseudocount, cluster_rows = cluster_rows, cluster_cols = cluster_cols,
                   hclust_method = hclust_method, show_dendrogram = show_dendrogram, fix_taxonomy = fix_taxonomy,
                   tax_fix_ranks = rank)
}
