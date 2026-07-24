# Functional profile heatmap

Clustered heatmap for functional (not taxonomic) data. Thin wrapper
around [`mp_taxa_heatmap()`](mp_taxa_heatmap.md) – identical engine,
different input names and defaults suited to functional-annotation tools
(eggNOG-mapper/KofamScan, see
[`mp_function_barplot()`](mp_function_barplot.md)).

## Usage

``` r
mp_function_heatmap(
  feature_table,
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
  fix_taxonomy = TRUE
)
```

## Arguments

- feature_table:

  Per-sample gene/KO/pathway abundance table (`Feature_ID` + one column
  per sample).

- function_annotation:

  `Feature_ID` + functional hierarchy columns (e.g. eggNOG-mapper's
  `COG_category`, `KEGG_ko`, `Description`).

- rank:

  Functional level to plot (rows). Default `"KEGG_ko"`.

- top_n, min_rel_abund, min_prevalence, detection, transform,
  pseudocount, cluster_rows, cluster_cols, hclust_method,
  show_dendrogram, fix_taxonomy:

  See [`mp_taxa_heatmap()`](mp_taxa_heatmap.md) – identical parameters
  and defaults.

## Value

A ggplot2 object (`show_dendrogram = FALSE`) or a patchwork object.
