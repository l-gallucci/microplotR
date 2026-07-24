# Functional treemap

Thin wrapper around [`mp_taxa_treemap()`](mp_taxa_treemap.md) for
functional (not taxonomic) data – see
[`mp_function_barplot()`](mp_function_barplot.md) for the input shape
(eggNOG-mapper / KofamScan).

## Usage

``` r
mp_function_treemap(
  feature_table,
  function_annotation,
  rank = "KEGG_ko",
  group_rank = "COG_category",
  top_n = NULL,
  min_rel_abund = 0,
  min_prevalence = 0,
  detection = 0,
  fix_taxonomy = TRUE,
  other_label = "Other"
)
```

## Arguments

- feature_table:

  Per-sample gene/KO/pathway abundance table.

- function_annotation:

  `Feature_ID` + functional hierarchy columns.

- rank, group_rank, top_n, min_rel_abund, min_prevalence, detection,
  fix_taxonomy, other_label:

  See [`mp_taxa_treemap()`](mp_taxa_treemap.md) – identical parameters
  and defaults.

## Value

A ggplot2 object.
