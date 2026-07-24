# Functional profile barplot

Stacked relative-abundance barplot for functional (not taxonomic) data –
gene family / KO / pathway abundance per sample, grouped/colored by a
broader functional category. Thin wrapper around
[`mp_taxa_barplot()`](mp_taxa_barplot.md): identical engine, different
input names and defaults suited to functional-annotation tools.

## Usage

``` r
mp_function_barplot(
  feature_table,
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
  other_label = "Other"
)
```

## Arguments

- feature_table:

  Per-sample gene/KO/pathway abundance table (`Feature_ID` + one column
  per sample).

- function_annotation:

  `Feature_ID` + functional hierarchy columns (e.g. eggNOG-mapper's
  `COG_category`, `KEGG_ko`, `Description`).

- metadata:

  Optional sample metadata (only needed for `facet_var`/ `sample_order`
  by metadata column); auto-stubbed from `feature_table`'s sample
  columns if `NULL`.

- rank:

  Functional level to plot (bar segments). Default `"KEGG_ko"`.

- group_rank:

  Broader functional category used to group/color/order `rank`. Default
  `"COG_category"`. Set `NULL` to disable (implies
  `nested_legend = FALSE`).

- top_n, min_rel_abund, min_prevalence, detection, facet_var,
  sample_order, fix_taxonomy, nested_legend, other_label:

  See [`mp_taxa_barplot()`](mp_taxa_barplot.md) – identical parameters
  and defaults.

## Value

A ggplot2 object.

## Details

Realistic input shape: a `feature_table` of per-sample gene/KO counts
(from any read quantifier, e.g. featureCounts/Salmon on predicted ORFs)
joined by `Feature_ID` to a `function_annotation` table from a gene
functional-annotation tool:

- **eggNOG-mapper** (`*.emapper.annotations`): use `rank = "KEGG_ko"`,
  `group_rank = "COG_category"` (defaults) – eggNOG-mapper's own column
  names, unchanged.

- **KofamScan** (KO-only, no broader category): pass a
  `function_annotation` with just `Feature_ID`/`KO`/`KO_definition`, set
  `rank = "KO"`, `group_rank = NULL`, `nested_legend = FALSE` (flat
  legend, no grouping tier to nest under).

Missing/blank category values (e.g. a gene with no COG hit) are resolved
by [`mp_tax_fix()`](mp_tax_fix.md) the same way an unclassified genus is
on the taxonomic side – see `fix_taxonomy`.
