# Taxa (or functional) treemap

Hierarchical composition as nested rectangles – a static, publication-
ready alternative to an interactive Krona chart: area = relative
abundance (mean across samples, so areas sum to ~100%), grouped/colored
by `group_rank` with a subgroup border, each rectangle labeled directly
(no external legend needed). Works generically on taxonomic or
functional data (see [`mp_function_barplot()`](mp_function_barplot.md)
for the same generalization on the barplot/heatmap).

## Usage

``` r
mp_taxa_treemap(
  data,
  rank = "Genus",
  group_rank = "Phylum",
  top_n = NULL,
  min_rel_abund = 0,
  min_prevalence = 0,
  detection = 0,
  fix_taxonomy = TRUE,
  tax_fix_ranks = NULL,
  other_label = "Other"
)
```

## Arguments

- data:

  List with `feature_table`, `taxonomy`, `metadata`.

- rank:

  Taxonomy/functional level to plot (rectangles). Default `"Genus"`.

- group_rank:

  Upper level used to group/color `rank` and draw subgroup borders.
  Default `"Phylum"`. Set `NULL` to color by `rank` directly instead (no
  subgroup borders).

- top_n:

  Number of `rank` taxa to show individually (ranked by mean relative
  abundance among taxa passing `min_rel_abund`/`min_prevalence`); the
  rest are pooled into `other_label`. `NULL` (default) shows every taxon
  that passes the filters, uncapped.

- min_rel_abund, min_prevalence, detection:

  See [`mp_taxa_barplot()`](mp_taxa_barplot.md).

- fix_taxonomy:

  Run [`mp_tax_fix()`](mp_tax_fix.md) on the taxonomy table first.

- tax_fix_ranks:

  Passed through to [`mp_tax_fix()`](mp_tax_fix.md)'s `ranks` argument;
  see [`mp_taxa_barplot()`](mp_taxa_barplot.md).

- other_label:

  Label for the pooled "everything else" bucket.

## Value

A ggplot2 object.

## Details

By default (`fix_taxonomy = TRUE`), unknown/missing values at `rank` are
first resolved via [`mp_tax_fix()`](mp_tax_fix.md).
