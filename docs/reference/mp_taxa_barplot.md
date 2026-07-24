# Stacked relative-abundance taxa barplot, optionally with a nested legend

Standard 16S/marker-gene stacked barplot: one bar per sample, segments
are relative abundance (%) of the top-N taxa at `rank` (everything else
pooled into `other_label`). With `nested_legend = TRUE` (default), taxa
are ordered and colored grouped by `group_rank` – each group gets a base
hue, its member taxa get shades of that hue, and the legend shows a bold
`**group_rank**` header above the first taxon of each group (rendered
via
[`ggtext::element_markdown()`](https://wilkelab.org/ggtext/reference/element_markdown.html)),
matching the grouped-legend convention used e.g. in the `ggnested`
package for microbiome barplots. Bars/legend are ordered by abundance at
`rank` (default Genus), grouped by `group_rank` (default Phylum)
abundance.

## Usage

``` r
mp_taxa_barplot(
  data,
  rank = "Genus",
  group_rank = "Phylum",
  top_n = 10,
  min_rel_abund = 0,
  min_prevalence = 0,
  detection = 0,
  facet_var = NULL,
  sample_order = NULL,
  fix_taxonomy = TRUE,
  tax_fix_ranks = NULL,
  nested_legend = TRUE,
  other_label = "Other"
)
```

## Arguments

- data:

  List with `feature_table`, `taxonomy`, `metadata` (as from
  [`mp_read_data()`](mp_read_data.md)).

- rank:

  Taxonomy rank to plot (bar segments). Default `"Genus"`.

- group_rank:

  Upper taxonomy rank used to group/color/order `rank` and (if
  `nested_legend = TRUE`) to build the bold legend headers. Default
  `"Phylum"`. Set to `NULL` to disable grouping entirely (implies
  `nested_legend = FALSE`).

- top_n:

  Number of `rank` taxa to show individually (ranked by total relative
  abundance among taxa passing `min_rel_abund`/`min_prevalence`); the
  rest are pooled into `other_label`. Set `NULL` to show every taxon
  that passes the filters, uncapped.

- min_rel_abund:

  Minimum mean relative abundance (%) across all samples for a taxon to
  be eligible; taxa below this are pooled into `other_label` regardless
  of `top_n`. `0` (default) disables this filter.

- min_prevalence:

  Minimum fraction of samples (`0-1`) – or, if `> 1`, an absolute sample
  count – in which a taxon must be detected (relative abundance
  `> detection`) to be eligible. `0` (default) disables this filter.

- detection:

  Relative-abundance (%) threshold above which a taxon counts as
  "detected" in a sample, for `min_prevalence`. Default `0`.

- facet_var:

  Optional metadata column to facet by (e.g. a treatment group), one
  panel per level.

- sample_order:

  Either `NULL` (keep the sample order from `feature_table`), a single
  metadata column name to sort samples by, or an explicit character
  vector of `Sample_ID` in the desired order.

- fix_taxonomy:

  Run [`mp_tax_fix()`](mp_tax_fix.md) on the taxonomy table first.

- tax_fix_ranks:

  Passed through to [`mp_tax_fix()`](mp_tax_fix.md)'s `ranks` argument.
  `NULL` (default) auto-detects standard 16S rank columns; pass
  explicitly (e.g. `c(group_rank, rank)`) when `data$taxonomy` is a
  non-taxonomic hierarchy (e.g. a functional annotation table).

- nested_legend:

  Use the grouped/bold-header legend. Requires `group_rank`.

- other_label:

  Label for the pooled "everything else" bucket (taxa failing
  `min_rel_abund`/`min_prevalence`, or outside `top_n`).

## Value

A ggplot2 object.

## Details

By default (`fix_taxonomy = TRUE`), unknown/missing values at `rank` are
first resolved via [`mp_tax_fix()`](mp_tax_fix.md) so no bar segment is
left unlabeled.

The plot carries no title/subtitle – only axis titles and the legend –
so it's ready to drop into a figure panel as-is.
