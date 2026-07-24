# Taxa bubble (dot) plot

One row per taxon at `rank` (top-N by total relative abundance), one
column per sample; dot size = relative abundance (%), dot color = the
upper taxonomy level `group_rank` – the
[`phyloseq::psmelt()`](https://rdrr.io/pkg/phyloseq/man/psmelt.html) +
ggplot bubble / MicrobiomeAnalyst dot-plot convention. Every sample x
taxon combination is drawn (size 0 where absent) so presence/absence
patterns read directly off the grid. Taxa (rows) are ordered grouped by
`group_rank` block (blocks and within-block taxa both ordered by
descending total abundance), matching the ordering convention used by
[`mp_taxa_barplot()`](mp_taxa_barplot.md).

## Usage

``` r
mp_taxa_bubbleplot(
  data,
  rank = "Genus",
  group_rank = "Phylum",
  top_n = 25,
  min_rel_abund = 0,
  min_prevalence = 0,
  detection = 0,
  facet_var = NULL,
  sample_order = NULL,
  fix_taxonomy = TRUE,
  tax_fix_ranks = NULL,
  max_size = 10
)
```

## Arguments

- data:

  List with `feature_table`, `taxonomy`, `metadata`.

- rank:

  Taxonomy rank to plot (rows). Default `"Genus"`.

- group_rank:

  Upper taxonomy rank used to color and order `rank`. Default
  `"Phylum"`. Set to `NULL` for a single-color plot.

- top_n:

  Number of `rank` taxa to show (ranked by total relative abundance
  among taxa passing `min_rel_abund`/`min_prevalence`). Set `NULL` to
  show every taxon that passes the filters, uncapped.

- min_rel_abund:

  Minimum mean relative abundance (%) across all samples for a taxon to
  be shown. `0` (default) disables this filter.

- min_prevalence:

  Minimum fraction of samples (`0-1`) – or, if `> 1`, an absolute sample
  count – in which a taxon must be detected (relative abundance
  `> detection`) to be shown. `0` (default) disables this filter.

- detection:

  Relative-abundance (%) threshold above which a taxon counts as
  "detected" in a sample, for `min_prevalence`. Default `0`.

- facet_var:

  Optional metadata column to facet by.

- sample_order:

  Either `NULL`, a single metadata column name to sort samples by, or an
  explicit character vector of `Sample_ID`.

- fix_taxonomy:

  Run [`mp_tax_fix()`](mp_tax_fix.md) on the taxonomy table first.

- tax_fix_ranks:

  Passed through to [`mp_tax_fix()`](mp_tax_fix.md)'s `ranks` argument.
  `NULL` (default) auto-detects standard 16S rank columns; pass
  explicitly when `data$taxonomy` is a non-taxonomic hierarchy.

- max_size:

  Maximum dot size (passed to
  [`ggplot2::scale_size_area()`](https://ggplot2.tidyverse.org/reference/scale_size.html)).

## Value

A ggplot2 object.

## Details

By default (`fix_taxonomy = TRUE`), unknown/missing values at `rank` are
first resolved via [`mp_tax_fix()`](mp_tax_fix.md).
