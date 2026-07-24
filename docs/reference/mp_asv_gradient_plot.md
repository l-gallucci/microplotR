# ASV/taxon abundance over a continuous gradient

Relative abundance (%) of the top-N taxa at `rank` plotted against a
continuous metadata gradient (depth, pH, time...) – the standard way to
show which taxa respond to an environmental gradient (e.g.
[`phyloseq::psmelt()`](https://rdrr.io/pkg/phyloseq/man/psmelt.html) +
`ggplot(aes(x = gradient, y = Abundance, color = Taxon))`, common in
depth-profile / time-series microbiome papers).

## Usage

``` r
mp_asv_gradient_plot(
  data,
  gradient_var,
  rank = "Genus",
  group_rank = "Phylum",
  top_n = 10,
  min_rel_abund = 0,
  min_prevalence = 0,
  detection = 0,
  facet = FALSE,
  smooth = TRUE,
  smooth_method = "loess",
  fix_taxonomy = TRUE,
  tax_fix_ranks = NULL,
  other_label = "Other"
)
```

## Arguments

- data:

  List with `feature_table`, `taxonomy`, `metadata`.

- gradient_var:

  Metadata column holding the continuous gradient (e.g. `"Depth_m"`,
  `"pH"`). Must be numeric.

- rank:

  Taxonomy rank to plot. Default `"Genus"`.

- group_rank:

  Upper taxonomy rank used to group/color/order `rank` (nested-legend
  convention, see [`mp_taxa_barplot()`](mp_taxa_barplot.md)). Default
  `"Phylum"`. Set `NULL` for a flat legend.

- top_n:

  Number of `rank` taxa to show (ranked by total relative abundance
  among taxa passing `min_rel_abund`/`min_prevalence`); the rest are
  pooled into `other_label`. `NULL` shows every taxon that passes the
  filters, uncapped.

- min_rel_abund, min_prevalence, detection:

  See [`mp_taxa_barplot()`](mp_taxa_barplot.md).

- facet:

  Facet into one panel per taxon (small multiples) instead of overlaying
  all taxa in one panel with a color legend.

- smooth:

  Add a per-taxon trend line
  ([`ggplot2::geom_smooth()`](https://ggplot2.tidyverse.org/reference/geom_smooth.html)).

- smooth_method:

  Smoothing method passed to `geom_smooth()`. Default `"loess"`.

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

Same taxa filtering
(`top_n`/`min_rel_abund`/`min_prevalence`/`detection`) and nested-legend
grouping (`group_rank`) conventions as
[`mp_taxa_barplot()`](mp_taxa_barplot.md), so a barplot and gradient
plot of the same data select/color/order taxa consistently.

By default (`fix_taxonomy = TRUE`), unknown/missing values at `rank` are
first resolved via [`mp_tax_fix()`](mp_tax_fix.md).
