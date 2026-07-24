# Filter taxa by minimum abundance and/or prevalence

Shared by [`mp_taxa_barplot()`](mp_taxa_barplot.md),
[`mp_taxa_heatmap()`](mp_taxa_heatmap.md),
[`mp_taxa_bubbleplot()`](mp_taxa_bubbleplot.md). Narrows the candidate
taxa *before* `top_n` is applied: `top_n` then picks the most abundant
among survivors, so a taxon is excluded from the "top" set either
because it failed a threshold here or because it wasn't in the top N of
what passed – both cases are treated identically by the calling plot
function (pooled into `other_label` for the barplot, simply dropped for
the heatmap/bubbleplot).

## Usage

``` r
mp_filter_taxa(
  agg,
  taxon_col,
  min_rel_abund = 0,
  min_prevalence = 0,
  detection = 0
)
```

## Arguments

- agg:

  A per-(Sample_ID, taxon) aggregated relative-abundance table (as built
  internally by the plot functions), with columns `Sample_ID`,
  `taxon_col`, `rel_abund`. Every taxon must have one row per sample (0
  where absent) for prevalence/mean to be computed correctly.

- taxon_col:

  Name of the taxon column in `agg`.

- min_rel_abund:

  Minimum mean relative abundance (%) across all samples. `0` (default)
  disables this filter.

- min_prevalence:

  Minimum fraction of samples (`0-1`) – or, if a value `> 1` is given,
  an absolute sample count – in which the taxon must be "detected" (see
  `detection`). `0` (default) disables this filter.

- detection:

  Relative-abundance (%) threshold above which a taxon counts as
  "detected" in a sample, for `min_prevalence` purposes. Default `0`
  (any nonzero abundance counts as present).

## Value

Character vector of taxa (from `taxon_col`) passing both filters.
