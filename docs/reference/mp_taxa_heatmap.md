# Clustered taxa heatmap

Top-N taxa at `rank` x samples, CLR- or log10-transformed relative
abundance, with hierarchical clustering (UPGMA) on both axes – samples
clustered by Bray-Curtis dissimilarity on relative abundance (the
community-ecology standard, e.g.
[`vegan::vegdist`](https://vegandevs.github.io/vegan/reference/vegdist.html)),
taxa clustered by Euclidean distance on the transformed values shown –
matching the `ampvis2::amp_heatmap` / `pheatmap` / `ComplexHeatmap`
convention for marker-gene heatmaps. Dendrograms are drawn on both
margins by default (`show_dendrogram = TRUE`); set to `FALSE` for a
plain ordered tile grid.

## Usage

``` r
mp_taxa_heatmap(
  data,
  rank = "Genus",
  top_n = 25,
  min_rel_abund = 0,
  min_prevalence = 0,
  detection = 0,
  transform = c("clr", "log10"),
  pseudocount = NULL,
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  hclust_method = "average",
  show_dendrogram = TRUE,
  fix_taxonomy = TRUE,
  tax_fix_ranks = NULL
)
```

## Arguments

- data:

  List with `feature_table`, `taxonomy`, `metadata`.

- rank:

  Taxonomy rank to plot (rows). Default `"Genus"`.

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

- transform:

  `"clr"` (centered log-ratio, computed per sample across the shown
  taxa) or `"log10"` (`log10(x + pseudocount)`).

- pseudocount:

  Added before log/CLR transform to handle zeros.

- cluster_rows, cluster_cols:

  Cluster taxa / samples.

- hclust_method:

  Agglomeration method passed to
  [`stats::hclust()`](https://rdrr.io/r/stats/hclust.html).

- show_dendrogram:

  Draw dendrograms on the margins (requires at least one of
  `cluster_rows`/`cluster_cols`).

- fix_taxonomy:

  Run [`mp_tax_fix()`](mp_tax_fix.md) on the taxonomy table first.

- tax_fix_ranks:

  Passed through to [`mp_tax_fix()`](mp_tax_fix.md)'s `ranks` argument.
  `NULL` (default) auto-detects standard 16S rank columns; pass
  explicitly when `data$taxonomy` is a non-taxonomic hierarchy (e.g. a
  functional annotation table).

## Value

A ggplot2 object (`show_dendrogram = FALSE`) or a patchwork object
(`show_dendrogram = TRUE`).

## Details

By default (`fix_taxonomy = TRUE`), unknown/missing values at `rank` are
first resolved via [`mp_tax_fix()`](mp_tax_fix.md).
