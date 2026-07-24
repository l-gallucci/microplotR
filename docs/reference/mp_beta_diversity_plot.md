# Beta diversity ordination plot

Bray-Curtis or Jaccard dissimilarity (on relative abundance), PCoA
([`stats::cmdscale()`](https://rdrr.io/r/stats/cmdscale.html)) or NMDS
([`vegan::metaMDS()`](https://vegandevs.github.io/vegan/reference/metaMDS.html))
ordination, colored by a grouping metadata variable with 95% confidence
ellipses, PERMANOVA
([`vegan::adonis2()`](https://vegandevs.github.io/vegan/reference/adonis.html))
annotated as a plot subtitle – the standard beta- diversity ordination
figure in amplicon/marker-gene studies. No phylogenetic (UniFrac)
distances – Bray-Curtis/Jaccard only.

## Usage

``` r
mp_beta_diversity_plot(
  data,
  group_var,
  method = c("bray", "jaccard"),
  ordination = c("pcoa", "nmds"),
  show_ellipse = TRUE,
  permanova = TRUE,
  permutations = 999
)
```

## Arguments

- data:

  List with `feature_table`, `taxonomy`, `metadata`.

- group_var:

  Metadata column (categorical) to color/test by.

- method:

  Dissimilarity method: `"bray"` (default) or `"jaccard"`.

- ordination:

  `"pcoa"` (default) or `"nmds"`.

- show_ellipse:

  Draw 95% confidence ellipses per group
  ([`ggplot2::stat_ellipse()`](https://ggplot2.tidyverse.org/reference/stat_ellipse.html)).

- permanova:

  Run and annotate PERMANOVA
  ([`vegan::adonis2()`](https://vegandevs.github.io/vegan/reference/adonis.html)).

- permutations:

  Number of permutations for PERMANOVA. Default 999.

## Value

A ggplot2 object.
