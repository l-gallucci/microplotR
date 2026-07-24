# Alpha diversity boxplot

Observed richness, Shannon, and Simpson diversity per sample (via
[`vegan::specnumber()`](https://vegandevs.github.io/vegan/reference/diversity.html)/[`vegan::diversity()`](https://vegandevs.github.io/vegan/reference/diversity.html)),
boxplot + jitter by a grouping metadata variable, faceted one panel per
metric, with a Wilcoxon (2 groups) or Kruskal-Wallis (\>2 groups) test
annotated via
[`ggpubr::stat_compare_means()`](https://rpkgs.datanovia.com/ggpubr/reference/stat_compare_means.html)
– the standard alpha-diversity comparison figure in amplicon/marker-gene
studies.

## Usage

``` r
mp_alpha_diversity_plot(
  data,
  group_var,
  metrics = c("Observed", "Shannon", "Simpson"),
  test = "auto"
)
```

## Arguments

- data:

  List with `feature_table`, `taxonomy`, `metadata`.

- group_var:

  Metadata column (categorical, \>=2 levels) to compare.

- metrics:

  Diversity metrics to compute/plot, any of `"Observed"`, `"Shannon"`,
  `"Simpson"`. Default all three.

- test:

  `"auto"` (default: Wilcoxon for 2 groups, Kruskal-Wallis for more than
  2), `"wilcoxon"`, `"kruskal"`, or `NULL` to disable the test
  annotation.

## Value

A ggplot2 object, faceted by metric.
