# Assembly summary barplot

A single summary statistic (default `N50`) per assembly, from a
QUAST-`report.tsv`-shaped table.

## Usage

``` r
mp_assembly_summary_barplot(assembly_summary, stat_col = "N50")
```

## Arguments

- assembly_summary:

  A data frame with an `Assembly_ID` column and the column named by
  `stat_col` (as validated by
  [`mp_validate_assembly_summary()`](mp_validate_assembly_summary.md)).

- stat_col:

  Which column to plot. Default `"N50"`.

## Value

A ggplot2 object.
