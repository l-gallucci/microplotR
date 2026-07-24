# MAG completeness / contamination distribution

Faceted histogram of Completeness and Contamination across all MAGs in
`mag_table` – the standard supplementary figure showing the overall
quality spread of a recovered MAG set. Same MIMAG threshold lines as
[`mp_mag_quality_plot()`](mp_mag_quality_plot.md), drawn per facet.
Accepts CheckM2 or CheckM (v1) column names (normalized internally).

## Usage

``` r
mp_mag_quality_distribution(mag_table, bins = 20, show_mimag_thresholds = TRUE)
```

## Arguments

- mag_table:

  A data frame with a bin-identifier column, `Completeness`,
  `Contamination`.

- bins:

  Number of histogram bins per facet.

- show_mimag_thresholds:

  Draw the MIMAG quality-tier reference lines.

## Value

A ggplot2 object.
