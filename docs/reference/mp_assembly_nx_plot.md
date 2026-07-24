# Assembly Nx curve

Cumulative contig-length distribution per assembly – the standard
QUAST-style plot: x = Nx percentile (0-100, cumulative % of total
assembly length), y = contig length (log10 scale), one line per
assembly.

## Usage

``` r
mp_assembly_nx_plot(contig_lengths)
```

## Arguments

- contig_lengths:

  A data frame with `Assembly_ID`, `Contig_ID`, `Length` columns (as
  validated by
  [`mp_validate_contig_lengths()`](mp_validate_contig_lengths.md)).

## Value

A ggplot2 object.
