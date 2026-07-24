# Validate a per-contig length table

Checks a tidy `Assembly_ID`, `Contig_ID`, `Length` table (one row per
contig) – the input [`mp_assembly_nx_plot()`](mp_assembly_nx_plot.md)
needs. Not the QUAST `report.tsv` summary shape (see
[`mp_validate_assembly_summary()`](mp_validate_assembly_summary.md)).

## Usage

``` r
mp_validate_contig_lengths(contig_lengths)
```

## Arguments

- contig_lengths:

  A data frame with `Assembly_ID`, `Contig_ID`, `Length` columns.

## Value

An `mp_validation_report` (see [`mp_validate()`](mp_validate.md)).
