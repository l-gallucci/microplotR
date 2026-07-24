# Validate an assembly summary table

Checks a QUAST `report.tsv`-shaped table (one row per assembly). Column
names from QUAST's own report (e.g. `"# contigs"`, `"Total length"`) are
accepted and normalized internally.

## Usage

``` r
mp_validate_assembly_summary(assembly_summary)
```

## Arguments

- assembly_summary:

  A data frame with an assembly-identifier column (`Assembly_ID` or
  QUAST's `"Assembly"`), and `N50`, `Total_length` (or QUAST's
  `"Total length"`).

## Value

An `mp_validation_report` (see [`mp_validate()`](mp_validate.md)).
