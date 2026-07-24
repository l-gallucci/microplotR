# Validate a MAG (metagenome-assembled genome) quality table

Checks a single flat table shaped like CheckM2's `quality_report.tsv` or
CheckM (v1)'s `qa`/`bin_stats` table (one row per bin/MAG) – not the
3-file taxonomy/metadata schema used by
[`mp_validate()`](mp_validate.md). Column names from either tool are
accepted (e.g. CheckM1's `"Bin Id"` or CheckM2's `"Name"`) and
normalized internally. Mirrors the same `{level, file, field, message}`
finding shape.

## Usage

``` r
mp_validate_mag(mag_table)
```

## Arguments

- mag_table:

  A data frame with at least a bin-identifier column (`Name` or
  CheckM1's `"Bin Id"`), `Completeness`, `Contamination`, and a
  genome-size column (`Genome_Size` or CheckM1's `"Genome size (bp)"`).

## Value

An `mp_validation_report` (see [`mp_validate()`](mp_validate.md)).
