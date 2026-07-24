# Validate tidy microbiome input (feature table + taxonomy + metadata)

Mirrors the checks documented in data-format.md and implemented on the
Python side by `validate()` in microplotpy – keep both in sync.

## Usage

``` r
mp_validate(
  data,
  gradient_column = NULL,
  group_column = NULL,
  required_ranks = c("Phylum", "Genus")
)
```

## Arguments

- data:

  List as returned by [`mp_read_data()`](mp_read_data.md), with
  `feature_table`, `taxonomy`, `metadata` tibbles.

- gradient_column:

  Optional metadata column to check is numeric (required by the
  ASV-level gradient plot).

- group_column:

  Optional metadata column to check has \>= 2 levels (required by
  alpha/beta diversity group comparisons).

- required_ranks:

  Column names that must be present in `data$taxonomy`. Default
  `c("Phylum", "Genus")` for the standard 16S taxonomy shape. Pass
  `character(0)` when `data$taxonomy` is actually a non-taxonomic
  hierarchy (e.g. a functional annotation table from eggNOG-mapper/
  KofamScan, see [`mp_function_barplot()`](mp_function_barplot.md))
  whose column names aren't ranks.

## Value

An `mp_validation_report`: a list with `findings`, a data frame with
columns `level`, `file`, `field`, `message`.
