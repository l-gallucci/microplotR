# Fix unknown/missing taxonomy labels by falling back to the last known rank

Genus-level (or any rank) plots otherwise show blank/NA/"uncultured"
bars, which is common in 16S taxonomy tables and misleading in a stacked
barplot. This fills an unknown value at a given rank with the last known
ancestor rank so every feature still gets a meaningful label – the same
idea as `tax_fix()` in the microViz package (Barnett et al. 2021,
Bioinformatics). Cascading unknowns (e.g. both Genus and Species
missing) all fall back to the same true ancestor rather than chaining
onto each other's fabricated labels.

## Usage

``` r
mp_tax_fix(
  taxonomy,
  ranks = NULL,
  unknown_strings = c("", "na", "unknown", "unclassified", "uncultured", "unidentified",
    "metagenome"),
  label_format = "Unclassified_{last}"
)
```

## Arguments

- taxonomy:

  A taxonomy data frame/tibble with `Feature_ID` and rank columns.

- ranks:

  Rank columns, high-to-low order. Default: auto-detected from
  `c("Domain","Kingdom","Phylum","Class","Order","Family","Genus","Species")`
  intersected with `names(taxonomy)`.

- unknown_strings:

  Values (case-insensitive, whitespace-trimmed) treated as "unknown".

- label_format:

  Template for the fallback label; `{last}` is replaced with the last
  known ancestor value.

## Value

`taxonomy` with unknown values replaced.
