# Build ordering, nested color palette, and nested legend labels for a taxon-colored plot

Orders taxa grouped by an upper taxonomy level (each group's taxa kept
contiguous, groups and within-group taxa both ordered by descending
abundance), assigns each group a base hue and shades its member taxa
from that hue (so related taxa read as a family of colors), and builds
legend labels where the first taxon in each group carries a bold
`**group**` markdown header (rendered via
[`ggtext::element_markdown()`](https://wilkelab.org/ggtext/reference/element_markdown.html)).
`other_label` (e.g. "Other") is always ordered last and colored grey,
with a plain (non-bold) label.

## Usage

``` r
mp_nested_order_palette(
  taxa,
  groups,
  totals,
  other_label = "Other",
  other_color = "grey80"
)
```

## Arguments

- taxa:

  Character vector, one entry per unique taxon (e.g. genus).

- groups:

  Character vector, same length as `taxa`, the upper taxonomy level
  (e.g. phylum) each taxon belongs to.

- totals:

  Numeric vector, same length as `taxa`, overall abundance used to rank
  groups and within-group taxa.

- other_label:

  The sentinel value used for the aggregated "everything else" bucket,
  if present in `taxa`.

- other_color:

  Color for `other_label`.

## Value

A list with `order` (character vector, plotting/legend order), `palette`
(named character vector, taxon -\> hex color), `labels` (named character
vector, taxon -\> markdown legend label).
