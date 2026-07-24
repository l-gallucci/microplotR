# Minimal publication-figure ggplot2 theme

No title/subtitle/caption chrome, no gridlines, no gratuitous background
– just axis titles, tick labels, and the legend. Intended so plot
functions never bake in a title (callers add one via
`patchwork`/`cowplot`/knitr chunk captions if they want one at all).

## Usage

``` r
mp_theme_pub(base_size = 11)
```

## Arguments

- base_size:

  Base font size.

## Value

A ggplot2 theme.
