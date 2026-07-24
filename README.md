# microplotr

Literature-grounded, publication-ready microbial ecology plots from plain
tidy input (feature table + taxonomy + metadata), with an upload-time
validator that reports exactly what is wrong with your data, and a Shiny
app for interactive use. Built on ggplot2. Python counterpart:
[`microplotpy`](https://github.com/l-gallucci/microplotPy) — same input
format, same validator rules, same plot catalog.

See [`data-format.md`](data-format.md) for the required input files. Full
documentation (one article per plot: input spec, parameters, literature
grounding) is a [pkgdown](https://pkgdown.r-lib.org/) site built from
`vignettes/` into `docs/`, served via GitHub Pages once enabled for this
repo (or open `docs/index.html` locally).

## Plot catalog

- **Taxonomy**: `mp_taxa_barplot()` (nested-legend stacked barplot),
  `mp_taxa_heatmap()` (CLR/log10, clustered), `mp_taxa_bubbleplot()`,
  `mp_asv_gradient_plot()` (abundance over a continuous gradient),
  `mp_alpha_diversity_plot()`, `mp_beta_diversity_plot()` (PCoA/NMDS +
  PERMANOVA), `mp_taxa_treemap()`.
- **Functional profiles** (eggNOG-mapper/KofamScan-shaped):
  `mp_function_barplot()`, `mp_function_heatmap()`, `mp_function_treemap()`
  — thin wrappers around the taxa engine.
- **Metagenomics**: `mp_mag_quality_plot()`/`mp_mag_quality_distribution()`
  (CheckM/CheckM2), `mp_assembly_nx_plot()`/`mp_assembly_summary_barplot()`
  (QUAST-shaped).

Every plot resolves missing/unclassified annotation via `mp_tax_fix()`,
supports abundance/prevalence filtering via `mp_filter_taxa()`, and returns
a plain `ggplot2` object you can extend further (`p + ggplot2::theme(...)`
etc.) — nothing here is a static image.

## Install

```r
install.packages("remotes")
remotes::install_github("l-gallucci/microplotR")
```

Not yet on CRAN — once published there, this becomes a plain
`install.packages("microplotr")`.

## Shiny app

```r
library(microplotr)
shiny::runApp(system.file("shiny", package = "microplotr"))
```

Four tabs (Taxonomy, Functional profile, MAG quality, Assembly QC), each:
upload the required file(s) → validator shows errors/warnings inline →
pick a plot type and parameters → view/download the plot as PNG or SVG at
a width/height you choose.

## Quick usage

```r
library(microplotr)

example_dir <- system.file("extdata", "example_valid", package = "microplotr")
data <- mp_read_data(
  file.path(example_dir, "feature_table.tsv"),
  file.path(example_dir, "taxonomy.tsv"),
  file.path(example_dir, "metadata.tsv")
)
report <- mp_validate(data, gradient_column = "Depth_m", group_column = "Group")
mp_is_valid(report)

mp_taxa_barplot(data, rank = "Genus", group_rank = "Phylum", top_n = 10)
```

The package ships several other bundled example datasets under
`system.file("extdata", package = "microplotr")` — including
`example_broken_*` sets that each trip a different validator check (ID
mismatches, negative counts, missing taxonomy columns, duplicate feature
IDs) — useful for seeing exactly what `mp_validate()` catches.

## Contributing / development

```r
git clone https://github.com/l-gallucci/microplotR.git && cd microplotR
```
```r
devtools::install_deps()
devtools::load_all()
devtools::test()
```

To rebuild the documentation site after a change:

```r
devtools::load_all()
pkgdown::build_site()
```
