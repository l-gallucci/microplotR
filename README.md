# microplotr

Literature-grounded, publication-ready microbial ecology plots from plain
tidy input (feature table + taxonomy + metadata), with an upload-time
validator that reports exactly what is wrong with your data, and a Shiny
app for interactive use. Built on ggplot2. Python counterpart:
[`microplotpy`](../microplotpy) — same input format, same validator rules,
same plot catalog.

See [`data-format.md`](data-format.md) for the required input files. Full
documentation (one article per plot: input spec, parameters, literature
grounding) is a [pkgdown](https://pkgdown.r-lib.org/) site built from
`vignettes/` into `docs/` — open `docs/index.html` locally, or serve it via
GitHub Pages once this repo is published (`devtools::load_all();
pkgdown::build_site()` to rebuild after any change).

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

## Install (dev)

```r
devtools::install_deps()
devtools::load_all()
devtools::test()
```

## Shiny app

```r
devtools::load_all()
shiny::runApp("inst/shiny")
```

Four tabs (Taxonomy, Functional profile, MAG quality, Assembly QC), each:
upload the required file(s) → validator shows errors/warnings inline →
pick a plot type and parameters → view/download the plot as PNG or SVG at
a width/height you choose.

## Quick usage (library)

```r
library(microplotr)

data <- mp_read_data(
  "inst/extdata/example_valid/feature_table.tsv",
  "inst/extdata/example_valid/taxonomy.tsv",
  "inst/extdata/example_valid/metadata.tsv"
)
report <- mp_validate(data, gradient_column = "Depth_m", group_column = "Group")
mp_is_valid(report)

mp_taxa_barplot(data, rank = "Genus", group_rank = "Phylum", top_n = 10)
```

Try validation against `inst/extdata/example_broken_*` to see what it
catches (ID mismatches, negative counts, missing taxonomy columns,
duplicate feature IDs).
