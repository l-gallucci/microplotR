# microplotr

Literature-grounded, publication-ready microbial ecology plots from
plain tidy input (feature table + taxonomy + metadata), with an
upload-time validator that reports exactly what is wrong with your data,
and a Shiny app for interactive use. Built on ggplot2. Python
counterpart: [`microplotpy`](../microplotpy) — same input format, same
validator rules, same plot catalog.

See [`data-format.md`](data-format.md) for the required input files, and
[`docs/`](docs) for one page per plot (input spec, parameters,
literature grounding).

## Plot catalog

- **Taxonomy**: [`mp_taxa_barplot()`](reference/mp_taxa_barplot.md)
  (nested-legend stacked barplot),
  [`mp_taxa_heatmap()`](reference/mp_taxa_heatmap.md) (CLR/log10,
  clustered), [`mp_taxa_bubbleplot()`](reference/mp_taxa_bubbleplot.md),
  [`mp_asv_gradient_plot()`](reference/mp_asv_gradient_plot.md)
  (abundance over a continuous gradient),
  [`mp_alpha_diversity_plot()`](reference/mp_alpha_diversity_plot.md),
  [`mp_beta_diversity_plot()`](reference/mp_beta_diversity_plot.md)
  (PCoA/NMDS + PERMANOVA),
  [`mp_taxa_treemap()`](reference/mp_taxa_treemap.md).
- **Functional profiles** (eggNOG-mapper/KofamScan-shaped):
  [`mp_function_barplot()`](reference/mp_function_barplot.md),
  [`mp_function_heatmap()`](reference/mp_function_heatmap.md),
  [`mp_function_treemap()`](reference/mp_function_treemap.md) — thin
  wrappers around the taxa engine.
- **Metagenomics**:
  [`mp_mag_quality_plot()`](reference/mp_mag_quality_plot.md)/[`mp_mag_quality_distribution()`](reference/mp_mag_quality_distribution.md)
  (CheckM/CheckM2),
  [`mp_assembly_nx_plot()`](reference/mp_assembly_nx_plot.md)/[`mp_assembly_summary_barplot()`](reference/mp_assembly_summary_barplot.md)
  (QUAST-shaped).

Every plot resolves missing/unclassified annotation via
[`mp_tax_fix()`](reference/mp_tax_fix.md), supports abundance/prevalence
filtering via [`mp_filter_taxa()`](reference/mp_filter_taxa.md), and
returns a plain `ggplot2` object you can extend further
(`p + ggplot2::theme(...)` etc.) — nothing here is a static image.

## Install (dev)

``` r
devtools::install_deps()
devtools::load_all()
devtools::test()
```

## Shiny app

``` r
devtools::load_all()
shiny::runApp("inst/shiny")
```

Four tabs (Taxonomy, Functional profile, MAG quality, Assembly QC),
each: upload the required file(s) → validator shows errors/warnings
inline → pick a plot type and parameters → view/download the plot as
PNG.

## Quick usage (library)

``` r
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
