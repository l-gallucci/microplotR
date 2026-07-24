# Package index

## Import & validate

Read tidy input and check it before plotting.

- [`mp_read_data()`](mp_read_data.md) : Read tidy microbiome input files
- [`mp_validate()`](mp_validate.md) : Validate tidy microbiome input
  (feature table + taxonomy + metadata)
- [`mp_validate_mag()`](mp_validate_mag.md) : Validate a MAG
  (metagenome-assembled genome) quality table
- [`mp_validate_contig_lengths()`](mp_validate_contig_lengths.md) :
  Validate a per-contig length table
- [`mp_validate_assembly_summary()`](mp_validate_assembly_summary.md) :
  Validate an assembly summary table
- [`mp_errors()`](mp_errors.md) : Errors from a validation report
- [`mp_warnings()`](mp_warnings.md) : Warnings from a validation report
- [`mp_is_valid()`](mp_is_valid.md) : Whether a validation report has no
  errors

## Taxonomy fallback & filtering

Shared machinery used by every taxa/function plot.

- [`mp_tax_fix()`](mp_tax_fix.md) : Fix unknown/missing taxonomy labels
  by falling back to the last known rank
- [`mp_filter_taxa()`](mp_filter_taxa.md) : Filter taxa by minimum
  abundance and/or prevalence
- [`mp_nested_order_palette()`](mp_nested_order_palette.md) : Build
  ordering, nested color palette, and nested legend labels for a
  taxon-colored plot
- [`mp_theme_pub()`](mp_theme_pub.md) : Minimal publication-figure
  ggplot2 theme

## Taxonomy plots

- [`mp_taxa_barplot()`](mp_taxa_barplot.md) : Stacked relative-abundance
  taxa barplot, optionally with a nested legend
- [`mp_taxa_heatmap()`](mp_taxa_heatmap.md) : Clustered taxa heatmap
- [`mp_taxa_bubbleplot()`](mp_taxa_bubbleplot.md) : Taxa bubble (dot)
  plot
- [`mp_asv_gradient_plot()`](mp_asv_gradient_plot.md) : ASV/taxon
  abundance over a continuous gradient
- [`mp_alpha_diversity_plot()`](mp_alpha_diversity_plot.md) : Alpha
  diversity boxplot
- [`mp_beta_diversity_plot()`](mp_beta_diversity_plot.md) : Beta
  diversity ordination plot
- [`mp_taxa_treemap()`](mp_taxa_treemap.md) : Taxa (or functional)
  treemap

## Functional profiles

- [`mp_function_barplot()`](mp_function_barplot.md) : Functional profile
  barplot
- [`mp_function_heatmap()`](mp_function_heatmap.md) : Functional profile
  heatmap
- [`mp_function_treemap()`](mp_function_treemap.md) : Functional treemap

## Metagenomics

- [`mp_mag_quality_plot()`](mp_mag_quality_plot.md) : MAG
  (metagenome-assembled genome) quality plot
- [`mp_mag_quality_distribution()`](mp_mag_quality_distribution.md) :
  MAG completeness / contamination distribution
- [`mp_assembly_nx_plot()`](mp_assembly_nx_plot.md) : Assembly Nx curve
- [`mp_assembly_summary_barplot()`](mp_assembly_summary_barplot.md) :
  Assembly summary barplot
