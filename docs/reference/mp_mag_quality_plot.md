# MAG (metagenome-assembled genome) quality plot

One dot per MAG/bin, x = Contamination (%), y = Completeness (%) – the
standard first figure in any MAG-recovery paper. Accepts CheckM2's
`quality_report.tsv` or CheckM (v1)'s `qa`/`bin_stats` table directly
(column-name differences between the two tools are normalized
internally). Draws dashed reference lines and a shaded region at the
MIMAG quality-tier thresholds by default (Bowers et al. 2017, *Nature
Biotechnology* – "Minimum information about a single amplified genome
and a metagenome-assembled genome"): high quality \>=90% completeness &
\<5% contamination; medium quality \>=50% completeness & \<10%
contamination.

## Usage

``` r
mp_mag_quality_plot(
  mag_table,
  size_col = "Genome_Size",
  color_col = NULL,
  show_mimag_thresholds = TRUE
)
```

## Arguments

- mag_table:

  A data frame with a bin-identifier column (`Name` or CheckM1's
  `"Bin Id"`), `Completeness`, `Contamination`, and optionally a
  genome-size column (`Genome_Size` or CheckM1's `"Genome size (bp)"`).

- size_col:

  Column to map to dot size (e.g. genome size). `NULL` disables sizing
  (uniform dots). Default `"Genome_Size"` (silently ignored if not
  present in `mag_table`).

- color_col:

  Optional column to map to dot color (e.g. a taxon or group column).

- show_mimag_thresholds:

  Draw the MIMAG quality-tier reference lines and shaded high-quality
  region.

## Value

A ggplot2 object.
