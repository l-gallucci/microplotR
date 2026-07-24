# Read tidy microbiome input files

Reads the three tidy input files as character-typed tibbles (no
validation, no type coercion) so [`mp_validate()`](mp_validate.md) can
distinguish "not numeric" from "numeric" without prior silent coercion.

## Usage

``` r
mp_read_data(feature_table_path, taxonomy_path, metadata_path)
```

## Arguments

- feature_table_path:

  Path to feature_table.tsv/.csv.

- taxonomy_path:

  Path to taxonomy.tsv/.csv.

- metadata_path:

  Path to metadata.tsv/.csv.

## Value

A list with `feature_table`, `taxonomy`, `metadata` tibbles.
