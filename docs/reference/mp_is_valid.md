# Whether a validation report has no errors

Whether a validation report has no errors

## Usage

``` r
mp_is_valid(report)
```

## Arguments

- report:

  An `mp_validation_report` from [`mp_validate()`](mp_validate.md).

## Value

`TRUE` if there are no `error`-level findings (warnings are fine).
