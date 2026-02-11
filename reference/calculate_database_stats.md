# Calculate Database Statistics

Computes statistics about the compound database coverage, similar to
Figure 3A in the paper.

## Usage

``` r
calculate_database_stats(report_data)
```

## Arguments

- report_data:

  Output from
  [`load_report_data()`](https://stanstrup.github.io/rePredRet/reference/load_report_data.md).

## Value

A tibble with per-system statistics.
