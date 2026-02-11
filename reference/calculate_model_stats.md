# Calculate Model Performance Statistics

Computes comprehensive performance metrics for all models, including
leave-one-out cross-validation errors when calibration data is
available.

## Usage

``` r
calculate_model_stats(model_results, datasets = NULL)
```

## Arguments

- model_results:

  Output from
  [`build_all_models()`](https://stanstrup.github.io/rePredRet/reference/build_all_models.md).

- datasets:

  Named list of datasets from
  [`load_report_data()`](https://stanstrup.github.io/rePredRet/reference/load_report_data.md).

## Value

A list with:

- model_stats:

  Tibble with per-model statistics

- overall_stats:

  Named list with aggregate statistics

- by_method:

  Tibble with stats grouped by method type
