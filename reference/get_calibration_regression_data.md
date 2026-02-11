# Aggregate Calibration Data for Regression Plot

Collects all calibration point predictions vs actual values for creating
the aggregate regression plot (Figure 3B in paper).

## Usage

``` r
get_calibration_regression_data(model_results)
```

## Arguments

- model_results:

  Output from
  [`build_all_models()`](https://stanstrup.github.io/rePredRet/reference/build_all_models.md).

## Value

A tibble with predicted and actual RT values for all calibration points
across all models.
