# Generate All Predictions (Optimized Batch Processing)

Generates predictions for all compounds in all systems where predictions
can be made. This version is optimized for batch processing.

## Usage

``` r
generate_all_predictions(
  models,
  datasets,
  studies,
  report_version = NA,
  ci_width_limit = .prediction_defaults$ci_width_limit,
  ci_width_limit_rel = .prediction_defaults$ci_width_limit_rel,
  density_bw_mult = .prediction_defaults$density_bw_mult,
  density_limit = .prediction_defaults$density_limit
)
```

## Arguments

- models:

  Named list of models from
  [`build_all_models()`](https://stanstrup.github.io/rePredRet/reference/build_all_models.md).

- datasets:

  Named list of datasets from
  [`load_report_data()`](https://stanstrup.github.io/rePredRet/reference/load_report_data.md).

- studies:

  Studies tibble from
  [`load_report_data()`](https://stanstrup.github.io/rePredRet/reference/load_report_data.md).

- report_version:

  RepoRT version string for provenance tracking.

- ci_width_limit:

  Maximum CI width in minutes (default 2.0).

- ci_width_limit_rel:

  Maximum relative CI width (default 0.2 = 20%).

- density_bw_mult:

  Bandwidth multiplier for density check (default 0.03).

- density_limit:

  Minimum density threshold (default 0.01).

## Value

A list with predictions organized by target system.
