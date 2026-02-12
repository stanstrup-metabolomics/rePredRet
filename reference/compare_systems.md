# Compare User's System to Repository Systems

Ranks repository systems by how well their models fit the user's data,
helping identify the most similar systems.

## Usage

``` r
compare_systems(
  models,
  metric = c("r_squared", "rmse", "median_ci_width", "n_calibration"),
  top_n = NULL
)
```

## Arguments

- models:

  A UserRTModels object from
  [`build_user_models`](https://stanstrup.github.io/rePredRet/reference/build_user_models.md).

- metric:

  Metric to rank by: "r_squared" (default), "rmse", "median_ci_width",
  or "n_calibration".

- top_n:

  Number of top systems to return (default: all).

## Value

A data.frame ranked by the specified metric, with columns for all
diagnostic metrics.

## Examples

``` r
if (FALSE) { # \dontrun{
models <- build_user_models(my_inchis, my_rts, "RP")

# Find most similar systems (best R-squared)
compare_systems(models)

# Rank by lowest RMSE
compare_systems(models, metric = "rmse")

# Top 5 systems
compare_systems(models, top_n = 5)
} # }
```
