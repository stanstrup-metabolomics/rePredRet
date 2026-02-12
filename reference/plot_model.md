# Plot a Pre-built Model from the Data Repository

Downloads and plots a calibration model showing the relationship between
retention times in two chromatographic systems.

## Usage

``` r
plot_model(
  from_id,
  to_id,
  interactive = FALSE,
  show_ci = TRUE,
  show_calibration = TRUE,
  data_path = NULL
)
```

## Arguments

- from_id:

  Source system ID (e.g., "0001").

- to_id:

  Target system ID (e.g., "0002").

- interactive:

  If TRUE, returns a plotly interactive plot. If FALSE (default),
  returns a ggplot2 static plot.

- show_ci:

  If TRUE (default), shows confidence interval bands.

- show_calibration:

  If TRUE (default), shows calibration points.

- data_path:

  Optional path to local rePredRet-data directory. If NULL,
  downloads/uses cached data.

## Value

A ggplot2 object (if interactive = FALSE) or a plotly object (if
interactive = TRUE).

## Details

This function fetches model data from the rePredRet-models repository
and creates a publication-ready plot showing:

- The calibration curve (monotonic GAM fit)

- Confidence interval bands

- Calibration points used to build the model

## See also

[`rePredRet_model`](https://stanstrup.github.io/rePredRet/reference/rePredRet_model.md)
for getting model data,
[`rePredRet_systems`](https://stanstrup.github.io/rePredRet/reference/rePredRet_systems.md)
for listing available systems

## Examples

``` r
if (FALSE) { # \dontrun{
# Static ggplot
plot_model("0001", "0002")

# Interactive plotly
plot_model("0001", "0002", interactive = TRUE)

# Without CI bands
plot_model("0001", "0002", show_ci = FALSE)
} # }
```
