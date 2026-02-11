# Plot System Network

Creates a network graph showing connections between the user's
chromatographic system and repository systems. Nodes represent systems,
edges represent calibration models, and edge opacity/width reflect model
quality (R-squared).

## Usage

``` r
plot_system_network(models, min_compounds = 10, interactive = FALSE)
```

## Arguments

- models:

  A `UserRTModels` object from
  [`build_user_models`](https://stanstrup.github.io/rePredRet/reference/build_user_models.md).

- min_compounds:

  Minimum number of calibration compounds for a model to be included in
  the network (default 10).

- interactive:

  If TRUE, returns an interactive plotly plot. If FALSE (default),
  returns a static ggplot2 plot.

## Value

A ggplot2 object (if interactive = FALSE) or a plotly object (if
interactive = TRUE).

## Details

The network uses a circular layout with the user's system ("USER")
placed at the center. Repository systems are arranged around the
perimeter. Edge width and opacity are proportional to the model's
R-squared value, making it easy to identify the strongest calibration
relationships at a glance.

Systems with fewer than `min_compounds` calibration compounds are
excluded from the visualization.

## See also

[`build_user_models`](https://stanstrup.github.io/rePredRet/reference/build_user_models.md),
[`compare_systems`](https://stanstrup.github.io/rePredRet/reference/compare_systems.md)

## Examples

``` r
if (FALSE) { # \dontrun{
models <- build_user_models(my_inchis, my_rts, "RP")

# Static network plot
plot_system_network(models)

# Only show models with >= 20 calibration compounds
plot_system_network(models, min_compounds = 20)

# Interactive version (requires plotly)
plot_system_network(models, interactive = TRUE)
} # }
```
