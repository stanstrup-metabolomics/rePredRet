# Plot Compound Coverage Heatmap

Creates a heatmap showing which compounds appear in which repository
systems, colored by retention time value. This helps identify compounds
with broad coverage across systems and systems with extensive compound
libraries.

## Usage

``` r
plot_compound_coverage(models, top_n = 30)
```

## Arguments

- models:

  A `UserRTModels` object from
  [`build_user_models`](https://stanstrup.github.io/rePredRet/reference/build_user_models.md).

- top_n:

  Number of most-covered compounds to display (default 30). Compounds
  are ranked by the number of systems in which they appear.

## Value

A ggplot2 object.

## Details

The heatmap displays compounds on the y-axis (labeled by truncated InChI
strings) and repository systems on the x-axis. Each tile is colored by
the retention time value in that system, using a viridis color scale.
Missing values (compound not measured in that system) are left blank.

Only systems that have at least one of the top `top_n` compounds are
included on the x-axis. Compound labels are truncated InChI strings
showing the molecular formula layer for readability.

## See also

[`build_user_models`](https://stanstrup.github.io/rePredRet/reference/build_user_models.md),
[`UserRTModels-class`](https://stanstrup.github.io/rePredRet/reference/UserRTModels-class.md)

## Examples

``` r
if (FALSE) { # \dontrun{
models <- build_user_models(my_inchis, my_rts, "RP")

# Top 30 most-covered compounds
plot_compound_coverage(models)

# Top 50 compounds
plot_compound_coverage(models, top_n = 50)
} # }
```
