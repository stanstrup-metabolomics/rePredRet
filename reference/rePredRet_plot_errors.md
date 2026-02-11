# Plot Prediction Error Distribution

Creates a histogram or density plot of prediction errors for a system.

## Usage

``` r
rePredRet_plot_errors(
  system_id,
  type = c("histogram", "density"),
  data_path = NULL
)
```

## Arguments

- system_id:

  Target system identifier.

- type:

  Plot type: "histogram" or "density".

- data_path:

  Path to rePredRet-data (or NULL to use cache).

## Value

A ggplot object.
