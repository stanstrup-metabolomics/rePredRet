# Plot Model Fit

Creates a visualization of a model showing the fit curve, confidence
interval, and calibration points.

## Usage

``` r
rePredRet_plot_model(from_id, to_id, data_path = NULL)
```

## Arguments

- from_id:

  Source system identifier.

- to_id:

  Target system identifier.

- data_path:

  Path to rePredRet-data (or NULL to use cache).

## Value

A ggplot object.

## Examples

``` r
if (FALSE) { # \dontrun{
p <- rePredRet_plot_model("0001", "0002")
print(p)
} # }
```
