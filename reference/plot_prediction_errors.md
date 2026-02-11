# Plot Prediction Errors

Visualizes the distribution of confidence interval widths (and
optionally model errors) from a set of RT predictions.

## Usage

``` r
plot_prediction_errors(
  predictions,
  type = c("density", "histogram", "boxplot"),
  by_system = FALSE
)
```

## Arguments

- predictions:

  A `UserRTPredictions` object from
  [`predict_rt_user`](https://stanstrup.github.io/rePredRet/reference/predict_rt_user.md).

- type:

  Type of plot to create. One of `"density"` (default), `"histogram"`,
  or `"boxplot"`.

- by_system:

  If TRUE, facets the plot by source system (showing the top 12 systems
  by number of predictions). Default FALSE.

## Value

A ggplot2 object.

## Details

This function creates distribution plots of the confidence interval
widths from predictions. CI width is the primary measure of prediction
uncertainty: narrower intervals indicate more confident predictions.

When `by_system = TRUE`, the plot is faceted by source system to reveal
which repository systems produce the tightest (or widest) predictions.
Only the top 12 systems (by prediction count) are shown to keep the plot
readable.

## See also

[`predict_rt_user`](https://stanstrup.github.io/rePredRet/reference/predict_rt_user.md),
[`UserRTPredictions-class`](https://stanstrup.github.io/rePredRet/reference/UserRTPredictions-class.md)

## Examples

``` r
if (FALSE) { # \dontrun{
predictions <- predict_rt_user(my_models)

# Density plot of CI widths
plot_prediction_errors(predictions)

# Histogram
plot_prediction_errors(predictions, type = "histogram")

# Boxplot faceted by source system
plot_prediction_errors(predictions, type = "boxplot", by_system = TRUE)
} # }
```
