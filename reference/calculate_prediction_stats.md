# Calculate Prediction Statistics

Computes statistics about predictions, if predictions have been
generated.

## Usage

``` r
calculate_prediction_stats(predictions, model_results = NULL)
```

## Arguments

- predictions:

  Output from
  [`generate_all_predictions()`](https://stanstrup-metabolomics.github.io/rePredRet/reference/generate_all_predictions.md).

- model_results:

  Output from
  [`build_all_models()`](https://stanstrup-metabolomics.github.io/rePredRet/reference/build_all_models.md).

## Value

A list with prediction statistics.
