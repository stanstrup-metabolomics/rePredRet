# Validate Models with Cross-Validation

Performs leave-one-out cross-validation on the calibration compounds to
assess prediction accuracy.

## Usage

``` r
validate_models(models, method = c("loo"))
```

## Arguments

- models:

  A UserRTModels object from
  [`build_user_models`](https://stanstrup.github.io/rePredRet/reference/build_user_models.md).

- method:

  Validation method: "loo" (leave-one-out, default).

## Value

A data.frame with cross-validation results for each model:

- source_system:

  System ID

- n_compounds:

  Number of calibration compounds

- cv_rmse:

  Cross-validated RMSE

- cv_mae:

  Cross-validated mean absolute error

- cv_median_ae:

  Cross-validated median absolute error

- cv_coverage:

  Proportion of actual values within CI

## Details

Leave-one-out cross-validation removes each calibration compound in
turn, refits the model, and predicts the held-out compound. This
provides an unbiased estimate of prediction error.

## Examples

``` r
if (FALSE) { # \dontrun{
models <- build_user_models(my_inchis, my_rts, "RP")
cv_results <- validate_models(models)
print(cv_results)
} # }
```
