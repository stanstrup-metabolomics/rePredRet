# Predict from a Monotone GAM Fit

Predict from a Monotone GAM Fit

## Usage

``` r
predict_gam_mono(fit, newdata)
```

## Arguments

- fit:

  Output of
  [`fit_gam_mono()`](https://stanstrup-metabolomics.github.io/rePredRet/reference/fit_gam_mono.md).

- newdata:

  Numeric vector of source RTs.

## Value

Predicted target RTs, truncated at 0.
