# Cross-Validated Residuals of a Monotone GAM

Cross-Validated Residuals of a Monotone GAM

## Usage

``` r
cv_residuals_gam_mono(fit, n_folds = 5)
```

## Arguments

- fit:

  Output of
  [`fit_gam_mono()`](https://stanstrup-metabolomics.github.io/rePredRet/reference/fit_gam_mono.md).

- n_folds:

  Number of folds, interleaved along x.

## Value

Out-of-fold residuals, one per compound.
