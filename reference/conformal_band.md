# Conformal Band with Local Scale

Turns out-of-sample residuals into a prediction interval: a local scale
model of the residuals
([`local_residual_scale()`](https://stanstrup-metabolomics.github.io/rePredRet/reference/local_residual_scale.md)),
a split-conformal multiplier and monotone bounds. Shared by
[`gam_adaptive_pi()`](https://stanstrup-metabolomics.github.io/rePredRet/reference/gam_adaptive_pi.md)
and the bootstrap method
([`boot2ci_PI()`](https://stanstrup-metabolomics.github.io/rePredRet/reference/boot2ci_PI.md)).

## Usage

``` r
conformal_band(fit, e, mu, newdata, alpha)
```

## Arguments

- fit:

  Output of
  [`fit_gam_mono()`](https://stanstrup-metabolomics.github.io/rePredRet/reference/fit_gam_mono.md)
  (gives x, weights and k).

- e:

  Out-of-sample residuals at `fit$x`.

- mu:

  Predictions at `newdata`.

- newdata:

  Numeric vector of source RTs.

- alpha:

  Significance level.

## Value

A matrix with 3 columns: `[predicted, lower, upper]`.
