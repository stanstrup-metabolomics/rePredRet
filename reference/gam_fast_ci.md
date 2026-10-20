# Fast Prediction Intervals with a Global Residual Scale

Fast prediction intervals from one monotone GAM fit and a single
residual variance.

Computes prediction intervals of nearly constant width from a single
monotone GAM fit. About four times faster than
[`gam_adaptive_pi()`](https://stanstrup-metabolomics.github.io/rePredRet/reference/gam_adaptive_pi.md),
whose width follows the local scatter.

## Usage

``` r
gam_fast_ci(rt_matrix, newdata, alpha = 0.05)
```

## Arguments

- rt_matrix:

  A 2-column matrix `[rt_source, rt_target]`.

- newdata:

  Numeric vector of x-values for predictions.

- alpha:

  Significance level (default 0.05 for 95% PI).

## Value

A matrix with 3 columns: `[predicted, lower, upper]`.

## Details

Algorithm:

1.  Fit the monotone GAM
    ([`fit_gam_mono()`](https://stanstrup-metabolomics.github.io/rePredRet/reference/fit_gam_mono.md)).

2.  SE: pointwise standard errors of the unconstrained GAM.

3.  sigma^2: mean squared leave-one-out residual. The residuals of the
    constrained fit are divided by `1 - h`, with `h` the leverages of
    the unconstrained GAM. In-sample residuals underestimate the error
    for new compounds, most of all with few compounds.

4.  Interval: prediction +/- z \* sqrt(SE^2 + sigma^2), with bounds made
    non-decreasing by isotonic regression and the lower bound truncated
    at 0
    ([`monotone_interval()`](https://stanstrup-metabolomics.github.io/rePredRet/reference/monotone_interval.md)).
