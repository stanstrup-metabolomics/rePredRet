# Locally Adaptive Prediction Intervals

Prediction intervals whose width follows the local scatter around the
monotone GAM, calibrated by cross-validated conformal prediction.

Computes prediction intervals for a monotonically constrained GAM whose
width varies along the source RT with the local scatter of the
calibration compounds.

## Usage

``` r
gam_adaptive_pi(rt_matrix, newdata, alpha = 0.05, n_folds = 5)
```

## Arguments

- rt_matrix:

  A 2-column matrix `[rt_source, rt_target]`.

- newdata:

  Numeric vector of x-values for predictions.

- alpha:

  Significance level (default 0.05 for 95% PI).

- n_folds:

  Number of cross-validation folds (default 5).

## Value

A matrix with 3 columns: `[predicted, lower, upper]`.

## Details

Algorithm:

1.  Fit the monotone GAM
    ([`fit_gam_mono()`](https://stanstrup-metabolomics.github.io/rePredRet/reference/fit_gam_mono.md)).

2.  Cross-validated residuals: the compounds are split into `n_folds`
    folds, interleaved along the source RT. Each fold is predicted by
    the monotone spline refitted on the other folds with the smoothing
    parameter and weights of the full fit. If a fold cannot be refitted,
    its leave-one-out residuals of the unconstrained GAM are used
    instead.

3.  Local scale s(x): a Gamma GAM (log link, REML, at most 5 basis
    functions) of the absolute cross-validated residuals, weighted by
    the sigmoid outlier weights so that single outliers do not widen the
    band locally. Between the weighted compounds s(x) is interpolated
    linearly (see
    [`local_residual_scale()`](https://stanstrup-metabolomics.github.io/rePredRet/reference/local_residual_scale.md)).
    s(x) is floored at 30% of the global weighted scale and capped at
    the largest residual.

4.  Calibration: the multiplier q is the weighted split-conformal
    quantile of the normalized residuals \|e\| / s(x), where each
    compound counts by its sigmoid outlier weight
    ([`weighted_conformal_quantile()`](https://stanstrup-metabolomics.github.io/rePredRet/reference/weighted_conformal_quantile.md)).
    Outliers are taken to be wrong annotations, so the interval is meant
    to cover `1 - alpha` of the well-annotated compounds, not the
    outliers. With fewer than about `1 / alpha` compounds q is the
    largest score of the compounds that are not down-weighted.

5.  Interval: prediction +/- q \* s(x). Each bound is made
    non-decreasing by isotonic regression, which keeps lower \<=
    prediction \<= upper. The lower bound is truncated at 0.

With few compounds the interval covers less than its nominal level: in a
cross-validation benchmark on RepoRT pairs with 10-15 shared compounds a
95% interval covered about 85% (about 95% from 30 compounds up).

Because the residuals are out-of-fold, q already includes the
uncertainty of the fitted curve, so no separate standard error term is
added.

The idea for this method, and the problems with the earlier interval
methods that it fixes, came from Barbora Hamidová's evaluation of
rePredRet in her MSc thesis (2026).

## References

Lei, J., G'Sell, M., Rinaldo, A., Tibshirani, R. J. & Wasserman, L.
(2018). Distribution-free predictive inference for regression. *Journal
of the American Statistical Association*, 113(523), 1094-1111.
