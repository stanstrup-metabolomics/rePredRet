# Fast Prediction Intervals with Monotonicity Enforcement

Fast calculation of prediction intervals using single GAM fit with
post-hoc monotonicity enforcement.

Computes prediction intervals by fitting the constrained GAM once,
estimating pointwise SE, and enforcing monotonicity with cummax().

## Usage

``` r
gam_fast_pi(rt_matrix, newdata, alpha = 0.05, prediction_interval = TRUE)
```

## Arguments

- rt_matrix:

  A 2-column matrix rt_source, rt_target.

- newdata:

  Numeric vector of x-values for predictions.

- alpha:

  Significance level (default 0.05 for 95% intervals).

- prediction_interval:

  Logical, if TRUE (default) adds observation variance for prediction
  intervals. If FALSE, gives confidence intervals.

## Value

A matrix with 3 columns: predicted, lower, upper.

## Details

Algorithm:

1.  Fit monotonically constrained GAM to full dataset

2.  Estimate pointwise SE from model variance

3.  Build intervals: mean ± k\*SE (+ observation variance if PI)

4.  Enforce monotonicity with cummax() on bounds

This is ~25x faster than bootstrap because:

- Fits model only once (vs 200 times)

- Uses analytic variance estimates (vs resampling)

Monotonicity enforcement:

- Predictions are monotonic by construction (constrained fit)

- Lower/upper bounds may violate due to varying SE

- cummax() ensures bounds are also monotonic

- May slightly overestimate uncertainty at inflection points
