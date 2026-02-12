# Fast Prediction Intervals using predict() + cummax

Fast calculation of prediction intervals by fitting GAM once and
enforcing monotonicity with cummax. 72x faster than bootstrap.

Computes prediction intervals efficiently by:

1.  Fitting constrained GAM to get monotonic predictions

2.  Estimating observation variance from residuals

3.  Building PIs: mean +/- z \* sqrt(SE^2 + sigma^2)

4.  Enforcing monotonicity on bounds with cummax()

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

1.  Fit constrained GAM for monotonic predictions

2.  Calculate residual variance (sigma^2) from constrained fit

3.  Get SE estimates from unconstrained fit

4.  Build PIs: mean +/- z\*sqrt(SE^2 + sigma^2) (includes observation
    variance)

5.  Enforce monotonicity: lower = cummax(lower), upper = cummax(upper)

This is ~72x faster than bootstrap (0.02 sec vs 1.5 sec per model)
because:

- Fits models once instead of 200 times

- Uses analytic SE + variance instead of resampling

- Monotonicity enforcement is instant (cummax)

Trade-off:

- PIs are wider than CIs (but narrower than bootstrap-based PIs)

- More realistic for predicting individual observations

## References

"Fit monotone GAM -\> Calculate sigma^2 from residuals -\> Get fit +/-
z\*sqrt(SE^2 + sigma^2) -\> Enforce monotonicity afterwards" - Hybrid
approach for fast PIs
