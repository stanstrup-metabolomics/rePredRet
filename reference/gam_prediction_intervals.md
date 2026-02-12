# GAM Prediction Intervals via Posterior Simulation

Fast calculation of prediction intervals for monotonically constrained
GAMs using posterior simulation instead of bootstrap.

Uses posterior simulation to generate prediction intervals for a
monotonically constrained GAM. This is much faster than bootstrap
because the model is only fitted once.

## Usage

``` r
gam_prediction_intervals(
  rt_matrix,
  newdata,
  n_sim = 1000,
  alpha = 0.05,
  n_cores = 1
)
```

## Arguments

- rt_matrix:

  A 2-column matrix where column 1 is RT in source system and column 2
  is RT in target system.

- newdata:

  Numeric vector of x-values at which to predict.

- n_sim:

  Number of posterior samples (default 1000).

- alpha:

  Significance level for PI (default 0.05 for 95% PI).

- n_cores:

  Number of cores for parallel simulation (default 1).

## Value

A matrix with 3 columns: `[predicted, lower, upper]` at each newdata
point.

## Details

The algorithm:

1.  Fit constrained GAM once to full data

2.  Extract coefficient vector β and covariance matrix V

3.  Simulate parameter vectors: β_sim = β + L*ν where V = L*L'
    (Cholesky) and ν ~ N(0,1)

4.  Generate predictions: pred_sim = Xp %\*% β_sim

5.  Add observation noise: y_sim ~ N(pred_sim, sigma^2)

6.  Calculate quantiles across simulations

This gives TRUE prediction intervals (for future observations), not
confidence intervals (for the mean function). Prediction intervals are
wider because they include both parameter uncertainty AND observation
variance.

## References

Marra & Wood (2012) "Coverage properties of confidence intervals for
generalized additive model components". Scandinavian Journal of
Statistics.
