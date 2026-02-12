# Calculate Prediction Intervals via Posterior Simulation (Parallel)

Parallelized version of posterior simulation for very large newdata
grids.

## Usage

``` r
gam_prediction_intervals_parallel(
  rt_matrix,
  newdata,
  n_sim = 1000,
  alpha = 0.01,
  n_cores = 2
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
