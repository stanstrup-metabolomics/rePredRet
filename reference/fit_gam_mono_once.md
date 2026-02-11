# Fit Monotonically Constrained GAM (Single Fit)

Fits a monotonically constrained GAM once to the full dataset. This is
used as the basis for posterior simulation.

## Usage

``` r
fit_gam_mono_once(rt_matrix, newdata)
```

## Arguments

- rt_matrix:

  A 2-column matrix where column 1 is RT in source system and column 2
  is RT in target system.

- newdata:

  Numeric vector of x-values at which to predict.

## Value

A list containing:

- fit:

  The fitted GAM object

- coef:

  Coefficient vector

- vcov:

  Variance-covariance matrix

- sigma:

  Residual standard deviation

- Xp:

  Linear predictor matrix at newdata

- pred:

  Point predictions at newdata

- newdata:

  The newdata vector
