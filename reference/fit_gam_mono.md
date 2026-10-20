# Fit a Monotonically Constrained GAM

The monotone fit shared by all interval methods.

## Usage

``` r
fit_gam_mono(x, y, y_max = max(y))
```

## Arguments

- x:

  RTs in the source system.

- y:

  RTs in the target system.

- y_max:

  Value the absolute residuals are divided by before the sigmoid
  weighting (default `max(y)`).

## Value

A list with the (possibly jittered) data `x` and `y`, the basis
dimension `k`, the unconstrained fit `fit_unc`, the sigmoid weights `w`,
the constrained smooth `sm` and its coefficients `p`.

## Details

1.  If there are fewer than 4 unique source RTs (rounded to 0.01), the
    source RTs are jittered by up to 0.01 with a fixed seed, so
    identical data always give the same model.

2.  An unconstrained thin plate regression spline GAM gives the
    smoothing parameter.

3.  Observations are weighted by a sigmoid of their unconstrained
    residual relative to `y_max`, which down-weights outliers: weights
    fall to 0.5 at a residual of 10% of `y_max`, with steepness 30.
    These two values were set by eye for the original PredRet (2015) and
    have never been optimized.

4.  A cubic regression spline is fitted with monotonic increasing
    constraints by penalized constrained least squares. If the weighted
    fit fails it is refitted without weights.
