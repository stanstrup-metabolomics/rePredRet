# Penalized Monotone Spline Coefficients

Penalized Monotone Spline Coefficients

## Usage

``` r
fit_mono_spline(sm, y, w, sp)
```

## Arguments

- sm:

  A smooth from
  [`mgcv::smoothCon()`](https://rdrr.io/pkg/mgcv/man/smoothCon.html).

- y:

  Response.

- w:

  Observation weights. If the weighted fit fails, it is refitted without
  weights.

- sp:

  Smoothing parameter.

## Value

Coefficient vector.
