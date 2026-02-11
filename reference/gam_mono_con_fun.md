# Monotonically Constrained GAM Bootstrap Function

Fits a monotonically constrained GAM to bootstrap-resampled data. This
is the core fitting function used by
[`boot::boot()`](https://rdrr.io/pkg/boot/man/boot.html) for confidence
interval estimation.

## Usage

``` r
gam_mono_con_fun(in_data, inds, newdata)
```

## Arguments

- in_data:

  A 2-column matrix where column 1 is RT in source system and column 2
  is RT in target system.

- inds:

  Bootstrap resampling indices (provided by `boot()`).

- newdata:

  Numeric vector of x-values at which to predict.

## Value

Numeric vector of predicted y-values at `newdata` points.

## Details

The algorithm:

1.  Resample data using bootstrap indices

2.  Add jitter if fewer than 4 unique x-values (needed for spline
    fitting)

3.  Fit initial unconstrained GAM to get smoothing parameter

4.  Weight observations by residuals using sigmoid function (downweights
    outliers)

5.  Build constrained spline with monotonic constraints

6.  Fit using penalized constrained least squares

7.  Predict at newdata points

The monotonic constraint ensures that predicted RTs maintain the elution
order assumption (if compound A elutes before B in system 1, A should
also elute before B in system 2).

## References

Stanstrup et al. (2015) "PredRet: Prediction of Retention Time by Direct
Mapping between Multiple Chromatographic Systems" Anal. Chem.
87:9421-9428
