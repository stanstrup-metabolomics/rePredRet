# Calculate Bootstrap Prediction Intervals

Prediction intervals from bootstrap refits of the monotone GAM, using
the out-of-bag predictions as out-of-sample residuals.

## Usage

``` r
boot2ci_PI(loess_boot, fit, newdata, alpha = 0.05)
```

## Arguments

- loess_boot:

  A boot object from
  [`boot::boot()`](https://rdrr.io/pkg/boot/man/boot.html) with
  [`gam_mono_con_fun()`](https://stanstrup-metabolomics.github.io/rePredRet/reference/gam_mono_con_fun.md),
  run with `newdata = c(newdata, rt_source)`, i.e. predicting at the
  grid followed by the calibration compounds. It must not be pruned yet.

- fit:

  Output of
  [`fit_gam_mono()`](https://stanstrup-metabolomics.github.io/rePredRet/reference/fit_gam_mono.md)
  on the full data.

- newdata:

  Numeric vector of x-values for the prediction grid.

- alpha:

  Significance level for PI (default 0.05 for 95% PI).

## Value

A matrix with 3 columns: `[predicted, lower, upper]` at each newdata
point.

## Details

Each compound is predicted by the refits whose resample did not contain
it (out-of-bag, on average 37% of the refits). Its residual from the
mean out-of-bag prediction is an out-of-sample error that includes the
uncertainty of the fitted curve. A compound that is in every resample
gets the leave-one-out residual of the unconstrained GAM instead. These
residuals then go through the same local scale model and conformal
calibration as
[`gam_adaptive_pi()`](https://stanstrup-metabolomics.github.io/rePredRet/reference/gam_adaptive_pi.md)
([`conformal_band()`](https://stanstrup-metabolomics.github.io/rePredRet/reference/conformal_band.md)),
around the full-data fit. This is the jackknife+-after-bootstrap idea of
Kim, Xu & Barber (2020).

## References

Kim, B., Xu, C. & Barber, R. F. (2020). Predictive inference is free
with the jackknife+-after-bootstrap. *Advances in Neural Information
Processing Systems*, 33, 4138-4149.
