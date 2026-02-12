# Calculate Bootstrap Prediction Intervals

Computes prediction intervals that account for both model uncertainty
and residual variance.

## Usage

``` r
boot2ci_PI(loess_boot, newdata, alpha = 0.05)
```

## Arguments

- loess_boot:

  A boot object from
  [`boot::boot()`](https://rdrr.io/pkg/boot/man/boot.html).

- newdata:

  Numeric vector of x-values used for prediction.

- alpha:

  Significance level for PI (default 0.05 for 95% PI).

## Value

A matrix with 3 columns: `[predicted, lower, upper]` at each newdata
point.

## Details

Prediction intervals are wider than confidence intervals because they
account for:

1.  Uncertainty in the model fit (from bootstrap)

2.  Residual variance (deviation of individual points from the fit)
