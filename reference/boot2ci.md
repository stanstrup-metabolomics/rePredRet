# Calculate Bootstrap Confidence Intervals

Computes percentile bootstrap confidence intervals from a boot object.

## Usage

``` r
boot2ci(loess_boot, alpha = 0.01)
```

## Arguments

- loess_boot:

  A boot object from
  [`boot::boot()`](https://rdrr.io/pkg/boot/man/boot.html).

- alpha:

  Significance level for CI (default 0.01 for 99% CI).

## Value

A matrix with 3 columns: predicted, lower, upper at each newdata point.
