# Monotone Interval from a Fit and Half-Widths

Each bound (`mu -/+ half_width`) is made non-decreasing in `newdata` by
isotonic regression, which keeps lower \<= mu \<= upper when `mu` is
monotone. The lower bound is truncated at 0.

## Usage

``` r
monotone_interval(newdata, mu, half_width)
```

## Arguments

- newdata:

  Numeric vector of source RTs.

- mu:

  Monotone predictions at `newdata`.

- half_width:

  Half-widths at `newdata`.

## Value

A matrix with 3 columns: `[predicted, lower, upper]`.
