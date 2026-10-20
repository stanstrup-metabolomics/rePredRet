# Interval Method Version

Version of each interval method's implementation. Bump it when a
method's intervals change, so that cached models built with the old
version are rebuilt by
[`build_all_models_furrr()`](https://stanstrup-metabolomics.github.io/rePredRet/reference/build_all_models_furrr.md).

## Usage

``` r
pi_method_version(method)
```

## Arguments

- method:

  Interval method name.

## Value

Integer version.
