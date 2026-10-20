# Filter Model Pairs Based on Changes

Determines which model pairs need to be built based on dataset changes
and on the interval method and alpha recorded for each cached model.

## Usage

``` r
filter_model_pairs(
  model_pairs,
  changes,
  cache,
  export_dir,
  method = NULL,
  alpha = NULL
)
```

## Arguments

- model_pairs:

  List of all potential model pairs

- changes:

  Dataset change analysis

- cache:

  Build cache object

- export_dir:

  Directory for model outputs

- method:

  Interval method of this build. Cached models built with another
  method, or another version of it
  ([`pi_method_version()`](https://stanstrup-metabolomics.github.io/rePredRet/reference/pi_method_version.md)),
  are rebuilt. NULL skips the check.

- alpha:

  Significance level of this build. Cached models built with another
  alpha are rebuilt. NULL skips the check.

## Value

Filtered list of pairs that need building
