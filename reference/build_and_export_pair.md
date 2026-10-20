# Build and Export a Single Model Pair

Worker function for
[`build_all_models_furrr()`](https://stanstrup-metabolomics.github.io/rePredRet/reference/build_all_models_furrr.md).
Defined at top level so that only its arguments are sent to parallel
workers.

## Usage

``` r
build_and_export_pair(pair, alpha, method, export_dir)
```

## Arguments

- pair:

  List with sys1_id, sys2_id and rt_matrix

- alpha:

  Significance level

- method:

  CI method

- export_dir:

  Directory to save models, or NULL

## Value

List with status and summary statistics
