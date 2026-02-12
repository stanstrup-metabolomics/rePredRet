# Build All Pairwise Models

Builds models between all pairs of chromatographic systems that share
sufficient common compounds.

## Usage

``` r
build_all_models(
  report_data,
  min_compounds = 10,
  n_boot = 200,
  alpha = 0.05,
  n_cores = parallel::detectCores(),
  method_match = TRUE,
  export_dir = NULL,
  method = c("fast_ci", "bootstrap"),
  save_plots = TRUE
)
```

## Arguments

- report_data:

  Data from
  [`load_report_data()`](https://stanstrup.github.io/rePredRet/reference/load_report_data.md).

- min_compounds:

  Minimum common compounds required (default 10).

- n_boot:

  Number of bootstrap iterations (default 1000).

- alpha:

  Significance level for CI (default 0.01).

- n_cores:

  Number of CPU cores for parallel processing.

- method_match:

  If TRUE (default), only build models between same method types (RP-RP,
  HILIC-HILIC).

- export_dir:

  Optional directory to incrementally save models as they're built. If
  NULL (default), models are only stored in memory.

- method:

  Method for confidence intervals: "fast_ci" (default, fast) or
  "bootstrap" (accurate, slower).

- save_plots:

  Logical, whether to save interactive HTML plots for each model
  (default TRUE).

## Value

A list with:

- models:

  Named list of model objects

- index:

  Tibble with model index and statistics
