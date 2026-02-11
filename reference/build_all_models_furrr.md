# Build All Models with furrr - Optimized with Caching and Scheduling

Builds retention time prediction models using furrr for efficient
outer-level parallelization with support for incremental builds and
scheduling tuning.

## Usage

``` r
build_all_models_furrr(
  report_data,
  min_compounds = 10,
  method = c("fast_ci", "bootstrap"),
  alpha = 0.05,
  n_workers = parallel::detectCores(),
  save_json = TRUE,
  export_dir = NULL,
  batch_size = NULL,
  verbose = TRUE
)
```

## Arguments

- report_data:

  Report data from load_report_data()

- min_compounds:

  Minimum overlapping compounds for a pair (default 10)

- method:

  "fast_ci" (default, fast) or "bootstrap" (slow, accurate)

- alpha:

  Significance level (default 0.05 for 95% intervals)

- n_workers:

  Number of parallel workers (default = available cores)

- save_json:

  Save model data as JSON (default TRUE, much faster than HTML)

- export_dir:

  Directory to save models (required for caching)

- batch_size:

  Number of models to send to each worker at once (NULL = auto, 1 =
  dynamic, large = upfront)

- verbose:

  Print progress messages (default TRUE)

## Value

List with:

- models: list of all built models

- index: data.frame with model metadata

- stats: summary statistics

## Details

This function uses furrr::future_map for outer-level parallelization,
avoiding nested parallelization overhead. Each worker builds models
sequentially.

Caching: Detects changed datasets and skips unchanged models (requires
export_dir)
