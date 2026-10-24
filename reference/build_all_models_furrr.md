# Build All Models with furrr - Optimized with Caching and Scheduling

Builds retention time prediction models using furrr for efficient
outer-level parallelization with support for incremental builds and
scheduling tuning.

## Usage

``` r
build_all_models_furrr(
  report_data,
  min_compounds = 10,
  method = c("adaptive_pi", "fast_ci", "bootstrap"),
  alpha = 0.05,
  n_workers = parallelly::availableCores(),
  save_json = TRUE,
  export_dir = NULL,
  batch_size = NULL,
  report_version = NULL,
  verbose = TRUE
)
```

## Arguments

- report_data:

  Report data from load_report_data()

- min_compounds:

  Minimum overlapping compounds for a pair (default 10)

- method:

  "adaptive_pi" (default), "fast_ci" or "bootstrap" (slow). See
  [`build_model()`](https://stanstrup-metabolomics.github.io/rePredRet/reference/build_model.md).
  Cached models built with another method, method version or alpha are
  rebuilt.

- alpha:

  Significance level (default 0.05 for 95% intervals)

- n_workers:

  Number of parallel workers (default = cores available to this process,
  respecting container/cgroup CPU limits)

- save_json:

  Kept for backwards compatibility.
  [`export_model_fast()`](https://stanstrup-metabolomics.github.io/rePredRet/reference/export_model_fast.md)
  always writes model.json.

- export_dir:

  Directory to save models (required for caching)

- batch_size:

  Number of models to send to each worker at once (NULL = auto, 1 =
  dynamic, large = upfront)

- report_version:

  RepoRT version identifier (e.g. commit SHA) to record in the build
  cache. NULL leaves the recorded version unchanged.

- verbose:

  Print progress messages (default TRUE)

## Value

List with:

- models: list of results for the models built in this run

- index: data.frame with metadata for all current models

- stats: summary statistics

- changed: TRUE if any model was built, added or removed in this run

## Details

This function uses furrr::future_map for outer-level parallelization,
avoiding nested parallelization overhead. Each worker builds models
sequentially.

Models are only built between systems of the same method type
(`method.type` in `report_data$studies`, e.g. RP or HILIC), since the
elution order differs between types.

Caching (requires export_dir): models are only rebuilt when one of their
datasets changed or their files are missing. Models whose datasets were
removed, or that no longer have enough common compounds, are deleted.
The build cache, `model_index.csv` and `systems.csv` (one row per
system: name, method type, number of compounds and number of models) are
updated in export_dir.
