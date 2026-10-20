# Changelog

## rePredRet (development version)

### Prediction intervals

- New default interval method `method = "adaptive_pi"`
  ([`gam_adaptive_pi()`](https://stanstrup-metabolomics.github.io/rePredRet/reference/gam_adaptive_pi.md))
  in
  [`build_model()`](https://stanstrup-metabolomics.github.io/rePredRet/reference/build_model.md),
  [`build_all_models()`](https://stanstrup-metabolomics.github.io/rePredRet/reference/build_all_models.md),
  [`build_all_models_furrr()`](https://stanstrup-metabolomics.github.io/rePredRet/reference/build_all_models_furrr.md)
  and
  [`build_user_models()`](https://stanstrup-metabolomics.github.io/rePredRet/reference/build_user_models.md).
  The interval width follows the local scatter of the calibration
  compounds (a Gamma GAM of cross-validated absolute residuals, weighted
  to ignore outliers) and is calibrated by cross-validated conformal
  prediction. On held-out compounds from 135 RepoRT system pairs it
  covered 94.6% (95% nominal), equally in low- and high-scatter regions
  (94.8% / 94.5%). It takes 60-90 ms per model.
- `method = "fast_ci"`
  ([`gam_fast_ci()`](https://stanstrup-metabolomics.github.io/rePredRet/reference/gam_fast_ci.md))
  is kept as the fast option (15-20 ms per model), with fixes: the
  residual variance comes from leave-one-out residuals instead of
  in-sample residuals (coverage 93.9% -\> 95.5%; with 10-15 compounds
  83% -\> 90%), and the bounds are made monotone by isotonic regression
  instead of a running maximum, which narrowed the lower side. Its width
  is still nearly constant (97.4% / 93.4% coverage in low- /
  high-scatter regions).
- `method = "bootstrap"` was a confidence interval of the curve, not a
  prediction interval (59% coverage; 53% with 200+ compounds). It now
  uses the out-of-bag predictions as out-of-sample residuals with the
  same local scale and conformal calibration as `adaptive_pi` (94.8%
  coverage).
  [`boot2ci_PI()`](https://stanstrup-metabolomics.github.io/rePredRet/reference/boot2ci_PI.md)
  has a new signature and the unused `boot2ci()` was removed.
- With 10-15 compounds no method reaches the nominal coverage (86-90%);
  this is now documented.
- All methods truncate the lower bound at 0.
- Model statistics (`model.json` `stats` and `model$stats`) gain
  `r_squared`, `coverage` (fraction of calibration compounds inside the
  PI) and `rel_mean_ci_width` (mean PI width / target RT range).
- Removed the unused `gam_prediction_intervals()`,
  `gam_prediction_intervals_parallel()`, `fit_gam_mono_once()` and
  `gam_fast_pi()`.
- [`build_all_models_furrr()`](https://stanstrup-metabolomics.github.io/rePredRet/reference/build_all_models_furrr.md)
  and
  [`build_all_models()`](https://stanstrup-metabolomics.github.io/rePredRet/reference/build_all_models.md)
  rebuild cached models that were built with another interval method,
  method version
  ([`pi_method_version()`](https://stanstrup-metabolomics.github.io/rePredRet/reference/pi_method_version.md))
  or alpha, so the next pipeline run rebuilds all models.
- Small source-RT jitter (fewer than 4 distinct RTs) uses a fixed seed,
  so the same data always give the same model.
- [`build_model()`](https://stanstrup-metabolomics.github.io/rePredRet/reference/build_model.md)
  returns status `"fit_failed"` instead of a model with non-finite
  predictions or intervals.
- The monotone GAM fit is shared by all interval methods
  ([`fit_gam_mono()`](https://stanstrup-metabolomics.github.io/rePredRet/reference/fit_gam_mono.md)).
  Its sigmoid outlier weights (0.1, 30) were set by eye for the original
  PredRet and have never been optimized; this is now documented.
- New article “Prediction interval benchmark”, and the benchmark script
  `inst/benchmark/pi_benchmark.R`.
- `withr` moved from Suggests to Imports; `MASS` is no longer needed.

### Breaking changes

- Removed `rePredRet_download()`, `rePredRet_predict()`,
  `rePredRet_system_predictions()` and `rePredRet_plot_errors()`. They
  depended on a `rePredRet-data` repository with pre-computed
  predictions that does not exist. Use
  [`build_user_models()`](https://stanstrup-metabolomics.github.io/rePredRet/reference/build_user_models.md)
  and
  [`predict_rt_user()`](https://stanstrup-metabolomics.github.io/rePredRet/reference/predict_rt_user.md)
  to predict RTs.
- Removed `rePredRet_plot_model()`; use
  [`plot_model()`](https://stanstrup-metabolomics.github.io/rePredRet/reference/plot_model.md).
- [`rePredRet_systems()`](https://stanstrup-metabolomics.github.io/rePredRet/reference/rePredRet_systems.md),
  [`rePredRet_models()`](https://stanstrup-metabolomics.github.io/rePredRet/reference/rePredRet_models.md),
  [`rePredRet_model()`](https://stanstrup-metabolomics.github.io/rePredRet/reference/rePredRet_model.md),
  [`plot_model()`](https://stanstrup-metabolomics.github.io/rePredRet/reference/plot_model.md)
  and
  [`rePredRet_plot_network()`](https://stanstrup-metabolomics.github.io/rePredRet/reference/rePredRet_plot_network.md)
  now fetch files on demand from the rePredRet-models repository and
  cache them locally. Their `data_path` argument is replaced by
  `refresh`.
  [`rePredRet_model()`](https://stanstrup-metabolomics.github.io/rePredRet/reference/rePredRet_model.md)
  returns the contents of `model.json`.
- Requires R \>= 4.4.0.

### Bug fixes

- [`save_model_plot()`](https://stanstrup-metabolomics.github.io/rePredRet/reference/save_model_plot.md)
  (and `build_all_models(save_plots = TRUE)`) failed because the
  internal plot function was masked by the exported
  [`plot_model()`](https://stanstrup-metabolomics.github.io/rePredRet/reference/plot_model.md).
- [`plot_model()`](https://stanstrup-metabolomics.github.io/rePredRet/reference/plot_model.md)
  read model statistics from the wrong index columns.
- [`search_compound()`](https://stanstrup-metabolomics.github.io/rePredRet/reference/search_compound.md)
  and
  [`get_system_info()`](https://stanstrup-metabolomics.github.io/rePredRet/reference/get_system_info.md)
  looked for RepoRT data in the wrong directory when downloading it.
- Model viewer: loads the index from the right location, parses it
  correctly, plots the model data and draws the interval band between
  the bounds. It is now part of the website.
- [`build_user_models()`](https://stanstrup-metabolomics.github.io/rePredRet/reference/build_user_models.md),
  [`build_model()`](https://stanstrup-metabolomics.github.io/rePredRet/reference/build_model.md)
  and
  [`build_all_models()`](https://stanstrup-metabolomics.github.io/rePredRet/reference/build_all_models.md)
  default to
  [`parallelly::availableCores()`](https://parallelly.futureverse.org/reference/availableCores.html).

### Website

- The Overview, Model Statistics and Network articles are computed from
  the published model index when the site is built, instead of from
  statistics files committed in January (a 642-model test run). The site
  is rebuilt after each model update.
- Removed the Predictions article and the committed prediction and
  statistics files; pre-computed predictions are no longer produced.
- Removed the old standalone Quarto site in `website/`; the pkgdown site
  replaces it.

### Other

- [`rePredRet_plot_network()`](https://stanstrup-metabolomics.github.io/rePredRet/reference/rePredRet_plot_network.md)
  gains `show_labels`.
- `inst/scripts/run_pipeline.R` writes to `output/` by default (was
  `website/data/`).
- Runner setup (`github_runner.yml`, `RUNNER_SETUP.md`): host paths are
  set in `.env` and must be on local disks; adds the
  `/opt/repredret-cache` mount.

### Pipeline

- [`build_all_models_furrr()`](https://stanstrup-metabolomics.github.io/rePredRet/reference/build_all_models_furrr.md)
  now builds incrementally: only models whose datasets changed or whose
  files are missing are rebuilt. Previously the build cache was loaded
  but never applied or saved, so every run rebuilt all models.
- [`build_all_models_furrr()`](https://stanstrup-metabolomics.github.io/rePredRet/reference/build_all_models_furrr.md)
  removes models for removed datasets and pairs that no longer have
  enough common compounds, keeps `model_index.csv` complete across
  incremental runs, records `report_version` in the build cache and
  returns `changed`.
- Model pairs are found with a sparse overlap matrix and RT matrices
  come from a per-dataset lookup, instead of calling
  [`get_common_compounds()`](https://stanstrup-metabolomics.github.io/rePredRet/reference/get_common_compounds.md)
  for every ordered pair (minutes to seconds on the full RepoRT data).
- The number of workers defaults to
  [`parallelly::availableCores()`](https://parallelly.futureverse.org/reference/availableCores.html),
  which respects container CPU limits
  ([`parallel::detectCores()`](https://rdrr.io/r/parallel/detectCores.html)
  reports all host cores).
- A failing model no longer aborts the whole build.
- Models whose build settings are unknown (no build cache entry,
  e.g. from an older pipeline) are rebuilt when an interval method is
  given, instead of being recorded as built with the current method.
- CI: tracks RepoRT master commits instead of releases and only pushes
  to the data repository when models changed. `max_systems` test runs no
  longer touch the published models.
- CI: models are built directly in a persistent clone of
  rePredRet-models on the runner (`/opt/repredret-cache`), which is
  reset to the published state at the start of each run. Only changed
  model directories are committed and pushed, in chunks, instead of
  force-pushing a freshly created repository with all models.

## rePredRet 0.2.0

### New Features

#### User Model Building

- [`build_user_models()`](https://stanstrup-metabolomics.github.io/rePredRet/reference/build_user_models.md):
  Build calibration models from your own RT data to all compatible
  RepoRT systems
- [`predict_rt_user()`](https://stanstrup-metabolomics.github.io/rePredRet/reference/predict_rt_user.md):
  Generate RT predictions using user-built models
- `UserRTModels` S4 class: Stores models with diagnostics, supports
  `show()`, [`summary()`](https://rdrr.io/r/base/summary.html),
  [`plot()`](https://rdrr.io/r/graphics/plot.default.html), `[[]]`
  methods
- `UserRTPredictions` S4 class: Stores predictions, supports `show()`,
  [`summary()`](https://rdrr.io/r/base/summary.html),
  [`as.data.frame()`](https://rdrr.io/r/base/as.data.frame.html),
  [`plot()`](https://rdrr.io/r/graphics/plot.default.html) methods

#### Utility Functions

- [`plot_model()`](https://stanstrup-metabolomics.github.io/rePredRet/reference/plot_model.md):
  Plot pre-built models from the data repository with ggplot2/plotly
- [`search_compound()`](https://stanstrup-metabolomics.github.io/rePredRet/reference/search_compound.md):
  Search for compounds by InChI across all systems
- [`search_compounds()`](https://stanstrup-metabolomics.github.io/rePredRet/reference/search_compounds.md):
  Batch search for multiple InChIs
- [`compare_systems()`](https://stanstrup-metabolomics.github.io/rePredRet/reference/compare_systems.md):
  Rank repository systems by similarity to user’s data
- [`get_system_info()`](https://stanstrup-metabolomics.github.io/rePredRet/reference/get_system_info.md):
  Get detailed system metadata (column, gradient, etc.)
- [`export_predictions()`](https://stanstrup-metabolomics.github.io/rePredRet/reference/export_predictions.md):
  Export predictions to CSV or Excel format
- [`validate_models()`](https://stanstrup-metabolomics.github.io/rePredRet/reference/validate_models.md):
  Cross-validation for assessing model quality
- [`check_input_data()`](https://stanstrup-metabolomics.github.io/rePredRet/reference/check_input_data.md):
  Validate InChIs and RTs before model building

#### Visualizations

- [`plot_system_network()`](https://stanstrup-metabolomics.github.io/rePredRet/reference/plot_system_network.md):
  Visualize system connectivity
- [`plot_prediction_errors()`](https://stanstrup-metabolomics.github.io/rePredRet/reference/plot_prediction_errors.md):
  Error distribution plots
- [`plot_compound_coverage()`](https://stanstrup-metabolomics.github.io/rePredRet/reference/plot_compound_coverage.md):
  Heatmap of compound presence across systems

### Improvements

- All InChI matching now consistently strips stereochemistry layers (/t,
  /b, /m, /s)
- [`predict_rt_user()`](https://stanstrup-metabolomics.github.io/rePredRet/reference/predict_rt_user.md)
  returns only the best prediction per compound (narrowest CI)
- Added progress bars for long-running operations
- Added comprehensive unit tests
- Added pkgdown documentation website
- Added GitHub Actions for R CMD check and code coverage

### Documentation

- New tutorial: “Predict RT for Your Own Data”
- New function reference page
- Vignette for offline documentation

------------------------------------------------------------------------

## rePredRet 0.1.0

- Initial release
- Core model building functions:
  [`build_model()`](https://stanstrup-metabolomics.github.io/rePredRet/reference/build_model.md),
  [`build_all_models()`](https://stanstrup-metabolomics.github.io/rePredRet/reference/build_all_models.md),
  [`build_all_models_furrr()`](https://stanstrup-metabolomics.github.io/rePredRet/reference/build_all_models_furrr.md)
- Prediction functions:
  [`predict_rt()`](https://stanstrup-metabolomics.github.io/rePredRet/reference/predict_rt.md),
  [`generate_all_predictions()`](https://stanstrup-metabolomics.github.io/rePredRet/reference/generate_all_predictions.md)
- Data loading:
  [`download_report()`](https://stanstrup-metabolomics.github.io/rePredRet/reference/download_report.md),
  [`load_report_data()`](https://stanstrup-metabolomics.github.io/rePredRet/reference/load_report_data.md),
  [`get_common_compounds()`](https://stanstrup-metabolomics.github.io/rePredRet/reference/get_common_compounds.md)
- User query functions: `rePredRet_download()`, `rePredRet_predict()`,
  [`rePredRet_systems()`](https://stanstrup-metabolomics.github.io/rePredRet/reference/rePredRet_systems.md),
  [`rePredRet_models()`](https://stanstrup-metabolomics.github.io/rePredRet/reference/rePredRet_models.md)
- GAM fitting with monotonic constraints and confidence intervals
- Fast CI calculation method (72x faster than bootstrap)
