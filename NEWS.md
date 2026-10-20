# rePredRet (development version)

## Prediction intervals

* New default interval method `method = "adaptive_pi"` (`gam_adaptive_pi()`) in
  `build_model()`, `build_all_models()`, `build_all_models_furrr()` and
  `build_user_models()`. The interval width follows the local scatter of the
  calibration compounds (a Gamma GAM of cross-validated absolute residuals,
  weighted to ignore outliers) and is calibrated by cross-validated conformal
  prediction. On held-out compounds from 135 RepoRT system pairs it covered
  94.6% (95% nominal), equally in low- and high-scatter regions, where the old
  near-constant interval covered 96.5% and 91.2%. It takes about 60 ms per
  model instead of 15 ms. The bounds are made monotone by isotonic regression
  instead of a running maximum, and the lower bound is truncated at 0.
* Removed `method = "fast_ci"` and `gam_fast_ci()`, and the unused
  `gam_prediction_intervals()`, `gam_prediction_intervals_parallel()`,
  `fit_gam_mono_once()` and `gam_fast_pi()`.
* `build_all_models_furrr()` and `build_all_models()` rebuild cached models that
  were built with another interval method or alpha, so the next pipeline run
  rebuilds all models.
* Small source-RT jitter (fewer than 4 distinct RTs) uses a fixed seed, so the
  same data always give the same model.
* `build_model()` returns status `"fit_failed"` instead of a model with
  non-finite predictions or intervals.
* The monotone GAM fit is shared by all interval methods (`fit_gam_mono()`).
* `withr` moved from Suggests to Imports; `MASS` is no longer needed.

## Breaking changes

* Removed `rePredRet_download()`, `rePredRet_predict()`,
  `rePredRet_system_predictions()` and `rePredRet_plot_errors()`. They depended on a
  `rePredRet-data` repository with pre-computed predictions that does not exist.
  Use `build_user_models()` and `predict_rt_user()` to predict RTs.
* Removed `rePredRet_plot_model()`; use `plot_model()`.
* `rePredRet_systems()`, `rePredRet_models()`, `rePredRet_model()`, `plot_model()` and
  `rePredRet_plot_network()` now fetch files on demand from the rePredRet-models
  repository and cache them locally. Their `data_path` argument is replaced by
  `refresh`. `rePredRet_model()` returns the contents of `model.json`.
* Requires R >= 4.4.0.

## Bug fixes

* `save_model_plot()` (and `build_all_models(save_plots = TRUE)`) failed because the
  internal plot function was masked by the exported `plot_model()`.
* `plot_model()` read model statistics from the wrong index columns.
* `search_compound()` and `get_system_info()` looked for RepoRT data in the wrong
  directory when downloading it.
* Model viewer: loads the index from the right location, parses it correctly, plots the
  model data and draws the interval band between the bounds. It is now part of the
  website.
* `build_user_models()`, `build_model()` and `build_all_models()` default to
  `parallelly::availableCores()`.

## Website

* The Overview, Model Statistics and Network articles are computed from the published
  model index when the site is built, instead of from statistics files committed in
  January (a 642-model test run). The site is rebuilt after each model update.
* Removed the Predictions article and the committed prediction and statistics files;
  pre-computed predictions are no longer produced.
* Removed the old standalone Quarto site in `website/`; the pkgdown site replaces it.

## Other

* `rePredRet_plot_network()` gains `show_labels`.
* `inst/scripts/run_pipeline.R` writes to `output/` by default (was `website/data/`).
* Runner setup (`github_runner.yml`, `RUNNER_SETUP.md`): host paths are set in `.env`
  and must be on local disks; adds the `/opt/repredret-cache` mount.

## Pipeline

* `build_all_models_furrr()` now builds incrementally: only models whose datasets
  changed or whose files are missing are rebuilt. Previously the build cache was
  loaded but never applied or saved, so every run rebuilt all models.
* `build_all_models_furrr()` removes models for removed datasets and pairs that no
  longer have enough common compounds, keeps `model_index.csv` complete across
  incremental runs, records `report_version` in the build cache and returns `changed`.
* Model pairs are found with a sparse overlap matrix and RT matrices come from a
  per-dataset lookup, instead of calling `get_common_compounds()` for every
  ordered pair (minutes to seconds on the full RepoRT data).
* The number of workers defaults to `parallelly::availableCores()`, which respects
  container CPU limits (`parallel::detectCores()` reports all host cores).
* A failing model no longer aborts the whole build.
* CI: tracks RepoRT master commits instead of releases and only pushes to the data
  repository when models changed. `max_systems` test runs no longer touch the
  published models.
* CI: models are built directly in a persistent clone of rePredRet-models on the
  runner (`/opt/repredret-cache`), which is reset to the published state at the
  start of each run. Only changed model directories are committed and pushed, in
  chunks, instead of force-pushing a freshly created repository with all models.

# rePredRet 0.2.0

## New Features

### User Model Building
* `build_user_models()`: Build calibration models from your own RT data to all compatible RepoRT systems
* `predict_rt_user()`: Generate RT predictions using user-built models
* `UserRTModels` S4 class: Stores models with diagnostics, supports `show()`, `summary()`, `plot()`, `[[]]` methods
* `UserRTPredictions` S4 class: Stores predictions, supports `show()`, `summary()`, `as.data.frame()`, `plot()` methods

### Utility Functions
* `plot_model()`: Plot pre-built models from the data repository with ggplot2/plotly
* `search_compound()`: Search for compounds by InChI across all systems
* `search_compounds()`: Batch search for multiple InChIs
* `compare_systems()`: Rank repository systems by similarity to user's data
* `get_system_info()`: Get detailed system metadata (column, gradient, etc.)
* `export_predictions()`: Export predictions to CSV or Excel format
* `validate_models()`: Cross-validation for assessing model quality
* `check_input_data()`: Validate InChIs and RTs before model building

### Visualizations
* `plot_system_network()`: Visualize system connectivity
* `plot_prediction_errors()`: Error distribution plots
* `plot_compound_coverage()`: Heatmap of compound presence across systems

## Improvements

* All InChI matching now consistently strips stereochemistry layers (/t, /b, /m, /s)
* `predict_rt_user()` returns only the best prediction per compound (narrowest CI)
* Added progress bars for long-running operations
* Added comprehensive unit tests
* Added pkgdown documentation website
* Added GitHub Actions for R CMD check and code coverage

## Documentation

* New tutorial: "Predict RT for Your Own Data"
* New function reference page
* Vignette for offline documentation

---

# rePredRet 0.1.0

* Initial release
* Core model building functions: `build_model()`, `build_all_models()`, `build_all_models_furrr()`
* Prediction functions: `predict_rt()`, `generate_all_predictions()`
* Data loading: `download_report()`, `load_report_data()`, `get_common_compounds()`
* User query functions: `rePredRet_download()`, `rePredRet_predict()`, `rePredRet_systems()`, `rePredRet_models()`
* GAM fitting with monotonic constraints and confidence intervals
* Fast CI calculation method (72x faster than bootstrap)
