# rePredRet (development version)

## Prediction intervals

* New default interval method `method = "adaptive_pi"` (`gam_adaptive_pi()`) in
  `build_model()`, `build_all_models()`, `build_all_models_furrr()` and
  `build_user_models()`. The interval width follows the local scatter of the
  calibration compounds (a Gamma GAM of cross-validated absolute residuals,
  weighted to ignore outliers) and is calibrated by cross-validated conformal
  prediction. On held-out compounds from 135 RepoRT system pairs it covered
  94.6% (95% nominal), equally in low- and high-scatter regions (94.8% /
  94.5%). It takes 60-90 ms per model. Thanks to Barbora Hamidová for the
  idea and for finding the problems with the other interval methods (MSc
  thesis, 2026).
* `adaptive_pi` and `bootstrap`: the local scale is no longer extrapolated
  across gaps between compounds or past down-weighted outliers, where the
  log-link smooth could grow without bound. In the first full build 1% of
  the models had a median PI width above 79 min, up to 1e36 min. The scale
  is now interpolated between the compounds and capped at the largest
  residual (method versions bumped, so the next pipeline run rebuilds them).
* `adaptive_pi` and `bootstrap`: the conformal multiplier q is now a weighted
  quantile in which each compound counts by its sigmoid outlier weight
  (`weighted_conformal_quantile()`). Before, with fewer than about 40
  compounds a single outlier (a wrong annotation the fit already ignores)
  set q and widened the whole band. The interval now aims to cover 95% of
  the well-annotated compounds (weight >= 0.5). In cross-validation on 400
  random published models it covered 94.9% of them (97.3% before) and 35%
  of the outliers; the median width fell from 3.6 to 2.1 min, and the
  models whose median width exceeds the target RT range from 71 to 17.
* `method = "fast_ci"` (`gam_fast_ci()`) is kept as the fast option (15-20 ms per
  model), with fixes: the residual variance comes from leave-one-out residuals
  instead of in-sample residuals (coverage 93.9% -> 95.5%; with 10-15
  compounds 83% -> 90%), and the bounds are made monotone by isotonic
  regression instead of a running maximum, which narrowed the lower side. Its
  width is still nearly constant (97.4% / 93.4% coverage in low- / high-scatter
  regions).
* `method = "bootstrap"` was a confidence interval of the curve, not a
  prediction interval (59% coverage; 53% with 200+ compounds). It now uses the
  out-of-bag predictions as out-of-sample residuals with the same local scale
  and conformal calibration as `adaptive_pi` (94.8% coverage).
  `boot2ci_PI()` has a new signature and the unused `boot2ci()` was removed.
* With 10-15 compounds no method reaches the nominal coverage (86-90%); this is
  now documented.
* All methods truncate the lower bound at 0.
* Model statistics (`model.json` `stats` and `model$stats`) gain `r_squared`,
  `coverage` (fraction of calibration compounds inside the PI) and
  `rel_mean_ci_width` (mean PI width / target RT range).
* Removed the unused `gam_prediction_intervals()`,
  `gam_prediction_intervals_parallel()`, `fit_gam_mono_once()` and
  `gam_fast_pi()`.
* `build_all_models_furrr()` and `build_all_models()` rebuild cached models that
  were built with another interval method, method version
  (`pi_method_version()`) or alpha, so the next pipeline run rebuilds all
  models.
* Small source-RT jitter (fewer than 4 distinct RTs) uses a fixed seed, so the
  same data always give the same model.
* `build_model()` returns status `"fit_failed"` instead of a model with
  non-finite predictions or intervals.
* The monotone GAM fit is shared by all interval methods (`fit_gam_mono()`). Its
  sigmoid outlier weights (0.1, 30) were set by eye for the original PredRet
  and have never been optimized; this is now documented.
* New article "Prediction interval benchmark", and the benchmark script
  `inst/benchmark/pi_benchmark.R`.
* `withr` moved from Suggests to Imports; `MASS` is no longer needed.

## Breaking changes

* `rePredRet_plot_network()` draws a different plot (see Website), and its
  `show_labels` now defaults to FALSE.
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

* `rePredRet_plot_network()` was a circle that is unreadable with hundreds of
  systems. It is now a force-directed network in which the similarity of two
  systems is the Jaccard index of their compounds. Systems are placed by the
  links with similarity of at least `layout_similarity` (default 0.4), and all
  links are drawn with width and opacity growing with similarity, so weak
  links fade into the background. Point size is the number of compounds in
  the system, shape the method type and colour groups of systems that share
  many compounds. `interactive = TRUE` gives a plotly version with the system
  names on hover and a slider for the minimum similarity of the links drawn;
  the Network article uses it. `show_labels` now defaults to FALSE. It needs
  the igraph package (Suggests).
* Six groups of systems are named (A-F, e.g. "Plant and food natural
  products", "Endogenous metabolite standards, Folberth et al. 2020"), from
  the studies that contributed their systems and the chemical classes of
  their compounds. Their members are listed in `inst/extdata/system_groups.csv`.
  Since the groups are recomputed from the current data, a group is named only
  when it matches a named group by mutual majority of their systems. The
  network and matrix plots show the letters and name them in the caption; the
  interactive network names them in its legend and on hover.
* New `rePredRet_plot_system_matrix()`: all system pairs as a similarity
  matrix ordered by group. It replaces the heatmap ordered by number of models
  in the Network article.
* Model viewer: systems are shown with their RepoRT name and method type, the
  system boxes can be searched by typing any part of the id or name, and the
  target list only offers systems with a model from the chosen source.
* Model statistics article: log axes are labelled with plain numbers
  (0.03, 0.1, 0.3, 1, ...) instead of scientific notation.
* The Overview, Model Statistics and Network articles are computed from the published
  model index when the site is built, instead of from statistics files committed in
  January (a 642-model test run). The site is rebuilt after each model update.
* Removed the Predictions article and the committed prediction and statistics files;
  pre-computed predictions are no longer produced.
* Removed the old standalone Quarto site in `website/`; the pkgdown site replaces it.

## Other

* `build_all_models()` runs its batches with furrr instead of a parallel cluster
  per batch; purrr is now imported and parallel no longer is.
* `rePredRet_plot_network()` gains `show_labels`.
* `inst/scripts/run_pipeline.R` writes to `output/` by default (was `website/data/`).
* Runner setup (`github_runner.yml`, `RUNNER_SETUP.md`): host paths are set in `.env`
  and must be on local disks; adds the `/opt/repredret-cache` mount.

## Pipeline

* The pipeline now also builds models between HILIC systems (about 3,000 from
  RepoRT's 71 HILIC systems; it was reversed-phase only).
  `build_all_models_furrr()` and `find_model_pairs()` only pair systems of the
  same method type (`method.type` in the studies table), so no RP-HILIC models
  are built.
* `build_all_models_furrr()` writes `systems.csv` next to `model_index.csv`:
  name, method type, number of compounds and number of models per system. The
  pipeline publishes it, and `rePredRet_systems()` gains `n_compounds` and
  `n_models` columns from it.
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
* Models whose build settings are unknown (no build cache entry, e.g. from an older
  pipeline) are rebuilt when an interval method is given, instead of being recorded
  as built with the current method.
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
