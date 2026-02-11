# Changelog

## rePredRet 0.2.0

### New Features

#### User Model Building

- [`build_user_models()`](https://stanstrup.github.io/rePredRet/reference/build_user_models.md):
  Build calibration models from your own RT data to all compatible
  RepoRT systems
- [`predict_rt_user()`](https://stanstrup.github.io/rePredRet/reference/predict_rt_user.md):
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

- [`plot_model()`](https://stanstrup.github.io/rePredRet/reference/plot_model.md):
  Plot pre-built models from the data repository with ggplot2/plotly
- [`search_compound()`](https://stanstrup.github.io/rePredRet/reference/search_compound.md):
  Search for compounds by InChI across all systems
- [`search_compounds()`](https://stanstrup.github.io/rePredRet/reference/search_compounds.md):
  Batch search for multiple InChIs
- [`compare_systems()`](https://stanstrup.github.io/rePredRet/reference/compare_systems.md):
  Rank repository systems by similarity to user’s data
- [`get_system_info()`](https://stanstrup.github.io/rePredRet/reference/get_system_info.md):
  Get detailed system metadata (column, gradient, etc.)
- [`export_predictions()`](https://stanstrup.github.io/rePredRet/reference/export_predictions.md):
  Export predictions to CSV or Excel format
- [`validate_models()`](https://stanstrup.github.io/rePredRet/reference/validate_models.md):
  Cross-validation for assessing model quality
- [`check_input_data()`](https://stanstrup.github.io/rePredRet/reference/check_input_data.md):
  Validate InChIs and RTs before model building

#### Visualizations

- [`plot_system_network()`](https://stanstrup.github.io/rePredRet/reference/plot_system_network.md):
  Visualize system connectivity
- [`plot_prediction_errors()`](https://stanstrup.github.io/rePredRet/reference/plot_prediction_errors.md):
  Error distribution plots
- [`plot_compound_coverage()`](https://stanstrup.github.io/rePredRet/reference/plot_compound_coverage.md):
  Heatmap of compound presence across systems

### Improvements

- All InChI matching now consistently strips stereochemistry layers (/t,
  /b, /m, /s)
- [`predict_rt_user()`](https://stanstrup.github.io/rePredRet/reference/predict_rt_user.md)
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
  [`build_model()`](https://stanstrup.github.io/rePredRet/reference/build_model.md),
  [`build_all_models()`](https://stanstrup.github.io/rePredRet/reference/build_all_models.md),
  [`build_all_models_furrr()`](https://stanstrup.github.io/rePredRet/reference/build_all_models_furrr.md)
- Prediction functions:
  [`predict_rt()`](https://stanstrup.github.io/rePredRet/reference/predict_rt.md),
  [`generate_all_predictions()`](https://stanstrup.github.io/rePredRet/reference/generate_all_predictions.md)
- Data loading:
  [`download_report()`](https://stanstrup.github.io/rePredRet/reference/download_report.md),
  [`load_report_data()`](https://stanstrup.github.io/rePredRet/reference/load_report_data.md),
  [`get_common_compounds()`](https://stanstrup.github.io/rePredRet/reference/get_common_compounds.md)
- User query functions:
  [`rePredRet_download()`](https://stanstrup.github.io/rePredRet/reference/rePredRet_download.md),
  [`rePredRet_predict()`](https://stanstrup.github.io/rePredRet/reference/rePredRet_predict.md),
  [`rePredRet_systems()`](https://stanstrup.github.io/rePredRet/reference/rePredRet_systems.md),
  [`rePredRet_models()`](https://stanstrup.github.io/rePredRet/reference/rePredRet_models.md)
- GAM fitting with monotonic constraints and confidence intervals
- Fast CI calculation method (72x faster than bootstrap)
