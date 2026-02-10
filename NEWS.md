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
