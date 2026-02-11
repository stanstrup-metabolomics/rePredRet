# Package index

## Package

- [`rePredRet`](https://stanstrup.github.io/rePredRet/reference/rePredRet.md)
  [`rePredRet-package`](https://stanstrup.github.io/rePredRet/reference/rePredRet.md)
  : rePredRet: Retention Time Prediction by Direct System Mapping

## User Workflow

Build models from your own retention time data and generate predictions
for new systems. This is the primary interface for end users.

- [`build_user_models()`](https://stanstrup.github.io/rePredRet/reference/build_user_models.md)
  : Build Calibration Models from Repository Systems to User's System
- [`predict_rt_user()`](https://stanstrup.github.io/rePredRet/reference/predict_rt_user.md)
  : Generate RT Predictions from User Models
- [`show(`*`<UserRTModels>`*`)`](https://stanstrup.github.io/rePredRet/reference/UserRTModels-class.md)
  [`summary(`*`<UserRTModels>`*`)`](https://stanstrup.github.io/rePredRet/reference/UserRTModels-class.md)
  [`` `[[`( ``*`<UserRTModels>`*`,`*`<character>`*`)`](https://stanstrup.github.io/rePredRet/reference/UserRTModels-class.md)
  [`names(`*`<UserRTModels>`*`)`](https://stanstrup.github.io/rePredRet/reference/UserRTModels-class.md)
  [`length(`*`<UserRTModels>`*`)`](https://stanstrup.github.io/rePredRet/reference/UserRTModels-class.md)
  [`plot(`*`<UserRTModels>`*`,`*`<missing>`*`)`](https://stanstrup.github.io/rePredRet/reference/UserRTModels-class.md)
  : UserRTModels S4 Class
- [`show(`*`<UserRTPredictions>`*`)`](https://stanstrup.github.io/rePredRet/reference/UserRTPredictions-class.md)
  [`summary(`*`<UserRTPredictions>`*`)`](https://stanstrup.github.io/rePredRet/reference/UserRTPredictions-class.md)
  [`as.data.frame(`*`<UserRTPredictions>`*`)`](https://stanstrup.github.io/rePredRet/reference/UserRTPredictions-class.md)
  [`length(`*`<UserRTPredictions>`*`)`](https://stanstrup.github.io/rePredRet/reference/UserRTPredictions-class.md)
  [`` `[`( ``*`<UserRTPredictions>`*`,`*`<ANY>`*`,`*`<ANY>`*`,`*`<ANY>`*`)`](https://stanstrup.github.io/rePredRet/reference/UserRTPredictions-class.md)
  [`plot(`*`<UserRTPredictions>`*`,`*`<missing>`*`)`](https://stanstrup.github.io/rePredRet/reference/UserRTPredictions-class.md)
  : UserRTPredictions S4 Class
- [`validate_models()`](https://stanstrup.github.io/rePredRet/reference/validate_models.md)
  : Validate Models with Cross-Validation
- [`check_input_data()`](https://stanstrup.github.io/rePredRet/reference/check_input_data.md)
  : Check Input Data Quality
- [`export_predictions()`](https://stanstrup.github.io/rePredRet/reference/export_predictions.md)
  : Export Predictions to Files

## Exploration & Search

Search for compounds across systems, compare chromatographic systems,
and retrieve system metadata.

- [`search_compound()`](https://stanstrup.github.io/rePredRet/reference/search_compound.md)
  : Search for a Compound Across All Systems
- [`search_compounds()`](https://stanstrup.github.io/rePredRet/reference/search_compounds.md)
  : Batch Search for Multiple Compounds
- [`compare_systems()`](https://stanstrup.github.io/rePredRet/reference/compare_systems.md)
  : Compare User's System to Repository Systems
- [`get_system_info()`](https://stanstrup.github.io/rePredRet/reference/get_system_info.md)
  [`print(`*`<system_info>`*`)`](https://stanstrup.github.io/rePredRet/reference/get_system_info.md)
  : Get Detailed System Information

## Visualization

Plot models, prediction errors, system networks, and compound coverage.

- [`plot_model()`](https://stanstrup.github.io/rePredRet/reference/plot_model.md)
  : Create Interactive Model Plot
- [`plot_system_network()`](https://stanstrup.github.io/rePredRet/reference/plot_system_network.md)
  : Plot System Network
- [`plot_prediction_errors()`](https://stanstrup.github.io/rePredRet/reference/plot_prediction_errors.md)
  : Plot Prediction Errors
- [`plot_compound_coverage()`](https://stanstrup.github.io/rePredRet/reference/plot_compound_coverage.md)
  : Plot Compound Coverage Heatmap
- [`save_model_plot()`](https://stanstrup.github.io/rePredRet/reference/save_model_plot.md)
  : Save Model Plot as Self-Contained HTML

## Pre-built Data Access

Access pre-built models, predictions, and system information from the
latest rePredRet release on GitHub.

- [`rePredRet_download()`](https://stanstrup.github.io/rePredRet/reference/rePredRet_download.md)
  : Download rePredRet Data
- [`rePredRet_systems()`](https://stanstrup.github.io/rePredRet/reference/rePredRet_systems.md)
  : List Available Chromatographic Systems
- [`rePredRet_models()`](https://stanstrup.github.io/rePredRet/reference/rePredRet_models.md)
  : List Available Models
- [`rePredRet_model()`](https://stanstrup.github.io/rePredRet/reference/rePredRet_model.md)
  : Get Model Details
- [`rePredRet_predict()`](https://stanstrup.github.io/rePredRet/reference/rePredRet_predict.md)
  : Get Predictions for a Compound
- [`rePredRet_system_predictions()`](https://stanstrup.github.io/rePredRet/reference/rePredRet_system_predictions.md)
  : Get All Predictions for a System
- [`rePredRet_plot_model()`](https://stanstrup.github.io/rePredRet/reference/rePredRet_plot_model.md)
  : Plot Model Fit
- [`rePredRet_plot_errors()`](https://stanstrup.github.io/rePredRet/reference/rePredRet_plot_errors.md)
  : Plot Prediction Error Distribution
- [`rePredRet_plot_network()`](https://stanstrup.github.io/rePredRet/reference/rePredRet_plot_network.md)
  : Plot System Network

## Core Pipeline

Lower-level functions for the full model-building and prediction
pipeline. Most users should use the User Workflow or Pre-built Data
Access functions instead.

- [`build_model()`](https://stanstrup.github.io/rePredRet/reference/build_model.md)
  : Build a Single Model Between Two Systems
- [`build_all_models()`](https://stanstrup.github.io/rePredRet/reference/build_all_models.md)
  : Build All Pairwise Models
- [`build_all_models_furrr()`](https://stanstrup.github.io/rePredRet/reference/build_all_models_furrr.md)
  : Build All Models with furrr - Optimized with Caching and Scheduling
- [`download_report()`](https://stanstrup.github.io/rePredRet/reference/download_report.md)
  : Download RepoRT Repository
- [`load_report_data()`](https://stanstrup.github.io/rePredRet/reference/load_report_data.md)
  : Load RepoRT Data
- [`get_common_compounds()`](https://stanstrup.github.io/rePredRet/reference/get_common_compounds.md)
  : Get Common Compounds Between Two Systems
- [`predict_rt()`](https://stanstrup.github.io/rePredRet/reference/predict_rt.md)
  : Predict RT for a Single Compound in a Target System
- [`generate_all_predictions()`](https://stanstrup.github.io/rePredRet/reference/generate_all_predictions.md)
  : Generate All Predictions (Optimized Batch Processing)
- [`get_calibration_regression_data()`](https://stanstrup.github.io/rePredRet/reference/get_calibration_regression_data.md)
  : Aggregate Calibration Data for Regression Plot
- [`get_latest_report_release()`](https://stanstrup.github.io/rePredRet/reference/get_latest_report_release.md)
  : Get Latest RepoRT Release Tag

## Statistics & GAM Utilities

Functions for computing model statistics, prediction statistics, and GAM
confidence/prediction intervals.

- [`gam_fast_ci()`](https://stanstrup.github.io/rePredRet/reference/gam_fast_ci.md)
  : Fast Prediction Intervals using predict() + cummax
- [`gam_prediction_intervals()`](https://stanstrup.github.io/rePredRet/reference/gam_prediction_intervals.md)
  : GAM Prediction Intervals via Posterior Simulation
- [`calculate_database_stats()`](https://stanstrup.github.io/rePredRet/reference/calculate_database_stats.md)
  : Calculate Database Statistics
- [`calculate_model_stats()`](https://stanstrup.github.io/rePredRet/reference/calculate_model_stats.md)
  : Calculate Model Performance Statistics
- [`calculate_prediction_stats()`](https://stanstrup.github.io/rePredRet/reference/calculate_prediction_stats.md)
  : Calculate Prediction Statistics
- [`export_models()`](https://stanstrup.github.io/rePredRet/reference/export_models.md)
  : Export Models to Files
- [`export_statistics()`](https://stanstrup.github.io/rePredRet/reference/export_statistics.md)
  : Export Statistics to Files
