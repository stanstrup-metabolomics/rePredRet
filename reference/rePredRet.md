# rePredRet: Retention Time Prediction by Direct System Mapping

rePredRet predicts retention times (RTs) by directly mapping between
chromatographic systems using monotonically constrained GAMs. This is a
reimplementation of PredRet that uses RepoRT data as input and publishes
predictions to GitHub.

## User Functions

For accessing pre-computed predictions from rePredRet-data:

- [`rePredRet_download`](https://stanstrup.github.io/rePredRet/reference/rePredRet_download.md):
  Download latest predictions

- [`rePredRet_systems`](https://stanstrup.github.io/rePredRet/reference/rePredRet_systems.md):
  List available systems

- [`rePredRet_models`](https://stanstrup.github.io/rePredRet/reference/rePredRet_models.md):
  List available models

- [`rePredRet_model`](https://stanstrup.github.io/rePredRet/reference/rePredRet_model.md):
  Get specific model details

- [`rePredRet_predict`](https://stanstrup.github.io/rePredRet/reference/rePredRet_predict.md):
  Get predictions for a compound

- [`rePredRet_system_predictions`](https://stanstrup.github.io/rePredRet/reference/rePredRet_system_predictions.md):
  Get system predictions

- [`rePredRet_plot_model`](https://stanstrup.github.io/rePredRet/reference/rePredRet_plot_model.md):
  Visualize model fit

- [`rePredRet_plot_errors`](https://stanstrup.github.io/rePredRet/reference/rePredRet_plot_errors.md):
  Plot error distribution

- [`rePredRet_plot_network`](https://stanstrup.github.io/rePredRet/reference/rePredRet_plot_network.md):
  Visualize system network

## Developer Functions

For building your own models:

- [`download_report`](https://stanstrup.github.io/rePredRet/reference/download_report.md):
  Download RepoRT repository

- [`load_report_data`](https://stanstrup.github.io/rePredRet/reference/load_report_data.md):
  Load datasets from RepoRT

- [`get_common_compounds`](https://stanstrup.github.io/rePredRet/reference/get_common_compounds.md):
  Find common compounds

- [`build_model`](https://stanstrup.github.io/rePredRet/reference/build_model.md):
  Build a single model

- [`build_all_models`](https://stanstrup.github.io/rePredRet/reference/build_all_models.md):
  Build all pairwise models

- [`predict_rt`](https://stanstrup.github.io/rePredRet/reference/predict_rt.md):
  Predict RT for a compound

- [`generate_all_predictions`](https://stanstrup.github.io/rePredRet/reference/generate_all_predictions.md):
  Generate all predictions

- [`export_models`](https://stanstrup.github.io/rePredRet/reference/export_models.md):
  Export models to files

- [`export_predictions`](https://stanstrup.github.io/rePredRet/reference/export_predictions.md):
  Export predictions to files

## Scientific Background

The core assumption is that elution order is largely conserved between
similar chromatographic systems. If compound A elutes before compound B
in system 1, it will likely also elute before B in a similar system 2.

The algorithm fits monotonically constrained GAMs using:

1.  Thin plate regression splines for flexibility

2.  Monotonic constraints to preserve elution order

3.  Bootstrap resampling (1000 iterations) for CI estimation

4.  Sigmoid weighting to downweight outliers

## References

Stanstrup, J. et al. (2015). PredRet: Prediction of Retention Time by
Direct Mapping between Multiple Chromatographic Systems. *Analytical
Chemistry*, 87(18), 9421-9428.

Kretschmer, F. et al. (2024). RepoRT: a comprehensive repository for
small molecule retention times. *Nature Methods*, 21, 153-155.

## See also

- RepoRT: <https://github.com/michaelwitting/RepoRT>

- rePredRet-data: <https://github.com/stanstrup/rePredRet-data>

- Original PredRet: <http://predret.org>

## Author

**Maintainer**: Jan Stanstrup <stanstrup@gmail.com>
([ORCID](https://orcid.org/0000-0003-0541-7369))
