# rePredRet: Retention Time Prediction by Direct System Mapping

rePredRet predicts retention times (RTs) by directly mapping between
chromatographic systems using monotonically constrained GAMs. This is a
reimplementation of PredRet that uses RepoRT data as input and publishes
pre-built models to GitHub.

## User Functions

For predicting RTs from your own data:

- [`build_user_models`](https://stanstrup-metabolomics.github.io/rePredRet/reference/build_user_models.md):
  Build models from your RT data to RepoRT systems

- [`predict_rt_user`](https://stanstrup-metabolomics.github.io/rePredRet/reference/predict_rt_user.md):
  Predict RTs with those models

For exploring the pre-built models (fetched on demand from
rePredRet-models):

- [`rePredRet_systems`](https://stanstrup-metabolomics.github.io/rePredRet/reference/rePredRet_systems.md):
  List available systems

- [`rePredRet_models`](https://stanstrup-metabolomics.github.io/rePredRet/reference/rePredRet_models.md):
  List available models

- [`rePredRet_model`](https://stanstrup-metabolomics.github.io/rePredRet/reference/rePredRet_model.md):
  Get specific model details

- [`plot_model`](https://stanstrup-metabolomics.github.io/rePredRet/reference/plot_model.md):
  Visualize model fit

- [`rePredRet_plot_network`](https://stanstrup-metabolomics.github.io/rePredRet/reference/rePredRet_plot_network.md):
  Visualize system network

## Developer Functions

For building your own models:

- [`download_report`](https://stanstrup-metabolomics.github.io/rePredRet/reference/download_report.md):
  Download RepoRT repository

- [`load_report_data`](https://stanstrup-metabolomics.github.io/rePredRet/reference/load_report_data.md):
  Load datasets from RepoRT

- [`get_common_compounds`](https://stanstrup-metabolomics.github.io/rePredRet/reference/get_common_compounds.md):
  Find common compounds

- [`build_model`](https://stanstrup-metabolomics.github.io/rePredRet/reference/build_model.md):
  Build a single model

- [`build_all_models`](https://stanstrup-metabolomics.github.io/rePredRet/reference/build_all_models.md):
  Build all pairwise models

- [`predict_rt`](https://stanstrup-metabolomics.github.io/rePredRet/reference/predict_rt.md):
  Predict RT for a compound

- [`generate_all_predictions`](https://stanstrup-metabolomics.github.io/rePredRet/reference/generate_all_predictions.md):
  Generate all predictions

- [`export_models`](https://stanstrup-metabolomics.github.io/rePredRet/reference/export_models.md):
  Export models to files

- [`export_predictions`](https://stanstrup-metabolomics.github.io/rePredRet/reference/export_predictions.md):
  Export predictions to files

## Scientific Background

The core assumption is that elution order is largely conserved between
similar chromatographic systems. If compound A elutes before compound B
in system 1, it will likely also elute before B in a similar system 2.

The algorithm fits monotonically constrained GAMs using:

1.  Thin plate regression splines for flexibility

2.  Monotonic constraints to preserve elution order

3.  Locally adaptive prediction intervals calibrated by cross-validated
    conformal prediction (`method = "adaptive_pi"`, default), faster
    intervals of nearly constant width (`method = "fast_ci"`) or
    bootstrap resampling

4.  Sigmoid weighting to downweight outliers

## References

Stanstrup, J. et al. (2015). PredRet: Prediction of Retention Time by
Direct Mapping between Multiple Chromatographic Systems. *Analytical
Chemistry*, 87(18), 9421-9428.

Kretschmer, F. et al. (2024). RepoRT: a comprehensive repository for
small molecule retention times. *Nature Methods*, 21, 153-155.

## See also

- RepoRT: <https://github.com/michaelwitting/RepoRT>

- rePredRet-models:
  <https://github.com/stanstrup-metabolomics/rePredRet-models>

- Original PredRet: <http://predret.org>

## Author

**Maintainer**: Jan Stanstrup <stanstrup@gmail.com>
([ORCID](https://orcid.org/0000-0003-0541-7369))

Authors:

- Jan Stanstrup <stanstrup@gmail.com>
  ([ORCID](https://orcid.org/0000-0003-0541-7369))

Other contributors:

- Barbora Hamidová (Idea of the adaptive prediction intervals; found the
  problems with the earlier interval methods) \[contributor\]
