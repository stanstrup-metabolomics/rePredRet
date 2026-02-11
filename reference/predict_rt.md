# Predict RT for a Single Compound in a Target System

Predicts the retention time of a compound in a target chromatographic
system using available models.

## Usage

``` r
predict_rt(
  compound_id,
  target_system_id,
  models,
  datasets,
  ci_width_limit = .prediction_defaults$ci_width_limit,
  ci_width_limit_rel = .prediction_defaults$ci_width_limit_rel,
  density_bw_mult = .prediction_defaults$density_bw_mult,
  density_limit = .prediction_defaults$density_limit,
  id_column = "inchi.std"
)
```

## Arguments

- compound_id:

  Compound identifier (InChI or InChIKey).

- target_system_id:

  Target system identifier.

- models:

  Named list of models from
  [`build_all_models()`](https://stanstrup.github.io/rePredRet/reference/build_all_models.md).

- datasets:

  Named list of datasets from
  [`load_report_data()`](https://stanstrup.github.io/rePredRet/reference/load_report_data.md).

- ci_width_limit:

  Maximum CI width in minutes (default 2.0).

- ci_width_limit_rel:

  Maximum relative CI width (default 0.2 = 20%).

- density_bw_mult:

  Bandwidth multiplier for density check (default 0.03).

- density_limit:

  Minimum density threshold (default 0.01).

- id_column:

  Column name for compound matching (default "inchi.std").

## Value

A list with prediction details including provenance, or NULL if no valid
prediction can be made.

## Details

The function:

1.  Finds all models that predict TO the target system

2.  Checks if the compound exists in any source system

3.  Makes predictions and applies quality filters:

    - CI width must be \< ci_width_limit OR \< ci_width_limit_rel

    - Source RT must be in a dense region of calibration data

4.  Returns the best prediction (narrowest relative CI)
