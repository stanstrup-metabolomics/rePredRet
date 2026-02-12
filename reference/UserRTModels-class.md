# UserRTModels S4 Class

S4 class for storing calibration models built from repository systems to
a user's chromatographic system.

## Usage

``` r
# S4 method for class 'UserRTModels'
show(object)

# S4 method for class 'UserRTModels'
summary(object, ...)

# S4 method for class 'UserRTModels,character'
x[[i]]

# S4 method for class 'UserRTModels'
names(x)

# S4 method for class 'UserRTModels'
length(x)

# S4 method for class 'UserRTModels,missing'
plot(
  x,
  y,
  type = c("diagnostics", "calibration", "r_squared", "ci_widths"),
  max_plots = 9,
  ...
)
```

## Arguments

- object:

  A UserRTModels object

- ...:

  Additional arguments

- x:

  A UserRTModels object

- i:

  System ID (character)

- y:

  Ignored

- type:

  Plot type: "diagnostics" (default), "calibration", "r_squared",
  "ci_widths"

- max_plots:

  Maximum number of calibration plots to show (default 9)

## Value

A list with summary statistics

The model object for that system

Character vector of system IDs

Integer

A ggplot object or list of ggplot objects

## Slots

- `models`:

  A named list of model objects (from build_model) keyed by
  source_system_id. Each model contains:

  status

  :   "success" or "insufficient_data"

  sys1_id

  :   Source system ID (repository system)

  sys2_id

  :   "USER" (target is user's system)

  n_points

  :   Number of calibration compounds

  newdata

  :   Prediction grid (1000 points)

  ci

  :   Matrix `[pred, lower, upper]` at newdata points

  calibration

  :   Matrix `[rt_source, rt_user]` for calibration compounds

  stats

  :   Error statistics at calibration points

- `diagnostics`:

  A data.frame with per-model diagnostic metrics:

  source_system

  :   Repository system ID

  n_calibration

  :   Number of common compounds used

  r_squared

  :   Model R-squared

  rmse

  :   Root mean square error

  mean_error

  :   Mean absolute error

  median_error

  :   Median absolute error

  q95_error

  :   95th percentile error

  mean_ci_width

  :   Mean CI width at calibration points

  median_ci_width

  :   Median CI width

  rt_range_source

  :   RT range in source system

  rt_range_user

  :   RT range in user's system

- `input_data`:

  A data.frame with the user's original input:

  inchi

  :   Original InChI strings

  inchi_clean

  :   Standardized InChI (no stereochemistry)

  rt

  :   Retention times provided

- `report_data`:

  A list containing cached RepoRT data for predictions

- `metadata`:

  A list containing:

  system_type

  :   RP or HILIC

  n_input_compounds

  :   Number of user compounds

  n_models_built

  :   Number of successful models

  n_models_failed

  :   Number of failed models

  created_at

  :   Timestamp

  repredret_version

  :   Package version

  method

  :   CI method used
