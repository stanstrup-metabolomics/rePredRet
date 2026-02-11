# UserRTPredictions S4 Class

S4 class for storing RT predictions generated from UserRTModels.

## Usage

``` r
# S4 method for class 'UserRTPredictions'
show(object)

# S4 method for class 'UserRTPredictions'
summary(object, ...)

# S4 method for class 'UserRTPredictions'
as.data.frame(x, row.names = NULL, optional = FALSE, ...)

# S4 method for class 'UserRTPredictions'
length(x)

# S4 method for class 'UserRTPredictions,ANY,ANY,ANY'
x[i, j, ..., drop = TRUE]

# S4 method for class 'UserRTPredictions,missing'
plot(x, y, type = c("ci_widths", "rt_distribution", "by_system"), ...)
```

## Arguments

- object:

  A UserRTPredictions object

- ...:

  Additional arguments

- x:

  A UserRTPredictions object

- row.names:

  Ignored

- optional:

  Ignored

- i:

  Row indices or logical vector

- j:

  Column names (optional)

- drop:

  Passed to data.frame subsetting

- y:

  Ignored

- type:

  Plot type: "ci_widths" (default), "rt_distribution", "by_system"

## Value

A list with summary statistics

A data.frame of predictions

Integer

Subset of predictions data.frame

A ggplot object

## Slots

- `predictions`:

  A data.frame containing all predictions with columns:

  compound_id

  :   InChI identifier of the predicted compound

  compound_name

  :   Human-readable name (if available)

  source_system

  :   Repository system ID where compound was measured

  source_rt

  :   Original RT in source system (minutes)

  predicted_rt

  :   Predicted RT in user's system (minutes)

  ci_lower

  :   Lower bound of confidence interval

  ci_upper

  :   Upper bound of confidence interval

  ci_width

  :   Width of confidence interval

  model_r_squared

  :   R-squared of the model used

  model_rmse

  :   RMSE of the model used

  n_calibration

  :   Number of calibration compounds in model

- `source_models`:

  Character vector of system IDs for models used

- `metadata`:

  A list containing:

  n_predictions

  :   Total number of predictions

  n_unique_compounds

  :   Number of unique compounds predicted

  n_models_used

  :   Number of models used

  created_at

  :   Timestamp
