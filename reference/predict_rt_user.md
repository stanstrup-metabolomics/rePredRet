# Generate RT Predictions from User Models

Generates retention time predictions using models built with
[`build_user_models`](https://stanstrup.github.io/rePredRet/reference/build_user_models.md).

## Usage

``` r
predict_rt_user(models, compounds = NULL, verbose = TRUE)
```

## Arguments

- models:

  A UserRTModels object from
  [`build_user_models`](https://stanstrup.github.io/rePredRet/reference/build_user_models.md).

- compounds:

  Optional character vector of InChI identifiers to predict. If NULL
  (default), predicts all compounds from all source systems.

- verbose:

  Print progress messages (default TRUE).

## Value

A UserRTPredictions S4 object containing:

- predictions:

  Data frame of all predictions with CIs and model quality

- source_models:

  Character vector of model system IDs used

- metadata:

  Run information and summary statistics

## Details

The function:

1.  Takes each model from the UserRTModels object

2.  For each source system, retrieves all measured compounds

3.  Generates predictions by interpolating from the model's CI grid

4.  Selects the best prediction per compound (narrowest CI width)

When a compound is measured in multiple source systems, only the
prediction with the smallest confidence interval width is returned.

## See also

[`build_user_models`](https://stanstrup.github.io/rePredRet/reference/build_user_models.md)
for building models,
[`UserRTPredictions-class`](https://stanstrup.github.io/rePredRet/reference/UserRTPredictions-class.md)
for the return object class

## Examples

``` r
if (FALSE) { # \dontrun{
# First build models
my_models <- build_user_models(
  inchis = my_inchis,
  rts = my_rts,
  system_type = "RP"
)

# Generate all predictions
predictions <- predict_rt_user(my_models)

# View summary
show(predictions)

# Extract as data.frame
preds_df <- as.data.frame(predictions)

# Filter to high-quality predictions
good_preds <- preds_df[preds_df$ci_width < 1.0 & preds_df$model_r_squared > 0.95, ]

# Or predict specific compounds only
caffeine_preds <- predict_rt_user(
  my_models,
  compounds = c("InChI=1S/C8H10N4O2/c1-10-4-9-6-5(10)7(13)12(3)8(14)11(6)2/h4H,1-3H3")
)
} # }
```
