# Build Calibration Models from Repository Systems to User's System

Builds calibration models from repository systems to the user's
chromatographic system using the user's retention time data as
calibration.

## Usage

``` r
build_user_models(
  inchis,
  rts,
  system_type = c("RP", "HILIC"),
  report_data = NULL,
  min_compounds = 10,
  method = c("fast_ci", "bootstrap"),
  alpha = 0.05,
  n_boot = 200,
  n_workers = parallel::detectCores(),
  verbose = TRUE
)
```

## Arguments

- inchis:

  Character vector of InChI identifiers for user's compounds.

- rts:

  Numeric vector of retention times (in minutes) corresponding to
  inchis.

- system_type:

  Character, either "RP" (reversed-phase) or "HILIC". Only models from
  systems of the same type will be built.

- report_data:

  Optional. Output from load_report_data(). If NULL, downloads and loads
  the latest RepoRT data.

- min_compounds:

  Minimum number of common compounds required to build a model (default
  10).

- method:

  Method for CI calculation: "fast_ci" (default, fast) or "bootstrap"
  (slower, more accurate).

- alpha:

  Significance level for CI (default 0.05 for 95% CI).

- n_boot:

  Number of bootstrap iterations if method = "bootstrap" (default 200).

- n_workers:

  Number of parallel workers (default = available cores). Currently not
  used (models built sequentially).

- verbose:

  Print progress messages (default TRUE).

## Value

A UserRTModels S4 object containing:

- models:

  Named list of fitted model objects (keyed by source system ID)

- diagnostics:

  Per-model quality metrics (R-squared, RMSE, CI widths)

- input_data:

  The user's original input data

- report_data:

  Cached RepoRT data for later predictions

- metadata:

  Run information and summary statistics

## Details

The function:

1.  Validates inputs and standardizes InChIs (removes stereochemistry)

2.  Loads RepoRT data and filters to matching system type

3.  For each compatible repository system:

    - Finds common compounds between repo system and user's data

    - Builds a monotonic GAM model (repo RT -\> user RT)

    - Calculates model diagnostics

4.  Returns all models in a UserRTModels object

Direction: FROM repository systems TO user's system (user's system is
TARGET)

## See also

[`predict_rt_user`](https://stanstrup.github.io/rePredRet/reference/predict_rt_user.md)
for generating predictions from models,
[`UserRTModels-class`](https://stanstrup.github.io/rePredRet/reference/UserRTModels-class.md)
for the return object class

## Examples

``` r
if (FALSE) { # \dontrun{
# User has measured RTs for some compounds
my_inchis <- c("InChI=1S/C8H10N4O2/...", "InChI=1S/C10H14N2/...")
my_rts <- c(3.2, 5.7)

# Build models from all RP systems
my_models <- build_user_models(
  inchis = my_inchis,
  rts = my_rts,
  system_type = "RP"
)

# View summary
show(my_models)

# Plot diagnostics
plot(my_models, type = "diagnostics")

# Save for later use
saveRDS(my_models, "my_rp_models.rds")

# Generate predictions
predictions <- predict_rt_user(my_models)
} # }
```
