# Function Reference

## Overview

rePredRet provides functions for three main use cases:

1.  **Explore pre-built models** - Models between RepoRT systems from
    the published model repository
2.  **Build custom models** - Use your own RT data to build models
3.  **Explore the database** - Search compounds and systems

------------------------------------------------------------------------

## Exploring Pre-built Models

These functions fetch files on demand from the
[rePredRet-models](https://github.com/stanstrup-metabolomics/rePredRet-models)
repository and cache them locally. Use `refresh = TRUE` to bypass the
cache.

### `rePredRet_systems()`

Lists all chromatographic systems that have pre-built models, with their
RepoRT study information.

``` r
systems <- rePredRet_systems()
head(systems[, c("id", "name", "method.type")])
#   id    name      method.type
# 1 0001  FEM_short RP
# 2 0002  FEM_long  RP
```

### `rePredRet_models()`

Lists all available models with quality metrics.

``` r
models <- rePredRet_models()
head(models)
#   sys1_id sys2_id n_compounds median_ci_width median_error
# 1 0001    0002             67           15.99        0.336
```

### `rePredRet_model()`

Gets a single model with its prediction grid, calibration data and
statistics.

``` r
model <- rePredRet_model("0001", "0002")
head(model$predictions)   # x, pred, lower, upper
head(model$calibration)   # x, y, compound
model$stats
```

### `plot_model()`

Plots a pre-built model showing the calibration curve and confidence
bands.

``` r
# Static ggplot
plot_model("0001", "0002")

# Interactive plotly
plot_model("0001", "0002", interactive = TRUE)

# Without confidence bands
plot_model("0001", "0002", show_ci = FALSE)

# Without calibration points
plot_model("0001", "0002", show_calibration = FALSE)
```

------------------------------------------------------------------------

## Building Custom Models

Use your own RT data to build models and generate predictions.

### `build_user_models()`

Builds calibration models from RepoRT systems to your system.

``` r
models <- build_user_models(
  inchis = my_inchis,           # Character vector of InChI strings
  rts = my_rts,                 # Numeric vector of retention times
  system_type = "RP",           # "RP" or "HILIC"
  min_compounds = 10,           # Minimum compounds for model building
  method = "adaptive_pi",       # "adaptive_pi" (default), "fast_ci" or "bootstrap"
  verbose = TRUE
)
```

**Returns:** `UserRTModels` S4 object with:

- `models`: Named list of model objects
- `diagnostics`: Per-model quality metrics
- `input_data`: Your original input
- `metadata`: Run information

**Methods:**

``` r
show(models)                    # Print summary
summary(models)                 # Detailed statistics
plot(models)                    # Diagnostic plots
plot(models, type = "calibration")  # Calibration curves
names(models)                   # List system IDs
models[["0001"]]               # Extract specific model
saveRDS(models, "models.rds")  # Save for later
```

### `predict_rt_user()`

Generates predictions using built models.

``` r
predictions <- predict_rt_user(models)

# Predict specific compounds only
predictions <- predict_rt_user(
  models,
  compounds = c("InChI=1S/C8H10N4O2/...")
)
```

**Returns:** `UserRTPredictions` S4 object with:

- `predictions`: Data frame of predictions (one per compound, best CI
  selected)
- `source_models`: System IDs used
- `metadata`: Run information

**Methods:**

``` r
show(predictions)               # Print summary
preds_df <- as.data.frame(predictions)  # Extract data frame
plot(predictions)               # CI width distribution
plot(predictions, type = "by_system")   # Predictions per system
```

### `compare_systems()`

Ranks repository systems by similarity to your system.

``` r
# Rank by R-squared (default)
compare_systems(models)

# Rank by RMSE
compare_systems(models, metric = "rmse")

# Top 5 most similar systems
compare_systems(models, metric = "r_squared", top_n = 5)
```

------------------------------------------------------------------------

## Exploring the Database

### `search_compound()`

Searches for a compound by InChI across all systems. Stereochemistry is
automatically removed for matching.

``` r
# Search by InChI (caffeine)
search_compound("InChI=1S/C8H10N4O2/c1-10-4-9-6-5(10)7(13)12(3)8(14)11(6)2/h4H,1-3H3")

# Partial InChI also works
search_compound("InChI=1S/C8H10N4O2")

# Filter to RP systems only
search_compound("InChI=1S/C15H14O6", system_type = "RP")
```

**Returns:** Data frame with columns:

- `system_id`, `system_name`, `system_type`
- `compound_name`, `rt`, `inchi`

### `get_system_info()`

Gets detailed information about a chromatographic system.

``` r
info <- get_system_info("0001")
info$name           # System name
info$type           # RP or HILIC
info$n_compounds    # Number of compounds
info$metadata       # Column, gradient, etc.
info$gradient       # Gradient profile
```

------------------------------------------------------------------------

## Data Loading (Advanced)

### `download_report()`

Downloads the RepoRT repository from GitHub.

``` r
# Download latest
path <- download_report()

# Download specific release
path <- download_report(release = "v1.0.0")
```

### `load_report_data()`

Loads datasets from a local RepoRT directory.

``` r
data <- load_report_data(path, method_types = c("RP", "HILIC"))
data$studies    # Study metadata
data$datasets   # List of datasets with RT data
```

### `get_common_compounds()`

Finds compounds present in two datasets.

``` r
common <- get_common_compounds(dataset1, dataset2)
# Returns matrix: [rt_source, rt_target]
```

------------------------------------------------------------------------

## Summary Table

| Function | Purpose |
|----|----|
| [`rePredRet_systems()`](https://stanstrup-metabolomics.github.io/rePredRet/reference/rePredRet_systems.md) | List systems with pre-built models |
| [`rePredRet_models()`](https://stanstrup-metabolomics.github.io/rePredRet/reference/rePredRet_models.md) | List pre-built models |
| [`rePredRet_model()`](https://stanstrup-metabolomics.github.io/rePredRet/reference/rePredRet_model.md) | Get a pre-built model |
| [`plot_model()`](https://stanstrup-metabolomics.github.io/rePredRet/reference/plot_model.md) | Plot a pre-built model |
| [`rePredRet_plot_network()`](https://stanstrup-metabolomics.github.io/rePredRet/reference/rePredRet_plot_network.md) | Plot the network of pre-built models |
| [`build_user_models()`](https://stanstrup-metabolomics.github.io/rePredRet/reference/build_user_models.md) | Build models from your data |
| [`predict_rt_user()`](https://stanstrup-metabolomics.github.io/rePredRet/reference/predict_rt_user.md) | Generate predictions |
| [`compare_systems()`](https://stanstrup-metabolomics.github.io/rePredRet/reference/compare_systems.md) | Find similar systems |
| [`search_compound()`](https://stanstrup-metabolomics.github.io/rePredRet/reference/search_compound.md) | Search for compounds by InChI |
| [`get_system_info()`](https://stanstrup-metabolomics.github.io/rePredRet/reference/get_system_info.md) | Get system details |
