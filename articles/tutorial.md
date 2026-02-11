# Tutorial: Predict RT for Your Own Data

## Overview

This tutorial shows how to use rePredRet to predict retention times for
compounds you haven’t measured, using your own chromatographic data as
input. Instead of uploading your data to RepoRT, you can build
calibration models locally and get predictions from all compatible
systems in the database.

**Workflow:**

1.  Provide your measured retention times (InChI + RT pairs)
2.  Build models from RepoRT systems to your system
3.  Generate predictions for all compounds in those systems

## Prerequisites

``` r
# Install the package (if not already installed)
devtools::install_github("stanstrup/rePredRet")
```

``` r
library(rePredRet)
```

## Example: Using Wine Polyphenol Data

For this tutorial, we’ll use retention time data from a reversed-phase
(RP) system that measured wine polyphenols. This is based on dataset
0001 (FEM_short) from RepoRT.

### Step 1: Prepare Your Data

Your data should be two vectors of equal length:

- `inchis`: InChI identifiers for your compounds
- `rts`: Retention times in minutes

``` r
# Example: Wine polyphenol retention times from an RP system
# These are real measurements from a C18 column with water/acetonitrile gradient

my_inchis <- c(
  # Epicatechin
  "InChI=1S/C15H14O6/c16-8-4-11(18)9-6-13(20)15(21-14(9)5-8)7-1-2-10(17)12(19)3-7/h1-5,13,15-20H,6H2",
  # Gallocatechin
  "InChI=1S/C15H14O7/c16-7-3-9(17)8-5-12(20)15(22-13(8)4-7)6-1-10(18)14(21)11(19)2-6/h1-4,12,15-21H,5H2",
  # 4-hydroxybenzoic acid
  "InChI=1S/C7H6O3/c8-6-3-1-5(2-4-6)7(9)10/h1-4,8H,(H,9,10)",
  # Caftaric acid
  "InChI=1S/C13H12O9/c14-7-3-1-6(5-8(7)15)2-4-9(16)22-11(13(20)21)10(17)12(18)19/h1-5,10-11,14-15,17H,(H,18,19)(H,20,21)",
  # Esculetin
  "InChI=1S/C9H6O4/c10-6-3-5-1-2-9(12)13-8(5)4-7(6)11/h1-4,10-11H",
  # Astringin
  "InChI=1S/C20H22O9/c21-9-16-17(25)18(26)19(27)20(29-16)28-13-6-11(5-12(22)8-13)2-1-10-3-4-14(23)15(24)7-10/h1-8,16-27H,9H2",
  # Myricetin
  "InChI=1S/C15H10O8/c16-6-3-7(17)11-10(4-6)23-15(14(22)13(11)21)5-1-8(18)12(20)9(19)2-5/h1-4,16-20,22H",
  # Naringenin
  "InChI=1S/C15H12O5/c16-9-3-1-8(2-4-9)13-7-12(19)15-11(18)5-10(17)6-14(15)20-13/h1-6,13,16-18H,7H2",
  # Quercetin 3-glucoside
  "InChI=1S/C21H20O12/c22-6-12(27)19-16(29)17(30)21(32-19)33-20-15(28)14-11(26)4-8(23)5-13(14)31-18(20)7-1-2-9(24)10(25)3-7/h1-5,12,16-17,19,21-27,29-30H,6H2",
  # Quercetin 3-rutinoside (Rutin)
  "InChI=1S/C27H30O16/c1-8-17(32)20(35)22(37)26(40-8)39-7-15-18(33)21(36)23(38)27(42-15)43-25-19(34)16-13(31)5-10(28)6-14(16)41-24(25)9-2-3-11(29)12(30)4-9/h2-6,8,15,17-18,20-23,26-33,35-38H,7H2,1H3",
  # Resveratrol
  "InChI=1S/C14H12O3/c15-12-5-3-10(4-6-12)1-2-11-7-13(16)9-14(17)8-11/h1-9,15-17H",
  # Piceatannol
  "InChI=1S/C14H12O4/c15-11-5-10(6-12(16)8-11)2-1-9-3-4-13(17)14(18)7-9/h1-8,15-18H",
  # Kaempferol 3-glucoside
  "InChI=1S/C21H20O11/c22-7-13-15(26)17(28)18(29)21(31-13)32-20-16(27)14-11(25)5-10(24)6-12(14)30-19(20)8-1-3-9(23)4-2-8/h1-6,13,15,17-18,21-26,28-29H,7H2"
)

my_rts <- c(
  12.07,  # Epicatechin
  8.675,  # Gallocatechin
  10.6,   # 4-hydroxybenzoic acid
  8.53,   # Caftaric acid
  12.19,  # Esculetin
  14.2,   # Astringin
  20.46,  # Myricetin
  21.04,  # Naringenin
  20.0,   # Quercetin 3-glucoside
  20.0,   # Quercetin 3-rutinoside
  20.37,  # Resveratrol
  17.88,  # Piceatannol
  20.285  # Kaempferol 3-glucoside
)

# Verify we have matching lengths
length(my_inchis) == length(my_rts)
```

### Step 2: Build Calibration Models

Use
[`build_user_models()`](https://stanstrup.github.io/rePredRet/reference/build_user_models.md)
to build models from all compatible RepoRT systems to your system.
Specify whether your system is RP (reversed-phase) or HILIC.

``` r
# Build models from all RP systems in RepoRT
my_models <- build_user_models(
  inchis = my_inchis,
  rts = my_rts,
  system_type = "RP",
  min_compounds = 10,  # Minimum compounds needed for a model
  method = "fast_ci",  # Use fast CI calculation
  verbose = TRUE
)
```

    ================================================
      build_user_models
      System type: RP
      User compounds: 13
      Method: fast_ci
    ================================================

    Step 1/3: Standardizing InChIs...
      Unique compounds after deduplication: 13

    Step 2/3: Loading RepoRT data...
      Found 45 RP systems in RepoRT

    Step 3/3: Building calibration models...
      [10/45] Building model from 0010 (12 compounds)...
      [20/45] Building model from 0019 (11 compounds)...
      [30/45] Building model from 0029 (10 compounds)...
      [45/45] Building model from 0045 (11 compounds)...

      Models built: 8
      Failed (insufficient compounds or error): 37

    ================================================
      Complete!
      Time elapsed: 12.3 seconds
      Models built: 8
    ================================================

### Step 3: Inspect the Models

The `UserRTModels` object contains all built models and their
diagnostics.

``` r
# View summary
my_models
```

    UserRTModels Object
    ===================

    System Type: RP
    Created: 2024-01-15 14:23:15

    Input:
      User compounds: 13
      Unique (after dedup): 13

    Models:
      Built successfully: 8
      Failed (insufficient data): 37

    Model Quality (median across models):
      R-squared: 0.987
      RMSE: 0.42 min
      Median CI width: 0.85 min

    Calibration compounds per model:
      Min: 10
      Median: 11
      Max: 13

``` r
# Plot model diagnostics
plot(my_models, type = "diagnostics")
```

This shows model quality (R²) versus number of calibration compounds,
colored by CI width.

``` r
# View calibration curves for top models
plot(my_models, type = "calibration")
```

### Step 4: Generate Predictions

Use
[`predict_rt_user()`](https://stanstrup.github.io/rePredRet/reference/predict_rt_user.md)
to generate predictions for all compounds in the source systems.

``` r
# Generate all predictions
predictions <- predict_rt_user(my_models, verbose = TRUE)
```

    ================================================
      predict_rt_user
      Models available: 8
    ================================================

    Generating predictions...
      [8/8] Processed 0045

    ================================================
      Complete!
      Time elapsed: 2.1 seconds
      Total predictions: 1847
      Unique compounds: 423
    ================================================

### Step 5: Explore Predictions

``` r
# View summary
predictions
```

    UserRTPredictions Object
    ========================

    Created: 2024-01-15 14:25:32

    Predictions:
      Total predictions: 1847
      Unique compounds: 423
      Models used: 8

    Predicted RT range: 2.15 - 28.43 min
    CI width (median): 0.92 min

    Preview (first 5 predictions):
                                     compound_id source_system predicted_rt ci_width
     InChI=1S/C15H14O6/c16-8-4-11(18)9-6-13...         0002        11.82     0.65
     InChI=1S/C22H18O10/c23-11-6-14(25)12-8-...         0002        15.43     0.78
     InChI=1S/C15H14O7/c16-7-3-9(17)8-5-12(...         0002         8.21     0.54

### Step 6: Extract and Filter Predictions

The predictions automatically contain only the **best prediction per
compound** (the one with the narrowest confidence interval). Convert to
a data.frame for further analysis:

``` r
# Extract as data.frame
preds_df <- as.data.frame(predictions)

# View structure
head(preds_df)
```

                  compound_id     compound_name source_system source_rt predicted_rt
    1 InChI=1S/C15H14O6/c16...       epicatechin          0002     11.45        11.82
    2 InChI=1S/C22H18O10/c2... epicatechin gall          0002     15.12        15.43
    3 InChI=1S/C15H14O7/c16...     gallocatechin          0002      8.05         8.21
      ci_lower ci_upper ci_width model_r_squared model_rmse n_calibration
    1    11.15    12.49     1.34           0.992       0.38            12
    2    14.62    16.24     1.62           0.992       0.38            12
    3     7.65     8.77     1.12           0.992       0.38            12

Optionally filter to high-quality predictions:

``` r
# Keep only predictions with narrow CI and good model fit
good_preds <- preds_df[
  preds_df$ci_width < 1.5 &           # CI width < 1.5 minutes
  preds_df$model_r_squared > 0.95 &   # Model R² > 0.95
  preds_df$n_calibration >= 10,       # At least 10 calibration compounds
]

nrow(good_preds)
# [1] 380

# Summary of filtered predictions
summary(good_preds$ci_width)
#    Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
#   0.210   0.580   0.820   0.785   0.980   1.490
```

### Step 7: Save and Reuse Models

Models can be saved and reused later:

``` r
# Save models for later use
saveRDS(my_models, "my_rp_models.rds")

# Load later
my_models <- readRDS("my_rp_models.rds")

# Generate new predictions (e.g., for specific compounds)
caffeine_inchi <- "InChI=1S/C8H10N4O2/c1-10-4-9-6-5(10)7(13)12(3)8(14)11(6)2/h4H,1-3H3"
specific_preds <- predict_rt_user(my_models, compounds = caffeine_inchi)
```

## Tips

### Choosing System Type

- **RP (Reversed-Phase)**: Most common; uses C18, C8, or similar columns
  with water/organic gradients
- **HILIC**: Hydrophilic interaction; uses polar columns with high
  organic mobile phase

Models are only built between systems of the same type.

### Minimum Compounds

The `min_compounds` parameter (default: 10) sets the minimum number of
compounds that must be shared between your data and a repository system
to build a model. Lower values allow more models but with potentially
lower quality.

### CI Calculation Methods

- **`fast_ci`** (default): Uses analytical confidence intervals. Fast
  (~0.02 sec/model).
- **`bootstrap`**: Uses bootstrap resampling. More accurate but slower
  (~1.5 sec/model).

### Filtering Predictions

After generating predictions, filter based on:

- **CI width**: Narrower intervals indicate more confident predictions
- **Model R²**: Higher values indicate better calibration fit
- **Number of calibration compounds**: More compounds = more reliable
  model

## Summary

| Function                                                                                      | Purpose                                    |
|-----------------------------------------------------------------------------------------------|--------------------------------------------|
| [`build_user_models()`](https://stanstrup.github.io/rePredRet/reference/build_user_models.md) | Build calibration models from your RT data |
| [`predict_rt_user()`](https://stanstrup.github.io/rePredRet/reference/predict_rt_user.md)     | Generate predictions using built models    |
| [`plot()`](https://rdrr.io/r/graphics/plot.default.html)                                      | Visualize model diagnostics                |
| [`as.data.frame()`](https://rdrr.io/r/base/as.data.frame.html)                                | Extract predictions as data.frame          |

This workflow allows you to get RT predictions for compounds you haven’t
measured, without needing to submit your data to RepoRT.
