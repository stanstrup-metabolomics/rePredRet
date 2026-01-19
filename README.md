# rePredRet

**Retention Time Prediction by Direct System Mapping**

rePredRet is a reimplementation of [PredRet](http://predret.org) that uses data from the [RepoRT](https://github.com/michaelwitting/RepoRT) community repository to build retention time prediction models.

**Features:**
- ~50,000 pre-computed RT prediction models between LC-MS systems
- Interactive web viewer for model exploration
- Incremental model building with change detection
- Bootstrap confidence intervals for predictions
- Monotonically constrained GAMs for robust predictions

**Links:**
- **Website & Model Viewer:** https://stanstrup.github.io/rePredRet/
- **Model Repository:** https://github.com/stanstrup/rePredRet-models (~3.7 GB, 49,570 models)

## Installation

```r
# Install from GitHub
remotes::install_github("stanstrup/rePredRet")
```

## Quick Start

### For End Users: Access Pre-computed Predictions

```r
library(rePredRet)

# Download the latest predictions
data_path <- rePredRet_download()

# List available chromatographic systems
systems <- rePredRet_systems()

# Get predictions for a compound (by InChI or InChIKey)
preds <- rePredRet_predict("RYYVLZVUVIJVGH-UHFFFAOYSA-N")

# Get all predictions for a specific target system
system_preds <- rePredRet_system_predictions("0001")

# Visualize a model
rePredRet_plot_model("0001", "0002")
```

### For Developers: Build Your Own Models

```r
library(rePredRet)

# Download RepoRT data
report_path <- download_report(dest_dir = tempdir())

# Load datasets
report_data <- load_report_data(report_path)

# Build all pairwise models
model_results <- build_all_models(report_data)

# Generate predictions
predictions <- generate_all_predictions(
  model_results$models,
  report_data$datasets,
  report_data$studies
)

# Export results
export_models(model_results, "my_models")
export_predictions(predictions, "my_predictions")
```

## How It Works

rePredRet predicts retention times (RTs) by directly mapping between chromatographic systems using monotonically constrained Generalized Additive Models (GAMs).

**Key assumptions:**
- Elution order is largely conserved between similar chromatographic systems
- Compounds that elute early in one system tend to elute early in similar systems

**The algorithm:**
1. Find compounds measured in both source and target systems
2. Fit a monotonically constrained GAM to the RT pairs
3. Use bootstrap resampling (1000 iterations) to estimate confidence intervals
4. Predictions are filtered by CI width and calibration data density

## Functions

### User Functions (for accessing pre-computed data)

| Function | Description |
|----------|-------------|
| `rePredRet_download()` | Download latest predictions from GitHub |
| `rePredRet_systems()` | List available chromatographic systems |
| `rePredRet_models()` | List all available models |
| `rePredRet_model(from, to)` | Get a specific model with details |
| `rePredRet_predict(compound)` | Get predictions for a compound |
| `rePredRet_system_predictions(id)` | Get all predictions for a system |
| `rePredRet_plot_model(from, to)` | Visualize a model fit |
| `rePredRet_plot_errors(id)` | Plot prediction error distribution |
| `rePredRet_plot_network()` | Visualize system connections |

### Developer Functions (for building models)

| Function | Description |
|----------|-------------|
| `download_report()` | Download RepoRT repository |
| `load_report_data()` | Load datasets from RepoRT |
| `get_common_compounds()` | Find compounds in both systems |
| `build_model()` | Build a single model |
| `build_all_models()` | Build all pairwise models |
| `predict_rt()` | Predict RT for a compound |
| `generate_all_predictions()` | Generate all predictions |
| `export_models()` | Save models to files |
| `export_predictions()` | Save predictions to files |

## Dependencies

- mgcv (GAM fitting)
- boot (bootstrap resampling)
- pracma (sigmoid function)
- tidyverse (data manipulation)
- ggplot2 (visualization)

## References

1. **Original method**: Stanstrup, J. et al. (2015). PredRet: Prediction of Retention Time by Direct Mapping between Multiple Chromatographic Systems. *Analytical Chemistry*, 87(18), 9421-9428.

2. **Data source**: Kretschmer, F. et al. (2024). RepoRT: a comprehensive repository for small molecule retention times. *Nature Methods*, 21, 153-155.

## License

GPL (>= 3)

## Repository Structure

This project is split into two repositories:

- **rePredRet** (this repo): R package code, website, and documentation (~15 files)
- **[rePredRet-models](https://github.com/stanstrup/rePredRet-models)**: Pre-computed models (~3.7 GB, 49,570 models)

This separation keeps the main repository lightweight for fast cloning, while the model data is stored separately and accessed on-demand by the website.

## Related

- **[rePredRet-models](https://github.com/stanstrup/rePredRet-models)** - Pre-computed model repository
- **[RepoRT](https://github.com/michaelwitting/RepoRT)** - Community retention time data source
- **[Original PredRet](http://predret.org)** - Original web application
- **[Website](https://stanstrup.github.io/rePredRet/)** - Interactive model viewer
