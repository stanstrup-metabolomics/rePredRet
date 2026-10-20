# rePredRet

**Retention Time Prediction by Direct System Mapping**

rePredRet is a reimplementation of [PredRet](http://predret.org) that uses data from the [RepoRT](https://github.com/michaelwitting/RepoRT) community repository to build retention time prediction models.

**Features:**
- Predict retention times in your own chromatographic system from your measured RTs
- ~50,000 pre-built models between RepoRT systems, with an interactive model viewer
- Monotonically constrained GAMs with prediction intervals
- Incremental, parallel model building with change detection

**Links:**
- **Website, tutorial & model viewer:** https://stanstrup-metabolomics.github.io/rePredRet/
- **Pre-built models:** https://github.com/stanstrup-metabolomics/rePredRet-models

## Installation

```r
# Install from GitHub
remotes::install_github("stanstrup-metabolomics/rePredRet")
```

## Quick Start

### Predict RTs for Your Own System

Provide compounds you have measured (InChI + RT). rePredRet builds models from all compatible RepoRT systems to your system and predicts RTs for the compounds measured in those systems.

```r
library(rePredRet)

# Build models from RepoRT systems to your system
models <- build_user_models(
  inchis = my_inchis,      # InChI strings of your measured compounds
  rts = my_rts,            # Their retention times (minutes)
  system_type = "RP"       # "RP" or "HILIC"
)

# Inspect model quality and find the most similar systems
summary(models)
compare_systems(models, top_n = 5)

# Predict RTs (one prediction per compound, narrowest interval)
predictions <- predict_rt_user(models)
export_predictions(predictions, "my_predictions.csv")
```

See the [tutorial](https://stanstrup-metabolomics.github.io/rePredRet/articles/tutorial.html) for a full example.

### Explore the Pre-built Models

Files are fetched on demand from the rePredRet-models repository and cached locally.

```r
# Systems and models
systems <- rePredRet_systems()
models <- rePredRet_models()

# A single model: prediction grid, calibration data and statistics
model <- rePredRet_model("0001", "0002")

# Visualize
plot_model("0001", "0002")
plot_model("0001", "0002", interactive = TRUE)
rePredRet_plot_network()
```

### Build All Models from RepoRT

```r
# Download and load RepoRT data
report_path <- download_report(dest_dir = "RepoRT_data")
report_data <- load_report_data(report_path, method_types = "RP")

# Build all pairwise models (incremental: rerunning only rebuilds what changed)
results <- build_all_models_furrr(report_data, export_dir = "models")
```

The production pipeline used by CI is `inst/scripts/run_pipeline.R`.

## How It Works

rePredRet predicts retention times (RTs) by directly mapping between chromatographic systems using monotonically constrained Generalized Additive Models (GAMs).

**Key assumptions:**
- Elution order is largely conserved between similar chromatographic systems
- Compounds that elute early in one system tend to elute early in similar systems

**The algorithm:**
1. Find compounds measured in both source and target systems (InChI without stereochemistry)
2. Fit a monotonically constrained GAM to the RT pairs, down-weighting outliers
3. Compute 95% prediction intervals whose width follows the local scatter of the compounds, calibrated by cross-validation (`method = "adaptive_pi"`, default), or by bootstrap resampling (`method = "bootstrap"`)

See the [methods article](https://stanstrup-metabolomics.github.io/rePredRet/articles/methods.html) for details.

## Functions

### Your Own Data

| Function | Description |
|----------|-------------|
| `check_input_data()` | Validate InChIs and RTs |
| `build_user_models()` | Build models from RepoRT systems to your system |
| `predict_rt_user()` | Predict RTs with those models |
| `compare_systems()` | Rank RepoRT systems by similarity to yours |
| `validate_models()` | Cross-validate model quality |
| `export_predictions()` | Export predictions to CSV or Excel |

### Pre-built Models

| Function | Description |
|----------|-------------|
| `rePredRet_systems()` | List systems with pre-built models |
| `rePredRet_models()` | List pre-built models with statistics |
| `rePredRet_model(from, to)` | Get a pre-built model |
| `plot_model(from, to)` | Plot a pre-built model |
| `rePredRet_plot_network()` | Plot the network of pre-built models |

### RepoRT Data

| Function | Description |
|----------|-------------|
| `search_compound()` | Find all systems where a compound was measured |
| `get_system_info()` | Get system metadata (column, gradient, ...) |
| `download_report()` | Download RepoRT |
| `load_report_data()` | Load datasets from RepoRT |

### Model Building

| Function | Description |
|----------|-------------|
| `build_model()` | Build a single model |
| `build_all_models_furrr()` | Build all pairwise models (incremental, parallel) |
| `get_common_compounds()` | Find compounds in both systems |
| `export_models()` | Save models to files |

## References

1. **Original method**: Stanstrup, J. et al. (2015). PredRet: Prediction of Retention Time by Direct Mapping between Multiple Chromatographic Systems. *Analytical Chemistry*, 87(18), 9421-9428.

2. **Data source**: Kretschmer, F. et al. (2024). RepoRT: a comprehensive repository for small molecule retention times. *Nature Methods*, 21, 153-155.

## License

GPL (>= 3)

## Related

- **[rePredRet-models](https://github.com/stanstrup-metabolomics/rePredRet-models)** - Pre-built models, kept in a separate repository so this one stays small
- **[RepoRT](https://github.com/michaelwitting/RepoRT)** - Community retention time data source
- **[Original PredRet](http://predret.org)** - Original web application
