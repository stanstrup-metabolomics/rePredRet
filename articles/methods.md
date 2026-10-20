# Methods

## Overview

rePredRet predicts retention times by directly mapping between
chromatographic systems using monotonically constrained Generalized
Additive Models (GAMs). This approach is based on the observation that
elution order is largely conserved between similar chromatographic
systems.

## Algorithm

### Core Principle

If compound A elutes before compound B in chromatographic system 1, it
will likely also elute before B in a similar system 2. This relationship
can be modeled as a monotonic function:

$$RT_{target} = f(RT_{source})$$

where $f$ is a monotonically increasing function.

### Model Fitting

1.  **Identify Common Compounds**: Find compounds present in both source
    and target systems using standardized InChI identifiers (with
    stereochemistry stripped for better matching).

2.  **Fit Monotonic GAM** (using `mgcv`):

    - If there are fewer than 4 distinct source RTs, they are jittered
      by up to 0.01 with a fixed seed, so the same data always give the
      same model
    - An unconstrained thin plate regression spline GAM gives the
      smoothing parameter
    - Compounds with large residuals are down-weighted with a sigmoid
      weight (0.5 at a residual of 10% of the maximum target RT,
      steepness 30). These two values were set by eye for the original
      PredRet and have never been optimized
    - A cubic regression spline is refitted with monotonic increasing
      constraints
      ([`mgcv::mono.con()`](https://rdrr.io/pkg/mgcv/man/mono.con.html)
      and [`mgcv::pcls()`](https://rdrr.io/pkg/mgcv/man/pcls.html))

3.  **Prediction Intervals**: three methods share the fit above. See the
    [benchmark
    article](https://stanstrup-metabolomics.github.io/rePredRet/articles/benchmark.md)
    for their coverage and speed.

    **`method = "adaptive_pi"`** (default): the width of the interval
    follows the local scatter of the calibration compounds, and its
    level is calibrated by cross-validated conformal prediction.

    - *Cross-validated residuals*: the compounds are split into 5 folds,
      interleaved along the source RT. Each fold is predicted by the
      monotone spline refitted on the other folds (with the smoothing
      parameter and weights of the full fit). These residuals include
      the uncertainty of the fitted curve.
    - *Local scale* s(x): a Gamma GAM (log link, REML, at most 5 basis
      functions) of the absolute residuals, weighted by the sigmoid
      outlier weights so that a single outlier does not widen the band
      around it. s(x) is floored at 30% of the global scale.
    - *Calibration*: the multiplier q is the conformal quantile of the
      normalized residuals \|e\| / s(x): the ⌈(n + 1) · 0.95⌉-th
      smallest. With fewer than 19 compounds this is the largest one.
    - *Interval*: fit ± q · s(x). Each bound is made non-decreasing by
      isotonic regression, which keeps lower ≤ fit ≤ upper, and the
      lower bound is truncated at 0.

    **`method = "fast_ci"`**: about 4 times faster, with an interval of
    nearly constant width.

    - fit ± z · sqrt(SE² + σ²), where SE is the pointwise standard error
      of the unconstrained GAM and σ² the mean squared leave-one-out
      residual (the constrained fit’s residuals divided by 1 − leverage)
    - Bounds are made monotone and truncated at 0 as above

    **`method = "bootstrap"`**: about 40 times slower than
    `adaptive_pi`.

    - The monotone fit is repeated on 200 bootstrap resamples. Each
      compound is predicted by the refits whose resample did not include
      it (out-of-bag), which gives out-of-sample residuals
    - These residuals go through the same local scale and conformal
      calibration as `adaptive_pi` (jackknife+-after-bootstrap)

    On held-out compounds from 135 RepoRT system pairs the 95% intervals
    covered 94.6% (`adaptive_pi`), 95.5% (`fast_ci`) and 94.8%
    (`bootstrap`) of the compounds. In the regions with the least / most
    scatter `adaptive_pi` covered 94.8% / 94.5%, while `fast_ci` covered
    97.4% / 93.4%: a constant width is too wide where compounds agree
    and too narrow where they scatter.

    **Few compounds**: with only 10–15 shared compounds no method
    reaches the nominal level. Coverage was 86% (`adaptive_pi`), 90%
    (`fast_ci`) and 89% (`bootstrap`). With so few compounds, a 95%
    interval would have to be very wide to hold. Treat intervals from
    models with fewer than about 20 compounds as roughly 85–90%
    intervals.

4.  **Quality Filters**: Predictions are only made where:

    - CI width \< 2 minutes OR relative CI width \< 20%
    - Source RT is in a dense region of calibration data (density check)

### GAM Specification

``` r
# Unconstrained fit for the smoothing parameter
fit <- mgcv::gam(y ~ s(x, k = k, bs = "tp"), data = calibration_data)

# Monotonic constrained refit
sm <- mgcv::smoothCon(mgcv::s(x, k = k, bs = "cr"), calibration_data)[[1]]
con <- mgcv::mono.con(sm$xp)
p <- mgcv::pcls(list(X = sm$X, y = y, w = w, sp = fit$sp, S = sm$S,
                     Ain = con$A, bin = con$b, p = sm$xp, C = matrix(0, 0, 0), off = 0))
```

- `k`: number of unique source RTs, at most 10
- `w`: sigmoid weights from the unconstrained residuals

## Data Source

### RepoRT Repository

Retention time data comes from
[RepoRT](https://github.com/michaelwitting/RepoRT) (Kretschmer et
al. 2024), a comprehensive repository containing:

- Multiple chromatographic method types (RP, HILIC, etc.)
- Standardized compound identifiers (InChI, InChIKey, SMILES)
- Rich metadata for each chromatographic system

### Compound Matching

Compounds are matched between systems using standardized InChI:

1.  Strip stereochemistry layers (`/t`, `/b`, `/m`, `/s`) from InChI
2.  Match on the remaining structural connectivity
3.  For duplicate measurements, use median RT value

## Implementation

### R Package Structure

    rePredRet/
    ├── R/
    │   ├── data_loading.R           # RepoRT data loading
    │   ├── gam_fitting.R            # Monotone GAM fit and bootstrap intervals
    │   ├── gam_adaptive_pi.R        # Default (locally adaptive) prediction intervals
    │   ├── gam_fast_ci.R            # Faster prediction intervals of nearly constant width
    │   ├── model_building.R         # Model construction
    │   ├── build_all_models_furrr.R # Incremental parallel pipeline
    │   ├── build_user_models.R      # Models from user RT data
    │   ├── prediction.R             # RT predictions
    │   ├── statistics.R             # Performance metrics
    │   └── user_functions.R         # Access to pre-built models
    ├── inst/scripts/run_pipeline.R  # Pipeline used by CI
    ├── tests/
    │   └── testthat/                # Unit tests
    └── vignettes/                   # This documentation

### Key Functions

| Function                                                                                                                 | Purpose                                           |
|--------------------------------------------------------------------------------------------------------------------------|---------------------------------------------------|
| [`download_report()`](https://stanstrup-metabolomics.github.io/rePredRet/reference/download_report.md)                   | Download RepoRT data                              |
| [`load_report_data()`](https://stanstrup-metabolomics.github.io/rePredRet/reference/load_report_data.md)                 | Load datasets from RepoRT                         |
| [`get_common_compounds()`](https://stanstrup-metabolomics.github.io/rePredRet/reference/get_common_compounds.md)         | Find shared compounds                             |
| [`build_model()`](https://stanstrup-metabolomics.github.io/rePredRet/reference/build_model.md)                           | Build single pairwise model                       |
| [`build_all_models_furrr()`](https://stanstrup-metabolomics.github.io/rePredRet/reference/build_all_models_furrr.md)     | Build all pairwise models (incremental, parallel) |
| [`predict_rt()`](https://stanstrup-metabolomics.github.io/rePredRet/reference/predict_rt.md)                             | Predict RT for a compound                         |
| [`generate_all_predictions()`](https://stanstrup-metabolomics.github.io/rePredRet/reference/generate_all_predictions.md) | Generate all predictions                          |

### Computational Considerations

- **Parallel Processing**: Models are built in parallel with `furrr`,
  one worker per available CPU core
- **Incremental Builds**: Only models whose datasets changed are rebuilt
- **Model Pruning**: Boot objects are pruned to reduce storage by ~95%

## Validation

### Cross-Validation

Model quality is assessed using leave-one-out style evaluation at
calibration points:

- **R²**: Coefficient of determination for predicted vs actual RT
- **Absolute Error**: \|predicted - actual\| in minutes
- **Relative Error**: \|predicted - actual\| / actual × 100%
- **CI Coverage**: Proportion of actual values within CI

### Expected Performance

Based on the original PredRet paper (Stanstrup et al. 2015):

| Metric                | Expected Range  |
|-----------------------|-----------------|
| R²                    | \> 0.99         |
| Median Absolute Error | 0.01 - 0.28 min |
| Median Relative Error | 1.8 - 3.7%      |
| Median CI Width       | 0.08 - 1.86 min |

## References

1.  **PredRet Algorithm**: Stanstrup, J., Neumann, S., & Ivanova, I. V.
    (2015). PredRet: Prediction of Retention Time by Direct Mapping
    between Multiple Chromatographic Systems. *Analytical Chemistry*,
    87(18), 9421-9428. doi:10.1021/acs.analchem.5b02287

2.  **RepoRT Database**: Kretschmer, F., Witting, M., & al. (2024).
    RepoRT: a comprehensive repository for small molecule retention
    times. *Nature Methods*, 21, 153-155. doi:10.1038/s41592-023-02143-z

3.  **Generalized Additive Models**: Wood, S. N. (2017). *Generalized
    Additive Models: An Introduction with R* (2nd ed.). Chapman and
    Hall/CRC.

4.  **Conformal Prediction Intervals**: Lei, J., G’Sell, M., Rinaldo,
    A., Tibshirani, R. J., & Wasserman, L. (2018). Distribution-free
    predictive inference for regression. *Journal of the American
    Statistical Association*, 113(523), 1094-1111.

## Source Code

The complete source code is available on GitHub:

- **rePredRet package**:
  [github.com/stanstrup-metabolomics/rePredRet](https://github.com/stanstrup-metabolomics/rePredRet)
- **Pre-built models**:
  [github.com/stanstrup-metabolomics/rePredRet-models](https://github.com/stanstrup-metabolomics/rePredRet-models)
- **RepoRT**:
  [github.com/michaelwitting/RepoRT](https://github.com/michaelwitting/RepoRT)
