# rePredRet: Retention Time Prediction Statistics

## Overview

rePredRet is a reimplementation of [PredRet](http://predret.org) that
predicts retention times (RTs) by directly mapping between
chromatographic systems using monotonically constrained GAMs. It uses
data from the [RepoRT
repository](https://github.com/michaelwitting/RepoRT) and publishes
predictions to GitHub.

### Key Statistics

| Metric                   |     Value |
|:-------------------------|----------:|
| Number of Models         |       642 |
| Number of Systems        |        30 |
| Total Calibration Points |   124,996 |
| Median R²                |    0.8805 |
| Median Absolute Error    | 0.199 min |
| Median Relative Error    |     4.61% |
| Median CI Width          | 0.643 min |
| Median Relative CI Width |    13.90% |

### Comparison with Original PredRet Paper

The original PredRet paper (Stanstrup et al. 2015) reported:

| Metric                   | Paper Value   | rePredRet Value |
|--------------------------|---------------|-----------------|
| R² (predicted vs actual) | 0.9992        | 0.8805          |
| Median Error             | 0.01-0.28 min | 0.199 min       |
| Median Relative Error    | 1.8-3.7%      | 4.61%           |

### Data Source

    **RepoRT Version:** test

    **Generated:** 2026-01-15 17:24:11.700818 

## Quick Links

- [Model
  Statistics](https://stanstrup.github.io/rePredRet/articles/models.md) -
  Detailed model performance metrics
- [Predictions](https://stanstrup.github.io/rePredRet/articles/predictions.md) -
  Prediction coverage and CI widths
- [Network](https://stanstrup.github.io/rePredRet/articles/network.md) -
  System connectivity visualization
- [Methods](https://stanstrup.github.io/rePredRet/articles/methods.md) -
  Algorithm documentation

## References

1.  Stanstrup, J. et al. (2015). PredRet: Prediction of Retention Time
    by Direct Mapping between Multiple Chromatographic Systems.
    *Analytical Chemistry*, 87(18), 9421-9428.

2.  Kretschmer, F. et al. (2024). RepoRT: a comprehensive repository for
    small molecule retention times. *Nature Methods*, 21, 153-155.
