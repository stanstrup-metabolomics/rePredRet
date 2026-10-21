# rePredRet: Retention Time Prediction by Direct System Mapping

## Overview

rePredRet is a reimplementation of [PredRet](http://predret.org) that
predicts retention times (RTs) by directly mapping between
chromatographic systems using monotonically constrained GAMs. It uses
data from the [RepoRT
repository](https://github.com/michaelwitting/RepoRT) and publishes
pre-built models between all RepoRT systems that share enough compounds
to
[rePredRet-models](https://github.com/stanstrup-metabolomics/rePredRet-models).

To predict RTs in your own chromatographic system, see the
[tutorial](https://stanstrup-metabolomics.github.io/rePredRet/articles/tutorial.md).
To explore individual models, use the [model
viewer](https://stanstrup-metabolomics.github.io/rePredRet/model_viewer.md).

### Key Statistics

| Metric                           |     Value |
|:---------------------------------|----------:|
| Number of Models                 |    49,570 |
| Number of Systems                |       328 |
| Total Calibration Points         | 3,356,376 |
| Median Compounds per Model       |        37 |
| Median Absolute Error            | 0.142 min |
| Median Prediction Interval Width | 3.822 min |

Errors are the median absolute differences between fitted and measured
RTs of the calibration compounds of each model; the table shows the
median over all models.

### Comparison with Original PredRet Paper

The original PredRet paper (Stanstrup et al. 2015) reported:

| Metric       | Paper Value   | rePredRet Value |
|--------------|---------------|-----------------|
| Median Error | 0.01-0.28 min | 0.142 min       |

### Data Source

**RepoRT version (commit):**
[9de8d60](https://github.com/michaelwitting/RepoRT/tree/9de8d603377bbb6cc0f74e250eb7533fd8874df1)

## Quick Links

- [Tutorial](https://stanstrup-metabolomics.github.io/rePredRet/articles/tutorial.md) -
  Predict RTs for your own data
- [Model
  Statistics](https://stanstrup-metabolomics.github.io/rePredRet/articles/models.md) -
  Model performance metrics
- [Network](https://stanstrup-metabolomics.github.io/rePredRet/articles/network.md) -
  System connectivity
- [Methods](https://stanstrup-metabolomics.github.io/rePredRet/articles/methods.md) -
  Algorithm documentation

## References

1.  Stanstrup, J. et al. (2015). PredRet: Prediction of Retention Time
    by Direct Mapping between Multiple Chromatographic Systems.
    *Analytical Chemistry*, 87(18), 9421-9428.

2.  Kretschmer, F. et al. (2024). RepoRT: a comprehensive repository for
    small molecule retention times. *Nature Methods*, 21, 153-155.
