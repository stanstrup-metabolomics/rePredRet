# Build a Single Model Between Two Systems

Builds a monotonically constrained GAM model to predict retention times
from one chromatographic system to another.

## Usage

``` r
build_model(
  rt_matrix,
  sys1_id,
  sys2_id,
  n_boot = 200,
  alpha = 0.05,
  n_cores = parallelly::availableCores(),
  method = c("adaptive_pi", "fast_ci", "bootstrap"),
  save_plot = FALSE,
  plot_dir = NULL,
  dataset1 = NULL,
  dataset2 = NULL
)
```

## Arguments

- rt_matrix:

  A 2-column matrix `[rt_source, rt_target]` of common compound
  retention times.

- sys1_id:

  Source system identifier.

- sys2_id:

  Target system identifier.

- n_boot:

  Number of bootstrap iterations (default 1000). Only used if method =
  "bootstrap".

- alpha:

  Significance level for PI (default 0.05 for 95% PI).

- n_cores:

  Number of CPU cores for parallel computation. Default uses all
  available cores.

- method:

  Method for prediction intervals: "adaptive_pi" (default, width follows
  the local scatter, see
  [`gam_adaptive_pi()`](https://stanstrup-metabolomics.github.io/rePredRet/reference/gam_adaptive_pi.md)),
  "fast_ci" (nearly constant width, about 4x faster, see
  [`gam_fast_ci()`](https://stanstrup-metabolomics.github.io/rePredRet/reference/gam_fast_ci.md))
  or "bootstrap" (width follows the local scatter, about 40x slower, see
  [`boot2ci_PI()`](https://stanstrup-metabolomics.github.io/rePredRet/reference/boot2ci_PI.md)).

- save_plot:

  Logical, whether to save an interactive plot (default FALSE).

- plot_dir:

  Directory to save plot HTML files (required if save_plot = TRUE).

- dataset1:

  Optional dataset object for compound names in plots.

- dataset2:

  Optional dataset object for compound names in plots.

## Value

A list with model information:

- status:

  "success", "insufficient_data" or "fit_failed"

- sys1_id:

  Source system ID

- sys2_id:

  Target system ID

- n_points:

  Number of calibration compounds

- newdata:

  Vector of x-values for prediction grid (1000 points)

- ci:

  Matrix `[pred, lower, upper]` at newdata points

- calibration:

  Original RT matrix used for training

- compounds:

  InChI strings of calibration compounds

- build_time:

  Timestamp when model was built

- stats:

  Error statistics (mean, median, q95, max error and PI width) and fit
  quality: `r_squared` (of the fit at the calibration compounds),
  `coverage` (fraction of calibration compounds inside the PI) and
  `rel_mean_ci_width` (mean PI width divided by the target RT range)

- method:

  Method used for PI calculation

- plot_path:

  Path to saved plot HTML (if save_plot = TRUE)
