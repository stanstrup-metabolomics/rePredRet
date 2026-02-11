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
  n_cores = parallel::detectCores(),
  method = c("fast_ci", "bootstrap"),
  save_plot = FALSE,
  plot_dir = NULL,
  dataset1 = NULL,
  dataset2 = NULL
)
```

## Arguments

- rt_matrix:

  A 2-column matrix rt_source, rt_target of common compound retention
  times.

- sys1_id:

  Source system identifier.

- sys2_id:

  Target system identifier.

- n_boot:

  Number of bootstrap iterations (default 200). Only used if method =
  "bootstrap".

- alpha:

  Significance level for PI (default 0.05 for 95% PI).

- n_cores:

  Number of CPU cores for parallel computation. Default uses all
  available cores.

- method:

  Method for confidence intervals: "fast_ci" (72x faster, default) or
  "bootstrap" (accurate, slower).

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

  "success" or "insufficient_data"

- sys1_id:

  Source system ID

- sys2_id:

  Target system ID

- n_points:

  Number of calibration compounds

- newdata:

  Vector of x-values for prediction grid (1000 points)

- ci:

  Matrix pred, lower, upper at newdata points

- calibration:

  Original RT matrix used for training

- compounds:

  InChI strings of calibration compounds

- build_time:

  Timestamp when model was built

- stats:

  Error statistics (mean, median, q95, max error and CI width)

- method:

  Method used for PI calculation

- plot_path:

  Path to saved plot HTML (if save_plot = TRUE)
