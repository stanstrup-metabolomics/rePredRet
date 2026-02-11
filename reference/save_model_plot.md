# Save Model Plot as Self-Contained HTML

Saves an interactive model plot as a self-contained HTML file.

## Usage

``` r
save_model_plot(
  model,
  output_dir,
  filename = NULL,
  dataset1 = NULL,
  dataset2 = NULL
)
```

## Arguments

- model:

  A model object from
  [`build_model()`](https://stanstrup.github.io/rePredRet/reference/build_model.md).

- output_dir:

  Directory to save the HTML file.

- filename:

  Optional filename (default: auto-generated from model IDs).

- dataset1:

  Optional dataset object for compound names (source system).

- dataset2:

  Optional dataset object for compound names (target system).

## Value

Path to the saved HTML file (invisibly).
