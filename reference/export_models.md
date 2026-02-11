# Export Models to Files

Saves models as individual RDS files with accompanying CSV files for
transparency.

## Usage

``` r
export_models(model_results, output_dir = "models")
```

## Arguments

- model_results:

  Output from
  [`build_all_models()`](https://stanstrup.github.io/rePredRet/reference/build_all_models.md).

- output_dir:

  Directory to save models.

## Value

Invisibly returns the output directory path.
