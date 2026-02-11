# Fast Model Export (No Plots)

Exports model to RDS + CSV + JSON without generating HTML plots. Much
faster for large-scale pipelines.

## Usage

``` r
export_model_fast(model, model_dir)
```

## Arguments

- model:

  A model object

- model_dir:

  Directory for output files

## Value

Invisibly returns paths to created files
