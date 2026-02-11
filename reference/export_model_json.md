# Export Model Data as JSON for Web Viewer

Exports model predictions and calibration data as JSON for efficient
web-based visualization with Plotly.

## Usage

``` r
export_model_json(model, output_file)
```

## Arguments

- model:

  A model object from build_model()

- output_file:

  Path to save JSON file

## Value

Invisibly returns the JSON string
