# Get Model Details

Retrieves a specific model with CI grid and calibration data.

## Usage

``` r
rePredRet_model(from_id, to_id, data_path = NULL)
```

## Arguments

- from_id:

  Source system identifier.

- to_id:

  Target system identifier.

- data_path:

  Path to rePredRet-data (or NULL to use cache).

## Value

A list with:

- ci_grid:

  Tibble with x, pred, lower, upper

- calibration_data:

  Tibble with rt_source, rt_target, compound

- model:

  RDS model object if available

## Examples

``` r
if (FALSE) { # \dontrun{
model_info <- rePredRet_model("0001", "0002")
} # }
```
