# Get All Predictions for a System

Retrieves all predictions for a specific chromatographic system.

## Usage

``` r
rePredRet_system_predictions(system_id, data_path = NULL)
```

## Arguments

- system_id:

  Target system identifier.

- data_path:

  Path to rePredRet-data (or NULL to use cache).

## Value

Tibble with predictions.

## Examples

``` r
if (FALSE) { # \dontrun{
preds <- rePredRet_system_predictions("0001")
} # }
```
