# Get Predictions for a Compound

Retrieves all predictions for a compound across all systems.

## Usage

``` r
rePredRet_predict(compound_id, target_system = NULL, data_path = NULL)
```

## Arguments

- compound_id:

  Compound identifier (InChI, InChIKey, or partial match).

- target_system:

  Optional: filter to specific target system.

- data_path:

  Path to rePredRet-data (or NULL to use cache).

## Value

Tibble with predictions.

## Examples

``` r
if (FALSE) { # \dontrun{
preds <- rePredRet_predict("RYYVLZVUVIJVGH-UHFFFAOYSA-N")
} # }
```
