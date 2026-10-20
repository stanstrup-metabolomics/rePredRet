# List Available Models

Returns the model index with statistics for all pre-built models.

## Usage

``` r
rePredRet_models(refresh = FALSE)
```

## Arguments

- refresh:

  If TRUE, download a fresh copy instead of using the cache.

## Value

Tibble with one row per model: sys1_id (source system), sys2_id (target
system), n_compounds, median_ci_width and median_error.

## Examples

``` r
if (FALSE) { # \dontrun{
models <- rePredRet_models()
} # }
```
