# List Available Chromatographic Systems

Returns RepoRT study information for all chromatographic systems that
have pre-built models.

## Usage

``` r
rePredRet_systems(refresh = FALSE)
```

## Arguments

- refresh:

  If TRUE, download fresh copies instead of using the cache.

## Value

Tibble with system information from RepoRT (id, name, url, method.type,
authors, ...), plus `n_compounds` (compounds with an RT in the system)
and `n_models` (models from or to it) from the published `systems.csv`.
Model sets published before `systems.csv` was added give NA for these
two columns.

## Examples

``` r
if (FALSE) { # \dontrun{
systems <- rePredRet_systems()
} # }
```
