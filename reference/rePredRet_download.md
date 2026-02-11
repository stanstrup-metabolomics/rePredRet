# Download rePredRet Data

Downloads the latest predictions and models from the rePredRet-data
GitHub repository.

## Usage

``` r
rePredRet_download(cache_dir = NULL, force_download = FALSE)
```

## Arguments

- cache_dir:

  Directory to cache downloaded data. Default uses a user-specific cache
  directory.

- force_download:

  If TRUE, re-download even if cached data exists.

## Value

Path to the cached data directory.

## Examples

``` r
if (FALSE) { # \dontrun{
data_path <- rePredRet_download()
} # }
```
