# Fetch a File into the Local Cache

Downloads a file once and reuses the cached copy until it is older than
`max_age_days`.

## Usage

``` r
.fetch_file(url, cache_path, refresh = FALSE, max_age_days = 1)
```

## Arguments

- url:

  URL to download.

- cache_path:

  Path of the file relative to the cache directory.

- refresh:

  If TRUE, download even if a recent cached copy exists.

- max_age_days:

  Age after which the cached copy is downloaded again.

## Value

Local path to the cached file.
