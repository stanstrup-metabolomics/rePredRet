# Get Detailed System Information

Retrieves detailed information about a chromatographic system including
column, gradient, and metadata.

## Usage

``` r
get_system_info(system_id, data_path = NULL)

# S3 method for class 'system_info'
print(x, ...)
```

## Arguments

- system_id:

  System identifier (e.g., "0001").

- data_path:

  Optional path to local RepoRT data directory.

## Value

A list with components: id, name, type, n_compounds, info, metadata,
gradient.

## Examples

``` r
if (FALSE) { # \dontrun{
# Get info about system 0001
info <- get_system_info("0001")
info$name
info$n_compounds
info$metadata
} # }
```
