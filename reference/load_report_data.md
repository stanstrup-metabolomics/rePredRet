# Load RepoRT Data

Loads all datasets from a RepoRT repository into a structured list.

## Usage

``` r
load_report_data(report_path, method_types = c("RP", "HILIC"))
```

## Arguments

- report_path:

  Path to the extracted RepoRT repository.

- method_types:

  Character vector of method types to include. Default is c("RP",
  "HILIC"). Use NULL to include all.

## Value

A list with two elements:

- studies:

  Tibble with all study metadata from studies.tsv

- datasets:

  Named list of datasets, each containing: id, info, metadata, gradient,
  rtdata

## Examples

``` r
if (FALSE) { # \dontrun{
report_data <- load_report_data("path/to/RepoRT-master")
} # }
```
