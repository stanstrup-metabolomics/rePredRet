# Analyze Dataset Changes

Compares current datasets against cache to detect changes

## Usage

``` r
analyze_dataset_changes(datasets, cache)
```

## Arguments

- datasets:

  Named list of current datasets

- cache:

  Build cache object

## Value

List with new_ids, changed_ids, unchanged_ids, removed_ids
