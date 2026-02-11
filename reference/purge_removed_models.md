# Purge Models for Removed Datasets

Deletes model directories and cache entries for models referencing
removed datasets

## Usage

``` r
purge_removed_models(removed_ids, cache, export_dir)
```

## Arguments

- removed_ids:

  Vector of removed dataset IDs

- cache:

  Build cache object

- export_dir:

  Directory containing model subdirectories

## Value

Updated cache object
