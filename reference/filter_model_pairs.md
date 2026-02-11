# Filter Model Pairs Based on Changes

Determines which model pairs need to be built based on dataset changes

## Usage

``` r
filter_model_pairs(model_pairs, changes, cache, export_dir)
```

## Arguments

- model_pairs:

  List of all potential model pairs

- changes:

  Dataset change analysis

- cache:

  Build cache object

- export_dir:

  Directory for model outputs

## Value

Filtered list of pairs that need building
