# Find Model Pairs by Compound Overlap

Counts common compounds for all ordered dataset pairs with a sparse
incidence matrix and returns the pairs with at least `min_compounds`.

## Usage

``` r
find_model_pairs(rt_lookup, min_compounds, groups = NULL)
```

## Arguments

- rt_lookup:

  Output of
  [`prepare_rt_lookup()`](https://stanstrup-metabolomics.github.io/rePredRet/reference/prepare_rt_lookup.md)

- min_compounds:

  Minimum common compounds

- groups:

  Optional group of each dataset (e.g. its method type), in the order of
  `rt_lookup`. Pairs are only formed within a group; datasets with a
  missing group form a group of their own.

## Value

List of pairs with id1 and id2
