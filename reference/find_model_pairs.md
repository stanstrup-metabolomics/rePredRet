# Find Model Pairs by Compound Overlap

Counts common compounds for all ordered dataset pairs with a sparse
incidence matrix and returns the pairs with at least `min_compounds`.

## Usage

``` r
find_model_pairs(rt_lookup, min_compounds)
```

## Arguments

- rt_lookup:

  Output of
  [`prepare_rt_lookup()`](https://stanstrup-metabolomics.github.io/rePredRet/reference/prepare_rt_lookup.md)

- min_compounds:

  Minimum common compounds

## Value

List of pairs with id1 and id2
