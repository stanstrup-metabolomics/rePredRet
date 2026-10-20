# Prepare Per-Dataset RT Lookup

Strips stereochemistry from InChIs and takes the median RT per compound,
once per dataset.
[`common_rt_matrix()`](https://stanstrup-metabolomics.github.io/rePredRet/reference/common_rt_matrix.md)
then matches two lookups, giving the same result as
[`get_common_compounds()`](https://stanstrup-metabolomics.github.io/rePredRet/reference/get_common_compounds.md)
much faster.

## Usage

``` r
prepare_rt_lookup(datasets)
```

## Arguments

- datasets:

  Named list of datasets with rtdata

## Value

Named list of lists with sorted `key` and matching median `rt`
