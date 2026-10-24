# Table of Systems

One row per dataset for `systems.csv`: id, name and method type from
RepoRT's studies table, the number of compounds (unique InChIs without
stereochemistry) and the number of models from or to the system.

## Usage

``` r
systems_table(rt_lookup, studies, index_df)
```

## Arguments

- rt_lookup:

  Output of
  [`prepare_rt_lookup()`](https://stanstrup-metabolomics.github.io/rePredRet/reference/prepare_rt_lookup.md)

- studies:

  RepoRT studies table, or NULL

- index_df:

  Model index with sys1_id and sys2_id

## Value

data.frame with sys_id, name, method_type, n_compounds, n_models
