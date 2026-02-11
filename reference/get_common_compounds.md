# Get Common Compounds Between Two Systems

Finds compounds that exist in both chromatographic systems and returns
their retention times as a matrix suitable for model building.

## Usage

``` r
get_common_compounds(dataset1, dataset2, by = "inchi.std")
```

## Arguments

- dataset1:

  First dataset (from load_report_data).

- dataset2:

  Second dataset (from load_report_data).

- by:

  Column name to match compounds. Default "inchi.std".

## Value

A 2-column matrix with retention times rt1, rt2 for compounds present in
both systems. If multiple measurements exist for the same compound, the
median RT is used.
