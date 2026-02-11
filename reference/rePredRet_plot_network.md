# Plot System Network

Creates a network visualization showing connections between
chromatographic systems based on model availability.

## Usage

``` r
rePredRet_plot_network(min_compounds = 20, data_path = NULL)
```

## Arguments

- min_compounds:

  Minimum compounds to show edge (default 20).

- data_path:

  Path to rePredRet-data (or NULL to use cache).

## Value

A ggplot object or network diagram.
