# Plot System Network

Creates a network visualization showing connections between
chromatographic systems based on pre-built model availability.

## Usage

``` r
rePredRet_plot_network(min_compounds = 20, show_labels = TRUE, refresh = FALSE)
```

## Arguments

- min_compounds:

  Minimum compounds to show edge (default 20).

- show_labels:

  If TRUE (default), label each system with its ID.

- refresh:

  If TRUE, download a fresh model index instead of using the cache.

## Value

A ggplot object.
