# Plot System Similarity Matrix

All pairs of systems with pre-built models as a matrix, coloured by
their similarity (Jaccard index of their compounds, see
[`rePredRet_plot_network()`](https://stanstrup-metabolomics.github.io/rePredRet/reference/rePredRet_plot_network.md)).
Systems are ordered by group (the colours of the network) and within
each group by average-linkage clustering of their similarities, so that
systems that share most of their compounds form blocks.

## Usage

``` r
rePredRet_plot_system_matrix(min_compounds = 20, refresh = FALSE)
```

## Arguments

- min_compounds:

  Minimum shared compounds for a pair to count (default 20).

- refresh:

  If TRUE, download fresh copies instead of using the cache.

## Value

A ggplot object.

## Examples

``` r
if (FALSE) { # \dontrun{
rePredRet_plot_system_matrix()
} # }
```
