# Network Data

Systems and links for
[`rePredRet_plot_network()`](https://stanstrup-metabolomics.github.io/rePredRet/reference/rePredRet_plot_network.md)
and
[`rePredRet_plot_system_matrix()`](https://stanstrup-metabolomics.github.io/rePredRet/reference/rePredRet_plot_system_matrix.md).

## Usage

``` r
network_data(min_compounds, layout_similarity, refresh = FALSE)
```

## Arguments

- min_compounds:

  Minimum shared compounds for a pair to count (default 20).

- layout_similarity:

  Minimum similarity of the links that place the systems, or NULL to
  skip the layout.

- refresh:

  If TRUE, download fresh copies instead of using the cache.

## Value

List with `nodes` (id, name, method_type, n_compounds, n_models,
community and, with a layout, x and y), `edges` (from, to, n_compounds,
similarity and, with a layout, the coordinates of both ends), `palette`
(one colour per group) and `size_name` (legend title for n_compounds).
