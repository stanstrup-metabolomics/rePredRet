# Interactive System Network

plotly version of
[`rePredRet_plot_network()`](https://stanstrup-metabolomics.github.io/rePredRet/reference/rePredRet_plot_network.md):
the links are drawn in ten similarity bins (one trace each, wider and
more opaque for higher similarity) so that a slider can hide the weaker
bins in the browser.

## Usage

``` r
plot_network_interactive(net, height = 750)
```

## Arguments

- net:

  Output of
  [`network_data()`](https://stanstrup-metabolomics.github.io/rePredRet/reference/network_data.md)
  with a layout.

- height:

  Height of the plot in pixels.

## Value

A plotly object.
