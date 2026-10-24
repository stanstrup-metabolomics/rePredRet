# Plot System Network

Network of the chromatographic systems with pre-built models. Each point
is a system and each line a pair of systems with models between them.

## Usage

``` r
rePredRet_plot_network(
  min_compounds = 20,
  layout_similarity = 0.4,
  interactive = FALSE,
  show_labels = FALSE,
  refresh = FALSE
)
```

## Arguments

- min_compounds:

  Minimum shared compounds for a pair to count (default 20).

- layout_similarity:

  Only links with at least this similarity (see Details) place the
  systems (default 0.4). All links are drawn.

- interactive:

  If TRUE, return an interactive plotly graph: hovering over a system
  shows its name, and a slider sets the minimum similarity of the links
  drawn.

- show_labels:

  If TRUE, label each system with its ID (static plot only, default
  FALSE).

- refresh:

  If TRUE, download fresh copies instead of using the cache.

## Value

A ggplot object, or a plotly object if `interactive = TRUE`.

## Details

The similarity of two systems is the Jaccard index of their compounds:
shared compounds / compounds in either system. Unlike the number of
shared compounds, it does not make every system look close to the few
very large ones. Model sets published before `systems.csv` existed lack
the number of compounds per system; it is then approximated by the most
compounds the system shares with any other system.

The systems are placed by a force-directed layout (Fruchterman-Reingold)
of the links with similarity of at least `layout_similarity`, plus each
system's most similar link so that no system is left unconnected. All
pairs with at least `min_compounds` shared compounds are then drawn,
with width and opacity growing with their similarity, so that weak links
fade into the background. Point size is the number of compounds in the
system and shape the method type. The colours are groups of systems that
share many compounds: Louvain communities of all pairs, weighted by
similarity. The layout and the communities use a fixed seed.

Requires the igraph package.

## See also

[`rePredRet_plot_system_matrix()`](https://stanstrup-metabolomics.github.io/rePredRet/reference/rePredRet_plot_system_matrix.md)
for all pairs as a matrix.

## Examples

``` r
if (FALSE) { # \dontrun{
rePredRet_plot_network()
rePredRet_plot_network(interactive = TRUE)
} # }
```
