#' @title System Network
#' @description Network and matrix views of how the chromatographic systems
#'   with pre-built models are connected through shared compounds.
#' @name system_network
#' @keywords internal
NULL


#' Plot System Network
#'
#' Network of the chromatographic systems with pre-built models. Each point
#' is a system and each line a pair of systems with models between them.
#'
#' @param min_compounds Minimum shared compounds for a pair to count
#'   (default 20).
#' @param layout_similarity Only links with at least this similarity (see
#'   Details) place the systems (default 0.4). All links are drawn.
#' @param interactive If TRUE, return an interactive plotly graph: hovering
#'   over a system shows its name, and a slider sets the minimum similarity
#'   of the links drawn.
#' @param show_labels If TRUE, label each system with its ID (static plot
#'   only, default FALSE).
#' @param refresh If TRUE, download fresh copies instead of using the cache.
#'
#' @return A ggplot object, or a plotly object if `interactive = TRUE`.
#'
#' @details
#' The similarity of two systems is the Jaccard index of their compounds:
#' shared compounds / compounds in either system. Unlike the number of shared
#' compounds, it does not make every system look close to the few very large
#' ones. Model sets published before `systems.csv` existed lack the number of
#' compounds per system; it is then approximated by the most compounds the
#' system shares with any other system.
#'
#' The systems are placed by a force-directed layout (Fruchterman-Reingold)
#' of the links with similarity of at least `layout_similarity`, plus each
#' system's most similar link so that no system is left unconnected. All
#' pairs with at least `min_compounds` shared compounds are then drawn, with
#' width and opacity growing with their similarity, so that weak links fade
#' into the background. Point size is the number of compounds in the system
#' and shape the method type. The colours are groups of systems that share
#' many compounds: Louvain communities of all pairs, weighted by similarity.
#' The layout and the communities use a fixed seed.
#'
#' Requires the igraph package.
#'
#' @seealso [rePredRet_plot_system_matrix()] for all pairs as a matrix.
#' @export
#' @examples
#' \dontrun{
#' rePredRet_plot_network()
#' rePredRet_plot_network(interactive = TRUE)
#' }
rePredRet_plot_network <- function(min_compounds = 20, layout_similarity = 0.4,
                                   interactive = FALSE, show_labels = FALSE,
                                   refresh = FALSE) {
  net <- network_data(min_compounds, layout_similarity, refresh = refresh)

  if (interactive) {
    return(plot_network_interactive(net))
  }

  nodes <- net$nodes
  edges <- net$edges[order(net$edges$similarity), ]  # strong links on top

  # Log legend for compounds: 10, 30, 100, ... 10,000
  log_breaks <- c(1, 3) * rep(10^(0:5), each = 2)
  plain_labels <- function(b) format(b, big.mark = ",", scientific = FALSE, trim = TRUE)

  p <- ggplot2::ggplot() +
    ggplot2::geom_segment(
      data = edges,
      ggplot2::aes(x = x_from, y = y_from, xend = x_to, yend = y_to,
                   linewidth = similarity, alpha = similarity),
      colour = "grey20"
    ) +
    ggplot2::geom_point(
      data = nodes,
      ggplot2::aes(x = x, y = y, size = n_compounds, colour = community, shape = method_type),
      alpha = 0.9
    ) +
    ggplot2::scale_linewidth(range = c(0.05, 1.2), limits = c(0, 1), name = "Similarity\n(Jaccard)") +
    ggplot2::scale_alpha(range = c(0.002, 0.3), limits = c(0, 1), name = "Similarity\n(Jaccard)") +
    ggplot2::scale_size(range = c(0.8, 7), trans = "log10", breaks = log_breaks,
                        labels = plain_labels, name = net$size_name) +
    ggplot2::scale_colour_manual(
      values = stats::setNames(net$palette, net$group_names$community),
      breaks = as.character(net$group_names$community[order(is.na(net$group_names$group), net$group_names$label)]),
      labels = stats::setNames(net$group_names$label, net$group_names$community),
      name = "Group"
    ) +
    ggplot2::scale_shape_discrete(name = "Method")

  # Letters of the named groups at the centre of their systems
  centres <- group_centres(nodes, net$group_names)
  if (nrow(centres) > 0) {
    p <- p + ggplot2::geom_label(
      data = centres, ggplot2::aes(x = x, y = y, label = group),
      fill = "white", alpha = 0.8, fontface = "bold", size = 5, linewidth = 0
    )
  }

  if (show_labels) {
    p <- p + ggplot2::geom_text(
      data = nodes,
      ggplot2::aes(x = x, y = y, label = id),
      size = 2.5, vjust = -1
    )
  }

  p +
    ggplot2::coord_equal() +
    ggplot2::theme_void() +
    ggplot2::theme(plot.caption = ggplot2::element_text(hjust = 0)) +
    ggplot2::labs(
      title = "Chromatographic System Network",
      subtitle = sprintf("%d systems, %s pairs sharing at least %d compounds",
                         nrow(nodes), format(nrow(edges), big.mark = ","), min_compounds),
      caption = group_caption(net$group_names)
    )
}


#' Plot System Similarity Matrix
#'
#' All pairs of systems with pre-built models as a matrix, coloured by their
#' similarity (Jaccard index of their compounds, see
#' [rePredRet_plot_network()]). Systems are ordered by group (the colours of
#' the network) and within each group by average-linkage clustering of their
#' similarities, so that systems that share most of their compounds form
#' blocks.
#'
#' @inheritParams rePredRet_plot_network
#'
#' @return A ggplot object.
#' @export
#' @examples
#' \dontrun{
#' rePredRet_plot_system_matrix()
#' }
rePredRet_plot_system_matrix <- function(min_compounds = 20, refresh = FALSE) {
  net <- network_data(min_compounds, layout_similarity = NULL, refresh = refresh)
  nodes <- net$nodes
  edges <- net$edges

  sim <- matrix(0, nrow(nodes), nrow(nodes), dimnames = list(nodes$id, nodes$id))
  sim[cbind(edges$from, edges$to)] <- edges$similarity
  sim[cbind(edges$to, edges$from)] <- edges$similarity

  ord <- unlist(purrr::map(levels(nodes$community), function(cm) {
    ids <- nodes$id[nodes$community == cm]
    if (length(ids) < 3) return(ids)
    ids[stats::hclust(stats::as.dist(1 - sim[ids, ids]), method = "average")$order]
  }))

  pairs <- dplyr::bind_rows(edges, dplyr::rename(edges, from = to, to = from)) |>
    dplyr::mutate(from = factor(from, levels = ord), to = factor(to, levels = rev(ord)))

  # Letters of the named groups next to their block on the diagonal
  blocks <- net$group_names[!is.na(net$group_names$group), ]
  blocks$centre <- purrr::map_dbl(as.character(blocks$community), function(cm) {
    stats::median(match(nodes$id[nodes$community == cm], ord))
  })

  ggplot2::ggplot(pairs, ggplot2::aes(x = from, y = to, fill = similarity)) +
    ggplot2::geom_tile() +
    ggplot2::geom_label(
      data = data.frame(x = blocks$centre, y = length(ord) - blocks$centre + 1, group = blocks$group),
      ggplot2::aes(x = x, y = y, label = group), inherit.aes = FALSE,
      fill = "white", alpha = 0.8, fontface = "bold", size = 5, linewidth = 0
    ) +
    ggplot2::scale_fill_viridis_c(option = "magma", direction = -1, limits = c(0, 1),
                                  name = "Similarity\n(Jaccard)") +
    ggplot2::coord_equal() +
    ggplot2::theme_void() +
    ggplot2::theme(plot.caption = ggplot2::element_text(hjust = 0)) +
    ggplot2::labs(
      title = "System Similarity Matrix",
      subtitle = sprintf("%d systems, ordered by group and similarity", nrow(nodes)),
      caption = group_caption(net$group_names)
    )
}


#' Caption Naming the Groups
#'
#' @param group_names Output of [name_communities()].
#' @return One line per named group ("A: name"), or a generic caption.
#' @keywords internal
group_caption <- function(group_names) {
  named <- group_names[!is.na(group_names$group), ]
  named <- named[order(named$group), ]
  if (nrow(named) == 0) {
    return("Colour: groups of systems that share many compounds")
  }
  paste(c("Colour: groups of systems that share many compounds. Named groups:",
          paste0(named$group, ": ", named$name)), collapse = "\n")
}


#' Centres of the Named Groups
#'
#' @param nodes Systems with x, y and community.
#' @param group_names Output of [name_communities()].
#' @return data.frame with group (letter), x and y (medians of its systems).
#' @keywords internal
group_centres <- function(nodes, group_names) {
  named <- group_names[!is.na(group_names$group), ]
  data.frame(
    group = named$group,
    x = purrr::map_dbl(as.character(named$community), function(cm) stats::median(nodes$x[nodes$community == cm])),
    y = purrr::map_dbl(as.character(named$community), function(cm) stats::median(nodes$y[nodes$community == cm]))
  )
}


#' Network Data
#'
#' Systems and links for [rePredRet_plot_network()] and
#' [rePredRet_plot_system_matrix()].
#'
#' @inheritParams rePredRet_plot_network
#' @param layout_similarity Minimum similarity of the links that place the
#'   systems, or NULL to skip the layout.
#'
#' @return List with `nodes` (id, name, method_type, n_compounds, n_models,
#'   community and, with a layout, x and y), `edges` (from, to, n_compounds,
#'   similarity and, with a layout, the coordinates of both ends), `palette`
#'   (one colour per group) and `size_name` (legend title for n_compounds).
#' @keywords internal
network_data <- function(min_compounds, layout_similarity, refresh = FALSE) {
  if (!requireNamespace("igraph", quietly = TRUE)) {
    stop("Package 'igraph' is required. Install with: install.packages('igraph')")
  }

  models <- rePredRet_models(refresh = refresh)

  # One undirected link per system pair
  edges <- models |>
    dplyr::filter(n_compounds >= min_compounds) |>
    dplyr::mutate(from = pmin(sys1_id, sys2_id), to = pmax(sys1_id, sys2_id)) |>
    dplyr::distinct(from, to, .keep_all = TRUE) |>
    dplyr::select(from, to, n_compounds)

  if (nrow(edges) == 0) {
    stop("No models meet the minimum compound threshold")
  }

  systems <- rePredRet_systems(refresh = refresh)
  nodes <- data.frame(id = sort(unique(c(edges$from, edges$to))))
  info <- match(nodes$id, systems$id)
  nodes$name <- systems$name[info]
  nodes$method_type <- dplyr::coalesce(systems$method.type[info], "unknown")
  nodes$n_models <- as.integer(table(factor(c(models$sys1_id, models$sys2_id), levels = nodes$id)))
  nodes$n_compounds <- systems$n_compounds[info]
  size_name <- "Compounds\nin system"
  if (anyNA(nodes$n_compounds)) {
    # Before systems.csv: the most compounds shared with any other system
    largest <- dplyr::bind_rows(dplyr::select(edges, id = from, n_compounds),
                                dplyr::select(edges, id = to, n_compounds)) |>
      dplyr::group_by(id) |>
      dplyr::summarize(n = max(n_compounds))
    nodes$n_compounds <- largest$n[match(nodes$id, largest$id)]
    size_name <- "Compounds in\nsystem (at least)"
  }

  size <- stats::setNames(nodes$n_compounds, nodes$id)
  edges$similarity <- unname(edges$n_compounds / (size[edges$from] + size[edges$to] - edges$n_compounds))
  edges$edge <- seq_len(nrow(edges))

  g <- igraph::graph_from_data_frame(edges, directed = FALSE, vertices = nodes["id"])

  # Groups of systems: Louvain communities of all links, numbered from the largest
  membership <- withr::with_seed(1, igraph::membership(
    igraph::cluster_louvain(g, weights = edges$similarity)))
  by_size <- as.integer(names(sort(table(membership), decreasing = TRUE)))
  nodes$community <- factor(match(membership, by_size))
  palette <- grDevices::hcl.colors(nlevels(nodes$community), "Dark 3")

  # Names of the groups that match a named reference group
  group_names <- name_communities(nodes$id, nodes$community)
  nodes$group <- group_names$label[match(nodes$community, group_names$community)]

  if (!is.null(layout_similarity)) {
    # Links that place the systems: the similar ones, plus each system's most
    # similar link so that no system is left unconnected
    most_similar <- dplyr::bind_rows(dplyr::select(edges, id = from, edge, similarity),
                                     dplyr::select(edges, id = to, edge, similarity)) |>
      dplyr::group_by(id) |>
      dplyr::slice_max(similarity, n = 1, with_ties = FALSE) |>
      dplyr::pull(edge)
    keep <- sort(unique(c(which(edges$similarity >= layout_similarity), most_similar)))

    g_layout <- igraph::subgraph_from_edges(g, keep, delete.vertices = FALSE)
    xy <- withr::with_seed(1, igraph::layout_with_fr(g_layout, weights = edges$similarity[keep],
                                                      niter = 2000))
    nodes$x <- xy[, 1]
    nodes$y <- xy[, 2]

    edges$x_from <- nodes$x[match(edges$from, nodes$id)]
    edges$y_from <- nodes$y[match(edges$from, nodes$id)]
    edges$x_to <- nodes$x[match(edges$to, nodes$id)]
    edges$y_to <- nodes$y[match(edges$to, nodes$id)]
  }

  list(nodes = nodes, edges = edges, palette = palette, size_name = size_name,
       group_names = group_names)
}


#' Name Groups of Systems
#'
#' Gives each group (community) of systems the letter and name of the named
#' reference group it matches. The reference groups
#' (`inst/extdata/system_groups.csv` and `system_group_names.csv`) were named
#' by hand from the studies that contributed their systems and the classes of
#' their compounds (RepoRT data of October 2026, see the Network article).
#' Since the groups are recomputed from the current data, a group gets a
#' reference name only when they match by mutual majority: at least half of
#' the group's systems are in the reference group and at least half of the
#' reference group's current systems are in the group. Each name is used at
#' most once.
#'
#' @param ids System ids.
#' @param community Group of each system (factor).
#' @param reference data.frame with sys_id and group (letter).
#' @param reference_names data.frame with group and name.
#'
#' @return data.frame with one row per community: community, group (letter or
#'   NA), name (or NA) and label (letter, or "unnamed" with a number).
#' @keywords internal
name_communities <- function(ids, community,
                             reference = utils::read.csv(
                               system.file("extdata", "system_groups.csv", package = "rePredRet"),
                               colClasses = "character"),
                             reference_names = utils::read.csv(
                               system.file("extdata", "system_group_names.csv", package = "rePredRet"),
                               colClasses = "character")) {
  reference <- reference[reference$sys_id %in% ids, ]
  community_size <- table(community)
  group_size <- table(reference$group)

  # Overlap of every community with every reference group
  overlap <- data.frame(community = community[match(reference$sys_id, ids)], group = reference$group) |>
    dplyr::count(community, group, name = "n") |>
    dplyr::mutate(
      share_of_community = n / as.integer(community_size[as.character(community)]),
      share_of_group = n / as.integer(group_size[group])
    ) |>
    dplyr::filter(share_of_community >= 0.5, share_of_group >= 0.5) |>
    dplyr::arrange(dplyr::desc(n))

  # Mutual majority makes the matches unique; drop duplicates defensively
  overlap <- overlap[!duplicated(overlap$community) & !duplicated(overlap$group), ]

  result <- data.frame(community = factor(levels(community), levels = levels(community)))
  result$group <- overlap$group[match(result$community, overlap$community)]
  result$name <- reference_names$name[match(result$group, reference_names$group)]
  unnamed <- is.na(result$group)
  result$label <- result$group
  result$label[unnamed] <- paste("unnamed", seq_len(sum(unnamed)))
  result
}


#' Interactive System Network
#'
#' plotly version of [rePredRet_plot_network()]: the links are drawn in ten
#' similarity bins (one trace each, wider and more opaque for higher
#' similarity) so that a slider can hide the weaker bins in the browser.
#'
#' @param net Output of [network_data()] with a layout.
#' @return A plotly object.
#' @keywords internal
plot_network_interactive <- function(net) {
  nodes <- net$nodes
  edges <- net$edges
  edges$bin <- pmin(floor(edges$similarity * 10), 9)  # 0.0-0.1 is bin 0, 0.9-1.0 bin 9

  # Line segments separated by NA, one trace per similarity bin
  edge_traces <- purrr::map(sort(unique(edges$bin)), function(b) {
    e <- edges[edges$bin == b, ]
    mid <- (b + 0.5) / 10
    list(
      x = round(as.vector(rbind(e$x_from, e$x_to, NA)), 2),
      y = round(as.vector(rbind(e$y_from, e$y_to, NA)), 2),
      lower = b / 10,
      width = 0.2 + 1.6 * mid,
      opacity = 0.005 + 0.4 * mid^2
    )
  })

  p <- purrr::reduce(edge_traces, function(p, tr) {
    plotly::add_trace(
      p, type = "scattergl", mode = "lines", x = tr$x, y = tr$y,
      line = list(width = tr$width, color = "rgb(40,40,40)"), opacity = tr$opacity,
      hoverinfo = "skip", showlegend = FALSE
    )
  }, .init = plotly::plot_ly())

  log_n <- log10(nodes$n_compounds)
  nodes$marker_size <- 5 + 15 * (log_n - min(log_n)) / max(diff(range(log_n)), 1e-9)
  # Legend and hover name of each group: "A: name" or "unnamed 1"
  gn <- net$group_names
  gn$title <- ifelse(is.na(gn$group), gn$label, paste0(gn$group, ": ", gn$name))
  nodes$hover <- sprintf(
    "<b>%s</b> %s (%s)<br>%s compounds, %s models<br>Group %s",
    nodes$id, dplyr::coalesce(nodes$name, ""), nodes$method_type,
    format(nodes$n_compounds, big.mark = ",", trim = TRUE),
    format(nodes$n_models, big.mark = ",", trim = TRUE),
    gn$title[match(nodes$community, gn$community)]
  )

  p <- purrr::reduce(levels(nodes$community), function(p, cm) {
    d <- nodes[nodes$community == cm, ]
    # WebGL like the links, so the systems are drawn on top of them
    plotly::add_trace(
      p, type = "scattergl", mode = "markers", x = d$x, y = d$y,
      marker = list(size = d$marker_size, color = net$palette[as.integer(cm)], opacity = 0.9,
                    symbol = ifelse(d$method_type == "HILIC", "diamond", "circle"),
                    line = list(color = "white", width = 0.5)),
      text = d$hover, hoverinfo = "text", name = gn$title[gn$community == cm]
    )
  }, .init = p)

  # Slider: hide the links below a minimum similarity; systems stay visible
  thresholds <- seq(0, 0.9, by = 0.1)
  lower <- purrr::map_dbl(edge_traces, "lower")
  steps <- purrr::map(thresholds, function(t) {
    list(label = format(t), method = "restyle",
         args = list("visible", as.list(c(lower >= t - 1e-9, rep(TRUE, nlevels(nodes$community))))))
  })

  axis <- list(visible = FALSE)
  plotly::layout(
    p,
    title = list(text = sprintf("Chromatographic System Network: %d systems, %s pairs",
                                nrow(nodes), format(nrow(edges), big.mark = ","))),
    xaxis = axis,
    yaxis = c(axis, list(scaleanchor = "x")),
    hovermode = "closest",
    legend = list(title = list(text = "Group")),
    annotations = purrr::pmap(group_centres(nodes, gn), function(group, x, y) {
      list(x = x, y = y, text = paste0("<b>", group, "</b>"), showarrow = FALSE,
           font = list(size = 18), bgcolor = "rgba(255,255,255,0.8)")
    }),
    sliders = list(list(
      active = 0, steps = steps,
      currentvalue = list(prefix = "Links drawn from similarity (Jaccard) ")
    ))
  )
}
