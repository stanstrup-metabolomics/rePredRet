#' @title User Convenience Functions
#' @description Easy-to-use functions for end users to access pre-built
#'   rePredRet models. Files are fetched on demand from the
#'   rePredRet-models GitHub repository and cached locally.
#' @name user_functions
#' @keywords internal
NULL

# GitHub repository with the published models
.models_repo <- "stanstrup-metabolomics/rePredRet-models"
.models_url_base <- paste0("https://raw.githubusercontent.com/", .models_repo, "/main")
.studies_url <- "https://raw.githubusercontent.com/michaelwitting/RepoRT/master/processed_data/studies.tsv"


#' Fetch a File into the Local Cache
#'
#' Downloads a file once and reuses the cached copy until it is older than
#' `max_age_days`.
#'
#' @param url URL to download.
#' @param cache_path Path of the file relative to the cache directory.
#' @param refresh If TRUE, download even if a recent cached copy exists.
#' @param max_age_days Age after which the cached copy is downloaded again.
#'
#' @return Local path to the cached file.
#' @keywords internal
.fetch_file <- function(url, cache_path, refresh = FALSE, max_age_days = 1) {
  local_file <- file.path(rappdirs::user_cache_dir("rePredRet"), cache_path)

  is_recent <- file.exists(local_file) &&
    difftime(Sys.time(), file.mtime(local_file), units = "days") < max_age_days

  if (refresh || !is_recent) {
    dir.create(dirname(local_file), recursive = TRUE, showWarnings = FALSE)
    tmp_file <- paste0(local_file, ".part")
    ok <- tryCatch(
      utils::download.file(url, tmp_file, mode = "wb", quiet = TRUE) == 0,
      error = function(e) FALSE,
      warning = function(w) FALSE
    )
    if (!ok) {
      unlink(tmp_file)
      stop("Could not download ", url, call. = FALSE)
    }
    file.rename(tmp_file, local_file)
  }

  local_file
}


#' List Available Chromatographic Systems
#'
#' Returns RepoRT study information for all chromatographic systems that
#' have pre-built models.
#'
#' @param refresh If TRUE, download fresh copies instead of using the cache.
#'
#' @return Tibble with system information from RepoRT (id, name, url,
#'   method.type, authors, ...), plus `n_compounds` (compounds with an RT in
#'   the system) and `n_models` (models from or to it) from the published
#'   `systems.csv`. Model sets published before `systems.csv` was added give
#'   NA for these two columns.
#'
#' @importFrom readr read_tsv
#' @export
#' @examples
#' \dontrun{
#' systems <- rePredRet_systems()
#' }
rePredRet_systems <- function(refresh = FALSE) {
  studies_file <- .fetch_file(.studies_url, "RepoRT/studies.tsv", refresh = refresh)
  studies <- readr::read_tsv(studies_file, col_types = readr::cols(id = "c", .default = "c"),
                             progress = FALSE)

  models <- rePredRet_models(refresh = refresh)
  system_ids <- unique(c(models$sys1_id, models$sys2_id))
  studies <- studies[studies$id %in% system_ids, ]

  systems_file <- tryCatch(
    .fetch_file(paste0(.models_url_base, "/systems.csv"), "models/systems.csv", refresh = refresh),
    error = function(e) NULL
  )
  counts <- if (!is.null(systems_file)) {
    readr::read_csv(systems_file, col_types = readr::cols(sys_id = "c"), progress = FALSE)
  } else {
    data.frame(sys_id = character(), n_compounds = integer(), n_models = integer())
  }
  studies$n_compounds <- counts$n_compounds[match(studies$id, counts$sys_id)]
  studies$n_models <- counts$n_models[match(studies$id, counts$sys_id)]

  studies
}


#' List Available Models
#'
#' Returns the model index with statistics for all pre-built models.
#'
#' @param refresh If TRUE, download a fresh copy instead of using the cache.
#'
#' @return Tibble with one row per model: sys1_id (source system), sys2_id
#'   (target system), n_compounds, median_ci_width and median_error.
#'
#' @importFrom readr read_csv
#' @export
#' @examples
#' \dontrun{
#' models <- rePredRet_models()
#' }
rePredRet_models <- function(refresh = FALSE) {
  index_file <- .fetch_file(paste0(.models_url_base, "/model_index.csv"),
                            "models/model_index.csv", refresh = refresh)
  readr::read_csv(index_file,
                  col_types = readr::cols(sys1_id = "c", sys2_id = "c"),
                  progress = FALSE)
}


#' Get Model Details
#'
#' Retrieves a pre-built model with its prediction grid and calibration data.
#'
#' @param from_id Source system identifier.
#' @param to_id Target system identifier.
#' @param refresh If TRUE, download a fresh copy instead of using the cache.
#'
#' @return A list with:
#'   \describe{
#'     \item{sys1_id, sys2_id}{Source and target system}
#'     \item{n_compounds}{Number of calibration compounds}
#'     \item{method}{Interval method}
#'     \item{predictions}{Data frame with x, pred, lower, upper}
#'     \item{calibration}{Data frame with x, y, compound}
#'     \item{stats}{List with median_ci_width, median_error, r_squared,
#'       coverage (fraction of calibration compounds inside the PI),
#'       rel_mean_ci_width (mean PI width / target RT range) and build_time.
#'       Models built before these fields were added lack them.}
#'   }
#'
#' @export
#' @examples
#' \dontrun{
#' model <- rePredRet_model("0001", "0002")
#' }
rePredRet_model <- function(from_id, to_id, refresh = FALSE) {
  model_key <- paste0(from_id, "_to_", to_id)

  json_file <- tryCatch(
    .fetch_file(paste0(.models_url_base, "/models/", model_key, "/model.json"),
                file.path("models", model_key, "model.json"),
                refresh = refresh, max_age_days = 7),
    error = function(e) {
      stop("Model not found: ", from_id, " -> ", to_id,
           "\nUse rePredRet_models() to see available models.", call. = FALSE)
    }
  )

  jsonlite::fromJSON(json_file)
}


#' Plot System Network
#'
#' Network of the chromatographic systems with pre-built models. Each point
#' is a system and each line a pair of systems with models between them.
#'
#' @param min_compounds Minimum shared compounds for a pair to count
#'   (default 20).
#' @param n_links Number of links drawn per system (default 3): the pairs
#'   with the most shared compounds. A line is drawn if it is among the
#'   `n_links` strongest of either of its systems. Drawing all pairs gives an
#'   unreadable tangle.
#' @param show_labels If TRUE, label each system with its ID (default FALSE).
#' @param refresh If TRUE, download a fresh model index instead of using the cache.
#'
#' @return A ggplot object.
#'
#' @details
#' Both directions of a pair share the same compounds, so each pair is one
#' line. Line width is the number of shared compounds and point size the
#' number of compounds in the system (the number of models for model sets
#' published before `systems.csv` existed), both on log scales. The shape
#' shows the method type. The colour shows groups of systems that share
#' many compounds: Louvain communities of all pairs (weighted by log10
#' shared compounds), not only of the lines drawn, so a group can appear in
#' more than one place. The layout (Fruchterman-Reingold on the lines drawn)
#' and the communities use a fixed seed.
#'
#' Requires the igraph package.
#'
#' @export
rePredRet_plot_network <- function(min_compounds = 20, n_links = 3, show_labels = FALSE,
                                   refresh = FALSE) {
  if (!requireNamespace("igraph", quietly = TRUE)) {
    stop("Package 'igraph' is required. Install with: install.packages('igraph')")
  }

  # One undirected edge per system pair
  edges <- rePredRet_models(refresh = refresh) |>
    dplyr::filter(n_compounds >= min_compounds) |>
    dplyr::mutate(from = pmin(sys1_id, sys2_id), to = pmax(sys1_id, sys2_id)) |>
    dplyr::distinct(from, to, .keep_all = TRUE) |>
    dplyr::select(from, to, n_compounds)

  if (nrow(edges) == 0) {
    stop("No models meet the minimum compound threshold")
  }

  g <- igraph::graph_from_data_frame(edges, directed = FALSE)
  weights <- log10(edges$n_compounds)

  # Communities of all pairs, numbered from the largest
  membership <- withr::with_seed(1, igraph::membership(igraph::cluster_louvain(g, weights = weights)))
  by_size <- as.integer(names(sort(table(membership), decreasing = TRUE)))

  systems <- rePredRet_systems(refresh = refresh)
  nodes <- data.frame(id = igraph::V(g)$name)
  info <- match(nodes$id, systems$id)
  nodes$method_type <- dplyr::coalesce(systems$method.type[info], "unknown")
  nodes$community <- factor(match(membership, by_size))
  nodes$node_size <- systems$n_compounds[info]
  size_name <- "Compounds\nin system"
  if (anyNA(nodes$node_size)) {
    nodes$node_size <- igraph::degree(g)
    size_name <- "Models"
  }

  # The n_links strongest links of each system
  edges$edge <- seq_len(nrow(edges))
  keep <- dplyr::bind_rows(dplyr::select(edges, node = from, edge, n_compounds),
                           dplyr::select(edges, node = to, edge, n_compounds)) |>
    dplyr::group_by(node) |>
    dplyr::slice_max(n_compounds, n = n_links, with_ties = FALSE) |>
    dplyr::pull(edge) |>
    unique() |>
    sort()
  g_shown <- igraph::subgraph_from_edges(g, keep, delete.vertices = FALSE)

  xy <- withr::with_seed(1, igraph::layout_with_fr(g_shown, weights = weights[keep]))
  nodes$x <- xy[, 1]
  nodes$y <- xy[, 2]

  shown <- edges[keep, ]
  shown$x_from <- nodes$x[match(shown$from, nodes$id)]
  shown$y_from <- nodes$y[match(shown$from, nodes$id)]
  shown$x_to <- nodes$x[match(shown$to, nodes$id)]
  shown$y_to <- nodes$y[match(shown$to, nodes$id)]

  # Log legends labelled 1, 3, 10, 30, ... 10,000
  log_breaks <- c(1, 3) * rep(10^(0:5), each = 2)
  plain_labels <- function(b) format(b, big.mark = ",", scientific = FALSE, trim = TRUE)

  p <- ggplot2::ggplot() +
    ggplot2::geom_segment(
      data = shown,
      ggplot2::aes(x = x_from, y = y_from, xend = x_to, yend = y_to, linewidth = n_compounds),
      colour = "grey50", alpha = 0.25
    ) +
    ggplot2::geom_point(
      data = nodes,
      ggplot2::aes(x = x, y = y, size = node_size, colour = community, shape = method_type),
      alpha = 0.85
    ) +
    ggplot2::scale_linewidth(range = c(0.15, 2.5), trans = "log10", breaks = log_breaks,
                             labels = plain_labels, name = "Shared\ncompounds") +
    ggplot2::scale_size(range = c(0.8, 7), trans = "log10", breaks = log_breaks,
                        labels = plain_labels, name = size_name) +
    ggplot2::scale_shape_discrete(name = "Method") +
    ggplot2::guides(colour = "none")

  if (nlevels(nodes$community) <= 8) {
    p <- p + ggplot2::scale_colour_brewer(palette = "Dark2")
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
    ggplot2::labs(
      title = "Chromatographic System Network",
      subtitle = sprintf("%d systems, %s of %s pairs drawn (the %d strongest links of each system)",
                         nrow(nodes), format(length(keep), big.mark = ","),
                         format(nrow(edges), big.mark = ","), n_links),
      caption = "Colour: groups of systems that share many compounds"
    )
}
