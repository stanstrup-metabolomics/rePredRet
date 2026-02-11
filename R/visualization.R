#' @title Visualization Functions for User-Built RT Models
#' @description Functions for creating network plots, error distributions, and
#'   compound coverage heatmaps from UserRTModels and UserRTPredictions objects.
#' @name visualization
#' @keywords internal
NULL


#' @title Plot System Network
#' @description Creates a network graph showing connections between the user's
#'   chromatographic system and repository systems. Nodes represent systems,
#'   edges represent calibration models, and edge opacity/width reflect model
#'   quality (R-squared).
#'
#' @param models A \code{\link{UserRTModels}} object from
#'   \code{\link{build_user_models}}.
#' @param min_compounds Minimum number of calibration compounds for a model to
#'   be included in the network (default 10).
#' @param interactive If TRUE, returns an interactive plotly plot. If FALSE
#'   (default), returns a static ggplot2 plot.
#'
#' @return A ggplot2 object (if interactive = FALSE) or a plotly object
#'   (if interactive = TRUE).
#'
#' @details
#' The network uses a circular layout with the user's system ("USER") placed
#' at the center. Repository systems are arranged around the perimeter. Edge
#' width and opacity are proportional to the model's R-squared value, making
#' it easy to identify the strongest calibration relationships at a glance.
#'
#' Systems with fewer than \code{min_compounds} calibration compounds are
#' excluded from the visualization.
#'
#' @examples
#' \dontrun{
#' models <- build_user_models(my_inchis, my_rts, "RP")
#'
#' # Static network plot
#' plot_system_network(models)
#'
#' # Only show models with >= 20 calibration compounds
#' plot_system_network(models, min_compounds = 20)
#'
#' # Interactive version (requires plotly)
#' plot_system_network(models, interactive = TRUE)
#' }
#'
#' @seealso \code{\link{build_user_models}}, \code{\link{compare_systems}}
#'
#' @export
#' @importFrom ggplot2 ggplot aes geom_segment geom_point geom_text
#'   scale_alpha_continuous scale_linewidth_continuous scale_color_manual
#'   scale_size_manual labs theme_void theme element_text coord_equal
plot_system_network <- function(models,
                                min_compounds = 10,
                                interactive = FALSE) {

  if (!methods::is(models, "UserRTModels")) {
    stop("models must be a UserRTModels object")
  }

  diag <- models@diagnostics

  if (nrow(diag) == 0) {
    stop("No model diagnostics available. Cannot create network plot.")
  }

  # Filter by minimum calibration compounds
  diag <- diag[diag$n_calibration >= min_compounds, ]

  if (nrow(diag) == 0) {
    stop("No models meet the minimum compound threshold (min_compounds = ",
         min_compounds, "). Try lowering min_compounds.")
  }

  # Build node positions: USER at center, repo systems in a circle
  n_systems <- nrow(diag)
  angles <- seq(0, 2 * pi, length.out = n_systems + 1)[-(n_systems + 1)]

  radius <- 1.0
  node_data <- data.frame(
    system = c("USER", diag$source_system),
    x = c(0, radius * cos(angles)),
    y = c(0, radius * sin(angles)),
    type = c("user", rep("repository", n_systems)),
    n_calibration = c(NA, diag$n_calibration),
    r_squared = c(NA, diag$r_squared),
    stringsAsFactors = FALSE
  )

  # Build edge data
  edge_data <- data.frame(
    x_start = rep(0, n_systems),
    y_start = rep(0, n_systems),
    x_end = radius * cos(angles),
    y_end = radius * sin(angles),
    source_system = diag$source_system,
    r_squared = diag$r_squared,
    n_calibration = diag$n_calibration,
    rmse = diag$rmse,
    stringsAsFactors = FALSE
  )

  # Clamp R-squared for visual mapping (avoid issues with very low values)
  edge_data$r_sq_display <- pmax(edge_data$r_squared, 0)

  # Create plot
  p <- ggplot2::ggplot() +
    # Edges
    ggplot2::geom_segment(
      data = edge_data,
      ggplot2::aes(
        x = x_start, y = y_start,
        xend = x_end, yend = y_end,
        alpha = r_sq_display,
        linewidth = r_sq_display
      ),
      color = "steelblue"
    ) +
    ggplot2::scale_alpha_continuous(
      name = "R-squared",
      range = c(0.15, 0.9),
      limits = c(0, 1)
    ) +
    ggplot2::scale_linewidth_continuous(
      name = "R-squared",
      range = c(0.3, 2.5),
      limits = c(0, 1)
    ) +
    # Repository system nodes
    ggplot2::geom_point(
      data = node_data[node_data$type == "repository", ],
      ggplot2::aes(x = x, y = y),
      color = "steelblue",
      size = 3,
      alpha = 0.8
    ) +
    # User node (larger, distinct color)
    ggplot2::geom_point(
      data = node_data[node_data$type == "user", ],
      ggplot2::aes(x = x, y = y),
      color = "firebrick",
      size = 6
    ) +
    # User label
    ggplot2::geom_text(
      data = node_data[node_data$type == "user", ],
      ggplot2::aes(x = x, y = y, label = system),
      vjust = -1.2,
      fontface = "bold",
      size = 4.5,
      color = "firebrick"
    ) +
    # Repo system labels (positioned outside the circle)
    ggplot2::geom_text(
      data = node_data[node_data$type == "repository", ],
      ggplot2::aes(
        x = x * 1.12,
        y = y * 1.12,
        label = system
      ),
      size = 2.5,
      color = "gray30"
    ) +
    ggplot2::coord_equal() +
    ggplot2::labs(
      title = "System Network: User vs. Repository Systems",
      subtitle = paste0(
        n_systems, " systems | ",
        models@metadata$system_type, " | ",
        "Median R\u00b2 = ", round(stats::median(diag$r_squared, na.rm = TRUE), 3)
      )
    ) +
    ggplot2::theme_void() +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold", hjust = 0.5, size = 14),
      plot.subtitle = ggplot2::element_text(hjust = 0.5, color = "gray40",
                                            size = 11),
      legend.position = "right"
    )

  if (interactive) {
    if (!requireNamespace("plotly", quietly = TRUE)) {
      warning("Package 'plotly' is required for interactive plots. ",
              "Returning static ggplot.")
      return(p)
    }
    return(plotly::ggplotly(p))
  }

  return(p)
}


#' @title Plot Prediction Errors
#' @description Visualizes the distribution of confidence interval widths (and
#'   optionally model errors) from a set of RT predictions.
#'
#' @param predictions A \code{\link{UserRTPredictions}} object from
#'   \code{\link{predict_rt_user}}.
#' @param type Type of plot to create. One of \code{"density"} (default),
#'   \code{"histogram"}, or \code{"boxplot"}.
#' @param by_system If TRUE, facets the plot by source system (showing the top
#'   12 systems by number of predictions). Default FALSE.
#'
#' @return A ggplot2 object.
#'
#' @details
#' This function creates distribution plots of the confidence interval widths
#' from predictions. CI width is the primary measure of prediction uncertainty:
#' narrower intervals indicate more confident predictions.
#'
#' When \code{by_system = TRUE}, the plot is faceted by source system to reveal
#' which repository systems produce the tightest (or widest) predictions. Only
#' the top 12 systems (by prediction count) are shown to keep the plot readable.
#'
#' @examples
#' \dontrun{
#' predictions <- predict_rt_user(my_models)
#'
#' # Density plot of CI widths
#' plot_prediction_errors(predictions)
#'
#' # Histogram
#' plot_prediction_errors(predictions, type = "histogram")
#'
#' # Boxplot faceted by source system
#' plot_prediction_errors(predictions, type = "boxplot", by_system = TRUE)
#' }
#'
#' @seealso \code{\link{predict_rt_user}}, \code{\link{UserRTPredictions-class}}
#'
#' @export
#' @importFrom ggplot2 ggplot aes geom_density geom_histogram geom_boxplot
#'   geom_vline labs theme_bw facet_wrap coord_flip
plot_prediction_errors <- function(predictions,
                                   type = c("density", "histogram", "boxplot"),
                                   by_system = FALSE) {

  if (!methods::is(predictions, "UserRTPredictions")) {
    stop("predictions must be a UserRTPredictions object")
  }

  type <- match.arg(type)

  preds <- predictions@predictions

  if (nrow(preds) == 0) {
    stop("No predictions available. Cannot create error plot.")
  }

  # If faceting by system, limit to top 12 systems by count
  if (by_system) {
    system_counts <- as.data.frame(table(preds$source_system),
                                   stringsAsFactors = FALSE)
    names(system_counts) <- c("source_system", "count")
    system_counts <- system_counts[order(-system_counts$count), ]

    top_systems <- utils::head(system_counts$source_system, 12)
    preds <- preds[preds$source_system %in% top_systems, ]

    # Factor to preserve ordering
    preds$source_system <- factor(preds$source_system, levels = top_systems)
  }

  median_ci <- stats::median(preds$ci_width, na.rm = TRUE)

  # Build the base plot depending on type
  if (type == "density") {
    p <- ggplot2::ggplot(preds, ggplot2::aes(x = ci_width)) +
      ggplot2::geom_density(fill = "steelblue", alpha = 0.5, color = "steelblue") +
      ggplot2::geom_vline(
        xintercept = median_ci,
        color = "red", linetype = "dashed", linewidth = 0.8
      ) +
      ggplot2::labs(
        x = "Confidence Interval Width (minutes)",
        y = "Density",
        title = "Distribution of Prediction CI Widths",
        subtitle = paste0(
          nrow(preds), " predictions | ",
          "Median CI width: ", round(median_ci, 3), " min"
        )
      )

  } else if (type == "histogram") {
    p <- ggplot2::ggplot(preds, ggplot2::aes(x = ci_width)) +
      ggplot2::geom_histogram(
        bins = 40, fill = "steelblue", color = "white", alpha = 0.8
      ) +
      ggplot2::geom_vline(
        xintercept = median_ci,
        color = "red", linetype = "dashed", linewidth = 0.8
      ) +
      ggplot2::labs(
        x = "Confidence Interval Width (minutes)",
        y = "Number of Predictions",
        title = "Distribution of Prediction CI Widths",
        subtitle = paste0(
          nrow(preds), " predictions | ",
          "Median CI width: ", round(median_ci, 3), " min"
        )
      )

  } else if (type == "boxplot") {
    if (by_system) {
      p <- ggplot2::ggplot(preds, ggplot2::aes(
        x = source_system, y = ci_width
      )) +
        ggplot2::geom_boxplot(fill = "steelblue", alpha = 0.5,
                              outlier.alpha = 0.4, outlier.size = 1) +
        ggplot2::coord_flip() +
        ggplot2::labs(
          x = "Source System",
          y = "Confidence Interval Width (minutes)",
          title = "Prediction CI Widths by Source System",
          subtitle = paste0(
            "Top ", length(unique(preds$source_system)), " systems | ",
            "Overall median: ", round(median_ci, 3), " min"
          )
        )
    } else {
      p <- ggplot2::ggplot(preds, ggplot2::aes(x = "", y = ci_width)) +
        ggplot2::geom_boxplot(fill = "steelblue", alpha = 0.5,
                              outlier.alpha = 0.4, outlier.size = 1) +
        ggplot2::labs(
          x = NULL,
          y = "Confidence Interval Width (minutes)",
          title = "Distribution of Prediction CI Widths",
          subtitle = paste0(
            nrow(preds), " predictions | ",
            "Median CI width: ", round(median_ci, 3), " min"
          )
        )
    }
  }

  # Add faceting if by_system and not boxplot (boxplot handles it differently)
  if (by_system && type != "boxplot") {
    p <- p + ggplot2::facet_wrap(~source_system, scales = "free_y")
  }

  p <- p + ggplot2::theme_bw()

  return(p)
}


#' @title Plot Compound Coverage Heatmap
#' @description Creates a heatmap showing which compounds appear in which
#'   repository systems, colored by retention time value. This helps identify
#'   compounds with broad coverage across systems and systems with extensive
#'   compound libraries.
#'
#' @param models A \code{\link{UserRTModels}} object from
#'   \code{\link{build_user_models}}.
#' @param top_n Number of most-covered compounds to display (default 30).
#'   Compounds are ranked by the number of systems in which they appear.
#'
#' @return A ggplot2 object.
#'
#' @details
#' The heatmap displays compounds on the y-axis (labeled by truncated InChI
#' strings) and repository systems on the x-axis. Each tile is colored by the
#' retention time value in that system, using a viridis color scale. Missing
#' values (compound not measured in that system) are left blank.
#'
#' Only systems that have at least one of the top \code{top_n} compounds are
#' included on the x-axis. Compound labels are truncated InChI strings showing
#' the molecular formula layer for readability.
#'
#' @examples
#' \dontrun{
#' models <- build_user_models(my_inchis, my_rts, "RP")
#'
#' # Top 30 most-covered compounds
#' plot_compound_coverage(models)
#'
#' # Top 50 compounds
#' plot_compound_coverage(models, top_n = 50)
#' }
#'
#' @seealso \code{\link{build_user_models}}, \code{\link{UserRTModels-class}}
#'
#' @export
#' @importFrom ggplot2 ggplot aes geom_tile scale_fill_viridis_c labs theme_bw
#'   theme element_text element_blank
plot_compound_coverage <- function(models,
                                   top_n = 30) {

  if (!methods::is(models, "UserRTModels")) {
    stop("models must be a UserRTModels object")
  }

  report_data <- models@report_data

  if (is.null(report_data) || length(report_data$datasets) == 0) {
    stop("No report_data available in models object. ",
         "Cannot create coverage heatmap.")
  }

  studies <- report_data$studies
  datasets <- report_data$datasets

  # Collect all (compound, system, rt) triples from the datasets
  records <- list()

  for (sys_id in names(datasets)) {
    ds <- datasets[[sys_id]]
    if (is.null(ds$rtdata) || nrow(ds$rtdata) == 0) next

    inchi_clean <- gsub("/(t|b|m|s)[^/]*.*$", "", ds$rtdata$inchi.std)

    records[[sys_id]] <- data.frame(
      compound = inchi_clean,
      system = sys_id,
      rt = ds$rtdata$rt,
      stringsAsFactors = FALSE
    )
  }

  if (length(records) == 0) {
    stop("No RT data found in report_data datasets.")
  }

  all_data <- do.call(rbind, records)
  rownames(all_data) <- NULL

  # Aggregate: take median RT per compound per system
  all_data <- stats::aggregate(
    rt ~ compound + system,
    data = all_data,
    FUN = stats::median
  )

  # Count number of systems each compound appears in
  compound_counts <- as.data.frame(
    table(all_data$compound),
    stringsAsFactors = FALSE
  )
  names(compound_counts) <- c("compound", "n_systems")
  compound_counts <- compound_counts[order(-compound_counts$n_systems), ]

  # Select top_n compounds
  top_compounds <- utils::head(compound_counts$compound, top_n)

  # Filter data to top compounds
  plot_data <- all_data[all_data$compound %in% top_compounds, ]

  # Filter systems to only those containing at least one top compound
  systems_with_data <- unique(plot_data$system)
  plot_data <- plot_data[plot_data$system %in% systems_with_data, ]

  # Create truncated compound labels from InChI
  # Extract the molecular formula layer (second layer after InChI=1S/)
  compound_labels <- vapply(top_compounds, function(inchi) {
    parts <- strsplit(inchi, "/")[[1]]
    if (length(parts) >= 2) {
      # Return formula layer (e.g., "C8H10N4O2")
      return(parts[2])
    }
    # Fallback: truncate to 30 characters
    return(substr(inchi, 1, 30))
  }, character(1))

  # Handle duplicate labels by appending a suffix
  if (any(duplicated(compound_labels))) {
    dup_mask <- duplicated(compound_labels) | duplicated(compound_labels, fromLast = TRUE)
    dup_labels <- compound_labels[dup_mask]
    dup_inchis <- top_compounds[dup_mask]
    for (i in seq_along(dup_inchis)) {
      parts <- strsplit(dup_inchis[i], "/")[[1]]
      suffix <- if (length(parts) >= 3) {
        paste0("/", parts[3])
      } else {
        paste0(" [", i, "]")
      }
      compound_labels[top_compounds == dup_inchis[i]] <- paste0(
        compound_labels[top_compounds == dup_inchis[i]], suffix
      )
    }
  }

  # Build a lookup for compound -> label
  label_lookup <- stats::setNames(compound_labels, top_compounds)
  plot_data$compound_label <- label_lookup[plot_data$compound]

  # Order compounds by coverage count (most covered at top)
  compound_order <- rev(compound_labels)  # rev so most-covered is at top
  plot_data$compound_label <- factor(plot_data$compound_label,
                                     levels = compound_order)

  # Create the heatmap
  p <- ggplot2::ggplot(
    plot_data,
    ggplot2::aes(x = system, y = compound_label, fill = rt)
  ) +
    ggplot2::geom_tile(color = "white", linewidth = 0.3) +
    ggplot2::scale_fill_viridis_c(
      name = "RT (min)",
      option = "viridis",
      na.value = "grey95"
    ) +
    ggplot2::labs(
      x = "Repository System",
      y = "Compound (Formula)",
      title = "Compound Coverage Across Systems",
      subtitle = paste0(
        "Top ", min(top_n, length(top_compounds)), " compounds by system coverage | ",
        length(systems_with_data), " systems"
      )
    ) +
    ggplot2::theme_bw() +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 90, hjust = 1, vjust = 0.5,
                                           size = 7),
      axis.text.y = ggplot2::element_text(size = 7),
      plot.title = ggplot2::element_text(face = "bold"),
      panel.grid = ggplot2::element_blank()
    )

  return(p)
}
