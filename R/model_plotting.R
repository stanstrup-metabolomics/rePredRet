#' @title Model Plotting Functions
#' @description Functions for creating interactive plots of RT prediction models.
#' @name model_plotting
NULL

#' Create Interactive Model Plot
#'
#' Creates an interactive scatter plot with prediction intervals using
#' ggplot2 and plotly.
#'
#' @param model A model object from `build_model()`.
#' @param dataset1 Optional dataset object for compound names (source system).
#' @param dataset2 Optional dataset object for compound names (target system).
#'
#' @return A plotly object (interactive plot).
#'
#' @importFrom ggplot2 ggplot aes geom_ribbon geom_line geom_point labs theme_minimal theme element_text
#' @importFrom plotly ggplotly
#' @export
plot_model <- function(model, dataset1 = NULL, dataset2 = NULL) {

  if (model$status != "success") {
    stop("Cannot plot model with status: ", model$status)
  }

  # Extract data
  calib <- model$calibration
  newdata <- model$newdata
  ci <- model$ci

  # Get compound names if available
  compound_names <- NULL
  if (!is.null(model$compounds)) {
    compound_names <- model$compounds
  } else if (!is.null(dataset1) && !is.null(dataset2)) {
    # Try to get names from datasets
    inchi1 <- gsub("/(t|b|m|s)[^/]*.*$", "", dataset1$rtdata$inchi.std)
    inchi2 <- gsub("/(t|b|m|s)[^/]*.*$", "", dataset2$rtdata$inchi.std)
    common <- intersect(inchi1, inchi2)

    if (length(common) == nrow(calib)) {
      compound_names <- sapply(common, function(ic) {
        idx <- which(inchi1 == ic)[1]
        name <- dataset1$rtdata$name[idx]
        if (is.null(name) || is.na(name) || name == "") {
          return(ic)
        }
        return(name)
      })
    }
  }

  # If still no names, use indices
  if (is.null(compound_names)) {
    compound_names <- paste0("Compound ", seq_len(nrow(calib)))
  }

  # Prepare data for plotting
  calib_df <- data.frame(
    x = calib[, 1],
    y = calib[, 2],
    name = compound_names,
    type = "Calibration"
  )

  pred_df <- data.frame(
    x = newdata,
    predicted = ci[, 1],
    lower = ci[, 2],
    upper = ci[, 3]
  )

  # Determine interval type based on width
  interval_type <- if (!is.null(model$method) && model$method == "posterior") {
    "95% Prediction Interval"
  } else {
    "95% Confidence Interval"
  }

  # Create plot
  p <- ggplot2::ggplot() +
    # Prediction interval (ribbon)
    ggplot2::geom_ribbon(
      data = pred_df,
      ggplot2::aes(x = x, ymin = lower, ymax = upper),
      fill = "lightblue",
      alpha = 0.3
    ) +
    # Predicted line
    ggplot2::geom_line(
      data = pred_df,
      ggplot2::aes(x = x, y = predicted),
      color = "blue",
      linewidth = 1
    ) +
    # Calibration points
    ggplot2::geom_point(
      data = calib_df,
      ggplot2::aes(x = x, y = y, label = name),
      color = "black",
      size = 2,
      alpha = 0.6
    ) +
    # Labels
    ggplot2::labs(
      title = paste0("RT Prediction Model: ", model$sys1_id, " → ", model$sys2_id),
      subtitle = paste0(model$n_points, " compounds, ", interval_type),
      x = paste0("RT in ", model$sys1_id, " (min)"),
      y = paste0("RT in ", model$sys2_id, " (min)")
    ) +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      plot.title = ggplot2::element_text(size = 14, face = "bold"),
      plot.subtitle = ggplot2::element_text(size = 11),
      axis.title = ggplot2::element_text(size = 12),
      axis.text = ggplot2::element_text(size = 10)
    )

  # Convert to interactive plotly
  p_interactive <- plotly::ggplotly(p, tooltip = c("x", "y", "label"))

  # Add custom hovertemplate for better tooltips
  p_interactive <- plotly::layout(
    p_interactive,
    hovermode = "closest",
    title = list(
      text = paste0("<b>", model$sys1_id, " → ", model$sys2_id, "</b><br>",
                   "<sub>", model$n_points, " compounds, ", interval_type, "</sub>"),
      x = 0.5,
      xanchor = "center"
    )
  )

  return(p_interactive)
}


#' Save Model Plot as Self-Contained HTML
#'
#' Saves an interactive model plot as a self-contained HTML file.
#'
#' @param model A model object from `build_model()`.
#' @param output_dir Directory to save the HTML file.
#' @param filename Optional filename (default: auto-generated from model IDs).
#' @param dataset1 Optional dataset object for compound names (source system).
#' @param dataset2 Optional dataset object for compound names (target system).
#'
#' @return Path to the saved HTML file (invisibly).
#'
#' @importFrom htmlwidgets saveWidget
#' @export
save_model_plot <- function(model,
                            output_dir,
                            filename = NULL,
                            dataset1 = NULL,
                            dataset2 = NULL) {

  if (model$status != "success") {
    warning("Skipping plot for model with status: ", model$status)
    return(invisible(NULL))
  }

  # Create output directory if needed
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }

  # Generate filename if not provided
  if (is.null(filename)) {
    filename <- paste0(model$sys1_id, "_to_", model$sys2_id, ".html")
  }

  output_path <- file.path(output_dir, filename)

  # Create plot
  p <- plot_model(model, dataset1, dataset2)

  # Save as self-contained HTML
  htmlwidgets::saveWidget(
    p,
    file = output_path,
    selfcontained = TRUE,
    title = paste0("RT Model: ", model$sys1_id, " → ", model$sys2_id)
  )

  return(invisible(output_path))
}


#' Create Model Plot Data for Batch Processing
#'
#' Prepares plot data that can be saved alongside model files.
#'
#' @param model A model object from `build_model()`.
#'
#' @return A list with plot data components.
#'
#' @keywords internal
extract_plot_data <- function(model) {
  if (model$status != "success") {
    return(NULL)
  }

  list(
    calibration = model$calibration,
    newdata = model$newdata,
    ci = model$ci,
    sys1_id = model$sys1_id,
    sys2_id = model$sys2_id,
    n_points = model$n_points,
    method = model$method
  )
}
