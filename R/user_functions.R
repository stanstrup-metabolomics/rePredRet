#' @title User Convenience Functions
#' @description Easy-to-use functions for end users to access rePredRet
#'   predictions and models from the published GitHub repository.
#' @name user_functions
NULL

# GitHub repository for published data
.data_repo <- "stanstrup/rePredRet-data"
.data_repo_url <- paste0("https://github.com/", .data_repo)
.raw_url_base <- paste0("https://raw.githubusercontent.com/", .data_repo, "/main")


#' Download rePredRet Data
#'
#' Downloads the latest predictions and models from the rePredRet-data
#' GitHub repository.
#'
#' @param cache_dir Directory to cache downloaded data. Default uses
#'   a user-specific cache directory.
#' @param force_download If TRUE, re-download even if cached data exists.
#'
#' @return Path to the cached data directory.
#'
#' @export
#' @examples
#' \dontrun{
#' data_path <- rePredRet_download()
#' }
rePredRet_download <- function(cache_dir = NULL,
                                force_download = FALSE) {

  if (is.null(cache_dir)) {
    cache_dir <- rappdirs::user_cache_dir("rePredRet")
  }

  if (!dir.exists(cache_dir)) {
    dir.create(cache_dir, recursive = TRUE)
  }

  # Download as zip
  zip_url <- paste0(.data_repo_url, "/archive/refs/heads/main.zip")
  zip_path <- file.path(cache_dir, "rePredRet-data.zip")

  if (!file.exists(zip_path) || force_download) {
    message("Downloading rePredRet data...")
    if (!requireNamespace("curl", quietly = TRUE)) {
      stop("Package 'curl' required. Install with install.packages('curl')")
    }
    curl::curl_download(zip_url, zip_path)

    # Extract
    utils::unzip(zip_path, exdir = cache_dir, overwrite = TRUE)
  }

  # Return path to extracted data
  data_path <- file.path(cache_dir, "rePredRet-data-main")

  if (!dir.exists(data_path)) {
    stop("Download failed or extraction failed.")
  }

  return(data_path)
}


#' Get Path to Cached Data
#'
#' Returns the path to cached rePredRet data, downloading if necessary.
#'
#' @param data_path Optional path to rePredRet-data. If NULL, uses
#'   cached data (downloading if needed).
#'
#' @return Path to data directory.
#' @keywords internal
.get_data_path <- function(data_path = NULL) {
  if (is.null(data_path)) {
    cache_dir <- rappdirs::user_cache_dir("rePredRet")
    data_path <- file.path(cache_dir, "rePredRet-data-main")

    if (!dir.exists(data_path)) {
      data_path <- rePredRet_download()
    }
  }

  if (!dir.exists(data_path)) {
    stop("Data path does not exist: ", data_path)
  }

  return(data_path)
}


#' List Available Chromatographic Systems
#'
#' Returns information about all chromatographic systems with available
#' predictions.
#'
#' @param data_path Path to rePredRet-data (or NULL to use cache).
#'
#' @return Tibble with system information.
#'
#' @importFrom readr read_csv
#' @export
#' @examples
#' \dontrun{
#' systems <- rePredRet_systems()
#' }
rePredRet_systems <- function(data_path = NULL) {
  data_path <- .get_data_path(data_path)

  systems_file <- file.path(data_path, "metadata", "systems_info.csv")

  if (!file.exists(systems_file)) {
    stop("systems_info.csv not found. Data may be incomplete.")
  }

  readr::read_csv(systems_file, show_col_types = FALSE)
}


#' List Available Models
#'
#' Returns the model index with statistics for all available models.
#'
#' @param data_path Path to rePredRet-data (or NULL to use cache).
#'
#' @return Tibble with model index.
#'
#' @importFrom readr read_csv
#' @export
#' @examples
#' \dontrun{
#' models <- rePredRet_models()
#' }
rePredRet_models <- function(data_path = NULL) {
  data_path <- .get_data_path(data_path)

  index_file <- file.path(data_path, "models", "model_index.csv")

  if (!file.exists(index_file)) {
    stop("model_index.csv not found. Data may be incomplete.")
  }

  readr::read_csv(index_file, show_col_types = FALSE)
}


#' Get Model Details
#'
#' Retrieves a specific model with CI grid and calibration data.
#'
#' @param from_id Source system identifier.
#' @param to_id Target system identifier.
#' @param data_path Path to rePredRet-data (or NULL to use cache).
#'
#' @return A list with:
#'   \describe{
#'     \item{ci_grid}{Tibble with x, pred, lower, upper}
#'     \item{calibration_data}{Tibble with rt_source, rt_target, compound}
#'     \item{model}{RDS model object if available}
#'   }
#'
#' @importFrom readr read_csv
#' @export
#' @examples
#' \dontrun{
#' model_info <- rePredRet_model("0001", "0002")
#' }
rePredRet_model <- function(from_id, to_id, data_path = NULL) {
  data_path <- .get_data_path(data_path)

  model_key <- paste0(from_id, "_to_", to_id)
  model_dir <- file.path(data_path, "models", model_key)

  if (!dir.exists(model_dir)) {
    stop("Model not found: ", model_key)
  }

  result <- list()

  # Load CI grid
  ci_file <- file.path(model_dir, "ci_grid.csv")
  if (file.exists(ci_file)) {
    result$ci_grid <- readr::read_csv(ci_file, show_col_types = FALSE)
  }

  # Load calibration data
  cal_file <- file.path(model_dir, "calibration_data.csv")
  if (file.exists(cal_file)) {
    result$calibration_data <- readr::read_csv(cal_file, show_col_types = FALSE)
  }

  # Load RDS model (optional)
  rds_file <- file.path(model_dir, "model.rds")
  if (file.exists(rds_file)) {
    result$model <- readRDS(rds_file)
  }

  return(result)
}


#' Get Predictions for a Compound
#'
#' Retrieves all predictions for a compound across all systems.
#'
#' @param compound_id Compound identifier (InChI, InChIKey, or partial match).
#' @param target_system Optional: filter to specific target system.
#' @param data_path Path to rePredRet-data (or NULL to use cache).
#'
#' @return Tibble with predictions.
#'
#' @importFrom readr read_csv
#' @importFrom dplyr filter
#' @export
#' @examples
#' \dontrun{
#' preds <- rePredRet_predict("RYYVLZVUVIJVGH-UHFFFAOYSA-N")
#' }
rePredRet_predict <- function(compound_id,
                               target_system = NULL,
                               data_path = NULL) {
  data_path <- .get_data_path(data_path)

  all_preds_file <- file.path(data_path, "predictions", "by_compound",
                               "all_predictions.csv")

  if (!file.exists(all_preds_file)) {
    stop("Predictions file not found. Data may be incomplete.")
  }

  preds <- readr::read_csv(all_preds_file, show_col_types = FALSE)

  # Filter by compound (supports partial matching)
  preds <- dplyr::filter(
    preds,
    grepl(compound_id, compound_id, fixed = TRUE) |
    grepl(compound_id, compound_name, ignore.case = TRUE, fixed = TRUE)
  )

  # Filter by target system if specified
  if (!is.null(target_system)) {
    preds <- dplyr::filter(preds, target_system == !!target_system)
  }

  return(preds)
}


#' Get All Predictions for a System
#'
#' Retrieves all predictions for a specific chromatographic system.
#'
#' @param system_id Target system identifier.
#' @param data_path Path to rePredRet-data (or NULL to use cache).
#'
#' @return Tibble with predictions.
#'
#' @importFrom readr read_csv
#' @export
#' @examples
#' \dontrun{
#' preds <- rePredRet_system_predictions("0001")
#' }
rePredRet_system_predictions <- function(system_id, data_path = NULL) {
  data_path <- .get_data_path(data_path)

  pred_file <- file.path(data_path, "predictions", "by_system",
                          paste0(system_id, "_predictions.csv"))

  if (!file.exists(pred_file)) {
    stop("No predictions found for system: ", system_id)
  }

  readr::read_csv(pred_file, show_col_types = FALSE)
}


#' Plot Model Fit
#'
#' Creates a visualization of a model showing the fit curve, confidence
#' interval, and calibration points.
#'
#' @param from_id Source system identifier.
#' @param to_id Target system identifier.
#' @param data_path Path to rePredRet-data (or NULL to use cache).
#'
#' @return A ggplot object.
#'
#' @importFrom ggplot2 ggplot aes geom_ribbon geom_line geom_point labs theme_bw
#' @export
#' @examples
#' \dontrun{
#' p <- rePredRet_plot_model("0001", "0002")
#' print(p)
#' }
rePredRet_plot_model <- function(from_id, to_id, data_path = NULL) {
  model_info <- rePredRet_model(from_id, to_id, data_path)

  if (is.null(model_info$ci_grid) || is.null(model_info$calibration_data)) {
    stop("Model data incomplete")
  }

  ci_grid <- model_info$ci_grid
  calibration <- model_info$calibration_data

  p <- ggplot2::ggplot() +
    ggplot2::geom_ribbon(
      data = ci_grid,
      ggplot2::aes(x = x, ymin = lower, ymax = upper),
      fill = "lightblue",
      alpha = 0.5
    ) +
    ggplot2::geom_line(
      data = ci_grid,
      ggplot2::aes(x = x, y = pred),
      color = "blue",
      linewidth = 1
    ) +
    ggplot2::geom_point(
      data = calibration,
      ggplot2::aes(x = rt_source, y = rt_target),
      color = "black",
      alpha = 0.6
    ) +
    ggplot2::labs(
      x = paste("RT in", from_id, "(min)"),
      y = paste("RT in", to_id, "(min)"),
      title = paste("Model:", from_id, "->", to_id),
      subtitle = paste(nrow(calibration), "calibration compounds")
    ) +
    ggplot2::theme_bw()

  return(p)
}


#' Plot Prediction Error Distribution
#'
#' Creates a histogram or density plot of prediction errors for a system.
#'
#' @param system_id Target system identifier.
#' @param type Plot type: "histogram" or "density".
#' @param data_path Path to rePredRet-data (or NULL to use cache).
#'
#' @return A ggplot object.
#'
#' @importFrom ggplot2 ggplot aes geom_histogram geom_density labs theme_bw
#' @export
rePredRet_plot_errors <- function(system_id,
                                   type = c("histogram", "density"),
                                   data_path = NULL) {
  type <- match.arg(type)

  preds <- rePredRet_system_predictions(system_id, data_path)

  p <- ggplot2::ggplot(preds, ggplot2::aes(x = ci_width))

  if (type == "histogram") {
    p <- p + ggplot2::geom_histogram(
      bins = 30,
      fill = "steelblue",
      color = "white"
    )
  } else {
    p <- p + ggplot2::geom_density(fill = "steelblue", alpha = 0.5)
  }

  p <- p +
    ggplot2::labs(
      x = "CI Width (min)",
      y = if (type == "histogram") "Count" else "Density",
      title = paste("Prediction CI Width Distribution for", system_id),
      subtitle = paste(nrow(preds), "predictions")
    ) +
    ggplot2::theme_bw()

  return(p)
}


#' Plot System Network
#'
#' Creates a network visualization showing connections between chromatographic
#' systems based on model availability.
#'
#' @param min_compounds Minimum compounds to show edge (default 20).
#' @param data_path Path to rePredRet-data (or NULL to use cache).
#'
#' @return A ggplot object or network diagram.
#'
#' @importFrom ggplot2 ggplot aes geom_segment geom_point geom_text
#'   coord_fixed theme_void labs
#' @export
rePredRet_plot_network <- function(min_compounds = 20, data_path = NULL) {
  models <- rePredRet_models(data_path)

  # Filter by minimum compounds
  models <- models[models$n_compounds >= min_compounds, ]

  if (nrow(models) == 0) {
    stop("No models meet the minimum compound threshold")
  }

  # Get unique nodes
  nodes <- unique(c(models$from_id, models$to_id))
  n_nodes <- length(nodes)

  # Create circular layout
  angles <- seq(0, 2 * pi, length.out = n_nodes + 1)[-(n_nodes + 1)]
  node_coords <- data.frame(
    id = nodes,
    x = cos(angles),
    y = sin(angles)
  )

  # Create edge coordinates
  edges <- merge(models, node_coords, by.x = "from_id", by.y = "id")
  names(edges)[names(edges) == "x"] <- "x_from"
  names(edges)[names(edges) == "y"] <- "y_from"
  edges <- merge(edges, node_coords, by.x = "to_id", by.y = "id")
  names(edges)[names(edges) == "x"] <- "x_to"
  names(edges)[names(edges) == "y"] <- "y_to"

  p <- ggplot2::ggplot() +
    ggplot2::geom_segment(
      data = edges,
      ggplot2::aes(
        x = x_from, y = y_from,
        xend = x_to, yend = y_to,
        alpha = n_compounds
      ),
      color = "gray50"
    ) +
    ggplot2::geom_point(
      data = node_coords,
      ggplot2::aes(x = x, y = y),
      size = 4,
      color = "steelblue"
    ) +
    ggplot2::geom_text(
      data = node_coords,
      ggplot2::aes(x = x * 1.15, y = y * 1.15, label = id),
      size = 3
    ) +
    ggplot2::coord_fixed() +
    ggplot2::theme_void() +
    ggplot2::labs(
      title = "Chromatographic System Network",
      subtitle = paste(n_nodes, "systems,", nrow(edges), "model pairs"),
      alpha = "Common\nCompounds"
    )

  return(p)
}
