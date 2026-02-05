#' @title Plot a Pre-built Model from the Data Repository
#' @description Downloads and plots a calibration model showing the relationship
#'   between retention times in two chromatographic systems.
#'
#' @param from_id Source system ID (e.g., "0001").
#' @param to_id Target system ID (e.g., "0002").
#' @param interactive If TRUE, returns a plotly interactive plot. If FALSE (default),
#'   returns a ggplot2 static plot.
#' @param show_ci If TRUE (default), shows confidence interval bands.
#' @param show_calibration If TRUE (default), shows calibration points.
#' @param data_path Optional path to local rePredRet-data directory. If NULL,
#'   downloads/uses cached data.
#'
#' @return A ggplot2 object (if interactive = FALSE) or a plotly object
#'   (if interactive = TRUE).
#'
#' @details
#' This function fetches model data from the rePredRet-models repository and
#' creates a publication-ready plot showing:
#' \itemize{
#'   \item The calibration curve (monotonic GAM fit)
#'   \item Confidence interval bands
#'   \item Calibration points used to build the model
#' }
#'
#' @examples
#' \dontrun{
#' # Static ggplot
#' plot_model("0001", "0002")
#'
#' # Interactive plotly
#' plot_model("0001", "0002", interactive = TRUE)
#'
#' # Without CI bands
#' plot_model("0001", "0002", show_ci = FALSE)
#' }
#'
#' @seealso \code{\link{rePredRet_model}} for getting model data,
#'   \code{\link{rePredRet_systems}} for listing available systems
#'
#' @export
#' @importFrom ggplot2 ggplot aes geom_ribbon geom_line geom_point labs theme_bw
#'   scale_x_continuous scale_y_continuous
plot_model <- function(from_id,
                       to_id,
                       interactive = FALSE,
                       show_ci = TRUE,
                       show_calibration = TRUE,
                       data_path = NULL) {

  # Get data path
  if (is.null(data_path)) {
    data_path <- .get_data_path()
  }

  model_dir <- file.path(data_path, "models", paste0(from_id, "_to_", to_id))

  if (!dir.exists(model_dir)) {
    stop("Model not found: ", from_id, " -> ", to_id,
         "\nUse rePredRet_models() to see available models.")
  }

  # Load CI grid
  ci_file <- file.path(model_dir, "ci_grid.csv")
  if (!file.exists(ci_file)) {
    stop("CI grid file not found for model ", from_id, " -> ", to_id)
  }
  ci_grid <- utils::read.csv(ci_file, stringsAsFactors = FALSE)

  # Load calibration data
  calib_file <- file.path(model_dir, "calibration_data.csv")
  calibration <- if (file.exists(calib_file)) {
    utils::read.csv(calib_file, stringsAsFactors = FALSE)
  } else {
    NULL
  }

  # Get system names for labels
  systems_file <- file.path(data_path, "metadata", "systems_info.csv")
  if (file.exists(systems_file)) {
    systems <- utils::read.csv(systems_file, stringsAsFactors = FALSE)
    from_name <- systems$name[systems$id == from_id]
    to_name <- systems$name[systems$id == to_id]
    if (length(from_name) == 0) from_name <- from_id
    if (length(to_name) == 0) to_name <- to_id
  } else {
    from_name <- from_id
    to_name <- to_id
  }

  # Load model index for statistics
  index_file <- file.path(data_path, "models", "model_index.csv")
  model_stats <- NULL
  if (file.exists(index_file)) {
    index <- utils::read.csv(index_file, stringsAsFactors = FALSE)
    model_stats <- index[index$from_id == from_id & index$to_id == to_id, ]
  }

  # Build subtitle with model stats
  subtitle <- NULL
  if (!is.null(model_stats) && nrow(model_stats) > 0) {
    subtitle <- paste0(
      "n = ", model_stats$n_compounds, " compounds | ",
      "Median error: ", round(model_stats$median_error, 2), " min | ",
      "Median CI: ", round(model_stats$median_ci_width, 2), " min"
    )
  }

  # Create ggplot
  p <- ggplot2::ggplot()

  # Add CI ribbon
  if (show_ci && "lower" %in% names(ci_grid) && "upper" %in% names(ci_grid)) {
    p <- p + ggplot2::geom_ribbon(
      data = ci_grid,
      ggplot2::aes(x = x, ymin = lower, ymax = upper),
      fill = "steelblue",
      alpha = 0.3
    )
  }

  # Add prediction line
  p <- p + ggplot2::geom_line(
    data = ci_grid,
    ggplot2::aes(x = x, y = pred),
    color = "steelblue",
    linewidth = 1
  )

  # Add calibration points
  if (show_calibration && !is.null(calibration)) {
    p <- p + ggplot2::geom_point(
      data = calibration,
      ggplot2::aes(x = rt_source, y = rt_target),
      color = "black",
      size = 2,
      alpha = 0.7
    )
  }

  # Add labels
  p <- p +
    ggplot2::labs(
      title = paste0("RT Model: ", from_name, " \u2192 ", to_name),
      subtitle = subtitle,
      x = paste0("RT in ", from_name, " (min)"),
      y = paste0("RT in ", to_name, " (min)")
    ) +
    ggplot2::theme_bw() +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold"),
      plot.subtitle = ggplot2::element_text(color = "gray40")
    )

  # Return interactive or static

  if (interactive) {
    if (!requireNamespace("plotly", quietly = TRUE)) {
      warning("Package 'plotly' is required for interactive plots. Returning ggplot.")
      return(p)
    }
    return(plotly::ggplotly(p))
  }

  return(p)
}


#' @title Search for a Compound Across All Systems
#' @description Searches for a compound by InChI and returns all systems where
#'   it has been measured. Stereochemistry is automatically stripped for matching.
#'
#' @param inchi InChI string to search for. Stereochemistry layers will be
#'   automatically removed for matching.
#' @param system_type Optional filter by system type ("RP" or "HILIC").
#' @param data_path Optional path to local RepoRT data directory.
#'
#' @return A data.frame with columns: system_id, system_name, system_type,
#'   compound_name, rt, inchi.
#'
#' @details
#' The search automatically sanitizes InChI strings by removing stereochemistry
#' layers (/t, /b, /m, /s) before matching, ensuring consistent results regardless
#' of whether the query or database entries include stereochemistry information.
#'
#' @examples
#' \dontrun{
#' # Search by InChI (caffeine)
#' search_compound("InChI=1S/C8H10N4O2/c1-10-4-9-6-5(10)7(13)12(3)8(14)11(6)2/h4H,1-3H3")
#'
#' # Partial InChI also works
#' search_compound("InChI=1S/C8H10N4O2")
#'
#' # Filter to RP systems only
#' search_compound("InChI=1S/C15H14O6", system_type = "RP")
#' }
#'
#' @export
search_compound <- function(inchi,
                            system_type = NULL,
                            data_path = NULL) {

  # Validate input
 if (!grepl("^InChI=", inchi, ignore.case = TRUE)) {
    stop("Query must be an InChI string (starting with 'InChI=')")
  }

  # Get data path
  if (is.null(data_path)) {
    # Try to use RepoRT data if available
    report_path <- file.path(getwd(), "RepoRT_data", "RepoRT-master", "processed_data")
    if (!dir.exists(report_path)) {
      report_path <- tryCatch(download_report(), error = function(e) NULL)
    }
    if (is.null(report_path) || !dir.exists(report_path)) {
      stop("RepoRT data not found. Use download_report() first.")
    }
    data_path <- report_path
  }

  # Load studies
  studies_file <- file.path(data_path, "studies.tsv")
  if (!file.exists(studies_file)) {
    stop("Studies file not found at: ", studies_file)
  }
  studies <- utils::read.delim(studies_file, stringsAsFactors = FALSE)

  # Filter by system type if specified
  if (!is.null(system_type)) {
    studies <- studies[studies$method.type == system_type, ]
  }

  if (nrow(studies) == 0) {
    return(data.frame())
  }

  # Sanitize query InChI (remove stereochemistry)
  query_clean <- gsub("/(t|b|m|s)[^/]*.*$", "", inchi)
  query_clean_lower <- tolower(query_clean)

  results <- list()

  for (sys_id in studies$id) {
    rt_file <- file.path(data_path, sys_id,
                         paste0(sys_id, "_rtdata_canonical_success.tsv"))
    if (!file.exists(rt_file)) next

    rtdata <- tryCatch(
      utils::read.delim(rt_file, stringsAsFactors = FALSE),
      error = function(e) NULL
    )
    if (is.null(rtdata) || nrow(rtdata) == 0) next

    # Sanitize database InChIs (remove stereochemistry)
    rtdata_inchi_clean <- gsub("/(t|b|m|s)[^/]*.*$", "", rtdata$inchi.std)

    # Match InChI (case-insensitive, partial match supported)
    matches <- grepl(query_clean_lower, tolower(rtdata_inchi_clean), fixed = TRUE)

    if (any(matches)) {
      matched_data <- rtdata[matches, ]
      sys_info <- studies[studies$id == sys_id, ]

      for (i in seq_len(nrow(matched_data))) {
        results[[length(results) + 1]] <- data.frame(
          system_id = sys_id,
          system_name = sys_info$name[1],
          system_type = sys_info$method.type[1],
          compound_name = if ("name" %in% names(matched_data)) matched_data$name[i] else NA,
          rt = matched_data$rt[i],
          inchi = matched_data$inchi.std[i],
          stringsAsFactors = FALSE
        )
      }
    }
  }

  if (length(results) == 0) {
    message("No compounds found matching '", inchi, "'")
    return(data.frame(
      system_id = character(),
      system_name = character(),
      system_type = character(),
      compound_name = character(),
      rt = numeric(),
      inchi = character(),
      stringsAsFactors = FALSE
    ))
  }

  result_df <- do.call(rbind, results)
  result_df <- result_df[order(result_df$compound_name, result_df$system_id), ]
  rownames(result_df) <- NULL

  return(result_df)
}


#' @title Compare User's System to Repository Systems
#' @description Ranks repository systems by how well their models fit the user's data,
#'   helping identify the most similar systems.
#'
#' @param models A UserRTModels object from \code{\link{build_user_models}}.
#' @param metric Metric to rank by: "r_squared" (default), "rmse", "median_ci_width",
#'   or "n_calibration".
#' @param top_n Number of top systems to return (default: all).
#'
#' @return A data.frame ranked by the specified metric, with columns for all
#'   diagnostic metrics.
#'
#' @examples
#' \dontrun{
#' models <- build_user_models(my_inchis, my_rts, "RP")
#'
#' # Find most similar systems (best R²)
#' compare_systems(models)
#'
#' # Rank by lowest RMSE
#' compare_systems(models, metric = "rmse")
#'
#' # Top 5 systems
#' compare_systems(models, top_n = 5)
#' }
#'
#' @export
compare_systems <- function(models,
                            metric = c("r_squared", "rmse", "median_ci_width", "n_calibration"),
                            top_n = NULL) {

  if (!methods::is(models, "UserRTModels")) {
    stop("models must be a UserRTModels object")
  }

  metric <- match.arg(metric)

  diag <- models@diagnostics

  if (nrow(diag) == 0) {
    warning("No model diagnostics available")
    return(data.frame())
  }

  # Sort by metric
  if (metric == "r_squared") {
    diag <- diag[order(-diag$r_squared), ]
  } else if (metric == "rmse") {
    diag <- diag[order(diag$rmse), ]
  } else if (metric == "median_ci_width") {
    diag <- diag[order(diag$median_ci_width), ]
  } else if (metric == "n_calibration") {
    diag <- diag[order(-diag$n_calibration), ]
  }

  # Limit to top_n

if (!is.null(top_n) && top_n < nrow(diag)) {
    diag <- diag[1:top_n, ]
  }

  # Add rank column
  diag$rank <- seq_len(nrow(diag))
  diag <- diag[, c("rank", names(diag)[names(diag) != "rank"])]

  rownames(diag) <- NULL
  return(diag)
}


#' @title Get Detailed System Information
#' @description Retrieves detailed information about a chromatographic system
#'   including column, gradient, and metadata.
#'
#' @param system_id System identifier (e.g., "0001").
#' @param data_path Optional path to local RepoRT data directory.
#'
#' @return A list with components: id, name, type, n_compounds, info, metadata, gradient.
#'
#' @examples
#' \dontrun{
#' # Get info about system 0001
#' info <- get_system_info("0001")
#' info$name
#' info$n_compounds
#' info$metadata
#' }
#'
#' @export
get_system_info <- function(system_id, data_path = NULL) {

  # Get data path
  if (is.null(data_path)) {
    report_path <- file.path(getwd(), "RepoRT_data", "RepoRT-master", "processed_data")
    if (!dir.exists(report_path)) {
      report_path <- tryCatch(download_report(), error = function(e) NULL)
    }
    if (is.null(report_path) || !dir.exists(report_path)) {
      stop("RepoRT data not found. Use download_report() first.")
    }
    data_path <- report_path
  }

  sys_dir <- file.path(data_path, system_id)
  if (!dir.exists(sys_dir)) {
    stop("System not found: ", system_id)
  }

  # Load studies for basic info
  studies_file <- file.path(data_path, "studies.tsv")
  studies <- utils::read.delim(studies_file, stringsAsFactors = FALSE)
  study_info <- studies[studies$id == system_id, ]

  if (nrow(study_info) == 0) {
    stop("System ", system_id, " not found in studies.tsv")
  }

  # Load info file
  info_file <- file.path(sys_dir, paste0(system_id, "_info.tsv"))
  info <- if (file.exists(info_file)) {
    utils::read.delim(info_file, stringsAsFactors = FALSE)
  } else {
    NULL
  }

  # Load metadata file
  meta_file <- file.path(sys_dir, paste0(system_id, "_metadata.tsv"))
  metadata <- if (file.exists(meta_file)) {
    utils::read.delim(meta_file, stringsAsFactors = FALSE)
  } else {
    NULL
  }

  # Load gradient file
  grad_file <- file.path(sys_dir, paste0(system_id, "_gradient.tsv"))
  gradient <- if (file.exists(grad_file)) {
    utils::read.delim(grad_file, stringsAsFactors = FALSE)
  } else {
    NULL
  }

  # Count compounds
  rt_file <- file.path(sys_dir, paste0(system_id, "_rtdata_canonical_success.tsv"))
  n_compounds <- if (file.exists(rt_file)) {
    rtdata <- utils::read.delim(rt_file, stringsAsFactors = FALSE)
    length(unique(rtdata$inchi.std))
  } else {
    NA
  }

  result <- list(
    id = system_id,
    name = study_info$name,
    type = study_info$method.type,
    url = study_info$url,
    pmid = study_info$pmid,
    authors = study_info$authors,
    n_compounds = n_compounds,
    info = info,
    metadata = metadata,
    gradient = gradient
  )

  class(result) <- c("system_info", "list")
  return(result)
}


#' @export
print.system_info <- function(x, ...) {
  cat("System:", x$id, "-", x$name, "\n")
  cat("Type:", x$type, "\n")
  cat("Compounds:", x$n_compounds, "\n")
  if (!is.null(x$url) && x$url != "") {
    cat("URL:", x$url, "\n")
  }
  if (!is.null(x$authors) && x$authors != "") {
    cat("Authors:", x$authors, "\n")
  }
  if (!is.null(x$metadata) && nrow(x$metadata) > 0) {
    cat("\nMetadata:\n")
    for (i in seq_len(nrow(x$metadata))) {
      cat("  ", x$metadata[i, 1], ": ", x$metadata[i, 2], "\n", sep = "")
    }
  }
  invisible(x)
}
