#' @title Statistics Generation Functions
#' @description Functions for computing performance statistics comparable to
#'   the original PredRet papers.
#' @name statistics
NULL

#' Calculate Model Performance Statistics
#'
#' Computes comprehensive performance metrics for all models, including
#' leave-one-out cross-validation errors when calibration data is available.
#'
#' @param model_results Output from `build_all_models()`.
#' @param datasets Named list of datasets from `load_report_data()`.
#'
#' @return A list with:
#'   \describe{
#'     \item{model_stats}{Tibble with per-model statistics}
#'     \item{overall_stats}{Named list with aggregate statistics}
#'     \item{by_method}{Tibble with stats grouped by method type}
#'   }
#'
#' @importFrom dplyr tibble bind_rows group_by summarize
#' @export
calculate_model_stats <- function(model_results, datasets = NULL) {

  models <- model_results$models

  if (length(models) == 0) {
    warning("No models provided")
    return(NULL)
  }

  # Calculate per-model statistics
  model_stats_list <- lapply(names(models), function(model_key) {
    model <- models[[model_key]]

    if (model$status != "success") {
      return(NULL)
    }

    # Get calibration points
    cal <- model$calibration
    n_points <- nrow(cal)

    # Find predicted values at calibration x-values
    pred_at_cal <- sapply(cal[, 1], function(x) {
      idx <- which.min(abs(x - model$newdata))
      model$ci[idx, 1]
    })

    # Calculate errors
    errors_abs <- abs(cal[, 2] - pred_at_cal)
    errors_rel <- errors_abs / cal[, 2] * 100  # As percentage

    # CI widths
    ci_lower_at_cal <- sapply(cal[, 1], function(x) {
      idx <- which.min(abs(x - model$newdata))
      model$ci[idx, 2]
    })
    ci_upper_at_cal <- sapply(cal[, 1], function(x) {
      idx <- which.min(abs(x - model$newdata))
      model$ci[idx, 3]
    })
    ci_widths <- ci_upper_at_cal - ci_lower_at_cal
    ci_widths_rel <- ci_widths / pred_at_cal * 100

    # R-squared
    ss_res <- sum((cal[, 2] - pred_at_cal)^2)
    ss_tot <- sum((cal[, 2] - mean(cal[, 2]))^2)
    r_squared <- 1 - (ss_res / ss_tot)

    # Correlation coefficient
    correlation <- cor(cal[, 1], cal[, 2])

    dplyr::tibble(
      model_key = model_key,
      from_id = model$sys1_id,
      to_id = model$sys2_id,
      n_compounds = n_points,
      r_squared = r_squared,
      correlation = correlation,
      # Absolute errors (minutes)
      mean_error_abs = mean(errors_abs),
      median_error_abs = median(errors_abs),
      sd_error_abs = sd(errors_abs),
      q95_error_abs = as.numeric(stats::quantile(errors_abs, 0.95)),
      max_error_abs = max(errors_abs),
      # Relative errors (%)
      mean_error_rel = mean(errors_rel),
      median_error_rel = median(errors_rel),
      sd_error_rel = sd(errors_rel),
      q95_error_rel = as.numeric(stats::quantile(errors_rel, 0.95)),
      max_error_rel = max(errors_rel),
      # CI widths (minutes)
      mean_ci_width = mean(ci_widths),
      median_ci_width = median(ci_widths),
      sd_ci_width = sd(ci_widths),
      q95_ci_width = as.numeric(stats::quantile(ci_widths, 0.95)),
      max_ci_width = max(ci_widths),
      # Relative CI widths (%)
      mean_ci_width_rel = mean(ci_widths_rel),
      median_ci_width_rel = median(ci_widths_rel),
      sd_ci_width_rel = sd(ci_widths_rel),
      q95_ci_width_rel = as.numeric(stats::quantile(ci_widths_rel, 0.95)),
      max_ci_width_rel = max(ci_widths_rel)
    )
  })

  model_stats <- dplyr::bind_rows(model_stats_list)

  # Aggregate statistics (as in papers)
  overall_stats <- list(
    n_models = nrow(model_stats),
    n_systems = length(unique(c(model_stats$from_id, model_stats$to_id))),
    total_calibration_points = sum(model_stats$n_compounds),
    # R-squared
    mean_r_squared = mean(model_stats$r_squared),
    median_r_squared = median(model_stats$r_squared),
    min_r_squared = min(model_stats$r_squared),
    # Absolute errors (weighted by n_compounds)
    weighted_mean_error_abs = weighted.mean(
      model_stats$mean_error_abs,
      model_stats$n_compounds
    ),
    overall_median_error_abs = median(model_stats$median_error_abs),
    overall_q95_error_abs = as.numeric(
      stats::quantile(model_stats$median_error_abs, 0.95)
    ),
    # Relative errors (weighted)
    weighted_mean_error_rel = weighted.mean(
      model_stats$mean_error_rel,
      model_stats$n_compounds
    ),
    overall_median_error_rel = median(model_stats$median_error_rel),
    overall_q95_error_rel = as.numeric(
      stats::quantile(model_stats$median_error_rel, 0.95)
    ),
    # CI widths
    weighted_mean_ci_width = weighted.mean(
      model_stats$mean_ci_width,
      model_stats$n_compounds
    ),
    overall_median_ci_width = median(model_stats$median_ci_width),
    overall_q95_ci_width = as.numeric(
      stats::quantile(model_stats$median_ci_width, 0.95)
    ),
    # Relative CI widths
    weighted_mean_ci_width_rel = weighted.mean(
      model_stats$mean_ci_width_rel,
      model_stats$n_compounds
    ),
    overall_median_ci_width_rel = median(model_stats$median_ci_width_rel),
    overall_q95_ci_width_rel = as.numeric(
      stats::quantile(model_stats$median_ci_width_rel, 0.95)
    )
  )

  # Add method type from model index if available
  if (!is.null(model_results$index) && "method_type" %in% names(model_results$index)) {
    model_stats <- dplyr::left_join(
      model_stats,
      model_results$index[, c("model_key", "method_type")],
      by = "model_key"
    )

    # Group by method type
    by_method <- model_stats |>
      dplyr::group_by(method_type) |>
      dplyr::summarize(
        n_models = dplyr::n(),
        n_systems = length(unique(c(from_id, to_id))),
        mean_n_compounds = mean(n_compounds),
        median_r_squared = median(r_squared),
        median_error_abs = median(median_error_abs),
        median_error_rel = median(median_error_rel),
        median_ci_width = median(median_ci_width),
        median_ci_width_rel = median(median_ci_width_rel),
        .groups = "drop"
      )
  } else {
    by_method <- NULL
  }

  list(
    model_stats = model_stats,
    overall_stats = overall_stats,
    by_method = by_method
  )
}


#' Aggregate Calibration Data for Regression Plot
#'
#' Collects all calibration point predictions vs actual values for
#' creating the aggregate regression plot (Figure 3B in paper).
#'
#' @param model_results Output from `build_all_models()`.
#'
#' @return A tibble with predicted and actual RT values for all
#'   calibration points across all models.
#'
#' @importFrom dplyr tibble bind_rows
#' @export
get_calibration_regression_data <- function(model_results) {

  models <- model_results$models

  regression_data_list <- lapply(names(models), function(model_key) {
    model <- models[[model_key]]

    if (model$status != "success") {
      return(NULL)
    }

    cal <- model$calibration

    # Find predicted values at calibration x-values
    pred_at_cal <- sapply(cal[, 1], function(x) {
      idx <- which.min(abs(x - model$newdata))
      model$ci[idx, 1]
    })

    dplyr::tibble(
      model_key = model_key,
      from_id = model$sys1_id,
      to_id = model$sys2_id,
      rt_source = cal[, 1],
      rt_target_actual = cal[, 2],
      rt_target_predicted = pred_at_cal,
      error_abs = abs(cal[, 2] - pred_at_cal),
      error_rel = abs(cal[, 2] - pred_at_cal) / cal[, 2] * 100
    )
  })

  dplyr::bind_rows(regression_data_list)
}


#' Calculate Database Statistics
#'
#' Computes statistics about the compound database coverage,
#' similar to Figure 3A in the paper.
#'
#' @param report_data Output from `load_report_data()`.
#'
#' @return A tibble with per-system statistics.
#'
#' @importFrom dplyr tibble
#' @export
calculate_database_stats <- function(report_data) {

  datasets <- report_data$datasets
  studies <- report_data$studies

  stats_list <- lapply(names(datasets), function(id) {
    ds <- datasets[[id]]

    n_compounds <- 0
    n_unique <- 0
    rt_min <- NA
    rt_max <- NA
    rt_range <- NA

    if (!is.null(ds$rtdata) && nrow(ds$rtdata) > 0) {
      n_compounds <- nrow(ds$rtdata)
      n_unique <- length(unique(ds$rtdata$inchi.std))
      rt_min <- min(ds$rtdata$rt, na.rm = TRUE)
      rt_max <- max(ds$rtdata$rt, na.rm = TRUE)
      rt_range <- rt_max - rt_min
    }

    # Get method type from studies
    method_type <- NA
    if (!is.null(studies) && id %in% studies$id) {
      method_type <- studies$method.type[studies$id == id]
    }

    dplyr::tibble(
      system_id = id,
      method_type = method_type,
      n_compounds = n_compounds,
      n_unique_compounds = n_unique,
      rt_min = rt_min,
      rt_max = rt_max,
      rt_range = rt_range
    )
  })

  dplyr::bind_rows(stats_list)
}


#' Calculate Prediction Statistics
#'
#' Computes statistics about predictions, if predictions have been generated.
#'
#' @param predictions Output from `generate_all_predictions()`.
#' @param model_results Output from `build_all_models()`.
#'
#' @return A list with prediction statistics.
#'
#' @importFrom dplyr tibble bind_rows
#' @export
calculate_prediction_stats <- function(predictions, model_results = NULL) {

  if (length(predictions) == 0) {
    return(NULL)
  }

  # Combine all predictions
  all_preds <- dplyr::bind_rows(predictions)

  if (nrow(all_preds) == 0) {
    return(NULL)
  }

  per_system <- lapply(names(predictions), function(target_id) {
    preds <- predictions[[target_id]]

    dplyr::tibble(
      target_system = target_id,
      n_predictions = nrow(preds),
      n_unique_compounds = length(unique(preds$compound_id)),
      mean_ci_width = mean(preds$ci_width),
      median_ci_width = median(preds$ci_width),
      q95_ci_width = as.numeric(stats::quantile(preds$ci_width, 0.95)),
      mean_n_calibration = mean(preds$n_calibration)
    )
  })

  per_system_df <- dplyr::bind_rows(per_system)

  overall <- list(
    total_predictions = nrow(all_preds),
    n_target_systems = length(predictions),
    n_unique_compounds = length(unique(all_preds$compound_id)),
    mean_ci_width = mean(all_preds$ci_width),
    median_ci_width = median(all_preds$ci_width),
    q95_ci_width = as.numeric(stats::quantile(all_preds$ci_width, 0.95))
  )

  list(
    per_system = per_system_df,
    overall = overall
  )
}


#' Export Statistics to Files
#'
#' Saves computed statistics to CSV and JSON files for the website.
#'
#' @param stats_list A list containing model_stats, database_stats, and
#'   prediction_stats from the calculate_* functions.
#' @param output_dir Directory to save statistics.
#' @param metadata Additional metadata (e.g., RepoRT version, run date).
#'
#' @return Invisibly returns the output directory path.
#'
#' @importFrom readr write_csv
#' @importFrom jsonlite write_json
#' @export
export_statistics <- function(stats_list, output_dir = "statistics",
                              metadata = NULL) {

  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }

  # Add timestamp if not provided
  if (is.null(metadata)) {
    metadata <- list()
  }
  metadata$generated_at <- as.character(Sys.time())

  # Export model statistics
  if (!is.null(stats_list$model_stats)) {
    ms <- stats_list$model_stats

    # Per-model CSV
    readr::write_csv(
      ms$model_stats,
      file.path(output_dir, "model_stats.csv")
    )

    # Overall stats as JSON
    overall <- ms$overall_stats
    overall$metadata <- metadata
    jsonlite::write_json(
      overall,
      file.path(output_dir, "overall_stats.json"),
      pretty = TRUE,
      auto_unbox = TRUE
    )

    # By method type if available
    if (!is.null(ms$by_method)) {
      readr::write_csv(
        ms$by_method,
        file.path(output_dir, "stats_by_method.csv")
      )
    }
  }

  # Export regression data
  if (!is.null(stats_list$regression_data)) {
    readr::write_csv(
      stats_list$regression_data,
      file.path(output_dir, "regression_data.csv")
    )
  }

  # Export database statistics
  if (!is.null(stats_list$database_stats)) {
    readr::write_csv(
      stats_list$database_stats,
      file.path(output_dir, "database_stats.csv")
    )
  }

  # Export prediction statistics
  if (!is.null(stats_list$prediction_stats)) {
    ps <- stats_list$prediction_stats

    if (!is.null(ps$per_system)) {
      readr::write_csv(
        ps$per_system,
        file.path(output_dir, "prediction_stats_per_system.csv")
      )
    }

    if (!is.null(ps$overall)) {
      jsonlite::write_json(
        ps$overall,
        file.path(output_dir, "prediction_stats_overall.json"),
        pretty = TRUE,
        auto_unbox = TRUE
      )
    }
  }

  message("Exported statistics to ", output_dir)
  invisible(output_dir)
}
