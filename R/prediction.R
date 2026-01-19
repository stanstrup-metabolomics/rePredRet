#' @title Prediction Functions
#' @description Functions for generating RT predictions from models.
#' @name prediction
NULL

# Default prediction parameters (matching original PredRet)
.prediction_defaults <- list(
  ci_width_limit = 2.0,       # Maximum CI width in minutes
  ci_width_limit_rel = 0.2,   # Maximum relative CI width (20%)
  density_bw_mult = 0.03,     # Bandwidth multiplier for density estimation
  density_limit = 0.01        # Minimum density required for prediction
)


#' Predict RT for a Single Compound in a Target System
#'
#' Predicts the retention time of a compound in a target chromatographic
#' system using available models.
#'
#' @param compound_id Compound identifier (InChI or InChIKey).
#' @param target_system_id Target system identifier.
#' @param models Named list of models from `build_all_models()`.
#' @param datasets Named list of datasets from `load_report_data()`.
#' @param ci_width_limit Maximum CI width in minutes (default 2.0).
#' @param ci_width_limit_rel Maximum relative CI width (default 0.2 = 20%).
#' @param density_bw_mult Bandwidth multiplier for density check (default 0.03).
#' @param density_limit Minimum density threshold (default 0.01).
#' @param id_column Column name for compound matching (default "inchi.std").
#'
#' @return A list with prediction details including provenance, or NULL
#'   if no valid prediction can be made.
#'
#' @details
#' The function:
#' 1. Finds all models that predict TO the target system
#' 2. Checks if the compound exists in any source system
#' 3. Makes predictions and applies quality filters:
#'    - CI width must be < ci_width_limit OR < ci_width_limit_rel
#'    - Source RT must be in a dense region of calibration data
#' 4. Returns the best prediction (narrowest relative CI)
#'
#' @importFrom stats density approx median
#' @export
predict_rt <- function(compound_id,
                       target_system_id,
                       models,
                       datasets,
                       ci_width_limit = .prediction_defaults$ci_width_limit,
                       ci_width_limit_rel = .prediction_defaults$ci_width_limit_rel,
                       density_bw_mult = .prediction_defaults$density_bw_mult,
                       density_limit = .prediction_defaults$density_limit,
                       id_column = "inchi.std") {

  # Find models that predict TO the target system
  available_models <- names(models)[
    grepl(paste0("_to_", target_system_id, "$"), names(models))
  ]

  if (length(available_models) == 0) {
    return(NULL)
  }

  # Clean compound ID once
  compound_id_clean <- gsub("/(t|b|m|s)[^/]*.*$", "", compound_id)

  predictions <- list()

  for (model_key in available_models) {
    model <- models[[model_key]]
    source_id <- model$sys1_id

    # Check if dataset exists
    if (is.null(datasets[[source_id]])) next

    # Check if compound exists in source system
    source_data <- datasets[[source_id]]$rtdata

    # Remove stereochemistry from InChI for matching
    source_inchi <- gsub("/(t|b|m|s)[^/]*.*$", "", source_data[[id_column]])

    compound_rows <- which(source_inchi == compound_id_clean)

    if (length(compound_rows) == 0) next

    # Get RT in source system (median if multiple measurements)
    source_rt <- stats::median(source_data$rt[compound_rows])

    # Check if source_rt is within model range
    if (source_rt < min(model$newdata) || source_rt > max(model$newdata)) {
      next
    }

    # Find closest x in model grid
    idx <- which.min(abs(source_rt - model$newdata))

    pred_rt <- model$ci[idx, 1]
    ci_lower <- model$ci[idx, 2]
    ci_upper <- model$ci[idx, 3]
    ci_width <- ci_upper - ci_lower

    # Apply CI width filters
    # Prediction is valid if EITHER absolute OR relative CI is acceptable
    if (ci_width >= ci_width_limit && (ci_width / pred_rt) >= ci_width_limit_rel) {
      next
    }

    # Density check: only predict where calibration data is dense
    dens <- stats::density(
      model$calibration[, 1],
      n = 4096,
      bw = density_bw_mult * max(model$calibration[, 1])
    )
    dens_at_source <- stats::approx(dens$x, dens$y, xout = source_rt)$y

    if (is.na(dens_at_source) || dens_at_source < density_limit) {
      next
    }

    # Get compound name if available
    compound_name <- NA
    if ("name" %in% names(source_data)) {
      compound_name <- source_data$name[compound_rows[1]]
    }

    # Store valid prediction
    predictions[[model_key]] <- list(
      compound_id = compound_id,
      compound_name = compound_name,
      source_system = source_id,
      source_rt = source_rt,
      predicted_rt = pred_rt,
      ci_lower = ci_lower,
      ci_upper = ci_upper,
      ci_width = ci_width,
      ci_width_rel = ci_width / pred_rt,
      model_key = model_key,
      n_calibration = model$n_points,
      model_median_error = model$stats$median_error
    )
  }

  if (length(predictions) == 0) {
    return(NULL)
  }

  # Select best prediction (narrowest relative CI)
  rel_widths <- sapply(predictions, function(p) p$ci_width_rel)
  best <- predictions[[which.min(rel_widths)]]

  best$target_system <- target_system_id

  return(best)
}


#' Generate All Predictions (Optimized Batch Processing)
#'
#' Generates predictions for all compounds in all systems where predictions
#' can be made. This version is optimized for batch processing.
#'
#' @param models Named list of models from `build_all_models()`.
#' @param datasets Named list of datasets from `load_report_data()`.
#' @param studies Studies tibble from `load_report_data()`.
#' @param report_version RepoRT version string for provenance tracking.
#' @param ci_width_limit Maximum CI width in minutes (default 2.0).
#' @param ci_width_limit_rel Maximum relative CI width (default 0.2 = 20%).
#' @param density_bw_mult Bandwidth multiplier for density check (default 0.03).
#' @param density_limit Minimum density threshold (default 0.01).
#'
#' @return A list with predictions organized by target system.
#'
#' @importFrom dplyr bind_rows tibble
#' @importFrom stats density approx median
#' @export
generate_all_predictions <- function(models,
                                     datasets,
                                     studies,
                                     report_version = NA,
                                     ci_width_limit = .prediction_defaults$ci_width_limit,
                                     ci_width_limit_rel = .prediction_defaults$ci_width_limit_rel,
                                     density_bw_mult = .prediction_defaults$density_bw_mult,
                                     density_limit = .prediction_defaults$density_limit) {

  message("Pre-processing datasets...")

  # Pre-process: strip stereochemistry from all datasets ONCE
  dataset_lookup <- list()
  for (ds_id in names(datasets)) {
    ds <- datasets[[ds_id]]
    if (is.null(ds$rtdata) || nrow(ds$rtdata) == 0) next

    # Strip stereochemistry and create lookup
    inchi_clean <- gsub("/(t|b|m|s)[^/]*.*$", "", ds$rtdata$inchi.std)

    # Aggregate RT by cleaned InChI (take median)
    unique_inchi <- unique(inchi_clean)
    rt_lookup <- vapply(unique_inchi, function(ic) {
      stats::median(ds$rtdata$rt[inchi_clean == ic])
    }, numeric(1))

    # Get names if available
    name_lookup <- NULL
    if ("name" %in% names(ds$rtdata) && !all(is.na(ds$rtdata$name))) {
      name_lookup <- vapply(unique_inchi, function(ic) {
        val <- ds$rtdata$name[which(inchi_clean == ic)[1]]
        if (is.null(val) || is.na(val)) NA_character_ else as.character(val)
      }, character(1))
    }

    dataset_lookup[[ds_id]] <- list(
      inchi = unique_inchi,
      rt = rt_lookup,
      name = name_lookup
    )
  }

  message("Pre-computing model density functions...")

  # Pre-compute density functions for each model
  model_density <- list()
  for (model_key in names(models)) {
    model <- models[[model_key]]
    if (model$status != "success") next

    dens <- stats::density(
      model$calibration[, 1],
      n = 4096,
      bw = density_bw_mult * max(model$calibration[, 1])
    )
    model_density[[model_key]] <- list(
      dens_fun = stats::approxfun(dens$x, dens$y, rule = 2),
      rt_min = min(model$newdata),
      rt_max = max(model$newdata)
    )
  }

  # Group models by target system
  target_systems <- unique(sapply(models, function(m) m$sys2_id))

  all_predictions <- list()

  for (target_id in target_systems) {
    message("Generating predictions for system ", target_id, "...")

    # Find models that predict TO this target
    target_models <- names(models)[
      sapply(models, function(m) m$sys2_id == target_id)
    ]

    if (length(target_models) == 0) next

    # Collect all predictions for this target
    pred_list <- list()

    for (model_key in target_models) {
      model <- models[[model_key]]
      source_id <- model$sys1_id

      if (is.null(dataset_lookup[[source_id]])) next

      src <- dataset_lookup[[source_id]]
      dens_info <- model_density[[model_key]]

      # Process all compounds in source system at once
      for (i in seq_along(src$inchi)) {
        inchi <- src$inchi[i]
        source_rt <- src$rt[i]

        # Skip if outside model range
        if (source_rt < dens_info$rt_min || source_rt > dens_info$rt_max) next

        # Find closest x in model grid
        idx <- which.min(abs(source_rt - model$newdata))

        pred_rt <- model$ci[idx, 1]
        ci_lower <- model$ci[idx, 2]
        ci_upper <- model$ci[idx, 3]
        ci_width <- ci_upper - ci_lower

        # Apply CI width filters
        if (ci_width >= ci_width_limit && (ci_width / pred_rt) >= ci_width_limit_rel) {
          next
        }

        # Density check
        dens_at_source <- dens_info$dens_fun(source_rt)
        if (is.na(dens_at_source) || dens_at_source < density_limit) {
          next
        }

        # Create prediction entry
        pred_list[[length(pred_list) + 1]] <- list(
          inchi = inchi,
          name = if (!is.null(src$name)) src$name[i] else NA_character_,
          source_id = source_id,
          source_rt = source_rt,
          pred_rt = pred_rt,
          ci_lower = ci_lower,
          ci_upper = ci_upper,
          ci_width = ci_width,
          ci_width_rel = ci_width / pred_rt,
          model_key = model_key,
          n_calibration = model$n_points,
          model_median_error = model$stats$median_error
        )
      }
    }

    if (length(pred_list) == 0) next

    # Convert to data frame
    pred_df <- data.frame(
      target_system = target_id,
      compound_id = vapply(pred_list, `[[`, character(1), "inchi"),
      compound_name = vapply(pred_list, `[[`, character(1), "name"),
      predicted_rt = vapply(pred_list, `[[`, numeric(1), "pred_rt"),
      ci_lower = vapply(pred_list, `[[`, numeric(1), "ci_lower"),
      ci_upper = vapply(pred_list, `[[`, numeric(1), "ci_upper"),
      ci_width = vapply(pred_list, `[[`, numeric(1), "ci_width"),
      source_system = vapply(pred_list, `[[`, character(1), "source_id"),
      source_rt = vapply(pred_list, `[[`, numeric(1), "source_rt"),
      model_key = vapply(pred_list, `[[`, character(1), "model_key"),
      n_calibration = vapply(pred_list, `[[`, numeric(1), "n_calibration"),
      model_median_error = vapply(pred_list, `[[`, numeric(1), "model_median_error"),
      report_version = report_version,
      generated_at = Sys.time(),
      stringsAsFactors = FALSE
    )

    # Select best prediction per compound (narrowest relative CI)
    pred_df$ci_width_rel <- pred_df$ci_width / pred_df$predicted_rt
    pred_df <- pred_df[order(pred_df$compound_id, pred_df$ci_width_rel), ]
    pred_df <- pred_df[!duplicated(pred_df$compound_id), ]
    pred_df$ci_width_rel <- NULL

    all_predictions[[target_id]] <- dplyr::as_tibble(pred_df)
    message("  -> ", nrow(pred_df), " predictions")
  }

  return(all_predictions)
}


#' Export Predictions to Files
#'
#' Saves predictions to CSV and JSON files with full provenance.
#'
#' @param predictions Output from `generate_all_predictions()`.
#' @param output_dir Directory to save predictions.
#' @param studies Studies tibble for system metadata.
#'
#' @return Invisibly returns the output directory path.
#'
#' @importFrom readr write_csv
#' @importFrom jsonlite write_json
#' @export
export_predictions <- function(predictions, output_dir = "predictions",
                               studies = NULL) {

  # Create directories
  by_system_dir <- file.path(output_dir, "by_system")
  by_compound_dir <- file.path(output_dir, "by_compound")
  summary_dir <- file.path(output_dir, "summary")

  for (dir in c(by_system_dir, by_compound_dir, summary_dir)) {
    if (!dir.exists(dir)) {
      dir.create(dir, recursive = TRUE)
    }
  }

  # Write per-system files
  for (target_id in names(predictions)) {
    pred_df <- predictions[[target_id]]

    # CSV
    readr::write_csv(
      pred_df,
      file.path(by_system_dir, paste0(target_id, "_predictions.csv"))
    )

    # JSON with more detail
    pred_json <- lapply(seq_len(nrow(pred_df)), function(i) {
      row <- pred_df[i, ]
      list(
        compound = list(
          id = row$compound_id,
          name = row$compound_name
        ),
        target_system = row$target_system,
        prediction = list(
          rt_minutes = row$predicted_rt,
          ci_lower = row$ci_lower,
          ci_upper = row$ci_upper,
          ci_width = row$ci_width
        ),
        provenance = list(
          source_system = row$source_system,
          source_rt = row$source_rt,
          model_key = row$model_key,
          n_calibration = row$n_calibration,
          model_median_error = row$model_median_error
        ),
        report_version = row$report_version,
        generated_at = as.character(row$generated_at)
      )
    })
    jsonlite::write_json(
      pred_json,
      file.path(by_system_dir, paste0(target_id, "_predictions.json")),
      pretty = TRUE,
      auto_unbox = TRUE
    )
  }

  # Combine all predictions for by-compound file
  all_preds <- dplyr::bind_rows(predictions)
  if (nrow(all_preds) > 0) {
    readr::write_csv(
      all_preds,
      file.path(by_compound_dir, "all_predictions.csv")
    )
  }

  # Write summary statistics
  if (nrow(all_preds) > 0) {
    summary_stats <- dplyr::tibble(
      target_system = names(predictions),
      n_predictions = sapply(predictions, nrow),
      mean_ci_width = sapply(predictions, function(df) mean(df$ci_width)),
      median_ci_width = sapply(predictions, function(df) median(df$ci_width))
    )
    readr::write_csv(
      summary_stats,
      file.path(summary_dir, "prediction_stats.csv")
    )
  }

  message("Exported predictions to ", output_dir)

  invisible(output_dir)
}
