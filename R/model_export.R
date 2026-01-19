#' Export Model Data as JSON for Web Viewer
#'
#' Exports model predictions and calibration data as JSON for efficient
#' web-based visualization with Plotly.
#'
#' @param model A model object from build_model()
#' @param output_file Path to save JSON file
#'
#' @return Invisibly returns the JSON string
#' @keywords internal
export_model_json <- function(model, output_file) {

  # Extract core data
  model_data <- list(
    sys1_id = model$sys1_id,
    sys2_id = model$sys2_id,
    n_compounds = model$n_points,
    method = model$method,

    # Prediction grid (1000 points)
    predictions = data.frame(
      x = model$newdata,
      pred = model$ci[, 1],
      lower = model$ci[, 2],
      upper = model$ci[, 3]
    ),

    # Calibration data
    calibration = data.frame(
      x = model$calibration[, 1],
      y = model$calibration[, 2],
      compound = model$compounds %||% NA_character_
    ),

    # Statistics
    stats = list(
      median_pi_width = model$stats$median_pi_width,
      median_error = model$stats$median_error,
      build_time = as.character(model$build_time)
    )
  )

  # Convert to JSON (compact format to minimize file size)
  json_str <- jsonlite::toJSON(model_data,
                               digits = 4,
                               na = "null",
                               auto_unbox = TRUE,
                               pretty = FALSE)

  # Write to file
  writeLines(json_str, output_file)

  # Return file size
  file_size_kb <- file.size(output_file) / 1024

  invisible(list(
    file = output_file,
    size_kb = file_size_kb,
    json = json_str
  ))
}


#' Fast Model Export (No Plots)
#'
#' Exports model to RDS + CSV + JSON without generating HTML plots.
#' Much faster for large-scale pipelines.
#'
#' @param model A model object
#' @param model_dir Directory for output files
#'
#' @return Invisibly returns paths to created files
#' @keywords internal
export_model_fast <- function(model, model_dir) {

  if (!dir.exists(model_dir)) {
    dir.create(model_dir, recursive = TRUE)
  }

  files_created <- list()

  # 1. Save model as RDS (stripped of large objects)
  model_for_save <- model
  model_for_save$boot_object <- NULL
  saveRDS(model_for_save, file.path(model_dir, "model.rds"))
  files_created$model_rds <- file.path(model_dir, "model.rds")

  # 2. Save prediction grid as CSV
  ci_df <- data.frame(
    rt_source = model$newdata,
    rt_target_predicted = model$ci[, 1],
    ci_lower = model$ci[, 2],
    ci_upper = model$ci[, 3]
  )
  readr::write_csv(ci_df, file.path(model_dir, "ci_grid.csv"))
  files_created$ci_csv <- file.path(model_dir, "ci_grid.csv")

  # 3. Save calibration data as CSV
  calib_df <- data.frame(
    rt_source = model$calibration[, 1],
    rt_target = model$calibration[, 2]
  )
  if (!is.null(model$compounds)) {
    calib_df$compound_id <- model$compounds
  }
  readr::write_csv(calib_df, file.path(model_dir, "calibration_data.csv"))
  files_created$calib_csv <- file.path(model_dir, "calibration_data.csv")

  # 4. Save as JSON for web viewer
  export_model_json(model, file.path(model_dir, "model.json"))
  files_created$json <- file.path(model_dir, "model.json")

  invisible(files_created)
}
