#' @title Generate RT Predictions from User Models
#' @description Generates retention time predictions using models built with
#'   \code{\link{build_user_models}}.
#'
#' @param models A UserRTModels object from \code{\link{build_user_models}}.
#' @param compounds Optional character vector of InChI identifiers to predict.
#'   If NULL (default), predicts all compounds from all source systems.
#' @param verbose Print progress messages (default TRUE).
#'
#' @return A UserRTPredictions S4 object containing:
#'   \describe{
#'     \item{predictions}{Data frame of all predictions with CIs and model quality}
#'     \item{source_models}{Character vector of model system IDs used}
#'     \item{metadata}{Run information and summary statistics}
#'   }
#'
#' @details
#' The function:
#' \enumerate{
#'   \item Takes each model from the UserRTModels object
#'   \item For each source system, retrieves all measured compounds
#'   \item Generates predictions by interpolating from the model's CI grid
#'   \item Selects the best prediction per compound (narrowest CI width)
#' }
#'
#' When a compound is measured in multiple source systems, only the prediction
#' with the smallest confidence interval width is returned.
#'
#' @examples
#' \dontrun{
#' # First build models
#' my_models <- build_user_models(
#'   inchis = my_inchis,
#'   rts = my_rts,
#'   system_type = "RP"
#' )
#'
#' # Generate all predictions
#' predictions <- predict_rt_user(my_models)
#'
#' # View summary
#' show(predictions)
#'
#' # Extract as data.frame
#' preds_df <- as.data.frame(predictions)
#'
#' # Filter to high-quality predictions
#' good_preds <- preds_df[preds_df$ci_width < 1.0 & preds_df$model_r_squared > 0.95, ]
#'
#' # Or predict specific compounds only
#' caffeine_preds <- predict_rt_user(
#'   my_models,
#'   compounds = c("InChI=1S/C8H10N4O2/c1-10-4-9-6-5(10)7(13)12(3)8(14)11(6)2/h4H,1-3H3")
#' )
#' }
#'
#' @seealso \code{\link{build_user_models}} for building models,
#'   \code{\link{UserRTPredictions-class}} for the return object class
#'
#' @export
#' @importFrom methods is new
#' @importFrom cli cli_rule cli_alert_info cli_alert_success cli_progress_bar
#'   cli_progress_update cli_progress_done
predict_rt_user <- function(models,
                             compounds = NULL,
                             verbose = TRUE) {

  # Validate input
 if (!methods::is(models, "UserRTModels")) {
    stop("models must be a UserRTModels object from build_user_models()")
  }

  if (length(models@models) == 0) {
    warning("No models available in the UserRTModels object")
    return(methods::new(
      "UserRTPredictions",
      predictions = data.frame(),
      source_models = character(),
      metadata = list(
        n_predictions = 0,
        n_unique_compounds = 0,
        n_models_used = 0,
        created_at = Sys.time()
      )
    ))
  }

  start_time <- Sys.time()

  if (verbose) {
    cli::cli_rule("predict_rt_user")
    cli::cli_alert_info("Models available: {length(models@models)}")
    if (!is.null(compounds)) {
      cli::cli_alert_info("Filtering to {length(compounds)} specific compounds")
    }
  }

  # Standardize compound filter if provided
  if (!is.null(compounds)) {
    compounds_clean <- gsub("/(t|b|m|s)[^/]*.*$", "", compounds)
  }

  # Get report data from the models object
  report_data <- models@report_data
  diagnostics <- models@diagnostics

  predictions_list <- list()
  n_models <- length(models@models)

  if (verbose) {
    cli::cli_progress_bar(
      "Generating predictions",
      total = n_models,
      format = "{cli::pb_spin} Generating predictions [{cli::pb_current}/{cli::pb_total}] | {cli::pb_bar} {cli::pb_percent} | ETA: {cli::pb_eta}"
    )
  }

  for (sys_id in names(models@models)) {
    model <- models@models[[sys_id]]

    # Skip failed models
    if (is.null(model) || model$status != "success") next

    # Get repo dataset
    repo_dataset <- report_data$datasets[[sys_id]]
    if (is.null(repo_dataset)) next

    # Get all compounds in repo system
    repo_inchi <- gsub("/(t|b|m|s)[^/]*.*$", "", repo_dataset$rtdata$inchi.std)
    repo_rt <- repo_dataset$rtdata$rt
    repo_name <- if ("name" %in% names(repo_dataset$rtdata)) {
      repo_dataset$rtdata$name
    } else {
      rep(NA_character_, length(repo_inchi))
    }

    # Get unique compounds
    unique_inchi <- unique(repo_inchi)

    # Filter to requested compounds if specified
    if (!is.null(compounds)) {
      unique_inchi <- unique_inchi[unique_inchi %in% compounds_clean]
    }

    if (length(unique_inchi) == 0) next

    # Get model diagnostics
    model_diag <- diagnostics[diagnostics$source_system == sys_id, ]
    model_r_squared <- if (nrow(model_diag) > 0) model_diag$r_squared else NA
    model_rmse <- if (nrow(model_diag) > 0) model_diag$rmse else NA

    # Generate predictions for each compound
    for (inchi in unique_inchi) {
      rows <- which(repo_inchi == inchi)
      source_rt <- stats::median(repo_rt[rows])
      compound_name <- repo_name[rows[1]]

      # Check if in model range
      if (source_rt < min(model$newdata) || source_rt > max(model$newdata)) {
        next
      }

      # Find prediction by interpolating from grid
      idx <- which.min(abs(source_rt - model$newdata))
      pred_rt <- model$ci[idx, 1]
      ci_lower <- model$ci[idx, 2]
      ci_upper <- model$ci[idx, 3]

      predictions_list[[length(predictions_list) + 1]] <- data.frame(
        compound_id = inchi,
        compound_name = compound_name,
        source_system = sys_id,
        source_rt = source_rt,
        predicted_rt = pred_rt,
        ci_lower = ci_lower,
        ci_upper = ci_upper,
        ci_width = ci_upper - ci_lower,
        model_r_squared = model_r_squared,
        model_rmse = model_rmse,
        n_calibration = model$n_points,
        stringsAsFactors = FALSE
      )
    }

    if (verbose) cli::cli_progress_update()
  }

  # Combine all predictions
  predictions <- if (length(predictions_list) > 0) {
    do.call(rbind, predictions_list)
  } else {
    data.frame(
      compound_id = character(),
      compound_name = character(),
      source_system = character(),
      source_rt = numeric(),
      predicted_rt = numeric(),
      ci_lower = numeric(),
      ci_upper = numeric(),
      ci_width = numeric(),
      model_r_squared = numeric(),
      model_rmse = numeric(),
      n_calibration = integer(),
      stringsAsFactors = FALSE
    )
  }

  # Select best prediction per compound (narrowest CI width)
  # A compound may have predictions from multiple source systems
  if (nrow(predictions) > 0) {
    predictions <- predictions[order(predictions$compound_id, predictions$ci_width), ]
    predictions <- predictions[!duplicated(predictions$compound_id), ]
    rownames(predictions) <- NULL
  }

  elapsed_time <- as.numeric(difftime(Sys.time(), start_time, units = "secs"))

  metadata <- list(
    n_predictions = nrow(predictions),
    n_unique_compounds = length(unique(predictions$compound_id)),
    n_models_used = length(unique(predictions$source_system)),
    created_at = Sys.time(),
    elapsed_time = elapsed_time,
    filtered_compounds = !is.null(compounds)
  )

  result <- methods::new(
    "UserRTPredictions",
    predictions = predictions,
    source_models = names(models@models),
    metadata = metadata
  )

  if (verbose) {
    cli::cli_progress_done()
    cli::cli_rule("Complete")
    cli::cli_alert_success("Time elapsed: {round(elapsed_time, 1)} seconds")
    cli::cli_alert_success("Total predictions: {nrow(predictions)}")
    cli::cli_alert_success("Unique compounds: {metadata$n_unique_compounds}")
  }

  return(result)
}
