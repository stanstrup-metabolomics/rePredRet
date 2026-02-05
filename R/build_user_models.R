#' @title Build Calibration Models from Repository Systems to User's System
#' @description Builds calibration models from repository systems to the user's
#'   chromatographic system using the user's retention time data as calibration.
#'
#' @param inchis Character vector of InChI identifiers for user's compounds.
#' @param rts Numeric vector of retention times (in minutes) corresponding to inchis.
#' @param system_type Character, either "RP" (reversed-phase) or "HILIC".
#'   Only models from systems of the same type will be built.
#' @param report_data Optional. Output from load_report_data(). If NULL,
#'   downloads and loads the latest RepoRT data.
#' @param min_compounds Minimum number of common compounds required to build
#'   a model (default 10).
#' @param method Method for CI calculation: "fast_ci" (default, fast) or
#'   "bootstrap" (slower, more accurate).
#' @param alpha Significance level for CI (default 0.05 for 95% CI).
#' @param n_boot Number of bootstrap iterations if method = "bootstrap" (default 200).
#' @param n_workers Number of parallel workers (default = available cores).
#'   Currently not used (models built sequentially).
#' @param verbose Print progress messages (default TRUE).
#'
#' @return A UserRTModels S4 object containing:
#'   \describe{
#'     \item{models}{Named list of fitted model objects (keyed by source system ID)}
#'     \item{diagnostics}{Per-model quality metrics (R-squared, RMSE, CI widths)}
#'     \item{input_data}{The user's original input data}
#'     \item{report_data}{Cached RepoRT data for later predictions}
#'     \item{metadata}{Run information and summary statistics}
#'   }
#'
#' @details
#' The function:
#' \enumerate{
#'   \item Validates inputs and standardizes InChIs (removes stereochemistry)
#'   \item Loads RepoRT data and filters to matching system type
#'   \item For each compatible repository system:
#'     \itemize{
#'       \item Finds common compounds between repo system and user's data
#'       \item Builds a monotonic GAM model (repo RT -> user RT)
#'       \item Calculates model diagnostics
#'     }
#'   \item Returns all models in a UserRTModels object
#' }
#'
#' Direction: FROM repository systems TO user's system (user's system is TARGET)
#'
#' @examples
#' \dontrun{
#' # User has measured RTs for some compounds
#' my_inchis <- c("InChI=1S/C8H10N4O2/...", "InChI=1S/C10H14N2/...")
#' my_rts <- c(3.2, 5.7)
#'
#' # Build models from all RP systems
#' my_models <- build_user_models(
#'   inchis = my_inchis,
#'   rts = my_rts,
#'   system_type = "RP"
#' )
#'
#' # View summary
#' show(my_models)
#'
#' # Plot diagnostics
#' plot(my_models, type = "diagnostics")
#'
#' # Save for later use
#' saveRDS(my_models, "my_rp_models.rds")
#'
#' # Generate predictions
#' predictions <- predict_rt_user(my_models)
#' }
#'
#' @seealso \code{\link{predict_rt_user}} for generating predictions from models,
#'   \code{\link{UserRTModels-class}} for the return object class
#'
#' @export
#' @importFrom methods new
build_user_models <- function(inchis,
                               rts,
                               system_type = c("RP", "HILIC"),
                               report_data = NULL,
                               min_compounds = 10,
                               method = c("fast_ci", "bootstrap"),
                               alpha = 0.05,
                               n_boot = 200,
                               n_workers = parallel::detectCores(),
                               verbose = TRUE) {

  # Validate inputs
  system_type <- match.arg(system_type)
  method <- match.arg(method)

  if (length(inchis) != length(rts)) {
    stop("inchis and rts must have the same length")
  }

  if (length(inchis) < min_compounds) {
    stop("Need at least ", min_compounds, " compounds. Provided: ", length(inchis))
  }

  if (any(is.na(rts)) || any(rts <= 0)) {
    stop("All retention times must be positive, non-NA values")
  }

  if (!all(grepl("^InChI=", inchis))) {
    warning("Some inputs don't look like InChI strings. Proceeding anyway.")
  }

  # Start timing
  start_time <- Sys.time()

  if (verbose) {
    message("\n================================================")
    message("  build_user_models")
    message("  System type: ", system_type)
    message("  User compounds: ", length(inchis))
    message("  Method: ", method)
    message("================================================\n")
  }

  # Step 1: Standardize user's InChIs
  if (verbose) message("Step 1/3: Standardizing InChIs...")

  inchis_clean <- gsub("/(t|b|m|s)[^/]*.*$", "", inchis)

  # Handle duplicates: take median RT per unique InChI
  user_data <- data.frame(
    inchi_orig = inchis,
    inchi_clean = inchis_clean,
    rt = rts,
    stringsAsFactors = FALSE
  )

  user_data_agg <- stats::aggregate(rt ~ inchi_clean, data = user_data, FUN = stats::median)
  names(user_data_agg) <- c("inchi_clean", "rt")

  # Store input data for the result object
  input_data <- data.frame(
    inchi = inchis,
    inchi_clean = inchis_clean,
    rt = rts,
    stringsAsFactors = FALSE
  )

  if (verbose) {
    message("  Unique compounds after deduplication: ", nrow(user_data_agg))
  }

  # Step 2: Load RepoRT data
  if (verbose) message("\nStep 2/3: Loading RepoRT data...")

  if (is.null(report_data)) {
    report_path <- download_report()
    report_data <- load_report_data(report_path, method_types = system_type)
  } else {
    # Filter to requested system type
    report_data$studies <- report_data$studies[
      report_data$studies$method.type == system_type,
    ]
    report_data$datasets <- report_data$datasets[
      names(report_data$datasets) %in% report_data$studies$id
    ]
  }

  n_systems <- length(report_data$datasets)
  if (verbose) {
    message("  Found ", n_systems, " ", system_type, " systems in RepoRT")
  }

  if (n_systems == 0) {
    stop("No ", system_type, " systems found in RepoRT data")
  }

  # Step 3: Build models
  if (verbose) message("\nStep 3/3: Building calibration models...")

  models <- list()
  diagnostics_list <- list()
  failed_systems <- character()

  # Create user "dataset" object compatible with get_common_compounds
  user_dataset <- list(
    rtdata = data.frame(
      inchi.std = user_data_agg$inchi_clean,
      rt = user_data_agg$rt,
      stringsAsFactors = FALSE
    )
  )

  # Progress tracking
  n_processed <- 0
  n_total <- length(report_data$datasets)

  # Build models from each repo system to user's system
  for (sys_id in names(report_data$datasets)) {
    n_processed <- n_processed + 1
    repo_dataset <- report_data$datasets[[sys_id]]

    # Get common compounds
    # Note: Direction is REPO -> USER (repo is source/sys1, user is target/sys2)
    rt_matrix <- get_common_compounds(repo_dataset, user_dataset)

    if (nrow(rt_matrix) < min_compounds) {
      failed_systems <- c(failed_systems, sys_id)
      next
    }

    if (verbose && (length(models) %% 10 == 0 || n_processed == n_total)) {
      message("  [", n_processed, "/", n_total, "] Building model from ",
              sys_id, " (", nrow(rt_matrix), " compounds)...")
    }

    # Build model: FROM repo TO user
    model <- tryCatch({
      build_model(
        rt_matrix = rt_matrix,
        sys1_id = sys_id,
        sys2_id = "USER",
        n_boot = n_boot,
        alpha = alpha,
        n_cores = 1,  # Single-threaded per model
        method = method,
        save_plot = FALSE
      )
    }, error = function(e) {
      if (verbose) message("    Error building model for ", sys_id, ": ", e$message)
      list(status = "error", error = e$message)
    })

    if (!is.null(model$status) && model$status == "success") {
      models[[sys_id]] <- model

      # Calculate diagnostics
      pred_at_cal <- sapply(model$calibration[, 1], function(x) {
        idx <- which.min(abs(x - model$newdata))
        model$ci[idx, 1]
      })

      errors <- abs(model$calibration[, 2] - pred_at_cal)
      ss_res <- sum((model$calibration[, 2] - pred_at_cal)^2)
      ss_tot <- sum((model$calibration[, 2] - mean(model$calibration[, 2]))^2)
      r_squared <- if (ss_tot > 0) 1 - (ss_res / ss_tot) else NA
      rmse <- sqrt(mean(errors^2))

      diagnostics_list[[sys_id]] <- data.frame(
        source_system = sys_id,
        n_calibration = model$n_points,
        r_squared = r_squared,
        rmse = rmse,
        mean_error = model$stats$mean_error,
        median_error = model$stats$median_error,
        q95_error = model$stats$q95_error,
        mean_ci_width = model$stats$mean_ci_width,
        median_ci_width = model$stats$median_ci_width,
        q95_ci_width = model$stats$q95_ci_width,
        rt_range_source = max(model$calibration[, 1]) - min(model$calibration[, 1]),
        rt_range_user = max(model$calibration[, 2]) - min(model$calibration[, 2]),
        stringsAsFactors = FALSE
      )
    } else {
      failed_systems <- c(failed_systems, sys_id)
    }
  }

  diagnostics <- if (length(diagnostics_list) > 0) {
    do.call(rbind, diagnostics_list)
  } else {
    data.frame()
  }

  if (verbose) {
    message("\n  Models built: ", length(models))
    message("  Failed (insufficient compounds or error): ", length(failed_systems))
  }

  if (length(models) == 0) {
    warning("No models could be built. Check that your compounds overlap with RepoRT.")
  }

  # Create result object
  elapsed_time <- as.numeric(difftime(Sys.time(), start_time, units = "secs"))

  metadata <- list(
    system_type = system_type,
    n_input_compounds = length(inchis),
    n_input_unique = nrow(user_data_agg),
    n_models_built = length(models),
    n_models_failed = length(failed_systems),
    created_at = Sys.time(),
    elapsed_time = elapsed_time,
    repredret_version = tryCatch(
      as.character(utils::packageVersion("rePredRet")),
      error = function(e) "unknown"
    ),
    method = method,
    alpha = alpha,
    min_compounds = min_compounds
  )

  result <- methods::new(
    "UserRTModels",
    models = models,
    diagnostics = diagnostics,
    input_data = input_data,
    report_data = report_data,
    metadata = metadata
  )

  if (verbose) {
    message("\n================================================")
    message("  Complete!")
    message("  Time elapsed: ", round(elapsed_time, 1), " seconds")
    message("  Models built: ", length(models))
    message("================================================\n")
  }

  return(result)
}
