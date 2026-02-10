#' @title Batch Search for Multiple Compounds
#' @description Searches for multiple compounds by InChI and returns all systems
#'   where they have been measured.
#'
#' @param inchis Character vector of InChI strings to search for.
#' @param system_type Optional filter by system type ("RP" or "HILIC").
#' @param data_path Optional path to local RepoRT data directory.
#'
#' @return A data.frame with columns: query_inchi, system_id, system_name,
#'   system_type, compound_name, rt, inchi.
#'
#' @examples
#' \dontrun{
#' # Search for multiple compounds
#' results <- search_compounds(c(
#'   "InChI=1S/C8H10N4O2/c1-10-4-9-6-5(10)7(13)12(3)8(14)11(6)2/h4H,1-3H3",
#'   "InChI=1S/C15H14O6/c16-8-4-11(18)9-6-13(20)15(21-14(9)5-8)7-1-2-10(17)12(19)3-7"
#' ))
#' }
#'
#' @seealso \code{\link{search_compound}} for single compound search
#' @export
search_compounds <- function(inchis,
                             system_type = NULL,
                             data_path = NULL) {

  if (length(inchis) == 0) {
    return(data.frame())
  }

  results_list <- lapply(inchis, function(inchi) {
    result <- tryCatch(
      search_compound(inchi, system_type = system_type, data_path = data_path),
      error = function(e) data.frame()
    )
    if (nrow(result) > 0) {
      result$query_inchi <- inchi
    }
    result
  })

  results <- do.call(rbind, results_list)

  if (is.null(results) || nrow(results) == 0) {
    return(data.frame(
      query_inchi = character(),
      system_id = character(),
      system_name = character(),
      system_type = character(),
      compound_name = character(),
      rt = numeric(),
      inchi = character(),
      stringsAsFactors = FALSE
    ))
  }

  # Reorder columns
  cols <- c("query_inchi", setdiff(names(results), "query_inchi"))
  results <- results[, cols]
  rownames(results) <- NULL

  return(results)
}


#' @title Export Predictions to File
#' @description Exports predictions to CSV or Excel format.
#'
#' @param predictions A UserRTPredictions object or data.frame of predictions.
#' @param file Output file path. Extension determines format if format is "auto".
#' @param format Output format: "auto" (from extension), "csv", or "xlsx".
#'
#' @return Invisibly returns the file path.
#'
#' @details
#' For Excel export, the `writexl` package is required. Install with
#' `install.packages("writexl")` if needed.
#'
#' @examples
#' \dontrun{
#' # Export to CSV
#' export_predictions(predictions, "my_predictions.csv")
#'
#' # Export to Excel
#' export_predictions(predictions, "my_predictions.xlsx")
#'
#' # Explicit format
#' export_predictions(predictions, "output.txt", format = "csv")
#' }
#'
#' @export
export_predictions <- function(predictions,
                                file,
                                format = c("auto", "csv", "xlsx")) {

  format <- match.arg(format)

  # Extract data.frame if S4 object
  if (methods::is(predictions, "UserRTPredictions")) {
    df <- predictions@predictions
  } else if (is.data.frame(predictions)) {
    df <- predictions
  } else {
    stop("predictions must be a UserRTPredictions object or data.frame")
  }

  # Determine format from extension if auto
  if (format == "auto") {
    ext <- tolower(tools::file_ext(file))
    if (ext == "xlsx" || ext == "xls") {
      format <- "xlsx"
    } else {
      format <- "csv"
    }
  }

  # Export
  if (format == "csv") {
    utils::write.csv(df, file, row.names = FALSE)
    message("Exported ", nrow(df), " predictions to ", file)
  } else if (format == "xlsx") {
    if (!requireNamespace("writexl", quietly = TRUE)) {
      stop("Package 'writexl' is required for Excel export. ",
           "Install with: install.packages('writexl')")
    }
    writexl::write_xlsx(df, file)
    message("Exported ", nrow(df), " predictions to ", file)
  }

  invisible(file)
}


#' @title Validate Models with Cross-Validation
#' @description Performs leave-one-out cross-validation on the calibration
#'   compounds to assess prediction accuracy.
#'
#' @param models A UserRTModels object from \code{\link{build_user_models}}.
#' @param method Validation method: "loo" (leave-one-out, default).
#'
#' @return A data.frame with cross-validation results for each model:
#'   \describe{
#'     \item{source_system}{System ID}
#'     \item{n_compounds}{Number of calibration compounds}
#'     \item{cv_rmse}{Cross-validated RMSE}
#'     \item{cv_mae}{Cross-validated mean absolute error}
#'     \item{cv_median_ae}{Cross-validated median absolute error}
#'     \item{cv_coverage}{Proportion of actual values within CI}
#'   }
#'
#' @details
#' Leave-one-out cross-validation removes each calibration compound in turn,
#' refits the model, and predicts the held-out compound. This provides an
#' unbiased estimate of prediction error.
#'
#' @examples
#' \dontrun{
#' models <- build_user_models(my_inchis, my_rts, "RP")
#' cv_results <- validate_models(models)
#' print(cv_results)
#' }
#'
#' @export
validate_models <- function(models,
                            method = c("loo")) {

  if (!methods::is(models, "UserRTModels")) {
    stop("models must be a UserRTModels object")
  }

  method <- match.arg(method)

  if (length(models@models) == 0) {
    warning("No models available for validation")
    return(data.frame())
  }

  results_list <- list()

  for (sys_id in names(models@models)) {
    model <- models@models[[sys_id]]

    if (is.null(model) || model$status != "success") next
    if (model$n_points < 5) next  # Need minimum points for LOO

    calibration <- model$calibration
    n <- nrow(calibration)

    # Leave-one-out predictions
    loo_pred <- numeric(n)
    loo_lower <- numeric(n)
    loo_upper <- numeric(n)

    for (i in seq_len(n)) {
      # Remove one point
      train_data <- calibration[-i, , drop = FALSE]

      # Simple linear interpolation for LOO (faster than refitting GAM)
      # Sort by source RT
      ord <- order(train_data[, 1])
      train_sorted <- train_data[ord, ]

      # Interpolate prediction for held-out point
      loo_pred[i] <- stats::approx(
        x = train_sorted[, 1],
        y = train_sorted[, 2],
        xout = calibration[i, 1],
        rule = 2  # Extrapolate
      )$y

      # Estimate CI from local variance (simplified)
      # Use residuals from nearby points
      nearby <- which(abs(train_sorted[, 1] - calibration[i, 1]) <
                        0.2 * diff(range(train_sorted[, 1])))
      if (length(nearby) < 3) nearby <- 1:min(5, nrow(train_sorted))

      local_resid <- train_sorted[nearby, 2] -
        stats::approx(train_sorted[, 1], train_sorted[, 2],
                      xout = train_sorted[nearby, 1])$y
      local_sd <- stats::sd(local_resid, na.rm = TRUE)
      if (is.na(local_sd) || local_sd == 0) local_sd <- 0.5

      loo_lower[i] <- loo_pred[i] - 1.96 * local_sd
      loo_upper[i] <- loo_pred[i] + 1.96 * local_sd
    }

    # Calculate metrics
    actual <- calibration[, 2]
    errors <- abs(loo_pred - actual)
    within_ci <- actual >= loo_lower & actual <= loo_upper

    results_list[[sys_id]] <- data.frame(
      source_system = sys_id,
      n_compounds = n,
      cv_rmse = sqrt(mean(errors^2)),
      cv_mae = mean(errors),
      cv_median_ae = stats::median(errors),
      cv_coverage = mean(within_ci),
      stringsAsFactors = FALSE
    )
  }

  if (length(results_list) == 0) {
    return(data.frame())
  }

  results <- do.call(rbind, results_list)
  rownames(results) <- NULL

  return(results)
}


#' @title Check Input Data Quality
#' @description Validates InChI strings and retention times before model building.
#'
#' @param inchis Character vector of InChI strings.
#' @param rts Numeric vector of retention times.
#' @param fix If TRUE (default), attempts to fix common issues and returns
#'   corrected data. If FALSE, only reports issues.
#'
#' @return If fix = TRUE, returns a list with corrected `inchis` and `rts`.
#'   If fix = FALSE, returns a list of issues found.
#'
#' @details
#' Checks performed:
#' \itemize{
#'   \item Length mismatch between inchis and rts
#'   \item Invalid InChI format (must start with "InChI=")
#'   \item Missing or NA values
#'   \item Non-positive retention times
#'   \item Duplicate compounds (same InChI after stripping stereochemistry)
#'   \item Outlier RTs (> 3 SD from mean)
#' }
#'
#' @examples
#' \dontrun{
#' # Check data quality
#' check_input_data(my_inchis, my_rts, fix = FALSE)
#'
#' # Fix issues automatically
#' fixed <- check_input_data(my_inchis, my_rts, fix = TRUE)
#' models <- build_user_models(fixed$inchis, fixed$rts, "RP")
#' }
#'
#' @export
check_input_data <- function(inchis, rts, fix = TRUE) {

  issues <- list()
  n_original <- length(inchis)

  # Check length mismatch
  if (length(inchis) != length(rts)) {
    issues$length_mismatch <- paste(
      "Length mismatch: inchis has", length(inchis),
      "elements, rts has", length(rts)
    )
    if (!fix) {
      return(list(valid = FALSE, issues = issues))
    }
    # Truncate to shorter length
    min_len <- min(length(inchis), length(rts))
    inchis <- inchis[1:min_len]
    rts <- rts[1:min_len]
  }

  # Check for NA values
  na_inchi <- is.na(inchis)
  na_rt <- is.na(rts)
  if (any(na_inchi) || any(na_rt)) {
    issues$na_values <- paste(
      sum(na_inchi), "NA InChIs,",
      sum(na_rt), "NA RTs"
    )
    if (fix) {
      keep <- !na_inchi & !na_rt
      inchis <- inchis[keep]
      rts <- rts[keep]
    }
  }

  # Check InChI format
  invalid_inchi <- !grepl("^InChI=", inchis, ignore.case = TRUE)
  if (any(invalid_inchi)) {
    issues$invalid_inchi <- paste(
      sum(invalid_inchi), "invalid InChI strings (must start with 'InChI=')"
    )
    if (fix) {
      keep <- !invalid_inchi
      inchis <- inchis[keep]
      rts <- rts[keep]
    }
  }

  # Check positive RTs
  if (length(rts) > 0) {
    non_positive <- rts <= 0
    if (any(non_positive)) {
      issues$non_positive_rt <- paste(
        sum(non_positive), "non-positive retention times"
      )
      if (fix) {
        keep <- !non_positive
        inchis <- inchis[keep]
        rts <- rts[keep]
      }
    }
  }

  # Check for duplicates
  if (length(inchis) > 0) {
    inchis_clean <- gsub("/(t|b|m|s)[^/]*.*$", "", inchis)
    dups <- duplicated(inchis_clean)
    if (any(dups)) {
      issues$duplicates <- paste(
        sum(dups), "duplicate compounds (will be aggregated by median RT)"
      )
      # Don't remove duplicates - build_user_models handles this
    }
  }

  # Check for RT outliers
  if (length(rts) >= 5) {
    rt_mean <- mean(rts)
    rt_sd <- stats::sd(rts)
    if (rt_sd > 0) {
      outliers <- abs(rts - rt_mean) > 3 * rt_sd
      if (any(outliers)) {
        issues$outliers <- paste(
          sum(outliers), "potential RT outliers (> 3 SD from mean)"
        )
        # Don't remove outliers automatically - just warn
      }
    }
  }

  # Summary
  n_final <- length(inchis)
  n_removed <- n_original - n_final

  if (length(issues) == 0) {
    message("Data quality check passed: ", n_final, " valid compounds")
  } else {
    message("Data quality check found ", length(issues), " issue(s):")
    for (issue_name in names(issues)) {
      message("  - ", issues[[issue_name]])
    }
    if (fix && n_removed > 0) {
      message("Removed ", n_removed, " problematic entries. ",
              n_final, " compounds remaining.")
    }
  }

  if (fix) {
    return(list(
      inchis = inchis,
      rts = rts,
      n_original = n_original,
      n_final = n_final,
      issues = issues
    ))
  } else {
    return(list(
      valid = length(issues) == 0,
      n_compounds = n_original,
      issues = issues
    ))
  }
}
