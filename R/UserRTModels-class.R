#' @title UserRTModels S4 Class
#' @description S4 class for storing calibration models built from repository
#'   systems to a user's chromatographic system.
#'
#' @slot models A named list of model objects (from build_model) keyed by
#'   source_system_id. Each model contains:
#'   \describe{
#'     \item{status}{"success" or "insufficient_data"}
#'     \item{sys1_id}{Source system ID (repository system)}
#'     \item{sys2_id}{"USER" (target is user's system)}
#'     \item{n_points}{Number of calibration compounds}
#'     \item{newdata}{Prediction grid (1000 points)}
#'     \item{ci}{Matrix `[pred, lower, upper]` at newdata points}
#'     \item{calibration}{Matrix `[rt_source, rt_user]` for calibration compounds}
#'     \item{stats}{Error statistics at calibration points}
#'   }
#' @slot diagnostics A data.frame with per-model diagnostic metrics:
#'   \describe{
#'     \item{source_system}{Repository system ID}
#'     \item{n_calibration}{Number of common compounds used}
#'     \item{r_squared}{Model R-squared}
#'     \item{rmse}{Root mean square error}
#'     \item{mean_error}{Mean absolute error}
#'     \item{median_error}{Median absolute error}
#'     \item{q95_error}{95th percentile error}
#'     \item{mean_ci_width}{Mean CI width at calibration points}
#'     \item{median_ci_width}{Median CI width}
#'     \item{rt_range_source}{RT range in source system}
#'     \item{rt_range_user}{RT range in user's system}
#'   }
#' @slot input_data A data.frame with the user's original input:
#'   \describe{
#'     \item{inchi}{Original InChI strings}
#'     \item{inchi_clean}{Standardized InChI (no stereochemistry)}
#'     \item{rt}{Retention times provided}
#'   }
#' @slot report_data A list containing cached RepoRT data for predictions
#' @slot metadata A list containing:
#'   \describe{
#'     \item{system_type}{RP or HILIC}
#'     \item{n_input_compounds}{Number of user compounds}
#'     \item{n_models_built}{Number of successful models}
#'     \item{n_models_failed}{Number of failed models}
#'     \item{created_at}{Timestamp}
#'     \item{repredret_version}{Package version}
#'     \item{method}{CI method used}
#'   }
#'
#' @name UserRTModels-class
#' @rdname UserRTModels-class
#' @exportClass UserRTModels
#' @importFrom methods new setClass setMethod setGeneric show
setClass(
 "UserRTModels",
 slots = c(
   models = "list",
   diagnostics = "data.frame",
   input_data = "data.frame",
   report_data = "list",
   metadata = "list"
 ),
 prototype = list(
   models = list(),
   diagnostics = data.frame(),
   input_data = data.frame(),
   report_data = list(),
   metadata = list()
 )
)

setValidity("UserRTModels", function(object) {
 errors <- character()

 # Check metadata has required fields
 required_meta <- c("system_type", "n_input_compounds", "created_at")
 missing_meta <- setdiff(required_meta, names(object@metadata))
 if (length(missing_meta) > 0) {
   errors <- c(errors, paste("Missing metadata fields:",
                             paste(missing_meta, collapse = ", ")))
 }

 # Check system_type is valid
 if ("system_type" %in% names(object@metadata)) {
   if (!object@metadata$system_type %in% c("RP", "HILIC")) {
     errors <- c(errors, "system_type must be 'RP' or 'HILIC'")
   }
 }

 if (length(errors) == 0) TRUE else errors
})


#' Show Method for UserRTModels
#'
#' @param object A UserRTModels object
#' @rdname UserRTModels-class
#' @export
setMethod("show", "UserRTModels", function(object) {
 cat("UserRTModels Object\n")
 cat("===================\n\n")

 cat("System Type:", object@metadata$system_type, "\n")
 cat("Created:", as.character(object@metadata$created_at), "\n\n")

 cat("Input:\n")
 cat("  User compounds:", object@metadata$n_input_compounds, "\n")
 if (!is.null(object@metadata$n_input_unique)) {
   cat("  Unique (after dedup):", object@metadata$n_input_unique, "\n")
 }
 cat("\n")

 cat("Models:\n")
 cat("  Built successfully:", object@metadata$n_models_built, "\n")
 cat("  Failed (insufficient data):", object@metadata$n_models_failed, "\n\n")

 if (nrow(object@diagnostics) > 0) {
   cat("Model Quality (median across models):\n")
   cat("  R-squared:", round(median(object@diagnostics$r_squared, na.rm = TRUE), 3), "\n")
   cat("  RMSE:", round(median(object@diagnostics$rmse, na.rm = TRUE), 3), "min\n")
   cat("  Median CI width:", round(median(object@diagnostics$median_ci_width, na.rm = TRUE), 3), "min\n")
   cat("\n")
   cat("Calibration compounds per model:\n")
   cat("  Min:", min(object@diagnostics$n_calibration), "\n")
   cat("  Median:", median(object@diagnostics$n_calibration), "\n")
   cat("  Max:", max(object@diagnostics$n_calibration), "\n")
 }

 invisible(object)
})


#' Summary Method for UserRTModels
#'
#' @param object A UserRTModels object
#' @param ... Ignored.
#' @return A list with summary statistics
#' @rdname UserRTModels-class
#' @export
setMethod("summary", "UserRTModels", function(object, ...) {

 model_summary <- if (nrow(object@diagnostics) > 0) {
   list(
     n_models = nrow(object@diagnostics),
     r_squared = list(
       min = min(object@diagnostics$r_squared, na.rm = TRUE),
       median = median(object@diagnostics$r_squared, na.rm = TRUE),
       max = max(object@diagnostics$r_squared, na.rm = TRUE)
     ),
     rmse = list(
       min = min(object@diagnostics$rmse, na.rm = TRUE),
       median = median(object@diagnostics$rmse, na.rm = TRUE),
       max = max(object@diagnostics$rmse, na.rm = TRUE)
     ),
     n_calibration = list(
       min = min(object@diagnostics$n_calibration),
       median = median(object@diagnostics$n_calibration),
       max = max(object@diagnostics$n_calibration)
     ),
     ci_width = list(
       min = min(object@diagnostics$median_ci_width, na.rm = TRUE),
       median = median(object@diagnostics$median_ci_width, na.rm = TRUE),
       max = max(object@diagnostics$median_ci_width, na.rm = TRUE)
     )
   )
 } else {
   list(n_models = 0)
 }

 list(
   metadata = object@metadata,
   models = model_summary,
   source_systems = names(object@models)
 )
})


#' Extract Model by System ID
#'
#' @param x A UserRTModels object
#' @param i System ID (character)
#' @return The model object for that system
#' @rdname UserRTModels-class
#' @export
setMethod("[[", signature(x = "UserRTModels", i = "character"),
         function(x, i) {
           if (!i %in% names(x@models)) {
             stop("System '", i, "' not found. Available: ",
                  paste(head(names(x@models), 10), collapse = ", "),
                  if (length(x@models) > 10) "..." else "")
           }
           x@models[[i]]
         })


#' List Available System IDs
#'
#' @param x A UserRTModels object
#' @return Character vector of system IDs
#' @rdname UserRTModels-class
#' @export
setMethod("names", "UserRTModels", function(x) {
 names(x@models)
})


#' Get Number of Models
#'
#' @param x A UserRTModels object
#' @return Integer
#' @rdname UserRTModels-class
#' @export
setMethod("length", "UserRTModels", function(x) {
 length(x@models)
})


#' Plot Method for UserRTModels
#'
#' @param x A UserRTModels object
#' @param y Ignored
#' @param type Plot type: "diagnostics" (default), "calibration", "r_squared", "ci_widths"
#' @param max_plots Maximum number of calibration plots to show (default 9)
#' @param ... Additional arguments
#' @return A ggplot object or list of ggplot objects
#' @importFrom ggplot2 ggplot aes geom_point geom_histogram geom_line geom_ribbon
#'   labs theme_bw facet_wrap scale_color_viridis_c
#' @rdname UserRTModels-class
#' @export
setMethod("plot", signature(x = "UserRTModels", y = "missing"),
         function(x, y, type = c("diagnostics", "calibration", "r_squared", "ci_widths"),
                  max_plots = 9, ...) {

           type <- match.arg(type)

           switch(type,
                  diagnostics = .plot_models_diagnostics(x),
                  calibration = .plot_models_calibration(x, max_plots),
                  r_squared = .plot_models_rsquared(x),
                  ci_widths = .plot_models_ci_widths(x)
           )
         })


# Internal plotting helpers
.plot_models_diagnostics <- function(object) {
 diag <- object@diagnostics

 if (nrow(diag) == 0) {
   stop("No diagnostics data available")
 }

 ggplot2::ggplot(diag, ggplot2::aes(x = n_calibration, y = r_squared)) +
   ggplot2::geom_point(ggplot2::aes(color = median_ci_width), alpha = 0.7, size = 3) +
   ggplot2::scale_color_viridis_c(name = "Median CI\nWidth (min)") +
   ggplot2::labs(
     x = "Number of Calibration Compounds",
     y = "R-squared",
     title = "Model Quality Overview",
     subtitle = paste(nrow(diag), "models from", object@metadata$system_type, "systems")
   ) +
   ggplot2::theme_bw()
}


.plot_models_rsquared <- function(object) {
 diag <- object@diagnostics

 ggplot2::ggplot(diag, ggplot2::aes(x = r_squared)) +
   ggplot2::geom_histogram(bins = 20, fill = "steelblue", color = "white") +
   ggplot2::geom_vline(xintercept = median(diag$r_squared), color = "red",
                       linetype = "dashed", linewidth = 1) +
   ggplot2::labs(
     x = "R-squared",
     y = "Number of Models",
     title = "Distribution of Model R-squared Values",
     subtitle = paste("Median R-squared:", round(median(diag$r_squared), 3))
   ) +
   ggplot2::theme_bw()
}


.plot_models_ci_widths <- function(object) {
 diag <- object@diagnostics

 ggplot2::ggplot(diag, ggplot2::aes(x = median_ci_width)) +
   ggplot2::geom_histogram(bins = 20, fill = "steelblue", color = "white") +
   ggplot2::geom_vline(xintercept = median(diag$median_ci_width), color = "red",
                       linetype = "dashed", linewidth = 1) +
   ggplot2::labs(
     x = "Median CI Width (minutes)",
     y = "Number of Models",
     title = "Distribution of Confidence Interval Widths",
     subtitle = paste("Median:", round(median(diag$median_ci_width), 3), "min")
   ) +
   ggplot2::theme_bw()
}


.plot_models_calibration <- function(object, max_plots = 9) {
 models <- object@models
 diag <- object@diagnostics

 if (length(models) == 0) {
   stop("No models available")
 }

 # Sort by R-squared and take top N
 diag_sorted <- diag[order(-diag$r_squared), ]
 top_ids <- head(diag_sorted$source_system, max_plots)

 # Create data for faceted plot
 plot_data_calib <- list()
 plot_data_pred <- list()

 for (sys_id in top_ids) {
   model <- models[[sys_id]]
   if (is.null(model) || model$status != "success") next

   r_sq <- diag[diag$source_system == sys_id, "r_squared"]
   label <- paste0(sys_id, " (R\u00b2=", round(r_sq, 3), ")")

   # Calibration points
   plot_data_calib[[sys_id]] <- data.frame(
     source_system = label,
     rt_source = model$calibration[, 1],
     rt_user = model$calibration[, 2],
     stringsAsFactors = FALSE
   )

   # Prediction line (subsample for speed)
   idx <- seq(1, length(model$newdata), length.out = 200)
   plot_data_pred[[sys_id]] <- data.frame(
     source_system = label,
     rt_source = model$newdata[idx],
     rt_user = model$ci[idx, 1],
     ci_lower = model$ci[idx, 2],
     ci_upper = model$ci[idx, 3],
     stringsAsFactors = FALSE
   )
 }

 calib_df <- do.call(rbind, plot_data_calib)
 pred_df <- do.call(rbind, plot_data_pred)

 ggplot2::ggplot() +
   ggplot2::geom_ribbon(
     data = pred_df,
     ggplot2::aes(x = rt_source, ymin = ci_lower, ymax = ci_upper),
     fill = "lightblue", alpha = 0.5
   ) +
   ggplot2::geom_line(
     data = pred_df,
     ggplot2::aes(x = rt_source, y = rt_user),
     color = "blue"
   ) +
   ggplot2::geom_point(
     data = calib_df,
     ggplot2::aes(x = rt_source, y = rt_user),
     color = "black", alpha = 0.6
   ) +
   ggplot2::facet_wrap(~source_system, scales = "free") +
   ggplot2::labs(
     x = "RT in Repository System (min)",
     y = "RT in User's System (min)",
     title = "Calibration Curves (Top Models by R-squared)"
   ) +
   ggplot2::theme_bw()
}
