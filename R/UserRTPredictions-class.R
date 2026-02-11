#' @title UserRTPredictions S4 Class
#' @description S4 class for storing RT predictions generated from UserRTModels.
#'
#' @slot predictions A data.frame containing all predictions with columns:
#'   \describe{
#'     \item{compound_id}{InChI identifier of the predicted compound}
#'     \item{compound_name}{Human-readable name (if available)}
#'     \item{source_system}{Repository system ID where compound was measured}
#'     \item{source_rt}{Original RT in source system (minutes)}
#'     \item{predicted_rt}{Predicted RT in user's system (minutes)}
#'     \item{ci_lower}{Lower bound of confidence interval}
#'     \item{ci_upper}{Upper bound of confidence interval}
#'     \item{ci_width}{Width of confidence interval}
#'     \item{model_r_squared}{R-squared of the model used}
#'     \item{model_rmse}{RMSE of the model used}
#'     \item{n_calibration}{Number of calibration compounds in model}
#'   }
#' @slot source_models Character vector of system IDs for models used
#' @slot metadata A list containing:
#'   \describe{
#'     \item{n_predictions}{Total number of predictions}
#'     \item{n_unique_compounds}{Number of unique compounds predicted}
#'     \item{n_models_used}{Number of models used}
#'     \item{created_at}{Timestamp}
#'   }
#'
#' @name UserRTPredictions-class
#' @rdname UserRTPredictions-class
#' @exportClass UserRTPredictions
#' @importFrom methods new setClass setMethod show
setClass(
  "UserRTPredictions",
  slots = c(
    predictions = "data.frame",
    source_models = "character",
    metadata = "list"
  ),
  prototype = list(
    predictions = data.frame(),
    source_models = character(),
    metadata = list()
  )
)


#' Show Method for UserRTPredictions
#'
#' @param object A UserRTPredictions object
#' @rdname UserRTPredictions-class
#' @export
setMethod("show", "UserRTPredictions", function(object) {
  cat("UserRTPredictions Object\n")
  cat("========================\n\n")

  cat("Created:", as.character(object@metadata$created_at), "\n\n")

  cat("Predictions:\n")
  cat("  Total predictions:", object@metadata$n_predictions, "\n")
  cat("  Unique compounds:", object@metadata$n_unique_compounds, "\n")
  cat("  Models used:", object@metadata$n_models_used, "\n\n")

  if (nrow(object@predictions) > 0) {
    cat("Predicted RT range:",
        round(min(object@predictions$predicted_rt), 2), "-",
        round(max(object@predictions$predicted_rt), 2), "min\n")
    cat("CI width (median):",
        round(median(object@predictions$ci_width), 3), "min\n\n")

    cat("Preview (first 5 predictions):\n")
    preview <- head(object@predictions[, c("compound_id", "source_system",
                                            "predicted_rt", "ci_width")], 5)
    # Truncate long InChIs for display
    preview$compound_id <- substr(preview$compound_id, 1, 40)
    preview$compound_id <- paste0(preview$compound_id, "...")
    print(preview, row.names = FALSE)
  }

  invisible(object)
})


#' Summary Method for UserRTPredictions
#'
#' @param object A UserRTPredictions object
#' @return A list with summary statistics
#' @rdname UserRTPredictions-class
#' @export
setMethod("summary", "UserRTPredictions", function(object, ...) {

  pred_summary <- if (nrow(object@predictions) > 0) {
    list(
      n_predictions = nrow(object@predictions),
      n_unique_compounds = length(unique(object@predictions$compound_id)),
      n_source_systems = length(unique(object@predictions$source_system)),
      predicted_rt = list(
        min = min(object@predictions$predicted_rt),
        median = median(object@predictions$predicted_rt),
        max = max(object@predictions$predicted_rt)
      ),
      ci_width = list(
        min = min(object@predictions$ci_width),
        median = median(object@predictions$ci_width),
        max = max(object@predictions$ci_width),
        q95 = quantile(object@predictions$ci_width, 0.95)
      ),
      model_quality = list(
        r_squared_median = median(object@predictions$model_r_squared),
        rmse_median = median(object@predictions$model_rmse)
      )
    )
  } else {
    list(n_predictions = 0)
  }

  list(
    metadata = object@metadata,
    predictions = pred_summary,
    source_models = object@source_models
  )
})


#' Convert UserRTPredictions to data.frame
#'
#' @param x A UserRTPredictions object
#' @param row.names Ignored
#' @param optional Ignored
#' @param ... Ignored
#' @return A data.frame of predictions
#' @rdname UserRTPredictions-class
#' @export
setMethod("as.data.frame", "UserRTPredictions",
          function(x, row.names = NULL, optional = FALSE, ...) {
            x@predictions
          })


#' Get Number of Predictions
#'
#' @param x A UserRTPredictions object
#' @return Integer
#' @rdname UserRTPredictions-class
#' @export
setMethod("length", "UserRTPredictions", function(x) {
  nrow(x@predictions)
})


#' Subset Predictions
#'
#' @param x A UserRTPredictions object
#' @param i Row indices or logical vector
#' @param j Column names (optional)
#' @param drop Passed to data.frame subsetting
#' @return Subset of predictions data.frame
#' @rdname UserRTPredictions-class
#' @export
setMethod("[", signature(x = "UserRTPredictions"),
          function(x, i, j, ..., drop = TRUE) {
            if (missing(j)) {
              x@predictions[i, , drop = drop]
            } else {
              x@predictions[i, j, drop = drop]
            }
          })


#' Plot Method for UserRTPredictions
#'
#' @param x A UserRTPredictions object
#' @param y Ignored
#' @param type Plot type: "ci_widths" (default), "rt_distribution", "by_system"
#' @param ... Additional arguments
#' @return A ggplot object
#' @importFrom ggplot2 ggplot aes geom_histogram geom_boxplot geom_density
#'   labs theme_bw facet_wrap coord_flip
#' @rdname UserRTPredictions-class
#' @export
setMethod("plot", signature(x = "UserRTPredictions", y = "missing"),
          function(x, y, type = c("ci_widths", "rt_distribution", "by_system"), ...) {

            type <- match.arg(type)

            switch(type,
                   ci_widths = .plot_preds_ci_widths(x),
                   rt_distribution = .plot_preds_rt_dist(x),
                   by_system = .plot_preds_by_system(x)
            )
          })


# Internal plotting helpers
.plot_preds_ci_widths <- function(object) {
  preds <- object@predictions

  ggplot2::ggplot(preds, ggplot2::aes(x = ci_width)) +
    ggplot2::geom_histogram(bins = 50, fill = "steelblue", color = "white") +
    ggplot2::geom_vline(xintercept = median(preds$ci_width), color = "red",
                        linetype = "dashed", linewidth = 1) +
    ggplot2::labs(
      x = "Confidence Interval Width (minutes)",
      y = "Number of Predictions",
      title = "Distribution of Prediction Confidence Intervals",
      subtitle = paste("Median CI width:", round(median(preds$ci_width), 3), "min")
    ) +
    ggplot2::theme_bw()
}


.plot_preds_rt_dist <- function(object) {
  preds <- object@predictions

  ggplot2::ggplot(preds, ggplot2::aes(x = predicted_rt)) +
    ggplot2::geom_histogram(bins = 50, fill = "steelblue", color = "white") +
    ggplot2::labs(
      x = "Predicted Retention Time (minutes)",
      y = "Number of Predictions",
      title = "Distribution of Predicted Retention Times"
    ) +
    ggplot2::theme_bw()
}


.plot_preds_by_system <- function(object) {
  preds <- object@predictions

  # Summarize by system
  n_by_system <- as.data.frame(table(preds$source_system))
  names(n_by_system) <- c("source_system", "n_predictions")
  n_by_system <- n_by_system[order(-n_by_system$n_predictions), ]

  # Take top 20 systems
  if (nrow(n_by_system) > 20) {
    n_by_system <- head(n_by_system, 20)
  }

  n_by_system$source_system <- factor(n_by_system$source_system,
                                       levels = rev(n_by_system$source_system))

  ggplot2::ggplot(n_by_system, ggplot2::aes(x = source_system, y = n_predictions)) +
    ggplot2::geom_col(fill = "steelblue") +
    ggplot2::coord_flip() +
    ggplot2::labs(
      x = "Source System",
      y = "Number of Predictions",
      title = "Predictions by Source System",
      subtitle = if (length(unique(preds$source_system)) > 20)
        "Showing top 20 systems" else NULL
    ) +
    ggplot2::theme_bw()
}
