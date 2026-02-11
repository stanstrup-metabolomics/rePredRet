#' @title Package Initialization
#' @description Package startup messages and initialization.
#' @name zzz
#' @keywords internal
NULL

# Suppress R CMD check NOTEs for ggplot2 NSE and data.table/dplyr column names
utils::globalVariables(c(
  "ci_lower", "ci_upper", "ci_width", "compound_label", "compound_name",
  "from_id", "id", "lower", "median_ci_width", "median_ci_width_rel",
  "median_error_abs", "median_error_rel", "method_type", "n_calibration",
  "n_compounds", "n_predictions", "pred", "predicted_rt", "r_sq_display",
  "r_squared", "rt", "rt.1", "rt.2", "rt_source", "rt_target", "rt_user",
  "source_system", "to_id", "upper", "x", "x_end", "x_from", "x_start",
  "x_to", "y", "y_end", "y_from", "y_start", "y_to"
))

.onAttach <- function(libname, pkgname) {
  packageStartupMessage(
    "rePredRet: Retention Time Prediction by Direct System Mapping\n",
    "  - Use rePredRet_download() to get the latest predictions\n",
    "  - Use rePredRet_predict() to query predictions for a compound\n",
    "  - See ?rePredRet for more information"
  )
}
