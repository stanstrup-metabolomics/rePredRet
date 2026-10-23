#' @title Package Initialization
#' @description Package startup messages and initialization.
#' @name zzz
#' @keywords internal
NULL

# Suppress R CMD check NOTEs for ggplot2 NSE and data.table/dplyr column names
utils::globalVariables(c(
  "ci_lower", "ci_upper", "ci_width", "community", "compound", "compound_label", "compound_name",
  "edge", "from", "from_id", "id", "inchi_clean", "key", "lower", "median_ci_width", "median_ci_width_rel",
  "median_error_abs", "median_error_rel", "method_type", "n_calibration",
  "n_compounds", "n_predictions", "name", "node", "node_size", "pred", "predicted", "predicted_rt", "r_sq_display",
  "r_squared", "rt", "rt.1", "rt.2", "rt_source", "rt_target", "rt_user",
  "source_system", "sys1_id", "sys2_id", "system", "to", "to_id", "upper", "x", "x_end", "x_from", "x_start",
  "x_to", "y", "y_end", "y_from", "y_start", "y_to"
))

.onAttach <- function(libname, pkgname) {
  packageStartupMessage(
    "rePredRet: Retention Time Prediction by Direct System Mapping\n",
    "  - Use build_user_models() and predict_rt_user() to predict RTs for your data\n",
    "  - Use rePredRet_models() and plot_model() to explore pre-built models\n",
    "  - See ?rePredRet for more information"
  )
}
