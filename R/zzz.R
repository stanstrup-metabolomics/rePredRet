#' @title Package Initialization
#' @description Package startup messages and initialization.
#' @name zzz
NULL

.onAttach <- function(libname, pkgname) {
  packageStartupMessage(
    "rePredRet: Retention Time Prediction by Direct System Mapping\n",
    "  - Use rePredRet_download() to get the latest predictions\n",
    "  - Use rePredRet_predict() to query predictions for a compound\n",
    "  - See ?rePredRet for more information"
  )
}
