#' @title Fast Prediction Intervals with a Global Residual Scale
#' @description Fast prediction intervals from one monotone GAM fit and a
#'   single residual variance.
#' @name gam_fast_ci
NULL

#' Fast Monotonic Prediction Intervals
#'
#' Computes prediction intervals of nearly constant width from a single
#' monotone GAM fit. About four times faster than [gam_adaptive_pi()], whose
#' width follows the local scatter.
#'
#' @param rt_matrix A 2-column matrix `[rt_source, rt_target]`.
#' @param newdata Numeric vector of x-values for predictions.
#' @param alpha Significance level (default 0.05 for 95% PI).
#'
#' @return A matrix with 3 columns: `[predicted, lower, upper]`.
#'
#' @details
#' Algorithm:
#' 1. Fit the monotone GAM ([fit_gam_mono()]).
#' 2. SE: pointwise standard errors of the unconstrained GAM.
#' 3. sigma^2: mean squared leave-one-out residual. The residuals of the
#'    constrained fit are divided by `1 - h`, with `h` the leverages of the
#'    unconstrained GAM. In-sample residuals underestimate the error for new
#'    compounds, most of all with few compounds.
#' 4. Interval: prediction +/- z * sqrt(SE^2 + sigma^2), with bounds made
#'    non-decreasing by isotonic regression and the lower bound truncated
#'    at 0 ([monotone_interval()]).
#'
#' @importFrom stats qnorm predict
#' @export
gam_fast_ci <- function(rt_matrix, newdata, alpha = 0.05) {

  fit <- fit_gam_mono(rt_matrix[, 1], rt_matrix[, 2])

  # Step 2: SE estimates from the unconstrained fit
  se_param <- stats::predict(fit$fit_unc, newdata = data.frame(x = newdata), se.fit = TRUE)$se.fit

  # Step 3: observation variance from leave-one-out residuals
  residuals_const <- fit$y - predict_gam_mono(fit, fit$x)
  residuals_loo <- residuals_const / pmax(1 - fit$fit_unc$hat, 0.05)
  sigma_sq <- mean(residuals_loo^2)

  # Step 4: interval with monotone bounds
  z <- stats::qnorm(1 - alpha / 2)
  monotone_interval(
    newdata,
    predict_gam_mono(fit, newdata),
    z * sqrt(se_param^2 + sigma_sq)
  )
}
