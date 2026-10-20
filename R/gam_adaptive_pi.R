#' @title Locally Adaptive Prediction Intervals
#' @description Prediction intervals whose width follows the local scatter
#'   around the monotone GAM, calibrated by cross-validated conformal
#'   prediction.
#' @name gam_adaptive_pi
NULL

#' Locally Adaptive Monotonic Prediction Intervals
#'
#' Computes prediction intervals for a monotonically constrained GAM whose
#' width varies along the source RT with the local scatter of the
#' calibration compounds.
#'
#' @param rt_matrix A 2-column matrix `[rt_source, rt_target]`.
#' @param newdata Numeric vector of x-values for predictions.
#' @param alpha Significance level (default 0.05 for 95% PI).
#' @param n_folds Number of cross-validation folds (default 5).
#'
#' @return A matrix with 3 columns: `[predicted, lower, upper]`.
#'
#' @details
#' Algorithm:
#' 1. Fit the monotone GAM ([fit_gam_mono()]).
#' 2. Cross-validated residuals: the compounds are split into `n_folds`
#'    folds, interleaved along the source RT. Each fold is predicted by the
#'    monotone spline refitted on the other folds with the smoothing
#'    parameter and weights of the full fit. If a fold cannot be refitted,
#'    its leave-one-out residuals of the unconstrained GAM are used instead.
#' 3. Local scale s(x): a Gamma GAM (log link, REML, at most 5 basis
#'    functions) of the absolute cross-validated residuals, weighted by the
#'    sigmoid outlier weights so that single outliers do not widen the band
#'    locally. s(x) is floored at 30% of the global weighted scale.
#' 4. Calibration: the multiplier q is the split-conformal quantile (the
#'    `ceiling((n + 1) * (1 - alpha))`-th smallest value) of the normalized
#'    residuals |e| / s(x). With fewer than about `1 / alpha` compounds this is
#'    their maximum.
#' 5. Interval: prediction +/- q * s(x). Each bound is made non-decreasing by
#'    isotonic regression, which keeps lower <= prediction <= upper. The
#'    lower bound is truncated at 0.
#'
#' Because the residuals are out-of-fold, q already includes the uncertainty
#' of the fitted curve, so no separate standard error term is added.
#'
#' @references
#' Lei, J., G'Sell, M., Rinaldo, A., Tibshirani, R. J. & Wasserman, L. (2018).
#' Distribution-free predictive inference for regression. *Journal of the
#' American Statistical Association*, 113(523), 1094-1111.
#'
#' @importFrom stats qnorm isoreg predict Gamma
#' @export
gam_adaptive_pi <- function(rt_matrix, newdata, alpha = 0.05, n_folds = 5) {

  fit <- fit_gam_mono(rt_matrix[, 1], rt_matrix[, 2])
  n <- length(fit$y)

  # Step 2: cross-validated residuals
  e <- cv_residuals_gam_mono(fit, n_folds)

  # Step 3: local scale of the residuals
  scale_fun <- local_residual_scale(fit$x, e, fit$w, fit$k)
  s_obs <- scale_fun(fit$x)

  # Step 4: conformal multiplier
  scores <- sort(abs(e) / s_obs)
  q <- scores[min(ceiling((n + 1) * (1 - alpha)), n)]

  # Step 5: interval with monotone bounds
  mu <- predict_gam_mono(fit, newdata)
  half_width <- q * scale_fun(newdata)
  ord <- order(newdata)
  lower <- upper <- numeric(length(newdata))
  lower[ord] <- stats::isoreg(newdata[ord], (mu - half_width)[ord])$yf
  upper[ord] <- stats::isoreg(newdata[ord], (mu + half_width)[ord])$yf
  lower <- pmax(lower, 0)  # RTs cannot be negative

  pi <- cbind(mu, lower, upper)
  colnames(pi) <- c("predicted", "lower", "upper")

  return(pi)
}


#' Cross-Validated Residuals of a Monotone GAM
#'
#' @param fit Output of [fit_gam_mono()].
#' @param n_folds Number of folds, interleaved along x.
#'
#' @return Out-of-fold residuals, one per compound.
#' @keywords internal
cv_residuals_gam_mono <- function(fit, n_folds = 5) {
  x <- NULL  # Avoid R CMD check note

  n <- length(fit$x)
  folds <- rep_len(seq_len(n_folds), n)[order(order(fit$x))]

  # Fallback: leave-one-out residuals of the unconstrained fit
  e <- fit$fit_unc$residuals / pmax(1 - fit$fit_unc$hat, 0.05)

  for (f in seq_len(n_folds)) {
    test <- folds == f
    train <- data.frame(x = fit$x[!test], y = fit$y[!test])
    k_train <- min(fit$k, length(unique(train$x)))
    if (k_train < 4) next

    pred <- tryCatch({
      sm <- mgcv::smoothCon(mgcv::s(x, k = k_train, bs = "cr"), train, knots = NULL)[[1]]
      p <- fit_mono_spline(sm, train$y, fit$w[!test], fit$fit_unc$sp)
      pmax(as.numeric(mgcv::Predict.matrix(sm, data.frame(x = fit$x[test])) %*% p), 0)
    }, error = function(err) NULL)

    if (!is.null(pred) && all(is.finite(pred))) {
      e[test] <- fit$y[test] - pred
    }
  }

  e
}


#' Local Scale of Residuals
#'
#' Models the absolute residuals as a smooth function of x.
#'
#' @param x Source RTs.
#' @param e Residuals.
#' @param w Observation weights (sigmoid outlier weights).
#' @param k Basis dimension of the mean fit; the scale smooth uses at most 5.
#' @param floor_frac The scale is floored at this fraction of the global
#'   weighted scale (default 0.3).
#'
#' @return A function returning the residual standard deviation at new x.
#'
#' @details
#' A Gamma GAM with log link is fitted to |e| with the weights `w`. Under
#' normal errors E|e| = sigma * sqrt(2 / pi), so the fitted values are scaled
#' by sqrt(pi / 2). Absolute residuals are floored at 5% of their mean, since
#' exact zeros (e.g. RTs reported with few decimals) are not allowed by the
#' Gamma family. If the GAM fails, the global weighted scale is used.
#'
#' @keywords internal
local_residual_scale <- function(x, e, w, k, floor_frac = 0.3) {
  a <- abs(e)
  s_global <- max(sqrt(pi / 2) * sum(w * a) / sum(w), 1e-4)

  scale_fun <- tryCatch({
    dat <- data.frame(x = x, a = pmax(a, 0.05 * mean(a) + 1e-6))
    fit <- mgcv::gam(
      a ~ s(x, k = max(3, min(k, 5)), bs = "tp"),
      data = dat,
      family = stats::Gamma(link = "log"),
      weights = w / mean(w),
      method = "REML"
    )
    function(newx) {
      sqrt(pi / 2) * as.numeric(stats::predict(fit, data.frame(x = newx), type = "response"))
    }
  }, error = function(err) {
    function(newx) rep(s_global, length(newx))
  })

  function(newx) pmax(scale_fun(newx), floor_frac * s_global)
}
