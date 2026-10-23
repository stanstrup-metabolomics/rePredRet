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
#'    locally. Between the weighted compounds s(x) is interpolated linearly
#'    (see [local_residual_scale()]). s(x) is floored at 30% of the global
#'    weighted scale and capped at the largest residual.
#' 4. Calibration: the multiplier q is the weighted split-conformal quantile
#'    of the normalized residuals |e| / s(x), where each compound counts by
#'    its sigmoid outlier weight ([weighted_conformal_quantile()]). Outliers
#'    are taken to be wrong annotations, so the interval is meant to cover
#'    `1 - alpha` of the well-annotated compounds, not the outliers. With
#'    fewer than about `1 / alpha` compounds q is the largest score of the
#'    compounds that are not down-weighted.
#' 5. Interval: prediction +/- q * s(x). Each bound is made non-decreasing by
#'    isotonic regression, which keeps lower <= prediction <= upper. The
#'    lower bound is truncated at 0.
#'
#' With few compounds the interval covers less than its nominal level: in a
#' cross-validation benchmark on RepoRT pairs with 10-15 shared compounds a
#' 95% interval covered about 85% (about 95% from 30 compounds up).
#'
#' Because the residuals are out-of-fold, q already includes the uncertainty
#' of the fitted curve, so no separate standard error term is added.
#'
#' The idea for this method, and the problems with the earlier interval
#' methods that it fixes, came from Barbora Hamidová's evaluation of
#' rePredRet in her MSc thesis (2026).
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

  # Step 2: cross-validated residuals
  e <- cv_residuals_gam_mono(fit, n_folds)

  # Steps 3-5: local scale, conformal multiplier, monotone bounds
  conformal_band(fit, e, predict_gam_mono(fit, newdata), newdata, alpha)
}


#' Conformal Band with Local Scale
#'
#' Turns out-of-sample residuals into a prediction interval: a local scale
#' model of the residuals ([local_residual_scale()]), a split-conformal
#' multiplier and monotone bounds. Shared by [gam_adaptive_pi()] and the
#' bootstrap method ([boot2ci_PI()]).
#'
#' @param fit Output of [fit_gam_mono()] (gives x, weights and k).
#' @param e Out-of-sample residuals at `fit$x`.
#' @param mu Predictions at `newdata`.
#' @param newdata Numeric vector of source RTs.
#' @param alpha Significance level.
#'
#' @return A matrix with 3 columns: `[predicted, lower, upper]`.
#' @keywords internal
conformal_band <- function(fit, e, mu, newdata, alpha) {
  # Step 3: local scale of the residuals
  scale_fun <- local_residual_scale(fit$x, e, fit$w, fit$k)

  # Step 4: conformal multiplier
  q <- weighted_conformal_quantile(abs(e) / scale_fun(fit$x), fit$w, alpha)

  # Step 5: interval with monotone bounds
  monotone_interval(newdata, mu, q * scale_fun(newdata))
}


#' Weighted Conformal Quantile
#'
#' The `1 - alpha` quantile of the scores, where each compound counts by its
#' outlier weight, so that compounds the fit treats as outliers (wrong
#' annotations) do not set the width of the interval.
#'
#' @param scores Normalized residuals |e| / s(x).
#' @param w Sigmoid outlier weights of the compounds.
#' @param alpha Significance level.
#'
#' @return The multiplier q.
#'
#' @details
#' The smallest score at which the cumulative weight of the sorted scores
#' reaches `(1 - alpha) * (sum(w) + max(w))`; the extra `max(w)` stands for
#' the new compound, assumed to be well annotated. With equal weights this is
#' the split-conformal `ceiling((n + 1) * (1 - alpha))`-th smallest score.
#' If the total weight is too small to reach that level (fewer than about
#' `1 / alpha` well-fitting compounds), the largest score of the compounds
#' with at least half the maximum weight is used.
#'
#' @references
#' Tibshirani, R. J., Barber, R. F., Candès, E. J. & Ramdas, A. (2019).
#' Conformal prediction under covariate shift. *Advances in Neural
#' Information Processing Systems* 32.
#'
#' @keywords internal
weighted_conformal_quantile <- function(scores, w, alpha) {
  o <- order(scores)
  j <- which(cumsum(w[o]) >= (1 - alpha) * (sum(w) + max(w)))[1]
  if (!is.na(j)) return(scores[o][j])
  max(scores[w >= 0.5 * max(w)])
}


#' Monotone Interval from a Fit and Half-Widths
#'
#' Each bound (`mu -/+ half_width`) is made non-decreasing in `newdata` by
#' isotonic regression, which keeps lower <= mu <= upper when `mu` is
#' monotone. The lower bound is truncated at 0.
#'
#' @param newdata Numeric vector of source RTs.
#' @param mu Monotone predictions at `newdata`.
#' @param half_width Half-widths at `newdata`.
#'
#' @return A matrix with 3 columns: `[predicted, lower, upper]`.
#' @importFrom stats isoreg
#' @keywords internal
monotone_interval <- function(newdata, mu, half_width) {
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
#' The GAM is only evaluated at the compounds that carry weight (`w` at least
#' 5% of its maximum) and interpolated linearly between them; beyond them the
#' scale is held constant. In gaps between compounds and next to
#' down-weighted outliers the log-link smooth is unconstrained and can grow
#' without bound (scales above 1e30 min were seen). The scale is capped at
#' sqrt(pi / 2) times the largest absolute residual.
#'
#' @keywords internal
local_residual_scale <- function(x, e, w, k, floor_frac = 0.3) {
  a <- abs(e)
  s_global <- max(sqrt(pi / 2) * sum(w * a) / sum(w), 1e-4)
  s_max <- max(sqrt(pi / 2) * max(a), s_global)

  scale_fun <- tryCatch({
    dat <- data.frame(x = x, a = pmax(a, 0.05 * mean(a) + 1e-6))
    fit <- mgcv::gam(
      a ~ s(x, k = max(3, min(k, 5)), bs = "tp"),
      data = dat,
      family = stats::Gamma(link = "log"),
      weights = w / mean(w),
      method = "REML"
    )
    xs <- sort(unique(x[w >= 0.05 * max(w)]))
    s <- sqrt(pi / 2) * as.numeric(stats::predict(fit, data.frame(x = xs), type = "response"))
    if (length(xs) < 2 || any(!is.finite(s))) stop("too few weighted compounds")
    function(newx) stats::approx(xs, s, xout = newx, rule = 2)$y
  }, error = function(err) {
    function(newx) rep(s_global, length(newx))
  })

  function(newx) pmin(pmax(scale_fun(newx), floor_frac * s_global), s_max)
}
