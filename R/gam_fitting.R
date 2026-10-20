#' @title Core GAM Fitting Functions
#' @description Internal functions for monotonically constrained GAM fitting
#'   used by all prediction interval methods.
#' @name gam_fitting
#' @keywords internal
NULL

#' Monotonically Constrained GAM Bootstrap Function
#'
#' Fits a monotonically constrained GAM to bootstrap-resampled data.
#' This is the core fitting function used by `boot::boot()` for confidence
#' interval estimation.
#'
#' @param in_data A 2-column matrix where column 1 is RT in source system
#'   and column 2 is RT in target system.
#' @param inds Bootstrap resampling indices (provided by `boot()`).
#' @param newdata Numeric vector of x-values at which to predict.
#'
#' @return Numeric vector of predicted y-values at `newdata` points.
#'
#' @details
#' The algorithm:
#' 1. Resample data using bootstrap indices
#' 2. Add jitter if fewer than 4 unique x-values (needed for spline fitting)
#' 3. Fit initial unconstrained GAM to get smoothing parameter
#' 4. Weight observations by residuals using sigmoid function
#'    (downweights outliers)
#' 5. Build constrained spline with monotonic constraints
#' 6. Fit using penalized constrained least squares
#' 7. Predict at newdata points
#'
#' The monotonic constraint ensures that predicted RTs maintain the
#' elution order assumption (if compound A elutes before B in system 1,
#' A should also elute before B in system 2).
#'
#' @references
#' Stanstrup et al. (2015) "PredRet: Prediction of Retention Time by
#' Direct Mapping between Multiple Chromatographic Systems"
#' Anal. Chem. 87:9421-9428
#'
#' @keywords internal
gam_mono_con_fun <- function(in_data, inds, newdata) {
  fit <- fit_gam_mono(in_data[inds, 1], in_data[inds, 2], y_max = max(in_data[, 2]))
  predict_gam_mono(fit, newdata)
}


#' Fit a Monotonically Constrained GAM
#'
#' The monotone fit shared by all interval methods.
#'
#' @param x RTs in the source system.
#' @param y RTs in the target system.
#' @param y_max Value the absolute residuals are divided by before the
#'   sigmoid weighting (default `max(y)`).
#'
#' @return A list with the (possibly jittered) data `x` and `y`, the basis
#'   dimension `k`, the unconstrained fit `fit_unc`, the sigmoid weights `w`,
#'   the constrained smooth `sm` and its coefficients `p`.
#'
#' @details
#' 1. If there are fewer than 4 unique source RTs (rounded to 0.01), the
#'    source RTs are jittered by up to 0.01 with a fixed seed, so identical
#'    data always give the same model.
#' 2. An unconstrained thin plate regression spline GAM gives the smoothing
#'    parameter.
#' 3. Observations are weighted by a sigmoid of their unconstrained residual
#'    relative to `y_max`, which down-weights outliers: weights fall to 0.5
#'    at a residual of 10% of `y_max`, with steepness 30. These two values
#'    were set by eye for the original PredRet (2015) and have never been
#'    optimized.
#' 4. A cubic regression spline is fitted with monotonic increasing
#'    constraints by penalized constrained least squares. If the weighted fit
#'    fails it is refitted without weights.
#'
#' @importFrom mgcv gam smoothCon mono.con pcls Predict.matrix s
#' @importFrom pracma sigmoid
#' @keywords internal
fit_gam_mono <- function(x, y, y_max = max(y)) {
  # We need at least 4 unique x-values to do the fit
  if (length(unique(round(x, 2))) < 4) {
    x <- withr::with_seed(42, jitter(x, amount = 0.01))
  }

  dat <- data.frame(x = x, y = y)

  # Fit initial unconstrained GAM to get smoothing parameter
  k_val <- min(length(unique(round(x, 2))), 10)
  fit_unc <- mgcv::gam(y ~ s(x, k = k_val, bs = "tp"), data = dat)

  # Weight by residuals using sigmoid function
  # This downweights observations with large residuals (outliers)
  w <- abs(fit_unc$residuals) / y_max
  w <- pracma::sigmoid(w, a = -30, b = 0.1)

  # Build constrained spline basis
  sm <- mgcv::smoothCon(mgcv::s(x, k = k_val, bs = "cr"), dat, knots = NULL)[[1]]

  list(
    x = x,
    y = y,
    k = k_val,
    fit_unc = fit_unc,
    w = w,
    sm = sm,
    p = fit_mono_spline(sm, y, w, fit_unc$sp)
  )
}


#' Penalized Monotone Spline Coefficients
#'
#' @param sm A smooth from `mgcv::smoothCon()`.
#' @param y Response.
#' @param w Observation weights. If the weighted fit fails, it is refitted
#'   without weights.
#' @param sp Smoothing parameter.
#'
#' @return Coefficient vector.
#' @keywords internal
fit_mono_spline <- function(sm, y, w, sp) {
  con <- mgcv::mono.con(sm$xp)

  G <- list(
    X = sm$X,
    C = matrix(0, 0, 0),
    sp = sp,
    p = sm$xp,
    y = y,
    w = w,
    Ain = con$A,
    bin = con$b,
    S = sm$S,
    off = 0
  )

  p <- mgcv::pcls(G)

  # Some models fail when weighted - retry without weights
  if (any(is.na(p))) {
    G$w <- rep(1, length(G$w))
    p <- mgcv::pcls(G)
  }

  p
}


#' Predict from a Monotone GAM Fit
#'
#' @param fit Output of [fit_gam_mono()].
#' @param newdata Numeric vector of source RTs.
#'
#' @return Predicted target RTs, truncated at 0.
#' @keywords internal
predict_gam_mono <- function(fit, newdata) {
  y_pred <- as.numeric(mgcv::Predict.matrix(fit$sm, data.frame(x = newdata)) %*% fit$p)
  y_pred[y_pred < 0] <- 0  # RTs cannot be negative
  y_pred
}


#' Calculate Bootstrap Prediction Intervals
#'
#' Prediction intervals from bootstrap refits of the monotone GAM, using the
#' out-of-bag predictions as out-of-sample residuals.
#'
#' @param loess_boot A boot object from `boot::boot()` with [gam_mono_con_fun()],
#'   run with `newdata = c(newdata, rt_source)`, i.e. predicting at the grid
#'   followed by the calibration compounds. It must not be pruned yet.
#' @param fit Output of [fit_gam_mono()] on the full data.
#' @param newdata Numeric vector of x-values for the prediction grid.
#' @param alpha Significance level for PI (default 0.05 for 95% PI).
#'
#' @return A matrix with 3 columns: `[predicted, lower, upper]` at each
#'   newdata point.
#'
#' @details
#' Each compound is predicted by the refits whose resample did not contain
#' it (out-of-bag, on average 37% of the refits). Its residual from the mean
#' out-of-bag prediction is an out-of-sample error that includes the
#' uncertainty of the fitted curve. A compound that is in every resample
#' gets the leave-one-out residual of the unconstrained GAM instead. These
#' residuals then go through the same local scale model and conformal
#' calibration as [gam_adaptive_pi()] ([conformal_band()]), around the
#' full-data fit. This is the jackknife+-after-bootstrap idea of Kim, Xu &
#' Barber (2020).
#'
#' @references
#' Kim, B., Xu, C. & Barber, R. F. (2020). Predictive inference is free with
#' the jackknife+-after-bootstrap. *Advances in Neural Information Processing
#' Systems*, 33, 4138-4149.
#'
#' @keywords internal
boot2ci_PI <- function(loess_boot, fit, newdata, alpha = 0.05) {
  n_new <- length(newdata)
  n <- nrow(loess_boot$data)

  # Out-of-bag predictions at the calibration compounds
  pred_obs <- loess_boot$t[, n_new + seq_len(n), drop = FALSE]
  pred_obs[boot::boot.array(loess_boot) > 0] <- NA
  e <- loess_boot$data[, 2] - colMeans(pred_obs, na.rm = TRUE)

  # Compounds in every resample: leave-one-out residual of the unconstrained fit
  no_oob <- !is.finite(e)
  e[no_oob] <- (fit$fit_unc$residuals / pmax(1 - fit$fit_unc$hat, 0.05))[no_oob]

  conformal_band(fit, e, loess_boot$t0[seq_len(n_new)], newdata, alpha)
}
