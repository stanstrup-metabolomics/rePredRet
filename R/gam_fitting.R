#' @title Core GAM Fitting Functions
#' @description Internal functions for monotonically constrained GAM fitting
#'   used in bootstrap confidence interval estimation.
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
#'    relative to `y_max`, which down-weights outliers.
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


#' Calculate Bootstrap Confidence Intervals
#'
#' Computes percentile bootstrap confidence intervals from a boot object.
#'
#' @param loess_boot A boot object from `boot::boot()`.
#' @param alpha Significance level for CI (default 0.01 for 99% CI).
#'
#' @return A matrix with 3 columns: `[predicted, lower, upper]` at each
#'   newdata point.
#'
#' @importFrom boot boot.ci
#' @keywords internal
boot2ci <- function(loess_boot, alpha = 0.01) {
  # For each prediction point, calculate percentile CI
  temp <- list()
  for (i2 in seq_len(length(loess_boot$t0))) {
    temp[[i2]] <- list()
    temp[[i2]][[1]] <- i2
    temp[[i2]][[2]] <- loess_boot
  }

  ci <- lapply(temp, function(x) {
    # Get bootstrap estimates for this point, removing NAs
    boot_vals <- x[[2]]$t[, x[[1]]]
    boot_vals_clean <- boot_vals[!is.na(boot_vals)]

    # If too few valid bootstrap values, return NA
    if (length(boot_vals_clean) < 10) {
      return(c(NA, NA, NA))
    }

    # Check if there's any variance in bootstrap estimates
    diff_to_mean <- (mean(boot_vals_clean) - boot_vals_clean) / boot_vals_clean

    # If all bootstrap estimates are identical, return mean with no CI width
    if (all(abs(diff_to_mean) < 1e-10, na.rm = TRUE)) {
      ci <- rep(mean(boot_vals_clean), 3)
      return(ci)
    }

    # Calculate percentile CI using boot.ci
    # boot.ci handles NAs internally
    temp2 <- tryCatch({
      boot::boot.ci(
        x[[2]],
        index = x[[1]],
        type = "perc",
        conf = 1 - alpha
      )
    }, error = function(e) NULL)

    # If boot.ci failed, calculate percentile CI manually
    if (is.null(temp2)) {
      ci <- c(
        x[[2]]$t0[x[[1]]],
        quantile(boot_vals_clean, alpha/2),
        quantile(boot_vals_clean, 1 - alpha/2)
      )
      return(ci)
    }

    ci <- vector(mode = "numeric", length = 3)
    ci[1] <- temp2$t0
    ci[c(2, 3)] <- temp2$percent[, c(4, 5)]
    return(ci)
  })

  ci <- do.call(rbind, ci)
  return(ci)
}


#' Calculate Bootstrap Prediction Intervals
#'
#' Computes prediction intervals that account for both model uncertainty
#' and residual variance.
#'
#' @param loess_boot A boot object from `boot::boot()`.
#' @param newdata Numeric vector of x-values used for prediction.
#' @param alpha Significance level for PI (default 0.05 for 95% PI).
#'
#' @return A matrix with 3 columns: `[predicted, lower, upper]` at each
#'   newdata point.
#'
#' @details
#' Prediction intervals are wider than confidence intervals because they
#' account for:
#' 1. Uncertainty in the model fit (from bootstrap)
#' 2. Residual variance (deviation of individual points from the fit)
#'
#' @importFrom stats qnorm approx
#' @keywords internal
boot2ci_PI <- function(loess_boot, newdata, alpha = 0.05) {
  # SD at each x values (newdata) between all bootstrap iterations
  loess_sd_newdata <- apply(loess_boot$t, 2, sd)

  # Interpolate from newdata x values to values used for model building
  loess_sd <- numeric(length = nrow(loess_boot$data))
  res_sd <- numeric(length = nrow(loess_boot$data))

  for (i in seq_len(nrow(loess_boot$data))) {
    # Check if predicted points match model building x values
    if (any(newdata == loess_boot$data[i, 1])) {
      loc <- which(newdata == loess_boot$data[i, 1])
      res <- loess_boot$t[, loc] - loess_boot$data[i, 2]
      res_sd[i] <- sd(res)
    } else {
      # Interpolate linearly between two nearest prediction x values
      lower_idx <- which(newdata < loess_boot$data[i, 1])
      lower_idx <- lower_idx[length(lower_idx)]
      higher_idx <- which(newdata > loess_boot$data[i, 1])[1]

      pred_y <- numeric(length = nrow(loess_boot$t))
      for (i2 in seq_len(nrow(loess_boot$t))) {
        pred_y[i2] <- stats::approx(
          x = c(newdata[lower_idx], newdata[higher_idx]),
          y = c(loess_boot$t[i2, lower_idx], loess_boot$t[i2, higher_idx]),
          xout = loess_boot$data[i, 1]
        )$y
      }

      res <- pred_y - loess_boot$data[i, 2]
      res_sd[i] <- sd(res)

      loess_sd[i] <- stats::approx(
        x = c(newdata[lower_idx], newdata[higher_idx]),
        y = c(loess_sd_newdata[lower_idx], loess_sd_newdata[higher_idx]),
        xout = loess_boot$data[i, 1]
      )$y
    }
  }

  # Combine SDs (model uncertainty + residual variance)
  SD_combined <- sqrt(loess_sd^2 + res_sd^2)

  # Interpolate back to prediction points (newdata)
  SD_combined_newdata <- numeric(length = length(newdata))

  for (i in seq_len(length(newdata))) {
    if (any(newdata[i] == loess_boot$data[, 1])) {
      loc <- which(newdata[i] == loess_boot$data[, 1])
      # Take first if multiple identical x values
      SD_combined_newdata[i] <- SD_combined[loc][1]
    } else {
      lower_idx <- which(newdata[i] > loess_boot$data[, 1])
      lower_idx <- lower_idx[length(lower_idx)]
      higher_idx <- which(newdata[i] < loess_boot$data[, 1])[1]

      SD_combined_newdata[i] <- stats::approx(
        x = c(loess_boot$data[lower_idx, 1], loess_boot$data[higher_idx, 1]),
        y = c(SD_combined[lower_idx], SD_combined[higher_idx]),
        xout = newdata[i]
      )$y
    }
  }

  # Build CI matrix: [predicted, lower, upper]
  ci <- cbind(
    loess_boot$t0,
    loess_boot$t0 - SD_combined_newdata * stats::qnorm(1 - alpha / 2),
    loess_boot$t0 + SD_combined_newdata * stats::qnorm(1 - alpha / 2)
  )

  return(ci)
}
