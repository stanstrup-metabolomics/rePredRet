#' @title Core GAM Fitting Functions
#' @description Internal functions for monotonically constrained GAM fitting
#'   used in bootstrap confidence interval estimation.
#' @name gam_fitting
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
#' @importFrom mgcv gam smoothCon mono.con pcls Predict.matrix s
#' @importFrom pracma sigmoid
#' @keywords internal
gam_mono_con_fun <- function(in_data, inds, newdata) {
  x <- NULL  # Avoid R CMD check note

  # Resample data for this bootstrap iteration
  x.star <- in_data[, 1][inds]
  y.star <- in_data[, 2][inds]

  # We need at least 4 unique x-values to do the fit
  # Add small jitter if not the case

  if (length(unique(round(x.star, 2))) < 4) {
    x.star <- jitter(x.star, amount = 0.01)
  }

  dat <- data.frame(x = x.star, y = y.star)

  # Fit initial unconstrained GAM to get smoothing parameter
  k_val <- min(length(unique(round(x.star, 2))), 10)
  f.ug <- mgcv::gam(y ~ s(x, k = k_val, bs = "tp"), data = dat)

  # Weight by residuals using sigmoid function
  # This downweights observations with large residuals (outliers)
  w <- abs(f.ug$residuals) / max(in_data[, 2])
  w <- pracma::sigmoid(w, a = -30, b = 0.1)

  # Build constrained spline basis
  sm <- mgcv::smoothCon(
    mgcv::s(x, k = k_val, bs = "cr"),
    dat,
    knots = NULL
  )[[1]]

  # Get monotonic constraints
  con <- mgcv::mono.con(sm$xp)

  # Set up penalized constrained least squares problem
  G <- list(
    X = sm$X,
    C = matrix(0, 0, 0),
    sp = f.ug$sp,
    p = sm$xp,
    y = y.star,
    w = w
  )
  G$Ain <- con$A
  G$bin <- con$b
  G$S <- sm$S
  G$off <- 0

  # Fit constrained spline
  p <- mgcv::pcls(G)

  # Some models fail when weighted - retry without weights

if (any(is.na(p))) {
    G$w <- rep(1, length(G$w))
    p <- mgcv::pcls(G)
  }

  # Predict at newdata points
  fv <- mgcv::Predict.matrix(sm, data.frame(x = newdata)) %*% p

  y_pred <- as.numeric(fv)
  y_pred[y_pred < 0] <- 0  # RTs cannot be negative

  return(y_pred)
}


#' Calculate Bootstrap Confidence Intervals
#'
#' Computes percentile bootstrap confidence intervals from a boot object.
#'
#' @param loess_boot A boot object from `boot::boot()`.
#' @param alpha Significance level for CI (default 0.01 for 99% CI).
#'
#' @return A matrix with 3 columns: [predicted, lower, upper] at each
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
#' @return A matrix with 3 columns: [predicted, lower, upper] at each
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
