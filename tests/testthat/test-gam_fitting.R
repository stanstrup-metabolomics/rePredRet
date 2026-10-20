# Tests for GAM fitting functions

test_that("gam_mono_con_fun returns correct output", {
  # Create simple test data
  set.seed(123)
  n <- 50
  x <- sort(runif(n, 1, 20))
  y <- 0.5 * x + rnorm(n, sd = 0.5)  # Simple linear relationship with noise

  rt_matrix <- cbind(x, y)
  newdata <- seq(min(x), max(x), length.out = 100)

  # Run single iteration (indices = 1:n means use all data)
  result <- gam_mono_con_fun(rt_matrix, 1:n, newdata)

  # Check output

  expect_length(result, 100)
  expect_true(all(result >= 0))  # No negative predictions
  expect_true(all(is.finite(result)))  # No NA/Inf values

  # Check monotonicity (should be increasing)
  expect_true(all(diff(result) >= -0.01))  # Allow tiny numerical error
})


test_that("boot2ci returns correct structure", {
  skip_if_not_installed("boot")

  # Create simple test data
  set.seed(456)
  n <- 30
  x <- sort(runif(n, 1, 15))
  y <- 0.8 * x + rnorm(n, sd = 0.3)

  rt_matrix <- cbind(x, y)
  newdata <- seq(min(x), max(x), length.out = 50)

  # Run bootstrap with few iterations for speed
  loess_boot <- boot::boot(
    rt_matrix,
    gam_mono_con_fun,
    R = 50,  # Few iterations for testing
    newdata = newdata
  )

  # Calculate CIs
  ci <- boot2ci(loess_boot, alpha = 0.05)

  # Check output structure
  expect_equal(nrow(ci), 50)
  expect_equal(ncol(ci), 3)

  # Check CI properties
  expect_true(all(ci[, 2] <= ci[, 1]))  # Lower <= predicted

  expect_true(all(ci[, 1] <= ci[, 3]))  # Predicted <= upper
})


test_that("prune_boot_object reduces size", {
  skip_if_not_installed("boot")

  # Create minimal boot object
  set.seed(789)
  n <- 20
  x <- sort(runif(n, 1, 10))
  y <- x + rnorm(n, sd = 0.2)
  rt_matrix <- cbind(x, y)
  newdata <- seq(min(x), max(x), length.out = 20)

  loess_boot <- boot::boot(
    rt_matrix,
    gam_mono_con_fun,
    R = 20,
    newdata = newdata
  )

  # Prune (use ::: to access internal function)
  pruned <- rePredRet:::prune_boot_object(loess_boot)

  # Check that large components are removed
  # Use [[ to avoid partial matching ($ would match "t" to "t0")
  expect_null(pruned[["statistic"]])
  expect_null(pruned[["t"]])
  expect_false("t" %in% names(pruned))
  expect_false("statistic" %in% names(pruned))

  # Check that essential components remain
  expect_false(is.null(pruned[["t0"]]))
  expect_false(is.null(pruned[["data"]]))
})


test_that("fit_gam_mono jitters few unique source RTs reproducibly", {
  x <- rep(c(1, 2, 3), length.out = 12)
  y <- x + seq(0, 0.5, length.out = 12)

  fit1 <- fit_gam_mono(x, y)
  fit2 <- fit_gam_mono(x, y)

  expect_identical(fit1$x, fit2$x)
  expect_identical(fit1$p, fit2$p)
  expect_gte(length(unique(round(fit1$x, 2))), 4)
})


# Heteroscedastic test data: scatter grows with RT
make_het_data <- function(n, seed = 1) {
  set.seed(seed)
  x <- sort(runif(n, 1, 20))
  cbind(x, 0.8 * x + 1 + rnorm(n, sd = 0.02 + 0.04 * x))
}


test_that("gam_adaptive_pi returns ordered, monotone intervals", {
  rt_matrix <- make_het_data(60)
  newdata <- seq(min(rt_matrix[, 1]), max(rt_matrix[, 1]), length.out = 200)

  pi <- gam_adaptive_pi(rt_matrix, newdata)

  expect_equal(dim(pi), c(200, 3))
  expect_equal(colnames(pi), c("predicted", "lower", "upper"))
  expect_true(all(is.finite(pi)))
  expect_true(all(pi[, "lower"] <= pi[, "predicted"] + 1e-8))
  expect_true(all(pi[, "predicted"] <= pi[, "upper"] + 1e-8))
  expect_true(all(diff(pi[, "predicted"]) >= -1e-8))
  expect_true(all(diff(pi[, "lower"]) >= -1e-8))
  expect_true(all(diff(pi[, "upper"]) >= -1e-8))

  # Same input gives the same intervals
  expect_identical(pi, gam_adaptive_pi(rt_matrix, newdata))
})


test_that("gam_adaptive_pi width follows the local scatter", {
  rt_matrix <- make_het_data(150)
  newdata <- c(3, 18)

  pi <- gam_adaptive_pi(rt_matrix, newdata)
  width <- pi[, "upper"] - pi[, "lower"]

  # True 95% widths are about 0.56 at RT 3 and 2.9 at RT 18
  expect_gt(width[2] / width[1], 2.5)
})


test_that("gam_adaptive_pi covers about 95% of new compounds", {
  train <- make_het_data(200, seed = 2)
  test <- make_het_data(2000, seed = 3)
  test <- test[test[, 1] >= min(train[, 1]) & test[, 1] <= max(train[, 1]), ]

  pi <- gam_adaptive_pi(train, test[, 1])
  covered <- test[, 2] >= pi[, "lower"] & test[, 2] <= pi[, "upper"]

  expect_gt(mean(covered), 0.92)
  expect_lt(mean(covered), 0.98)
})


test_that("gam_adaptive_pi handles few compounds and outliers", {
  set.seed(4)
  x <- sort(runif(10, 1, 15))
  y <- 0.9 * x + rnorm(10, sd = 0.2)
  y[5] <- y[5] + 8  # Gross outlier
  newdata <- seq(min(x), max(x), length.out = 100)

  pi <- gam_adaptive_pi(cbind(x, y), newdata)

  expect_true(all(is.finite(pi)))
  expect_true(all(diff(pi[, "lower"]) >= -1e-8))
  expect_true(all(diff(pi[, "upper"]) >= -1e-8))
  # The band must not explode to many times the RT range
  expect_lt(stats::median(pi[, "upper"] - pi[, "lower"]), 2 * diff(range(y)))
})


test_that("build_model uses adaptive_pi by default", {
  model <- suppressMessages(build_model(make_het_data(30), "a", "b"))

  expect_equal(model$status, "success")
  expect_equal(model$method, "adaptive_pi")
  expect_equal(dim(model$ci), c(1000, 3))
  expect_true(model$stats$median_ci_width > 0)
})
