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
