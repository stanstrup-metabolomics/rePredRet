# Tests for access to pre-built models (downloads replaced by local fixtures)

local_model_repo <- function(env = parent.frame()) {
  repo <- withr::local_tempdir(.local_envir = env)

  set.seed(1)
  rt1 <- sort(runif(20, 1, 20))
  rt_matrix <- cbind(rt_source = rt1, rt_target = 0.8 * rt1 + 1 + rnorm(20, sd = 0.1))
  attr(rt_matrix, "compounds") <- sprintf("InChI=1S/C%dH4", 1:20)
  model <- suppressMessages(build_model(rt_matrix, "0001", "0002", n_cores = 1))

  dir.create(file.path(repo, "models", "0001_to_0002"), recursive = TRUE)
  export_model_json(model, file.path(repo, "models", "0001_to_0002", "model.json"))
  utils::write.csv(data.frame(
    sys1_id = c("0001", "0002"), sys2_id = c("0002", "0001"),
    n_compounds = c(20, 20), median_ci_width = c(1, 1), median_error = c(0.1, 0.1)
  ), file.path(repo, "models", "model_index.csv"), row.names = FALSE)

  dir.create(file.path(repo, "RepoRT"))
  writeLines(c(
    "id\tname\tmethod.type",
    "0001\tFEM_short\tRP",
    "0002\tFEM_long\tRP",
    "0003\tOther\tHILIC"
  ), file.path(repo, "RepoRT", "studies.tsv"))

  testthat::local_mocked_bindings(
    .fetch_file = function(url, cache_path, ...) {
      path <- file.path(repo, cache_path)
      if (!file.exists(path)) stop("Could not download ", url, call. = FALSE)
      path
    },
    .env = env
  )
  repo
}


test_that("rePredRet_models reads the index with character system IDs", {
  local_model_repo()
  models <- rePredRet_models()

  expect_equal(nrow(models), 2)
  expect_equal(models$sys1_id, c("0001", "0002"))
})


test_that("rePredRet_systems returns only systems with models", {
  local_model_repo()
  systems <- rePredRet_systems()

  expect_equal(systems$id, c("0001", "0002"))
  expect_equal(systems$name, c("FEM_short", "FEM_long"))
})


test_that("rePredRet_model returns prediction grid, calibration and stats", {
  local_model_repo()
  model <- rePredRet_model("0001", "0002")

  expect_equal(model$sys1_id, "0001")
  expect_equal(model$n_compounds, 20)
  expect_named(model$predictions, c("x", "pred", "lower", "upper"))
  expect_equal(nrow(model$calibration), 20)
  expect_true(is.numeric(model$stats$median_ci_width))
})


test_that("rePredRet_model gives a clear error for unknown models", {
  local_model_repo()
  expect_error(rePredRet_model("0001", "9999"), "Model not found: 0001 -> 9999")
})


test_that("plot_model builds a ggplot with system names", {
  local_model_repo()
  p <- plot_model("0001", "0002")

  expect_s3_class(p, "ggplot")
  expect_match(p$labels$title, "FEM_short")
})


test_that("rePredRet_plot_network uses the model index", {
  local_model_repo()
  p <- rePredRet_plot_network(min_compounds = 10)
  expect_length(p$layers, 3)
  expect_length(rePredRet_plot_network(min_compounds = 10, show_labels = FALSE)$layers, 2)

  expect_s3_class(p, "ggplot")
  expect_error(rePredRet_plot_network(min_compounds = 100), "No models meet")
})


test_that(".fetch_file caches downloads and refreshes on request", {
  cache_dir <- withr::local_tempdir()
  local_mocked_bindings(user_cache_dir = function(...) cache_dir, .package = "rappdirs")

  src <- withr::local_tempfile(fileext = ".txt")
  writeLines("v1", src)
  url <- paste0("file:///", normalizePath(src, winslash = "/"))

  path <- .fetch_file(url, "sub/file.txt")
  expect_equal(readLines(path), "v1")

  writeLines("v2", src)
  expect_equal(readLines(.fetch_file(url, "sub/file.txt")), "v1")
  expect_equal(readLines(.fetch_file(url, "sub/file.txt", refresh = TRUE)), "v2")

  expect_error(.fetch_file("file:///does/not/exist", "missing.txt"), "Could not download")
  expect_false(file.exists(file.path(cache_dir, "missing.txt.part")))
})
