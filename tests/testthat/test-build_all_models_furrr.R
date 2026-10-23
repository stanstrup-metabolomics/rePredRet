# Tests for build_all_models_furrr() incremental builds

make_dataset <- function(inchis, rt) {
  list(rtdata = data.frame(inchi.std = inchis, rt = rt, stringsAsFactors = FALSE))
}

make_report_data <- function(n_compounds = 30, shift = 0) {
  set.seed(42)
  inchis <- sprintf("InChI=1S/C%dH%dO", seq_len(n_compounds), seq_len(n_compounds))
  rt1 <- sort(runif(n_compounds, 1, 20))
  list(datasets = list(
    "0001" = make_dataset(inchis, rt1),
    "0002" = make_dataset(inchis, 0.8 * rt1 + 1 + rnorm(n_compounds, sd = 0.1)),
    "0003" = make_dataset(inchis, 1.2 * rt1 + shift + rnorm(n_compounds, sd = 0.1))
  ))
}

build_quietly <- function(report_data, export_dir, ...) {
  suppressMessages(build_all_models_furrr(
    report_data, n_workers = 1, export_dir = export_dir, verbose = FALSE, ...
  ))
}


test_that("common_rt_matrix matches get_common_compounds", {
  d1 <- make_dataset(
    c("InChI=1S/B", "InChI=1S/A", "InChI=1S/A", "InChI=1S/C/t1-/m0/s1", "InChI=1S/X", NA),
    c(2.0, 1.0, 1.2, 3.0, 9.0, 4.0)
  )
  d2 <- make_dataset(
    c("InChI=1S/A", "InChI=1S/C", "InChI=1S/B", "InChI=1S/B", "InChI=1S/Y", NA),
    c(5.0, 6.0, 7.0, 7.4, 8.0, 1.0)
  )
  lookup <- prepare_rt_lookup(list(d1 = d1, d2 = d2))

  fast <- common_rt_matrix(lookup$d1, lookup$d2)
  ref <- get_common_compounds(d1, d2)

  expect_equal(unname(unclass(fast)[, 1:2]), unname(as.matrix(ref)[, 1:2]), ignore_attr = TRUE)
  expect_identical(colnames(fast), colnames(ref))
  expect_identical(attr(fast, "compounds"), attr(ref, "compounds"))
})


test_that("common_rt_matrix returns an empty matrix without overlap", {
  lookup <- prepare_rt_lookup(list(
    a = make_dataset(c("InChI=1S/A"), 1),
    b = make_dataset(c("InChI=1S/B"), 2)
  ))
  result <- common_rt_matrix(lookup$a, lookup$b)
  expect_equal(dim(result), c(0, 2))
})


test_that("find_model_pairs respects min_compounds", {
  rd <- make_report_data(n_compounds = 12)
  lookup <- prepare_rt_lookup(rd$datasets)

  expect_length(find_model_pairs(lookup, min_compounds = 10), 6)
  expect_length(find_model_pairs(lookup, min_compounds = 13), 0)
})


test_that("find_model_pairs only pairs datasets of the same group", {
  rd <- make_report_data()
  lookup <- prepare_rt_lookup(rd$datasets)

  pairs <- find_model_pairs(lookup, 10, groups = c("RP", "RP", "HILIC"))
  keys <- purrr::map_chr(pairs, function(p) paste(p$id1, p$id2))
  expect_setequal(keys, c("0001 0002", "0002 0001"))
  # A missing group forms a group of its own
  expect_length(find_model_pairs(lookup, 10, groups = c("RP", NA, NA)), 2)
})


test_that("build_all_models_furrr builds within method types and writes systems.csv", {
  export_dir <- withr::local_tempdir()
  rd <- make_report_data()
  rd$studies <- data.frame(id = c("0001", "0002", "0003"), name = c("a", "b", "c"),
                           method.type = c("RP", "RP", "HILIC"))

  run <- build_quietly(rd, export_dir)
  expect_setequal(paste(run$index$sys1_id, run$index$sys2_id), c("0001 0002", "0002 0001"))

  systems <- utils::read.csv(file.path(export_dir, "systems.csv"), colClasses = c(sys_id = "character"))
  expect_equal(systems$sys_id, c("0001", "0002", "0003"))
  expect_equal(systems$method_type, c("RP", "RP", "HILIC"))
  expect_equal(systems$n_compounds, c(30, 30, 30))
  expect_equal(systems$n_models, c(2, 2, 0))
})


test_that("build_all_models_furrr only rebuilds what changed", {
  export_dir <- withr::local_tempdir()
  rd <- make_report_data()

  # First run builds all ordered pairs
  run1 <- build_quietly(rd, export_dir, report_version = "abc123")
  expect_equal(run1$stats$built, 6)
  expect_true(run1$changed)
  expect_equal(nrow(run1$index), 6)
  expect_true(file.exists(file.path(export_dir, "0001_to_0002", "model.json")))
  expect_equal(readRDS(file.path(export_dir, "build_cache.rds"))$report_version, "abc123")

  # Second run with identical data builds nothing and keeps the index
  run2 <- build_quietly(rd, export_dir)
  expect_equal(run2$stats$built, 0)
  expect_false(run2$changed)
  expect_equal(run2$index, run1$index)
  expect_equal(readRDS(file.path(export_dir, "build_cache.rds"))$report_version, "abc123")

  # Changing one dataset rebuilds only the pairs involving it
  rd_changed <- make_report_data(shift = 2)
  run3 <- build_quietly(rd_changed, export_dir)
  expect_equal(run3$stats$built, 4)
  expect_true(run3$changed)
  expect_equal(nrow(run3$index), 6)

  # Removing a dataset deletes its models from disk and index
  rd_removed <- rd_changed
  rd_removed$datasets[["0003"]] <- NULL
  run4 <- build_quietly(rd_removed, export_dir)
  expect_equal(run4$stats$built, 0)
  expect_true(run4$changed)
  expect_equal(nrow(run4$index), 2)
  expect_false(dir.exists(file.path(export_dir, "0001_to_0003")))

  index_csv <- utils::read.csv(file.path(export_dir, "model_index.csv"),
                               colClasses = c(sys1_id = "character", sys2_id = "character"))
  expect_equal(index_csv$sys1_id, c("0001", "0002"))
})


test_that("missing model files are rebuilt and missing index rows recovered", {
  export_dir <- withr::local_tempdir()
  rd <- make_report_data()
  build_quietly(rd, export_dir)

  # Simulate an interrupted run: one model lost, index lost
  unlink(file.path(export_dir, "0002_to_0001"), recursive = TRUE)
  unlink(file.path(export_dir, "model_index.csv"))

  run <- build_quietly(rd, export_dir)
  expect_equal(run$stats$built, 1)
  expect_equal(nrow(run$index), 6)
  expect_true(all(!is.na(run$index$median_ci_width)))
})


test_that("read_model_index reads the older from_id/to_id format", {
  index_file <- withr::local_tempfile(fileext = ".csv")
  utils::write.csv(data.frame(
    model_key = "0001_to_0002", from_id = "0001", to_id = "0002", method_type = NA,
    n_compounds = 67, median_error = 0.3, median_ci_width = 16, mean_ci_width = 16,
    build_time = NA
  ), index_file, row.names = FALSE)

  idx <- read_model_index(index_file)
  expect_equal(idx$model_key, "0001_to_0002")
  expect_equal(idx$sys1_id, "0001")
  expect_equal(idx$n_compounds, 67)
})


test_that("models built with another interval method are rebuilt", {
  export_dir <- withr::local_tempdir()
  rd <- make_report_data()
  build_quietly(rd, export_dir)

  # Pretend the existing models were built with an older method
  cache_file <- file.path(export_dir, "build_cache.rds")
  cache <- readRDS(cache_file)
  cache$models <- purrr::map(cache$models, function(m) {
    m$method <- "fast_ci"
    m
  })
  saveRDS(cache, cache_file)

  run <- build_quietly(rd, export_dir)
  expect_equal(run$stats$built, 6)
  expect_true(all(purrr::map_chr(readRDS(cache_file)$models, "method") == "adaptive_pi"))

  # Same settings again: nothing to rebuild
  expect_equal(build_quietly(rd, export_dir)$stats$built, 0)

  # An older version of the same method is rebuilt too
  cache <- readRDS(cache_file)
  cache$models[[1]]$method_version <- 0L
  saveRDS(cache, cache_file)
  expect_equal(build_quietly(rd, export_dir)$stats$built, 1)
})


test_that("models of unknown provenance are rebuilt", {
  export_dir <- withr::local_tempdir()
  rd <- make_report_data()
  build_quietly(rd, export_dir)
  cache_file <- file.path(export_dir, "build_cache.rds")

  # Model files without cache entries, e.g. built by an older pipeline
  cache <- readRDS(cache_file)
  cache$models <- list()
  saveRDS(cache, cache_file)
  expect_equal(build_quietly(rd, export_dir)$stats$built, 6)

  # Entries recorded for unknown model files (no build time), as older
  # versions wrote them with the current method
  cache <- readRDS(cache_file)
  cache$models[1:2] <- purrr::map(cache$models[1:2], function(m) {
    m$built <- NA_character_
    m
  })
  saveRDS(cache, cache_file)
  expect_equal(build_quietly(rd, export_dir)$stats$built, 2)
  expect_equal(build_quietly(rd, export_dir)$stats$built, 0)
})
