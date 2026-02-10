# Tests for user model building and prediction functions

test_that("build_user_models validates input lengths", {
 expect_error(
   build_user_models(
     inchis = c("InChI=1S/A", "InChI=1S/B"),
     rts = c(1.0),
     system_type = "RP"
   ),
   "same length"
 )
})

test_that("build_user_models validates system_type", {
 expect_error(
   build_user_models(
     inchis = c("InChI=1S/A"),
     rts = c(1.0),
     system_type = "INVALID"
   ),
   "RP.*HILIC"
 )
})

test_that("build_user_models validates minimum compounds", {
 expect_error(
   build_user_models(
     inchis = rep("InChI=1S/A", 5),
     rts = rep(1.0, 5),
     system_type = "RP",
     min_compounds = 10
   ),
   "at least 10"
 )
})

test_that("build_user_models validates positive RTs", {
 expect_error(
   build_user_models(
     inchis = rep("InChI=1S/A", 10),
     rts = c(rep(1.0, 9), -1.0),
     system_type = "RP"
   ),
   "positive"
 )
})

test_that("build_user_models warns on non-InChI strings", {
 expect_warning(
   build_user_models(
     inchis = rep("not_an_inchi", 10),
     rts = rep(1.0, 10),
     system_type = "RP",
     verbose = FALSE
   ),
   "InChI"
 )
})

test_that("UserRTModels S4 class is valid", {
 obj <- methods::new(
   "UserRTModels",
   models = list(),
   diagnostics = data.frame(),
   input_data = data.frame(inchi = "InChI=1S/A", rt = 1.0),
   report_data = list(),
   metadata = list(
     system_type = "RP",
     n_input_compounds = 1,
     created_at = Sys.time()
   )
 )

 expect_s4_class(obj, "UserRTModels")
 expect_equal(obj@metadata$system_type, "RP")
})

test_that("UserRTModels validity check catches invalid system_type", {
 expect_error(
   methods::new(
     "UserRTModels",
     metadata = list(
       system_type = "INVALID",
       n_input_compounds = 1,
       created_at = Sys.time()
     )
   ),
   "RP.*HILIC"
 )
})

test_that("predict_rt_user validates input class", {
 expect_error(
   predict_rt_user(list()),
   "UserRTModels"
 )
})

test_that("predict_rt_user handles empty models", {
 empty_models <- methods::new(
   "UserRTModels",
   models = list(),
   diagnostics = data.frame(),
   input_data = data.frame(),
   report_data = list(),
   metadata = list(
     system_type = "RP",
     n_input_compounds = 0,
     created_at = Sys.time()
   )
 )

 expect_warning(result <- predict_rt_user(empty_models, verbose = FALSE))
 expect_s4_class(result, "UserRTPredictions")
 expect_equal(nrow(result@predictions), 0)
})

test_that("UserRTPredictions S4 class is valid", {
 obj <- methods::new(
   "UserRTPredictions",
   predictions = data.frame(
     compound_id = "InChI=1S/A",
     predicted_rt = 5.0,
     ci_lower = 4.5,
     ci_upper = 5.5,
     ci_width = 1.0
   ),
   source_models = "0001",
   metadata = list(
     n_predictions = 1,
     n_unique_compounds = 1,
     n_models_used = 1,
     created_at = Sys.time()
   )
 )

 expect_s4_class(obj, "UserRTPredictions")
 expect_equal(obj@metadata$n_predictions, 1)
})

test_that("as.data.frame works for UserRTPredictions", {
 obj <- methods::new(
   "UserRTPredictions",
   predictions = data.frame(
     compound_id = c("A", "B"),
     predicted_rt = c(1.0, 2.0)
   ),
   source_models = character(),
   metadata = list()
 )

 df <- as.data.frame(obj)
 expect_s3_class(df, "data.frame")
 expect_equal(nrow(df), 2)
})


# Tests for search_compound
test_that("search_compound validates InChI input", {
 expect_error(
   search_compound("not_an_inchi"),
   "InChI"
 )
})

test_that("search_compound accepts valid InChI", {
 # This will fail if no data, but should not error on validation
 skip_if_not(dir.exists(file.path(getwd(), "RepoRT_data")))

 # Should not error on a valid InChI format
 result <- search_compound("InChI=1S/C8H10N4O2")
 expect_s3_class(result, "data.frame")
})


# Tests for compare_systems
test_that("compare_systems validates input class", {
 expect_error(
   compare_systems(list()),
   "UserRTModels"
 )
})

test_that("compare_systems handles empty diagnostics", {
 empty_models <- methods::new(
   "UserRTModels",
   models = list(),
   diagnostics = data.frame(),
   input_data = data.frame(),
   report_data = list(),
   metadata = list(
     system_type = "RP",
     n_input_compounds = 0,
     created_at = Sys.time()
   )
 )

 expect_warning(result <- compare_systems(empty_models))
 expect_equal(nrow(result), 0)
})

test_that("compare_systems sorts by metric correctly", {
 models <- methods::new(
   "UserRTModels",
   models = list(),
   diagnostics = data.frame(
     source_system = c("0001", "0002", "0003"),
     r_squared = c(0.95, 0.99, 0.90),
     rmse = c(0.5, 0.2, 0.8),
     median_ci_width = c(1.0, 0.5, 1.5),
     n_calibration = c(10, 15, 8)
   ),
   input_data = data.frame(),
   report_data = list(),
   metadata = list(
     system_type = "RP",
     n_input_compounds = 0,
     created_at = Sys.time()
   )
 )

 # Sort by r_squared (descending)
 result <- compare_systems(models, metric = "r_squared")
 expect_equal(result$source_system[1], "0002")

 # Sort by rmse (ascending)
 result <- compare_systems(models, metric = "rmse")
 expect_equal(result$source_system[1], "0002")

 # Top N
 result <- compare_systems(models, top_n = 2)
 expect_equal(nrow(result), 2)
})


# Tests for get_system_info
test_that("get_system_info errors on missing system", {
 skip_if_not(dir.exists(file.path(getwd(), "RepoRT_data")))

 expect_error(
   get_system_info("NONEXISTENT"),
   "not found"
 )
})
