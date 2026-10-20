#!/usr/bin/env Rscript
# Optimized pipeline: furrr parallelization with incremental builds
#
# Environment variables:
#   MAX_SYSTEMS     Only use the first N datasets (0 = all). Test runs write to
#                   output/models_test so the real models are untouched.
#   MODELS_DIR      Where models, build_cache.rds and model_index.csv are read
#                   and written (default output/models). Ignored for test runs.
#   REPORT_VERSION  RepoRT version (commit SHA) recorded in the build cache
#                   and in report_version.txt
#   GITHUB_ENV      Set by GitHub Actions; models_changed=true/false is written to it

message("\n")
message("╔════════════════════════════════════════════════════════╗")
message("║  rePredRet: OPTIMIZED PIPELINE (furrr)                ║")
message("║  • Outer-level furrr parallelization (3-4x faster)    ║")
message("║  • Smart caching for incremental builds                ║")
message("║  • Batch size = 50 (optimized scheduling)              ║")
message("╚════════════════════════════════════════════════════════╝")

# Format timestamp helper
ts <- function() format(Sys.time(), "%Y-%m-%d %H:%M:%S")
message("\n[", ts(), "] Pipeline started\n")

# Attach the installed package
library(rePredRet)

max_systems <- as.integer(Sys.getenv("MAX_SYSTEMS", "0"))
report_version <- Sys.getenv("REPORT_VERSION", "")
if (report_version == "") report_version <- NULL

# Load data (cached locally)
message("[", ts(), "] Loading RepoRT data...")
report_cache_dir <- file.path(getwd(), "RepoRT_data")
if (!dir.exists(report_cache_dir)) {
  dir.create(report_cache_dir, recursive = TRUE)
}
report_path <- download_report(dest_dir = report_cache_dir, overwrite = FALSE)
report_data <- load_report_data(report_path, method_types = c("RP"))
message("  ✓ Loaded ", length(report_data$datasets), " datasets")

if (!is.na(max_systems) && max_systems > 0) {
  report_data$datasets <- utils::head(report_data$datasets, max_systems)
  message("  TEST RUN: limited to ", length(report_data$datasets), " datasets")
}
message("[", ts(), "] Data loading complete\n")

# Setup output directory
is_test_run <- !is.na(max_systems) && max_systems > 0
output_dir <- if (is_test_run) {
  file.path(getwd(), "output", "models_test")
} else {
  Sys.getenv("MODELS_DIR", file.path(getwd(), "output", "models"))
}
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

n_workers <- parallelly::availableCores()
message("[", ts(), "] Building models with furrr parallelization...")
message("  Workers: ", n_workers)
message("  Batch size: 50 (models per worker)")
message("")

start_time <- Sys.time()

# Use optimized build_all_models_furrr with tuned batch_size
# batch_size=50 has been benchmarked as optimal for this workload
# Also updates build_cache.rds and model_index.csv in output_dir
results <- build_all_models_furrr(
  report_data = report_data,
  min_compounds = 10,
  method = "adaptive_pi",
  alpha = 0.05,
  n_workers = n_workers,
  export_dir = output_dir,
  batch_size = 50,  # Optimal balance: minimize overhead while allowing load balancing
  report_version = report_version,
  verbose = TRUE
)

if (!is.null(report_version)) {
  writeLines(report_version, file.path(output_dir, "report_version.txt"))
}

elapsed <- difftime(Sys.time(), start_time, units = "secs")
end_time <- Sys.time()

stats <- results$stats

message("\n")
message("═══════════════════════════════════════════════════════════")
message("[", ts(), "] PIPELINE COMPLETE!")
message("═══════════════════════════════════════════════════════════\n")

message("Summary:")
message("  Qualifying model pairs: ", stats$total_pairs)
message("  Built this run: ", stats$built, " (", stats$failed, " failed)")
message("  Removed this run: ", stats$removed)
message("  Models in index: ", stats$successful)
message("  Models changed: ", results$changed)
message("  Total runtime: ", round(as.numeric(elapsed)/60, 1), " minutes")
message("  Start time: ", format(start_time, "%Y-%m-%d %H:%M:%S"))
message("  End time:   ", format(end_time, "%Y-%m-%d %H:%M:%S"), "\n")

if (stats$successful > 0) {
  message("Model Statistics:")
  message("  Median CI width: ", stats$median_ci_width)
  message("  Mean CI width: ", stats$mean_ci_width)
  message("  Median error: ", stats$median_error, "\n")
}

# Tell the workflow whether the data repository needs updating
github_env <- Sys.getenv("GITHUB_ENV", "")
if (github_env != "") {
  cat("models_changed=", tolower(results$changed), "\n", sep = "", file = github_env, append = TRUE)
}

message("═══════════════════════════════════════════════════════════")
message("Explore models: https://stanstrup-metabolomics.github.io/rePredRet/model_viewer.html")
message("═══════════════════════════════════════════════════════════\n")
