#!/usr/bin/env Rscript
# Optimized pipeline: furrr parallelization with caching
# Much faster than nested parallelization approach

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

# Load data (cached locally)
message("[", ts(), "] Loading RepoRT data...")
report_cache_dir <- file.path(getwd(), "RepoRT_data")
if (!dir.exists(report_cache_dir)) {
  dir.create(report_cache_dir, recursive = TRUE)
}
report_path <- download_report(dest_dir = report_cache_dir, overwrite = FALSE)
report_data <- load_report_data(report_path, method_types = c("RP"))
message("  ✓ Loaded ", length(report_data$datasets), " datasets")
message("[", ts(), "] Data loading complete\n")

# Setup output directory
output_dir <- file.path(getwd(), "website", "data", "models")
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

message("[", ts(), "] Building all models with furrr parallelization...")
message("  Workers: ", parallel::detectCores())
message("  Batch size: 50 (models per worker)")
message("")

start_time <- Sys.time()

# Use optimized build_all_models_furrr with tuned batch_size
# batch_size=50 has been benchmarked as optimal for this workload
results <- build_all_models_furrr(
  report_data = report_data,
  min_compounds = 10,
  method = "fast_ci",
  alpha = 0.05,
  n_workers = parallel::detectCores(),
  export_dir = output_dir,
  save_json = TRUE,
  batch_size = 50,  # Optimal balance: minimize overhead while allowing load balancing
  verbose = TRUE
)

elapsed <- difftime(Sys.time(), start_time, units = "secs")
end_time <- Sys.time()

# Results
models <- results$models
index <- results$index
stats <- results$stats

message("\n")
message("═══════════════════════════════════════════════════════════")
message("[", ts(), "] PIPELINE COMPLETE!")
message("═══════════════════════════════════════════════════════════\n")

message("Summary:")
message("  Total model pairs to build: ", stats$total_pairs)
message("  Successful models: ", stats$successful)
message("  Success rate: ", stats$success_rate, "%")
message("  Total runtime: ", round(as.numeric(elapsed)/60, 1), " minutes")
message("  Speed: ", round(stats$successful / as.numeric(elapsed), 2), " models/sec")
message("  Start time: ", format(start_time, "%Y-%m-%d %H:%M:%S"))
message("  End time:   ", format(end_time, "%Y-%m-%d %H:%M:%S"), "\n")

if (stats$successful > 0) {
  message("Model Statistics:")
  message("  Median CI width: ", stats$median_ci_width)
  message("  Mean CI width: ", stats$mean_ci_width)
  message("  Median error: ", stats$median_error)

  # Save index
  write.csv(index, file.path(output_dir, "model_index.csv"), row.names = FALSE)
  message("\n✓ Model index saved to: ", file.path(output_dir, "model_index.csv"), "\n")
}

message("═══════════════════════════════════════════════════════════")
message("Open website/model_viewer.html to explore models!")
message("═══════════════════════════════════════════════════════════\n")
