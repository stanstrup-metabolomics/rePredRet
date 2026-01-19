#' Build All Models with Proper furrr Parallelization
#'
#' Builds retention time prediction models using furrr for efficient outer-level
#' parallelization. Much faster than nested parallelization approach.
#'
#' @param report_data Report data from load_report_data()
#' @param min_compounds Minimum overlapping compounds for a pair (default 10)
#' @param method "fast_ci" (default, fast) or "bootstrap" (slow, accurate)
#' @param alpha Significance level (default 0.05 for 95% intervals)
#' @param n_workers Number of parallel workers (default = available cores)
#' @param save_json Save model data as JSON (default TRUE, much faster than HTML)
#' @param export_dir Directory to save models (required)
#' @param verbose Print progress messages (default TRUE)
#'
#' @return List with:
#'   - models: list of all built models
#'   - index: data.frame with model metadata
#'   - stats: summary statistics
#'
#' @details
#' This function uses furrr::future_map for outer-level parallelization,
#' avoiding nested parallelization overhead. Each worker builds models
#' sequentially.
#'
#' Expected speedup vs nested parallelization: 3-4x (less overhead)
#'
#' @importFrom furrr future_map
#' @importFrom future plan multisession
#' @export
build_all_models_furrr <- function(report_data,
                                  min_compounds = 10,
                                  method = c("fast_ci", "bootstrap"),
                                  alpha = 0.05,
                                  n_workers = parallel::detectCores(),
                                  save_json = TRUE,
                                  export_dir = NULL,
                                  verbose = TRUE) {

  if (!requireNamespace("furrr", quietly = TRUE)) {
    stop("Package 'furrr' required. Install with: install.packages('furrr')")
  }

  method <- match.arg(method)

  datasets <- report_data$datasets
  studies <- report_data$studies
  dataset_ids <- names(datasets)

  if (verbose) {
    message("\n╔════════════════════════════════════════════════════════╗")
    message("║  Building Models with furrr Parallelization           ║")
    message("║  Method: ", method, " | Workers: ", n_workers, "                   ║")
    message("╚════════════════════════════════════════════════════════╝\n")
  }

  # Set up parallelization plan
  if (verbose) message("Setting up ", n_workers, " parallel workers...\n")
  future::plan(future::multisession, workers = n_workers)

  # Collect all model pairs to build
  model_pairs <- list()

  for (i in seq_along(dataset_ids)) {
    for (j in seq_along(dataset_ids)) {
      if (i == j) next  # Skip self-comparison

      id1 <- dataset_ids[i]
      id2 <- dataset_ids[j]

      # Get common compounds
      rt_matrix <- get_common_compounds(datasets[[id1]], datasets[[id2]])

      if (nrow(rt_matrix) < min_compounds) {
        next
      }

      # Get compound identifiers if available
      inchi1 <- datasets[[id1]]$rtdata$inchi.std
      inchi2 <- datasets[[id2]]$rtdata$inchi.std
      inchi1_std <- gsub("/(t|b|m|s)[^/]*.*$", "", inchi1)
      inchi2_std <- gsub("/(t|b|m|s)[^/]*.*$", "", inchi2)

      common_inchi <- intersect(inchi1_std, inchi2_std)

      model_pairs[[length(model_pairs) + 1]] <- list(
        sys1_id = id1,
        sys2_id = id2,
        rt_matrix = rt_matrix,
        compounds = common_inchi
      )
    }
  }

  if (verbose) {
    message("Building ", length(model_pairs), " models with ", n_workers, " workers...\n")
  }

  # Build models in parallel using furrr
  results <- furrr::future_map(
    model_pairs,
    function(pair) {
      # Build single model
      model <- build_model(
        rt_matrix = pair$rt_matrix,
        sys1_id = pair$sys1_id,
        sys2_id = pair$sys2_id,
        alpha = alpha,
        n_cores = 1,  # No internal parallelization (already parallelized)
        method = method,
        save_plot = FALSE  # Don't generate HTML plots
      )

      # Export model data
      if (!is.null(export_dir) && model$status == "success") {
        model_dir <- file.path(
          export_dir,
          paste0(model$sys1_id, "_to_", model$sys2_id)
        )
        export_model_fast(model, model_dir)

        if (save_json) {
          export_model_json(model, file.path(model_dir, "model.json"))
        }

        # Return metadata instead of full model
        return(list(
          status = "success",
          sys1_id = model$sys1_id,
          sys2_id = model$sys2_id,
          n_compounds = model$n_points,
          median_pi_width = model$stats$median_pi_width,
          median_error = model$stats$median_error
        ))
      } else {
        return(list(
          status = model$status,
          sys1_id = pair$sys1_id,
          sys2_id = pair$sys2_id,
          message = model$message %||% NA_character_
        ))
      }
    },
    .progress = verbose
  )

  # Clean up
  future::plan(future::sequential)

  # Convert results to index data.frame
  successful <- Filter(function(r) r$status == "success", results)

  index_df <- do.call(rbind, lapply(successful, function(r) {
    data.frame(
      sys1_id = r$sys1_id,
      sys2_id = r$sys2_id,
      n_compounds = r$n_compounds,
      median_pi_width = r$median_pi_width,
      median_error = r$median_error,
      stringsAsFactors = FALSE
    )
  }))

  # Calculate summary statistics
  if (nrow(index_df) > 0) {
    stats <- list(
      total_pairs = length(model_pairs),
      successful = nrow(index_df),
      success_rate = round(nrow(index_df) / length(model_pairs) * 100, 1),
      median_pi_width = round(median(index_df$median_pi_width), 3),
      mean_pi_width = round(mean(index_df$median_pi_width), 3),
      median_error = round(median(index_df$median_error), 3),
      mean_error = round(mean(index_df$median_error), 3)
    )

    if (verbose) {
      message("\n╔════════════════════════════════════════════════════════╗")
      message("║  Pipeline Complete!                                    ║")
      message("╚════════════════════════════════════════════════════════╝\n")
      message("Summary:")
      message("  Total pairs considered: ", stats$total_pairs)
      message("  Successful models: ", stats$successful)
      message("  Success rate: ", stats$success_rate, "%")
      message("  Median PI width: ", stats$median_pi_width)
      message("  Mean PI width: ", stats$mean_pi_width)
      message("\n")
    }
  } else {
    stats <- list(
      total_pairs = length(model_pairs),
      successful = 0,
      success_rate = 0
    )
  }

  list(
    models = structure(results, class = "model_list"),
    index = index_df,
    stats = stats
  )
}
