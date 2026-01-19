#' @title Model Building Functions
#' @description Functions for building RT prediction models between
#'   chromatographic systems.
#' @name model_building
NULL

#' Build a Single Model Between Two Systems
#'
#' Builds a monotonically constrained GAM model to predict retention times
#' from one chromatographic system to another.
#'
#' @param rt_matrix A 2-column matrix [rt_source, rt_target] of common
#'   compound retention times.
#' @param sys1_id Source system identifier.
#' @param sys2_id Target system identifier.
#' @param n_boot Number of bootstrap iterations (default 1000). Only used if
#'   method = "bootstrap".
#' @param alpha Significance level for PI (default 0.05 for 95% PI).
#' @param n_cores Number of CPU cores for parallel computation.
#'   Default uses all available cores.
#' @param method Method for confidence intervals: "fast_ci" (72x faster, default)
#'   or "bootstrap" (accurate, slower).
#' @param n_boot Number of bootstrap iterations (default 200). Only used if
#'   method = "bootstrap".
#' @param save_plot Logical, whether to save an interactive plot (default FALSE).
#' @param plot_dir Directory to save plot HTML files (required if save_plot = TRUE).
#' @param dataset1 Optional dataset object for compound names in plots.
#' @param dataset2 Optional dataset object for compound names in plots.
#'
#' @return A list with model information:
#'   \describe{
#'     \item{status}{"success" or "insufficient_data"}
#'     \item{sys1_id}{Source system ID}
#'     \item{sys2_id}{Target system ID}
#'     \item{n_points}{Number of calibration compounds}
#'     \item{newdata}{Vector of x-values for prediction grid (1000 points)}
#'     \item{ci}{Matrix [pred, lower, upper] at newdata points}
#'     \item{calibration}{Original RT matrix used for training}
#'     \item{compounds}{InChI strings of calibration compounds}
#'     \item{build_time}{Timestamp when model was built}
#'     \item{stats}{Error statistics (mean, median, q95, max error and CI width)}
#'     \item{method}{Method used for PI calculation}
#'     \item{plot_path}{Path to saved plot HTML (if save_plot = TRUE)}
#'   }
#'
#' @importFrom boot boot
#' @importFrom parallel detectCores
#' @export
build_model <- function(rt_matrix,
                        sys1_id,
                        sys2_id,
                        n_boot = 200,
                        alpha = 0.05,
                        n_cores = parallel::detectCores(),
                        method = c("fast_ci", "bootstrap"),
                        save_plot = FALSE,
                        plot_dir = NULL,
                        dataset1 = NULL,
                        dataset2 = NULL) {

  method <- match.arg(method)

  # Check minimum data requirement
  if (nrow(rt_matrix) < 10) {
    return(list(
      status = "insufficient_data",
      sys1_id = sys1_id,
      sys2_id = sys2_id,
      n_points = nrow(rt_matrix),
      message = "Less than 10 common compounds"
    ))
  }

  # Create prediction grid spanning the RT range
  newdata <- seq(
    from = min(rt_matrix[, 1]),
    to = max(rt_matrix[, 1]),
    length.out = 1000
  )

  message("Building model ", sys1_id, " -> ", sys2_id,
          " (", nrow(rt_matrix), " compounds, method=", method, ")...")

  # Calculate confidence intervals using selected method
  if (method == "fast_ci") {
    # Fast: Single fit + predict SE + cummax (72x faster)
    ci <- gam_fast_ci(
      rt_matrix,
      newdata,
      alpha = alpha
    )
    boot_object <- NULL  # No bootstrap object for fast method

  } else {
    # Bootstrap: Fit model many times (more accurate)
    loess_boot <- boot::boot(
      rt_matrix,
      gam_mono_con_fun,
      R = n_boot,
      newdata = newdata,
      parallel = "multicore",
      ncpus = n_cores
    )

    # Calculate prediction intervals (includes both parameter + observation variance)
    ci <- boot2ci_PI(loess_boot, newdata = newdata, alpha = alpha)
    boot_object <- prune_boot_object(loess_boot)
  }

  # Calculate error statistics at original calibration points
  predicted_at_orig_x <- sapply(
    rt_matrix[, 1],
    function(x) which.min(abs(x - newdata))
  )
  errors <- abs(rt_matrix[, 2] - ci[predicted_at_orig_x, 1])
  ci_widths <- ci[, 3] - ci[, 2]

  # Store compound identifiers if available
  compounds <- attr(rt_matrix, "compounds")

  # Build result list
  result <- list(
    status = "success",
    sys1_id = sys1_id,
    sys2_id = sys2_id,
    n_points = nrow(rt_matrix),
    newdata = newdata,
    ci = ci,
    calibration = rt_matrix,
    compounds = compounds,
    build_time = Sys.time(),
    method = method,
    stats = list(
      mean_error = mean(errors, na.rm = TRUE),
      median_error = median(errors, na.rm = TRUE),
      q95_error = as.numeric(stats::quantile(errors, 0.95, na.rm = TRUE)),
      max_error = max(errors, na.rm = TRUE),
      mean_ci_width = mean(ci_widths, na.rm = TRUE),
      median_ci_width = median(ci_widths, na.rm = TRUE),
      q95_ci_width = as.numeric(stats::quantile(ci_widths, 0.95, na.rm = TRUE)),
      max_ci_width = max(ci_widths, na.rm = TRUE)
    ),
    boot_object = boot_object
  )

  # Save interactive plot if requested
  if (save_plot) {
    if (is.null(plot_dir)) {
      warning("save_plot = TRUE but plot_dir is NULL. Skipping plot.")
    } else {
      tryCatch({
        plot_path <- save_model_plot(
          result,
          output_dir = plot_dir,
          dataset1 = dataset1,
          dataset2 = dataset2
        )
        result$plot_path <- plot_path
      }, error = function(e) {
        warning("Failed to save plot: ", e$message)
      })
    }
  }

  return(result)
}


#' Prune Bootstrap Object for Storage
#'
#' Removes large components from a boot object to reduce storage size
#' while preserving essential information.
#'
#' @param boot_object A boot object from `boot::boot()`.
#'
#' @return Pruned boot object with reduced size.
#'
#' @details
#' The function removes:
#' - `statistic`: References the calling environment (huge)
#' - `t`: All bootstrap iteration results (we only need final CI)
#'
#' This typically reduces object size by ~95%.
#'
#' @keywords internal
prune_boot_object <- function(boot_object) {
  # Remove function reference (stores entire calling environment)
  # Remove all bootstrap iteration results (we only need final CI)
  # Use list subsetting to properly remove elements
  keep_elements <- setdiff(names(boot_object), c("statistic", "t"))
  pruned <- boot_object[keep_elements]
  class(pruned) <- class(boot_object)
  return(pruned)
}


#' Build All Pairwise Models
#'
#' Builds models between all pairs of chromatographic systems that share
#' sufficient common compounds.
#'
#' @param report_data Data from `load_report_data()`.
#' @param min_compounds Minimum common compounds required (default 10).
#' @param n_boot Number of bootstrap iterations (default 1000).
#' @param alpha Significance level for CI (default 0.01).
#' @param n_cores Number of CPU cores for parallel processing.
#' @param method_match If TRUE (default), only build models between
#'   same method types (RP-RP, HILIC-HILIC).
#' @param export_dir Optional directory to incrementally save models as
#'   they're built. If NULL (default), models are only stored in memory.
#' @param save_plots Logical, whether to save interactive HTML plots for each
#'   model (default TRUE).
#'
#' @return A list with:
#'   \describe{
#'     \item{models}{Named list of model objects}
#'     \item{index}{Tibble with model index and statistics}
#'   }
#'
#' @importFrom dplyr bind_rows tibble
#' @importFrom readr write_csv
#' @export
build_all_models <- function(report_data,
                             min_compounds = 10,
                             n_boot = 200,
                             alpha = 0.05,
                             n_cores = parallel::detectCores(),
                             method_match = TRUE,
                             export_dir = NULL,
                             method = c("fast_ci", "bootstrap"),
                             save_plots = TRUE) {

  method <- match.arg(method)

  studies <- report_data$studies
  datasets <- report_data$datasets
  dataset_ids <- names(datasets)

  models <- list()
  index_rows <- list()

  # Load build cache for incremental updates
  cache_file <- if (!is.null(export_dir)) {
    file.path(export_dir, "build_cache.rds")
  } else {
    NULL
  }

  # Load existing model index for incremental updates
  index_file <- if (!is.null(export_dir)) {
    file.path(export_dir, "model_index.csv")
  } else {
    NULL
  }

  existing_index <- if (!is.null(index_file) && file.exists(index_file)) {
    message("Loading existing model index...")
    idx <- readr::read_csv(index_file, show_col_types = FALSE)
    message("  Found ", nrow(idx), " existing models in index")
    idx
  } else {
    NULL
  }

  if (!is.null(cache_file)) {
    cache <- load_build_cache(cache_file)

    # Analyze dataset changes
    message("\nAnalyzing dataset changes...")
    changes <- analyze_dataset_changes(datasets, cache)

    message("  NEW: ", length(changes$new_ids), " datasets")
    message("  CHANGED: ", length(changes$changed_ids), " datasets")
    message("  UNCHANGED: ", length(changes$unchanged_ids), " datasets")
    message("  REMOVED: ", length(changes$removed_ids), " datasets")

    # Purge models for removed datasets
    if (length(changes$removed_ids) > 0 && !is.null(export_dir)) {
      cache <- purge_removed_models(changes$removed_ids, cache, export_dir)
    }

    # Update cache with new/changed dataset hashes
    cache <- update_dataset_cache(datasets, cache, changes$new_ids, changes$changed_ids)
  } else {
    # No caching if export_dir not specified
    cache <- NULL
    changes <- NULL
  }

  # Group by method type if matching required
  if (method_match) {
    method_groups <- split(dataset_ids, studies$method.type[
      match(dataset_ids, studies$id)
    ])
  } else {
    method_groups <- list(all = dataset_ids)
  }

  # Intelligent parallelization based on method:
  # - fast_ci: single-threaded, so use all cores for models (8×1)
  # - bootstrap: needs cores for resampling, so use balanced (4×2)
  if (method == "fast_ci") {
    n_parallel <- n_cores               # All cores for parallel models
    n_cores_per_model <- 1              # 1 core per model (fast_ci doesn't parallelize)
  } else {
    n_parallel <- min(4, n_cores)       # 4 models in parallel
    n_cores_per_model <- max(1, floor(n_cores / n_parallel))  # 2 cores per model for bootstrap
  }

  message("Parallelization: ", n_parallel, " models × ", n_cores_per_model,
          " cores = ", n_parallel * n_cores_per_model, " total cores")
  message("(Optimized for ", method, " method)\n")

  # Build models within each method group
  for (method_type in names(method_groups)) {
    ids <- method_groups[[method_type]]
    n_ids <- length(ids)

    message("\n=== Building ", method_type, " models (",
            n_ids, " systems) ===\n")

    # Fast compound overlap detection using sparse matrix multiplication
    message("Fast overlap detection for ", n_ids, " systems using sparse matrices...")

    # Extract standardized InChI for each system
    system_compounds <- lapply(ids, function(id) {
      inchi <- datasets[[id]]$rtdata$inchi.std
      gsub("/(t|b|m|s)[^/]*.*$", "", inchi)  # Remove stereochemistry
    })
    names(system_compounds) <- ids

    # Build sparse incidence matrix
    all_compounds <- unique(unlist(system_compounds))
    compound_id <- match(unlist(system_compounds), all_compounds)
    system_id <- rep(seq_along(system_compounds), lengths(system_compounds))

    M <- Matrix::sparseMatrix(
      i = compound_id,
      j = system_id,
      x = 1L,
      dims = c(length(all_compounds), length(system_compounds))
    )

    # Compute overlap matrix (systems × systems)
    overlap_matrix <- as.matrix(Matrix::crossprod(M))

    message("  Found ", sum(overlap_matrix >= min_compounds & row(overlap_matrix) != col(overlap_matrix)),
            " system pairs with ≥", min_compounds, " compounds\n")

    # Collect model pairs that meet minimum compound threshold (metadata only)
    # RT data will be fetched on-demand during model building
    model_pairs <- list()
    for (i in seq_len(n_ids)) {
      for (j in seq_len(n_ids)) {
        if (i == j) next  # Skip self-comparison

        n_overlap <- overlap_matrix[i, j]
        if (n_overlap < min_compounds) next

        id1 <- ids[i]
        id2 <- ids[j]

        # Store only pair metadata, not RT data yet
        model_pairs[[length(model_pairs) + 1]] <- list(
          id1 = id1,
          id2 = id2,
          method_type = method_type
        )
      }
    }

    if (length(model_pairs) == 0) next

    # Filter model pairs based on dataset changes (skip unchanged models)
    if (!is.null(cache) && !is.null(changes)) {
      model_pairs <- filter_model_pairs(model_pairs, changes, cache, export_dir)
    }

    if (length(model_pairs) == 0) {
      message("  ✓ All models up-to-date, nothing to build\n")
      next
    }

    # Set up plot directory if needed
    plot_dir <- NULL
    if (save_plots && !is.null(export_dir)) {
      plot_dir <- file.path(export_dir, "plots")
      if (!dir.exists(plot_dir)) {
        dir.create(plot_dir, recursive = TRUE)
      }
    }

    message("Building ", length(model_pairs), " models...")
    if (!is.null(export_dir)) {
      message("(saving incrementally every ", n_parallel * 10, " models)")
      if (save_plots) {
        message(" with interactive plots")
      }
      message("\n")
    } else {
      message("\n")
    }

    # Process in batches for incremental saving
    batch_size <- n_parallel * 10  # Process 40 models at a time (10 batches of 4)
    n_batches <- ceiling(length(model_pairs) / batch_size)
    all_results <- list()

    for (batch_idx in 1:n_batches) {
      start_idx <- (batch_idx - 1) * batch_size + 1
      end_idx <- min(batch_idx * batch_size, length(model_pairs))
      batch_pairs <- model_pairs[start_idx:end_idx]

      if (batch_idx %% 5 == 1 || batch_idx == n_batches) {
        message("  Batch ", batch_idx, "/", n_batches,
                " (models ", start_idx, "-", end_idx, ")")
      }

      # Build models in this batch
      if (n_parallel > 1) {
        cl <- parallel::makeCluster(n_parallel)
        parallel::clusterEvalQ(cl, library(rePredRet))
        parallel::clusterExport(cl, c("n_boot", "alpha", "n_cores_per_model", "method",
                                      "save_plots", "plot_dir", "datasets"),
                               envir = environment())

        results <- parallel::parLapply(cl, batch_pairs, function(pair) {
          # Fetch RT data on-demand for this model
          rt_matrix <- get_common_compounds(datasets[[pair$id1]], datasets[[pair$id2]])

          model <- build_model(
            rt_matrix,
            sys1_id = pair$id1,
            sys2_id = pair$id2,
            n_boot = n_boot,
            alpha = alpha,
            n_cores = n_cores_per_model,
            method = method,
            save_plot = save_plots,
            plot_dir = plot_dir,
            dataset1 = datasets[[pair$id1]],
            dataset2 = datasets[[pair$id2]]
          )
          model$model_key <- paste0(pair$id1, "_to_", pair$id2)
          model$method_type <- pair$method_type
          model
        })

        parallel::stopCluster(cl)
      } else {
        results <- lapply(batch_pairs, function(pair) {
          # Fetch RT data on-demand for this model
          rt_matrix <- get_common_compounds(datasets[[pair$id1]], datasets[[pair$id2]])

          model <- build_model(
            rt_matrix,
            sys1_id = pair$id1,
            sys2_id = pair$id2,
            n_boot = n_boot,
            alpha = alpha,
            n_cores = n_cores,
            method = method,
            save_plot = save_plots,
            plot_dir = plot_dir,
            dataset1 = datasets[[pair$id1]],
            dataset2 = datasets[[pair$id2]]
          )
          model$model_key <- paste0(pair$id1, "_to_", pair$id2)
          model$method_type <- pair$method_type
          model
        })
      }

      # Save models from this batch immediately
      for (model in results) {
        if (model$status == "success") {
          # Store in memory
          models[[model$model_key]] <- model
          all_results[[length(all_results) + 1]] <- model

          # Add to index
          index_rows[[model$model_key]] <- dplyr::tibble(
            model_key = model$model_key,
            from_id = model$sys1_id,
            to_id = model$sys2_id,
            method_type = model$method_type,
            n_compounds = model$n_points,
            mean_error = model$stats$mean_error,
            median_error = model$stats$median_error,
            q95_error = model$stats$q95_error,
            mean_ci_width = model$stats$mean_ci_width,
            median_ci_width = model$stats$median_ci_width,
            build_time = model$build_time
          )

          # Incremental export if directory specified
          if (!is.null(export_dir)) {
            model_dir <- file.path(export_dir, model$model_key)
            if (!dir.exists(model_dir)) {
              dir.create(model_dir, recursive = TRUE, showWarnings = FALSE)
            }

            # Save model
            saveRDS(model, file.path(model_dir, "model.rds"))

            # Save calibration data as CSV
            calib_df <- data.frame(
              rt_source = model$calibration[, 1],
              rt_target = model$calibration[, 2]
            )
            if (!is.null(model$compounds)) {
              calib_df$compound_id <- model$compounds
            }
            readr::write_csv(calib_df, file.path(model_dir, "calibration_data.csv"))

            # Save JSON for web viewer (includes prediction grid data)
            model_json <- list(
              model_key = model$model_key,
              from_id = model$sys1_id,
              to_id = model$sys2_id,
              n_compounds = model$n_points,
              stats = list(
                mean_error = model$stats$mean_error,
                median_error = model$stats$median_error,
                q95_error = model$stats$q95_error,
                mean_ci_width = model$stats$mean_ci_width,
                median_ci_width = model$stats$median_ci_width
              ),
              calibration = list(
                rt_source = as.numeric(model$calibration[, 1]),
                rt_target = as.numeric(model$calibration[, 2])
              ),
              prediction = list(
                rt_source = as.numeric(model$newdata),
                rt_target_predicted = as.numeric(model$ci[, 1]),
                ci_lower = as.numeric(model$ci[, 2]),
                ci_upper = as.numeric(model$ci[, 3])
              )
            )
            json_str <- jsonlite::toJSON(model_json, digits = 4, auto_unbox = TRUE, pretty = FALSE)
            writeLines(json_str, file.path(model_dir, "model.json"))

            # Update cache for this model
            if (!is.null(cache)) {
              cache$models[[model$model_key]] <- list(
                built = as.character(Sys.time()),
                method = method,
                alpha = alpha,
                n_compounds = model$n_points
              )
            }
          }
        }
      }
    }

    # Results already saved above in batch loop
  }

  # Build updated index incrementally
  new_index_entries <- dplyr::bind_rows(index_rows)

  # Merge with existing index
  if (!is.null(existing_index) && nrow(existing_index) > 0) {
    # Remove entries that were rebuilt (will be replaced with new entries)
    if (nrow(new_index_entries) > 0) {
      existing_index <- existing_index[
        !existing_index$model_key %in% new_index_entries$model_key,
      ]
    }

    # Remove entries for purged models (from removed datasets)
    if (!is.null(changes) && length(changes$removed_ids) > 0) {
      # Find models involving removed datasets
      removed_models <- existing_index$from_id %in% changes$removed_ids |
                       existing_index$to_id %in% changes$removed_ids
      if (any(removed_models)) {
        message("  Removing ", sum(removed_models), " purged models from index")
        existing_index <- existing_index[!removed_models, ]
      }
    }

    # Combine old (not rebuilt) + new (rebuilt)
    index <- dplyr::bind_rows(existing_index, new_index_entries)
  } else {
    # No existing index, use only new entries
    index <- new_index_entries
  }

  message("\n=== Built ", length(models), " models ===")
  message("=== Total models in index: ", nrow(index), " ===\n")

  # Save updated cache
  if (!is.null(cache) && !is.null(cache_file)) {
    save_build_cache(cache, cache_file)
  }

  # Save updated model index
  if (!is.null(index_file) && nrow(index) > 0) {
    message("Saving updated model index...")
    readr::write_csv(index, index_file)
    message("  Saved ", nrow(index), " models to model_index.csv")
  }

  list(
    models = models,
    index = index
  )
}


#' Export Models to Files
#'
#' Saves models as individual RDS files with accompanying CSV files
#' for transparency.
#'
#' @param model_results Output from `build_all_models()`.
#' @param output_dir Directory to save models.
#'
#' @return Invisibly returns the output directory path.
#'
#' @importFrom readr write_csv
#' @export
export_models <- function(model_results, output_dir = "models") {

  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }

  # Write model index
  readr::write_csv(
    model_results$index,
    file.path(output_dir, "model_index.csv")
  )

  # Write each model
  for (model_key in names(model_results$models)) {
    model <- model_results$models[[model_key]]

    # Create model directory
    model_dir <- file.path(output_dir, model_key)
    if (!dir.exists(model_dir)) {
      dir.create(model_dir)
    }

    # Save pruned RDS
    model_for_save <- model
    # Already pruned, but remove boot_object entirely for even smaller file
    model_for_save$boot_object <- NULL
    saveRDS(model_for_save, file.path(model_dir, "model.rds"))

    # Save CI grid as CSV
    ci_grid <- data.frame(
      x = model$newdata,
      pred = model$ci[, 1],
      lower = model$ci[, 2],
      upper = model$ci[, 3]
    )
    readr::write_csv(ci_grid, file.path(model_dir, "ci_grid.csv"))

    # Save calibration data as CSV
    calibration <- data.frame(
      rt_source = model$calibration[, 1],
      rt_target = model$calibration[, 2]
    )
    if (!is.null(model$compounds)) {
      calibration$compound = model$compounds
    }
    readr::write_csv(calibration, file.path(model_dir, "calibration_data.csv"))
  }

  message("Exported ", length(model_results$models),
          " models to ", output_dir)

  invisible(output_dir)
}
