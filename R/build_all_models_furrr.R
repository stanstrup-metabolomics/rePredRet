#' Build All Models with furrr - Optimized with Caching and Scheduling
#'
#' Builds retention time prediction models using furrr for efficient outer-level
#' parallelization with support for incremental builds and scheduling tuning.
#'
#' @param report_data Report data from load_report_data()
#' @param min_compounds Minimum overlapping compounds for a pair (default 10)
#' @param method "fast_ci" (default, fast) or "bootstrap" (slow, accurate)
#' @param alpha Significance level (default 0.05 for 95% intervals)
#' @param n_workers Number of parallel workers (default = cores available to
#'   this process, respecting container/cgroup CPU limits)
#' @param save_json Kept for backwards compatibility. `export_model_fast()`
#'   always writes model.json.
#' @param export_dir Directory to save models (required for caching)
#' @param batch_size Number of models to send to each worker at once (NULL = auto, 1 = dynamic, large = upfront)
#' @param report_version RepoRT version identifier (e.g. commit SHA) to record
#'   in the build cache. NULL leaves the recorded version unchanged.
#' @param verbose Print progress messages (default TRUE)
#'
#' @return List with:
#'   - models: list of results for the models built in this run
#'   - index: data.frame with metadata for all current models
#'   - stats: summary statistics
#'   - changed: TRUE if any model was built, added or removed in this run
#'
#' @details
#' This function uses furrr::future_map for outer-level parallelization,
#' avoiding nested parallelization overhead. Each worker builds models sequentially.
#'
#' Caching (requires export_dir): models are only rebuilt when one of their
#' datasets changed or their files are missing. Models whose datasets were
#' removed, or that no longer have enough common compounds, are deleted.
#' The build cache and `model_index.csv` are updated in export_dir.
#'
#' @importFrom furrr future_map
#' @importFrom future plan multisession
#' @export
build_all_models_furrr <- function(report_data,
                                  min_compounds = 10,
                                  method = c("fast_ci", "bootstrap"),
                                  alpha = 0.05,
                                  n_workers = parallelly::availableCores(),
                                  save_json = TRUE,
                                  export_dir = NULL,
                                  batch_size = NULL,
                                  report_version = NULL,
                                  verbose = TRUE) {

  if (!requireNamespace("furrr", quietly = TRUE)) {
    stop("Package 'furrr' required. Install with: install.packages('furrr')")
  }

  method <- match.arg(method)

  datasets <- report_data$datasets
  dataset_ids <- names(datasets)

  if (verbose) {
    message("\n========================================================")
    message("  Building Models with furrr + Caching")
    message("  Method: ", method, " | Workers: ", n_workers)
    message("========================================================\n")
  }

  # Load build cache for incremental updates
  cache_file <- if (!is.null(export_dir)) {
    file.path(export_dir, "build_cache.rds")
  } else {
    NULL
  }
  index_file <- if (!is.null(export_dir)) {
    file.path(export_dir, "model_index.csv")
  } else {
    NULL
  }

  if (!is.null(cache_file)) {
    cache <- load_build_cache(cache_file)
    message("Analyzing dataset changes for incremental build...")
    changes <- analyze_dataset_changes(datasets, cache)

    message("  NEW: ", length(changes$new_ids), " datasets")
    message("  CHANGED: ", length(changes$changed_ids), " datasets")
    message("  UNCHANGED: ", length(changes$unchanged_ids), " datasets")
    message("  REMOVED: ", length(changes$removed_ids), " datasets\n")

    n_cached_models <- length(cache$models)
    if (length(changes$removed_ids) > 0) {
      cache <- purge_removed_models(changes$removed_ids, cache, export_dir)
    }
    n_purged <- n_cached_models - length(cache$models)

    cache <- update_dataset_cache(datasets, cache, changes$new_ids, changes$changed_ids)
  } else {
    cache <- NULL
    changes <- NULL
    n_purged <- 0
  }

  # Find all pairs with enough common compounds (sparse overlap counts)
  rt_lookup <- prepare_rt_lookup(datasets)
  candidate_pairs <- find_model_pairs(rt_lookup, min_compounds)
  candidate_keys <- vapply(candidate_pairs, function(p) paste0(p$id1, "_to_", p$id2), character(1))

  if (verbose) {
    message("Found ", length(candidate_pairs), " system pairs with >=",
            min_compounds, " common compounds")
  }

  # Skip models that are up to date
  if (!is.null(cache)) {
    pairs_to_build <- filter_model_pairs(candidate_pairs, changes, cache, export_dir)
  } else {
    pairs_to_build <- candidate_pairs
  }
  build_keys <- vapply(pairs_to_build, function(p) paste0(p$id1, "_to_", p$id2), character(1))

  # Remove models that no longer qualify (e.g. a changed dataset lost overlap)
  n_stale <- 0
  existing_index <- read_model_index(index_file)
  if (!is.null(cache)) {
    if (length(candidate_pairs) == 0 && length(cache$models) > 0) {
      stop("No model pairs found but the build cache lists ", length(cache$models),
           " models. Refusing to delete them; check the input data.")
    }
    known_keys <- union(names(cache$models), existing_index$model_key)
    stale_keys <- setdiff(known_keys, candidate_keys)
    if (length(stale_keys) > 0) {
      message("Removing ", length(stale_keys), " models that no longer qualify")
      unlink(file.path(export_dir, stale_keys), recursive = TRUE)
      cache$models[stale_keys] <- NULL
      n_stale <- length(stale_keys)
    }
  }

  if (verbose) {
    message("Building ", length(pairs_to_build), " models with ", n_workers, " workers...")
  }

  # Set batch size (how many models to send to each worker at once)
  if (is.null(batch_size)) {
    # Auto: send groups of n_models / (n_workers * 4) to balance responsiveness
    batch_size <- max(1, ceiling(length(pairs_to_build) / (n_workers * 4)))
  }

  all_results <- list()

  if (length(pairs_to_build) > 0) {
    if (verbose) {
      message("Batch size: ", batch_size, " models per worker")
      message("Starting build at ", format(Sys.time(), "%H:%M:%S"), "\n")
    }

    build_start_time <- Sys.time()

    # Set up parallelization plan
    if (n_workers > 1) {
      future::plan(future::multisession, workers = n_workers)
    } else {
      future::plan(future::sequential)
    }
    on.exit(future::plan(future::sequential), add = TRUE)

    # Process in batches and track progress after each batch
    n_batches <- ceiling(length(pairs_to_build) / batch_size)
    n_successful <- 0

    for (batch_idx in 1:n_batches) {
      start_idx <- (batch_idx - 1) * batch_size + 1
      end_idx <- min(batch_idx * batch_size, length(pairs_to_build))

      if (verbose) {
        message("  Batch ", batch_idx, "/", n_batches,
                " (models ", start_idx, "-", end_idx, ")...")
      }

      # Only the RT matrices are sent to the workers, not the full datasets
      batch_pairs <- lapply(pairs_to_build[start_idx:end_idx], function(p) {
        list(
          sys1_id = p$id1,
          sys2_id = p$id2,
          rt_matrix = common_rt_matrix(rt_lookup[[p$id1]], rt_lookup[[p$id2]])
        )
      })

      # Build this batch in parallel using furrr
      batch_results <- furrr::future_map(
        batch_pairs,
        build_and_export_pair,
        alpha = alpha,
        method = method,
        export_dir = export_dir,
        .progress = FALSE  # Don't show progress for individual batches, we'll show batch progress instead
      )

      all_results <- c(all_results, batch_results)

      # Log progress after batch completes
      if (verbose) {
        n_successful <- n_successful +
          sum(vapply(batch_results, function(r) r$status == "success", logical(1)))
        elapsed <- as.numeric(difftime(Sys.time(), build_start_time, units = "secs"))
        speed <- n_successful / max(1, elapsed)

        message("    Progress: ", length(all_results), " / ", length(pairs_to_build),
                " | Success: ", n_successful,
                " | Speed: ", round(speed, 2), " models/sec")

        if (speed > 0 && length(all_results) < length(pairs_to_build)) {
          remaining <- length(pairs_to_build) - length(all_results)
          message("    ETA: ", round(remaining / speed / 60, 1), " minutes\n")
        }
      }
    }

    if (verbose) {
      elapsed_build <- as.numeric(difftime(Sys.time(), build_start_time, units = "secs"))
      message("\n--- Build Complete ---")
      message("  Built: ", length(all_results), " | Successful: ", n_successful)
      message("  Total time: ", round(elapsed_build / 60, 1), " minutes\n")
    }
  } else if (verbose) {
    message("  \u2713 All models up-to-date, nothing to build\n")
  }

  succeeded <- Filter(function(r) r$status == "success", all_results)
  failed <- Filter(function(r) r$status != "success", all_results)
  failed_keys <- vapply(failed, function(r) paste0(r$sys1_id, "_to_", r$sys2_id), character(1))

  # Failed rebuilds must not leave stale files behind
  if (!is.null(export_dir) && length(failed_keys) > 0) {
    unlink(file.path(export_dir, failed_keys), recursive = TRUE)
  }

  new_index <- do.call(rbind, lapply(succeeded, function(r) {
    data.frame(
      sys1_id = r$sys1_id,
      sys2_id = r$sys2_id,
      n_compounds = r$n_compounds,
      median_ci_width = r$median_ci_width,
      median_error = r$median_error,
      stringsAsFactors = FALSE
    )
  }))

  # Merge with the existing index: keep rows for models that were skipped
  n_added <- 0
  if (!is.null(cache)) {
    built_at <- as.character(Sys.time())
    for (r in all_results) {
      cache$models[[paste0(r$sys1_id, "_to_", r$sys2_id)]] <- list(
        built = built_at,
        method = method,
        alpha = alpha,
        n_compounds = r$n_compounds %||% NA_integer_,
        status = if (r$status == "success") "success" else "failed"
      )
    }

    skipped_keys <- setdiff(candidate_keys, build_keys)
    skipped_ok <- skipped_keys[vapply(skipped_keys, function(k) {
      !identical(cache$models[[k]]$status, "failed")
    }, logical(1))]

    # Models present on disk but unknown to the cache were built elsewhere
    # (e.g. restored from a runner cache) and count as additions
    unknown_keys <- setdiff(skipped_ok, names(cache$models))
    for (k in unknown_keys) {
      cache$models[[k]] <- list(built = NA_character_, method = method,
                                alpha = alpha, status = "success")
    }
    n_added <- length(unknown_keys)

    kept_index <- existing_index[existing_index$model_key %in% skipped_ok, , drop = FALSE]
    missing_keys <- setdiff(skipped_ok, kept_index$model_key)
    kept_index$model_key <- NULL
    if (length(missing_keys) > 0) {
      message("Recovering ", length(missing_keys), " index rows from model.json files")
      recovered <- do.call(rbind, lapply(missing_keys, function(k) {
        index_row_from_json(file.path(export_dir, k, "model.json"))
      }))
      kept_index <- rbind(kept_index, recovered)
    }

    index_df <- rbind(kept_index, new_index)
  } else {
    index_df <- new_index
  }

  if (is.null(index_df)) {
    index_df <- data.frame(
      sys1_id = character(),
      sys2_id = character(),
      n_compounds = integer(),
      median_ci_width = numeric(),
      median_error = numeric(),
      stringsAsFactors = FALSE
    )
  }
  index_df <- index_df[order(index_df$sys1_id, index_df$sys2_id), , drop = FALSE]
  rownames(index_df) <- NULL

  changed <- length(all_results) > 0 || n_purged > 0 || n_stale > 0 || n_added > 0

  if (!is.null(cache)) {
    if (!is.null(report_version)) {
      cache$report_version <- report_version
    }
    cache$build_params <- list(method = method, alpha = alpha, min_compounds = min_compounds)
    save_build_cache(cache, cache_file)
    utils::write.csv(index_df, index_file, row.names = FALSE)
    message("Model index saved: ", index_file, " (", nrow(index_df), " models)")
  }

  # Summary statistics
  stats <- list(
    total_pairs = length(candidate_pairs),
    built = length(all_results),
    failed = length(failed),
    removed = n_purged + n_stale,
    successful = nrow(index_df),
    success_rate = if (length(candidate_pairs) > 0) {
      round(nrow(index_df) / length(candidate_pairs) * 100, 1)
    } else {
      0
    }
  )
  if (nrow(index_df) > 0) {
    stats$median_ci_width <- round(stats::median(index_df$median_ci_width), 3)
    stats$mean_ci_width <- round(mean(index_df$median_ci_width), 3)
    stats$median_error <- round(stats::median(index_df$median_error), 3)
    stats$mean_error <- round(mean(index_df$median_error), 3)
  }

  if (verbose) {
    message("\n========================================================")
    message("  Pipeline Complete!")
    message("========================================================\n")
    message("Summary:")
    message("  Qualifying pairs: ", stats$total_pairs)
    message("  Built this run: ", stats$built, " (", stats$failed, " failed)")
    message("  Removed this run: ", stats$removed)
    message("  Models in index: ", stats$successful)
    message("\n")
  }

  list(
    models = structure(all_results, class = "model_list"),
    index = index_df,
    stats = stats,
    changed = changed
  )
}


#' Build and Export a Single Model Pair
#'
#' Worker function for [build_all_models_furrr()]. Defined at top level so
#' that only its arguments are sent to parallel workers.
#'
#' @param pair List with sys1_id, sys2_id and rt_matrix
#' @param alpha Significance level
#' @param method CI method
#' @param export_dir Directory to save models, or NULL
#' @return List with status and summary statistics
#' @keywords internal
build_and_export_pair <- function(pair, alpha, method, export_dir) {
  model <- tryCatch(
    suppressMessages(build_model(
      rt_matrix = pair$rt_matrix,
      sys1_id = pair$sys1_id,
      sys2_id = pair$sys2_id,
      alpha = alpha,
      n_cores = 1,
      method = method,
      save_plot = FALSE
    )),
    error = function(e) list(status = "error", message = conditionMessage(e))
  )

  if (model$status != "success") {
    return(list(
      status = model$status,
      sys1_id = pair$sys1_id,
      sys2_id = pair$sys2_id,
      message = model$message %||% NA_character_
    ))
  }

  if (!is.null(export_dir)) {
    export_model_fast(model, file.path(export_dir, paste0(model$sys1_id, "_to_", model$sys2_id)))
  }

  list(
    status = "success",
    sys1_id = model$sys1_id,
    sys2_id = model$sys2_id,
    n_compounds = model$n_points,
    median_ci_width = model$stats$median_ci_width,
    median_error = model$stats$median_error
  )
}


#' Prepare Per-Dataset RT Lookup
#'
#' Strips stereochemistry from InChIs and takes the median RT per compound,
#' once per dataset. [common_rt_matrix()] then matches two lookups, giving the
#' same result as [get_common_compounds()] much faster.
#'
#' @param datasets Named list of datasets with rtdata
#' @return Named list of lists with sorted `key` and matching median `rt`
#' @keywords internal
prepare_rt_lookup <- function(datasets) {
  lapply(datasets, function(d) {
    if (is.null(d$rtdata)) {
      return(list(key = character(), rt = numeric()))
    }
    key <- gsub("/(t|b|m|s)[^/]*.*$", "", d$rtdata$inchi.std)
    # Same order as dplyr::group_by(): C locale, NA last
    ukey <- sort(unique(key), method = "radix", na.last = TRUE)
    groups <- split(d$rtdata$rt, factor(key, levels = ukey, exclude = NULL))
    list(
      key = ukey,
      rt = vapply(groups, stats::median, numeric(1), na.rm = TRUE, USE.NAMES = FALSE)
    )
  })
}


#' Common Compound RT Matrix from Lookups
#'
#' @param lookup1,lookup2 Elements of [prepare_rt_lookup()] output
#' @return Matrix like [get_common_compounds()]
#' @keywords internal
common_rt_matrix <- function(lookup1, lookup2) {
  idx <- match(lookup1$key, lookup2$key, nomatch = 0L)
  keep <- idx > 0L

  if (!any(keep)) {
    return(matrix(nrow = 0, ncol = 2))
  }

  rt_matrix <- cbind(rt_source = lookup1$rt[keep], rt_target = lookup2$rt[idx[keep]])
  attr(rt_matrix, "compounds") <- lookup1$key[keep]
  rt_matrix
}


#' Find Model Pairs by Compound Overlap
#'
#' Counts common compounds for all ordered dataset pairs with a sparse
#' incidence matrix and returns the pairs with at least `min_compounds`.
#'
#' @param rt_lookup Output of [prepare_rt_lookup()]
#' @param min_compounds Minimum common compounds
#' @return List of pairs with id1 and id2
#' @keywords internal
find_model_pairs <- function(rt_lookup, min_compounds) {
  ids <- names(rt_lookup)
  keys <- lapply(rt_lookup, `[[`, "key")
  all_keys <- unique(unlist(keys))

  if (length(all_keys) == 0) {
    return(list())
  }

  M <- Matrix::sparseMatrix(
    i = match(unlist(keys), all_keys),
    j = rep(seq_along(keys), lengths(keys)),
    x = 1L,
    dims = c(length(all_keys), length(keys))
  )
  overlap <- as.matrix(Matrix::crossprod(M))
  diag(overlap) <- 0

  hits <- which(overlap >= min_compounds, arr.ind = TRUE)
  hits <- hits[order(hits[, 1], hits[, 2]), , drop = FALSE]

  lapply(seq_len(nrow(hits)), function(r) {
    list(id1 = ids[hits[r, 1]], id2 = ids[hits[r, 2]])
  })
}


#' Read Model Index
#'
#' Reads model_index.csv in the current (sys1_id/sys2_id) or the older
#' (from_id/to_id) format.
#'
#' @param index_file Path to model_index.csv, or NULL
#' @return data.frame with model_key and the index columns
#' @keywords internal
read_model_index <- function(index_file) {
  cols <- c("sys1_id", "sys2_id", "n_compounds", "median_ci_width", "median_error")
  empty <- data.frame(
    model_key = character(), sys1_id = character(), sys2_id = character(),
    n_compounds = integer(), median_ci_width = numeric(), median_error = numeric(),
    stringsAsFactors = FALSE
  )

  if (is.null(index_file) || !file.exists(index_file)) {
    return(empty)
  }

  # Read system IDs as character to keep leading zeros
  header <- names(utils::read.csv(index_file, nrows = 0))
  col_classes <- ifelse(header %in% c("sys1_id", "sys2_id", "from_id", "to_id", "model_key"),
                        "character", NA)
  idx <- utils::read.csv(index_file, stringsAsFactors = FALSE, colClasses = col_classes)
  if (all(c("from_id", "to_id") %in% names(idx))) {
    names(idx)[names(idx) == "from_id"] <- "sys1_id"
    names(idx)[names(idx) == "to_id"] <- "sys2_id"
  }
  if (!all(cols %in% names(idx))) {
    warning("Ignoring model index with unexpected columns: ", index_file)
    return(empty)
  }

  idx <- idx[, cols]
  idx$model_key <- paste0(idx$sys1_id, "_to_", idx$sys2_id)
  idx[, c("model_key", cols)]
}


#' Index Row from model.json
#'
#' @param json_file Path to a model.json written by [export_model_json()]
#' @return One-row data.frame with the index columns
#' @keywords internal
index_row_from_json <- function(json_file) {
  j <- jsonlite::fromJSON(json_file)
  data.frame(
    sys1_id = j$sys1_id,
    sys2_id = j$sys2_id,
    n_compounds = j$n_compounds,
    median_ci_width = j$stats$median_ci_width,
    median_error = j$stats$median_error,
    stringsAsFactors = FALSE
  )
}
