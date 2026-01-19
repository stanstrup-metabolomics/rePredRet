#' Build Cache Management
#'
#' Functions for managing incremental model building with change detection
#' @name build_cache
NULL


#' Compute Hash for Dataset RT Data
#'
#' Creates an MD5 hash of a dataset's RT data to detect changes
#'
#' @param dataset Dataset object with rtdata component
#' @return Character string with MD5 hash
compute_dataset_hash <- function(dataset) {
  if (is.null(dataset$rtdata)) {
    return(NA_character_)
  }

  # Hash the RT data (InChI and RT values)
  data_to_hash <- dataset$rtdata[, c("inchi.std", "rt")]
  digest::digest(data_to_hash, algo = "md5")
}


#' Load Build Cache
#'
#' Loads the build cache from disk, or creates empty cache if none exists
#'
#' @param cache_file Path to cache file
#' @return List with datasets, models, version, and build_params
load_build_cache <- function(cache_file) {
  if (file.exists(cache_file)) {
    cache <- readRDS(cache_file)
    message("Loaded build cache: ", length(cache$datasets), " datasets, ",
            length(cache$models), " models")
    return(cache)
  }

  # Create empty cache
  list(
    datasets = list(),
    models = list(),
    version = "0.1.0",
    build_params = list()
  )
}


#' Save Build Cache
#'
#' Saves the build cache to disk
#'
#' @param cache Cache object
#' @param cache_file Path to cache file
save_build_cache <- function(cache, cache_file) {
  dir.create(dirname(cache_file), recursive = TRUE, showWarnings = FALSE)
  saveRDS(cache, cache_file)
  message("Build cache saved: ", cache_file)
}


#' Analyze Dataset Changes
#'
#' Compares current datasets against cache to detect changes
#'
#' @param datasets Named list of current datasets
#' @param cache Build cache object
#' @return List with new_ids, changed_ids, unchanged_ids, removed_ids
analyze_dataset_changes <- function(datasets, cache) {
  current_ids <- names(datasets)
  cached_ids <- names(cache$datasets)

  new_ids <- setdiff(current_ids, cached_ids)
  removed_ids <- setdiff(cached_ids, current_ids)
  potential_unchanged <- intersect(current_ids, cached_ids)

  changed_ids <- character()
  unchanged_ids <- character()

  # Check hashes for datasets that exist in both
  for (id in potential_unchanged) {
    current_hash <- compute_dataset_hash(datasets[[id]])
    cached_hash <- cache$datasets[[id]]$hash

    if (is.na(current_hash) || is.na(cached_hash) || current_hash != cached_hash) {
      changed_ids <- c(changed_ids, id)
    } else {
      unchanged_ids <- c(unchanged_ids, id)
    }
  }

  list(
    new_ids = new_ids,
    changed_ids = changed_ids,
    unchanged_ids = unchanged_ids,
    removed_ids = removed_ids
  )
}


#' Purge Models for Removed Datasets
#'
#' Deletes model directories and cache entries for models referencing removed datasets
#'
#' @param removed_ids Vector of removed dataset IDs
#' @param cache Build cache object
#' @param export_dir Directory containing model subdirectories
#' @return Updated cache object
purge_removed_models <- function(removed_ids, cache, export_dir) {
  if (length(removed_ids) == 0) {
    return(cache)
  }

  message("\nPurging models for ", length(removed_ids), " removed datasets...")

  # Find all models referencing removed datasets
  models_to_purge <- character()

  for (model_key in names(cache$models)) {
    # Parse model key (e.g., "0001_to_0002")
    parts <- strsplit(model_key, "_to_")[[1]]
    if (length(parts) == 2) {
      from_id <- parts[1]
      to_id <- parts[2]

      if (from_id %in% removed_ids || to_id %in% removed_ids) {
        models_to_purge <- c(models_to_purge, model_key)
      }
    }
  }

  if (length(models_to_purge) > 0) {
    message("  Found ", length(models_to_purge), " models to purge")

    # Delete from disk
    for (model_key in models_to_purge) {
      model_dir <- file.path(export_dir, model_key)
      if (dir.exists(model_dir)) {
        unlink(model_dir, recursive = TRUE)
      }

      # Remove from cache
      cache$models[[model_key]] <- NULL
    }

    message("  ✓ Purged ", length(models_to_purge), " models from disk and cache")
  }

  # Remove datasets from cache
  for (id in removed_ids) {
    cache$datasets[[id]] <- NULL
  }

  cache
}


#' Update Dataset Cache
#'
#' Updates cache with current dataset hashes
#'
#' @param datasets Named list of datasets
#' @param cache Build cache object
#' @param new_ids Vector of new dataset IDs
#' @param changed_ids Vector of changed dataset IDs
#' @return Updated cache object
update_dataset_cache <- function(datasets, cache, new_ids, changed_ids) {
  updated_ids <- c(new_ids, changed_ids)

  if (length(updated_ids) > 0) {
    for (id in updated_ids) {
      cache$datasets[[id]] <- list(
        hash = compute_dataset_hash(datasets[[id]]),
        last_updated = as.character(Sys.Date())
      )
    }
  }

  cache
}


#' Check if Model is Complete
#'
#' Checks if all required files exist for a model
#'
#' @param model_key Model key (e.g., "0001_to_0002")
#' @param export_dir Base directory for models
#' @return TRUE if all files exist, FALSE otherwise
model_is_complete <- function(model_key, export_dir) {
  model_dir <- file.path(export_dir, model_key)

  required_files <- c("model.json", "model.rds", "calibration_data.csv")

  all(file.exists(file.path(model_dir, required_files)))
}


#' Filter Model Pairs Based on Changes
#'
#' Determines which model pairs need to be built based on dataset changes
#'
#' @param model_pairs List of all potential model pairs
#' @param changes Dataset change analysis
#' @param cache Build cache object
#' @param export_dir Directory for model outputs
#' @return Filtered list of pairs that need building
filter_model_pairs <- function(model_pairs, changes, cache, export_dir) {
  changed_ids <- changes$changed_ids

  pairs_to_build <- list()
  pairs_skipped <- 0

  for (pair in model_pairs) {
    model_key <- paste0(pair$id1, "_to_", pair$id2)

    # Must rebuild if either dataset CHANGED (data modified)
    data_changed <- pair$id1 %in% changed_ids || pair$id2 %in% changed_ids

    if (data_changed) {
      # Data changed, must rebuild
      pairs_to_build[[length(pairs_to_build) + 1]] <- pair
    } else {
      # Data unchanged or new - check if model files exist (enables resume)
      if (model_is_complete(model_key, export_dir)) {
        pairs_skipped <- pairs_skipped + 1
      } else {
        # Model incomplete or missing, need to build
        pairs_to_build[[length(pairs_to_build) + 1]] <- pair
      }
    }
  }

  if (pairs_skipped > 0) {
    message("  ✓ Skipping ", pairs_skipped, " models (files complete)")
  }

  pairs_to_build
}
