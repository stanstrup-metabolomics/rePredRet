#' @title User Convenience Functions
#' @description Easy-to-use functions for end users to access pre-built
#'   rePredRet models. Files are fetched on demand from the
#'   rePredRet-models GitHub repository and cached locally.
#' @name user_functions
#' @keywords internal
NULL

# GitHub repository with the published models
.models_repo <- "stanstrup-metabolomics/rePredRet-models"
.models_url_base <- paste0("https://raw.githubusercontent.com/", .models_repo, "/main")
.studies_url <- "https://raw.githubusercontent.com/michaelwitting/RepoRT/master/processed_data/studies.tsv"


#' Fetch a File into the Local Cache
#'
#' Downloads a file once and reuses the cached copy until it is older than
#' `max_age_days`.
#'
#' @param url URL to download.
#' @param cache_path Path of the file relative to the cache directory.
#' @param refresh If TRUE, download even if a recent cached copy exists.
#' @param max_age_days Age after which the cached copy is downloaded again.
#'
#' @return Local path to the cached file.
#' @keywords internal
.fetch_file <- function(url, cache_path, refresh = FALSE, max_age_days = 1) {
  local_file <- file.path(rappdirs::user_cache_dir("rePredRet"), cache_path)

  is_recent <- file.exists(local_file) &&
    difftime(Sys.time(), file.mtime(local_file), units = "days") < max_age_days

  if (refresh || !is_recent) {
    dir.create(dirname(local_file), recursive = TRUE, showWarnings = FALSE)
    tmp_file <- paste0(local_file, ".part")
    ok <- tryCatch(
      utils::download.file(url, tmp_file, mode = "wb", quiet = TRUE) == 0,
      error = function(e) FALSE,
      warning = function(w) FALSE
    )
    if (!ok) {
      unlink(tmp_file)
      stop("Could not download ", url, call. = FALSE)
    }
    file.rename(tmp_file, local_file)
  }

  local_file
}


#' List Available Chromatographic Systems
#'
#' Returns RepoRT study information for all chromatographic systems that
#' have pre-built models.
#'
#' @param refresh If TRUE, download fresh copies instead of using the cache.
#'
#' @return Tibble with system information from RepoRT (id, name, url,
#'   method.type, authors, ...), plus `n_compounds` (compounds with an RT in
#'   the system) and `n_models` (models from or to it) from the published
#'   `systems.csv`. Model sets published before `systems.csv` was added give
#'   NA for these two columns.
#'
#' @importFrom readr read_tsv
#' @export
#' @examples
#' \dontrun{
#' systems <- rePredRet_systems()
#' }
rePredRet_systems <- function(refresh = FALSE) {
  studies_file <- .fetch_file(.studies_url, "RepoRT/studies.tsv", refresh = refresh)
  studies <- readr::read_tsv(studies_file, col_types = readr::cols(id = "c", .default = "c"),
                             progress = FALSE)

  models <- rePredRet_models(refresh = refresh)
  system_ids <- unique(c(models$sys1_id, models$sys2_id))
  studies <- studies[studies$id %in% system_ids, ]

  systems_file <- tryCatch(
    .fetch_file(paste0(.models_url_base, "/systems.csv"), "models/systems.csv", refresh = refresh),
    error = function(e) NULL
  )
  counts <- if (!is.null(systems_file)) {
    readr::read_csv(systems_file, col_types = readr::cols(sys_id = "c"), progress = FALSE)
  } else {
    data.frame(sys_id = character(), n_compounds = integer(), n_models = integer())
  }
  studies$n_compounds <- counts$n_compounds[match(studies$id, counts$sys_id)]
  studies$n_models <- counts$n_models[match(studies$id, counts$sys_id)]

  studies
}


#' List Available Models
#'
#' Returns the model index with statistics for all pre-built models.
#'
#' @param refresh If TRUE, download a fresh copy instead of using the cache.
#'
#' @return Tibble with one row per model: sys1_id (source system), sys2_id
#'   (target system), n_compounds, median_ci_width and median_error.
#'
#' @importFrom readr read_csv
#' @export
#' @examples
#' \dontrun{
#' models <- rePredRet_models()
#' }
rePredRet_models <- function(refresh = FALSE) {
  index_file <- .fetch_file(paste0(.models_url_base, "/model_index.csv"),
                            "models/model_index.csv", refresh = refresh)
  readr::read_csv(index_file,
                  col_types = readr::cols(sys1_id = "c", sys2_id = "c"),
                  progress = FALSE)
}


#' Get Model Details
#'
#' Retrieves a pre-built model with its prediction grid and calibration data.
#'
#' @param from_id Source system identifier.
#' @param to_id Target system identifier.
#' @param refresh If TRUE, download a fresh copy instead of using the cache.
#'
#' @return A list with:
#'   \describe{
#'     \item{sys1_id, sys2_id}{Source and target system}
#'     \item{n_compounds}{Number of calibration compounds}
#'     \item{method}{Interval method}
#'     \item{predictions}{Data frame with x, pred, lower, upper}
#'     \item{calibration}{Data frame with x, y, compound}
#'     \item{stats}{List with median_ci_width, median_error, r_squared,
#'       coverage (fraction of calibration compounds inside the PI),
#'       rel_mean_ci_width (mean PI width / target RT range) and build_time.
#'       Models built before these fields were added lack them.}
#'   }
#'
#' @export
#' @examples
#' \dontrun{
#' model <- rePredRet_model("0001", "0002")
#' }
rePredRet_model <- function(from_id, to_id, refresh = FALSE) {
  model_key <- paste0(from_id, "_to_", to_id)

  json_file <- tryCatch(
    .fetch_file(paste0(.models_url_base, "/models/", model_key, "/model.json"),
                file.path("models", model_key, "model.json"),
                refresh = refresh, max_age_days = 7),
    error = function(e) {
      stop("Model not found: ", from_id, " -> ", to_id,
           "\nUse rePredRet_models() to see available models.", call. = FALSE)
    }
  )

  jsonlite::fromJSON(json_file)
}
