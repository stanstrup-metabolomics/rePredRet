#' @title RepoRT Data Loading Functions
#' @description Functions for loading and processing data from the RepoRT
#'   repository.
#' @name data_loading
NULL

#' Download RepoRT Repository
#'
#' Downloads the RepoRT repository from GitHub.
#'
#' @param dest_dir Directory to save the downloaded data.
#' @param release Optional release tag. If NULL (default), downloads the
#'   master branch.
#' @param overwrite If TRUE, overwrite existing files.
#'
#' @return Path to the extracted RepoRT directory.
#'
#' @export
#' @examples
#' \dontrun{
#' report_path <- download_report(dest_dir = tempdir())
#' }
download_report <- function(dest_dir = tempdir(),
                            release = NULL,
                            overwrite = FALSE) {

  if (!requireNamespace("curl", quietly = TRUE)) {
    stop("Package 'curl' is required. Please install it.", call. = FALSE)
  }

  # Construct download URL
  if (is.null(release)) {
    url <- "https://github.com/michaelwitting/RepoRT/archive/refs/heads/master.zip"
    extract_name <- "RepoRT-master"
  } else {
    url <- paste0(
      "https://github.com/michaelwitting/RepoRT/archive/refs/tags/",
      release, ".zip"
    )
    extract_name <- paste0("RepoRT-", gsub("^v", "", release))
  }

  # Download
  zip_path <- file.path(dest_dir, "RepoRT.zip")

  if (!file.exists(zip_path) || overwrite) {
    message("Downloading RepoRT from GitHub...")
    curl::curl_download(url, zip_path)
  }

  # Extract
  extract_path <- file.path(dest_dir, extract_name)
  if (!dir.exists(extract_path) || overwrite) {
    message("Extracting RepoRT...")
    utils::unzip(zip_path, exdir = dest_dir)
  }

  return(extract_path)
}


#' Load RepoRT Data
#'
#' Loads all datasets from a RepoRT repository into a structured list.
#'
#' @param report_path Path to the extracted RepoRT repository.
#' @param method_types Character vector of method types to include.
#'   Default is c("RP", "HILIC"). Use NULL to include all.
#'
#' @return A list with two elements:
#'   \describe{
#'     \item{studies}{Tibble with all study metadata from studies.tsv}
#'     \item{datasets}{Named list of datasets, each containing:
#'       id, info, metadata, gradient, rtdata}
#'   }
#'
#' @importFrom readr read_tsv
#' @importFrom data.table fread
#' @importFrom tibble as_tibble
#' @export
#' @examples
#' \dontrun{
#' report_data <- load_report_data("path/to/RepoRT-master")
#' }
load_report_data <- function(report_path,
                             method_types = c("RP", "HILIC")) {

  processed_path <- file.path(report_path, "processed_data")

  if (!dir.exists(processed_path)) {
    stop("processed_data directory not found in ", report_path, call. = FALSE)
  }

  # Load studies index
  studies_file <- file.path(processed_path, "studies.tsv")
  if (!file.exists(studies_file)) {
    stop("studies.tsv not found in ", processed_path, call. = FALSE)
  }

  studies <- readr::read_tsv(
    studies_file,
    col_types = readr::cols(.default = readr::col_character()),
    show_col_types = FALSE
  )

  # Filter by method type if specified
  if (!is.null(method_types)) {
    studies <- studies[studies$method.type %in% method_types, ]
  }

  # Load each dataset
  datasets <- list()

  for (id in studies$id) {
    dir_path <- file.path(processed_path, id)

    if (!dir.exists(dir_path)) {
      message("Directory not found for dataset ", id, ", skipping...")
      next
    }

    # Load dataset components
    dataset <- list(id = id)

    # Info file
    info_file <- file.path(dir_path, paste0(id, "_info.tsv"))
    if (file.exists(info_file)) {
      dataset$info <- data.table::fread(info_file, sep = "\t") |>
        tibble::as_tibble()
    }

    # Metadata file
    meta_file <- file.path(dir_path, paste0(id, "_metadata.tsv"))
    if (file.exists(meta_file)) {
      dataset$metadata <- data.table::fread(meta_file, sep = "\t") |>
        tibble::as_tibble()
    }

    # Gradient file
    grad_file <- file.path(dir_path, paste0(id, "_gradient.tsv"))
    if (file.exists(grad_file)) {
      dataset$gradient <- data.table::fread(grad_file, sep = "\t") |>
        tibble::as_tibble()
    }

    # RT data (canonical - without stereochemistry)
    rt_file <- file.path(dir_path, paste0(id, "_rtdata_canonical_success.tsv"))
    if (file.exists(rt_file)) {
      dataset$rtdata <- data.table::fread(rt_file, sep = "\t") |>
        tibble::as_tibble()
    } else {
      message("No RT data found for dataset ", id, ", skipping...")
      next
    }

    datasets[[id]] <- dataset
  }

  message("Loaded ", length(datasets), " datasets")

  list(
    studies = studies,
    datasets = datasets
  )
}


#' Get Common Compounds Between Two Systems
#'
#' Finds compounds that exist in both chromatographic systems and returns
#' their retention times as a matrix suitable for model building.
#'
#' @param dataset1 First dataset (from load_report_data).
#' @param dataset2 Second dataset (from load_report_data).
#' @param by Column name to match compounds. Default "inchi.std".
#'
#' @return A 2-column matrix with retention times [rt1, rt2] for compounds
#'   present in both systems. If multiple measurements exist for the same
#'   compound, the median RT is used.
#'
#' @importFrom dplyr inner_join select all_of group_by summarize across
#' @export
get_common_compounds <- function(dataset1, dataset2, by = "inchi.std") {

  if (is.null(dataset1$rtdata) || is.null(dataset2$rtdata)) {
    return(matrix(nrow = 0, ncol = 2))
  }

  # Get columns we need
  rt1 <- dataset1$rtdata[, c(by, "rt"), drop = FALSE]
  rt2 <- dataset2$rtdata[, c(by, "rt"), drop = FALSE]

  # Remove stereochemistry from InChI (keep only main layer)
  # Stereochemistry layers are: /b (double bond), /t (tetrahedral), /m (chirality), /s (stereo type)
  # We remove from the first stereo layer onwards
  if (by == "inchi.std") {
    rt1[[by]] <- gsub("/(t|b|m|s)[^/]*.*$", "", rt1[[by]])
    rt2[[by]] <- gsub("/(t|b|m|s)[^/]*.*$", "", rt2[[by]])
  }

  # Inner join on compound identifier
  merged <- dplyr::inner_join(
    rt1,
    rt2,
    by = by,
    suffix = c(".1", ".2"),
    relationship = "many-to-many"
  )

  if (nrow(merged) == 0) {
    return(matrix(nrow = 0, ncol = 2))
  }

  # Handle duplicates: take median RT for same compound
  merged <- merged |>
    dplyr::group_by(dplyr::across(dplyr::all_of(by))) |>
    dplyr::summarize(
      rt.1 = stats::median(rt.1, na.rm = TRUE),
      rt.2 = stats::median(rt.2, na.rm = TRUE),
      .groups = "drop"
    )

  # Return as matrix
  rt_matrix <- as.matrix(merged[, c("rt.1", "rt.2")])
  colnames(rt_matrix) <- c("rt_source", "rt_target")

  # Store compound identifiers as attribute for traceability
  attr(rt_matrix, "compounds") <- merged[[by]]

  return(rt_matrix)
}


#' Get Latest RepoRT Release Tag
#'
#' Queries the GitHub API to get the latest release tag for RepoRT.
#'
#' @return Character string with the release tag, or NULL if no releases.
#'
#' @export
get_latest_report_release <- function() {
  url <- "https://api.github.com/repos/michaelwitting/RepoRT/releases/latest"

  tryCatch({
    response <- jsonlite::fromJSON(url)
    return(response$tag_name)
  }, error = function(e) {
    message("Could not fetch latest release: ", e$message)
    return(NULL)
  })
}
