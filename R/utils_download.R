#------------------------------------------------------------------------------------------
#   Project             : Replicating the Blood vs. Epigenetic PhenoAge Divergence
#   Repository          : BloodVsMethylation
#   Release Version     : 0.1.0.0
#   Author              : Iris Ivy Gauran
#   Description         : Shared Functions: Download Validation for Raw Files
#------------------------------------------------------------------------------------------

# Download helpers for R/01_pull_data.R, kept separate so they can be tested
# without triggering any downloads.

# A CDC error page comes back as a small HTML file with a 200 status, so a
# download is only trusted if it is big enough and does not start with "<".
looks_like_data <- function(path, min_bytes = 10000) {
  if (!file.exists(path) || file.size(path) < min_bytes) return(FALSE)
  first <- readBin(path, "raw", n = 1)
  !identical(rawToChar(first), "<")
}

download_raw <- function(file, url, dir = "data_raw") {
  dir.create(dir, showWarnings = FALSE, recursive = TRUE)
  dest <- file.path(dir, file)
  if (looks_like_data(dest)) {
    message("cached: ", file)
    return(invisible(dest))
  }
  message("downloading: ", file)
  status <- utils::download.file(url, dest, mode = "wb", quiet = TRUE)
  if (status != 0 || !looks_like_data(dest)) {
    unlink(dest)
    stop("download failed or returned a non-data page: ", url)
  }
  invisible(dest)
}
