# Shared project setup for the yacht resistance analysis pipeline.
# Source this first, or run run_yacht_pipeline.R from the project root.

find_project_root <- function(start = getwd()) {
  path <- normalizePath(start, mustWork = TRUE)
  repeat {
    if (file.exists(file.path(path, "yacht_resistance.Rmd")) &&
        dir.exists(file.path(path, "data"))) {
      return(path)
    }
    parent <- dirname(path)
    if (identical(parent, path)) {
      stop("Could not find project root containing yacht_resistance.Rmd and data/.", call. = FALSE)
    }
    path <- parent
  }
}

PROJECT_ROOT <- find_project_root()
setwd(PROJECT_ROOT)
dir.create("figures", showWarnings = FALSE)
options(stringsAsFactors = FALSE)

require_package <- function(package) {
  if (!requireNamespace(package, quietly = TRUE)) {
    stop(sprintf("Package '%s' is required. Install it with install.packages('%s').", package, package),
         call. = FALSE)
  }
}
