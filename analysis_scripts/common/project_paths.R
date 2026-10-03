# Shared project path helpers for the reorganized production analysis tree.

is_absolute_path <- function(path) {
  grepl("^([A-Za-z]:[\\\\/]|/)", path)
}

find_project_root <- function(start = getwd()) {
  cur <- normalizePath(start, winslash = "/", mustWork = TRUE)
  if (file.exists(cur) && !dir.exists(cur)) {
    cur <- dirname(cur)
  }
  repeat {
    has_layout <- dir.exists(file.path(cur, "analysis_scripts", "common")) &&
      file.exists(file.path(cur, "analysis_scripts", "pipeline_registry.R"))
    if (has_layout) return(normalizePath(cur, winslash = "/", mustWork = TRUE))
    parent <- dirname(cur)
    if (identical(parent, cur)) {
      stop("Could not locate project root from: ", start, call. = FALSE)
    }
    cur <- parent
  }
}

project_path <- function(..., root = getOption("sepsis_project_root", find_project_root())) {
  file.path(root, ...)
}

source_csv_path <- function(..., root = getOption("sepsis_project_root", find_project_root())) {
  file.path(root, "source_csv", ...)
}

computed_data_path <- function(..., root = getOption("sepsis_project_root", find_project_root())) {
  file.path(root, "computed_data", ...)
}

plot_common_path <- function(..., root = getOption("sepsis_project_root", find_project_root())) {
  file.path(root, "analysis_scripts", "common", ...)
}

analysis_script_path <- function(..., root = getOption("sepsis_project_root", find_project_root())) {
  file.path(root, "analysis_scripts", ...)
}

figure_output_path <- function(..., root = getOption("sepsis_project_root", find_project_root())) {
  file.path(root, ...)
}

resolve_env_path <- function(env_var, default_path, root = getOption("sepsis_project_root", find_project_root())) {
  value <- Sys.getenv(env_var, unset = default_path)
  if (is_absolute_path(value)) value else file.path(root, value)
}

set_project_root_from_script <- function(script_dir = getwd()) {
  root <- find_project_root(script_dir)
  options(sepsis_project_root = root)
  root
}

