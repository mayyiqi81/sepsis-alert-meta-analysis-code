# Rebuild canonical mortality models and all Figure S3A plotting summaries.

script_arg <- c(
  grep("^--file=", commandArgs(FALSE), value = TRUE),
  grep("\\.R$", commandArgs(FALSE), value = TRUE, ignore.case = TRUE)
)
script_dir <- if (length(script_arg) > 0) {
  dirname(normalizePath(sub("^--file=", "", script_arg[1]), mustWork = FALSE))
} else {
  getwd()
}

source(file.path(script_dir, "fit_canonical_main.R"), local = new.env(parent = globalenv()))
source(file.path(script_dir, "build_mortality_forest_cache.R"), local = new.env(parent = globalenv()))

message("Completed mortality-main computed-data pipeline.")
