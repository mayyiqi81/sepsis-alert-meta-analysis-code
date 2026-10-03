# Rebuild all publication-trend, cumulative-mortality, and TSA caches for Figure 4.

script_arg <- c(
  grep("^--file=", commandArgs(FALSE), value = TRUE),
  grep("\\.R$", commandArgs(FALSE), value = TRUE, ignore.case = TRUE)
)
script_dir <- if (length(script_arg) > 0) {
  dirname(normalizePath(sub("^--file=", "", script_arg[1]), mustWork = FALSE))
} else {
  getwd()
}

source(file.path(script_dir, "build_publication_trends_and_tsa.R"), local = new.env(parent = globalenv()))
source(file.path(script_dir, "fit_publication_cumulative_cache.R"), local = new.env(parent = globalenv()))

message("Completed publication/cumulative-mortality computed-data pipeline.")
