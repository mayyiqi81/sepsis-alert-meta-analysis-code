# Rebuild every computed-data artifact owned by the Figure 3 evidence-map domain.

script_arg <- c(
  grep("^--file=", commandArgs(FALSE), value = TRUE),
  grep("\\.R$", commandArgs(FALSE), value = TRUE, ignore.case = TRUE)
)
script_dir <- if (length(script_arg) > 0) {
  dirname(normalizePath(sub("^--file=", "", script_arg[1]), mustWork = FALSE))
} else {
  getwd()
}

source(file.path(script_dir, "build_small_study_effect_diagnostics.R"), local = new.env(parent = globalenv()))
source(file.path(script_dir, "build_evidence_map_overview.R"), local = new.env(parent = globalenv()))

message("Completed evidence-map computed-data pipeline.")
