# Produce the study-level data consumed by the Figure 5 narrative analysis.
script_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
script_dir <- if (length(script_arg)) {
  dirname(normalizePath(sub("^--file=", "", script_arg[1])))
} else getwd()
source(file.path(script_dir, "..", "..", "analysis_scripts", "common", "project_paths.R"))

build_narrative_analysis <- function(project_root) {
  suppressPackageStartupMessages({library(dplyr); library(tidyr)})
  PROJECT_ROOT <- project_root
  INPUT_CSV <- source_csv_path("narrative_analysis_combined_conclusions.csv", root=PROJECT_ROOT)
  DATA_DIR <- computed_data_path("narrative_analysis", root=PROJECT_ROOT)
  MODULE_DIR <- analysis_script_path("narrative_analysis", root=PROJECT_ROOT)
  dir.create(DATA_DIR, recursive=TRUE, showWarnings=FALSE)
  save_csv <- function(x, name) write.csv(x, file.path(DATA_DIR, paste0(name,".csv")),
    row.names=FALSE, fileEncoding="UTF-8")
  for (module in c("definitions", "data")) {
    sys.source(file.path(MODULE_DIR,paste0(module,".R")), envir=environment())
  }
  fields <- c("dat", "alert", "population", "process", "mortality", "links",
              "outer", "schemes", "subgroup_orders")
  cache <- mget(fields, envir=environment(), inherits=FALSE)
  inputs <- c(INPUT_CSV, source_csv_path("narrative_analysis_antibiotic_time_audit.csv", root=PROJECT_ROOT))
  cache$input_md5 <- setNames(unname(tools::md5sum(inputs)), basename(inputs))
  saveRDS(cache, file.path(DATA_DIR,"figure5_data.rds"))
  message("Figure 5 computed data saved: ", DATA_DIR)
}
build_narrative_analysis(set_project_root_from_script(script_dir))
