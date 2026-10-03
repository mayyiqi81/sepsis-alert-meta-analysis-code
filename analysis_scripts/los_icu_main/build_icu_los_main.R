script_arg <- c(
  grep("^--file=", commandArgs(FALSE), value = TRUE),
  grep("\\.R$", commandArgs(FALSE), value = TRUE, ignore.case = TRUE)
)
script_dir <- if (length(script_arg) > 0) {
  dirname(normalizePath(sub("^--file=", "", script_arg[1]), mustWork = FALSE))
} else {
  getwd()
}
project_root <- normalizePath(file.path(script_dir, "..", ".."), winslash = "/", mustWork = TRUE)
source(file.path(project_root, "analysis_scripts", "common", "los_main_analysis_functions.R"))

build_los_main_analysis(
  input_file = file.path(project_root, "source_csv", "ICU_los.csv"),
  output_dir = file.path(project_root, "computed_data", "los_icu_main"),
  outcome_label = "ICU LOS",
  effect_file = "icu_los_effect_sizes.csv",
  summary_file = "icu_los_meta_summary.csv",
  excluded_file = "icu_los_excluded.csv"
)

message("Built ICU LOS main-analysis data.")
