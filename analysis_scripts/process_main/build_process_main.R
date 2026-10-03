library(dplyr)

script_arg <- c(grep("^--file=", commandArgs(FALSE), value = TRUE), grep("\\.R$", commandArgs(FALSE), value = TRUE, ignore.case = TRUE))
script_dir <- if (length(script_arg) > 0) dirname(normalizePath(sub("^--file=", "", script_arg[1]), mustWork = FALSE)) else getwd()
project_root <- normalizePath(file.path(script_dir, "..", ".."), winslash = "/", mustWork = TRUE)
options(sepsis_code_dir = project_root)
source(file.path(project_root, "analysis_scripts", "common", "process_meta_analysis_functions.R"))

out_dir <- file.path(project_root, "computed_data", "process_main")
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

specs <- list(
  list(1, "process_lactate_measurement.csv", "lactate_measurement_completion_proportion_meta_analysis", tolower(process_titles[[1]])),
  list(2, "process_blood_culture.csv", "blood_culture_collection_completion_proportion_meta_analysis", tolower(process_titles[[2]])),
  list(3, "process_antibiotic_proportion.csv", "timely_antibiotic_administration_proportion_meta_analysis", tolower(process_titles[[3]])),
  list(4, "process_fluid_resuscitation.csv", "fluid_resuscitation_completion_proportion_meta_analysis", tolower(process_titles[[4]])),
  list(5, "process_sepsis_bundle.csv", "overall_sepsis_bundle_completion_proportion_meta_analysis", tolower(process_titles[[5]]))
)

for (spec in specs) {
  build_binary_process_main(
    index = spec[[1]],
    csv_file = file.path(project_root, "source_csv", spec[[2]]),
    folder_slug = spec[[3]],
    title = spec[[4]],
    out_dir = out_dir
  )
}
build_antibiotic_time_main(out_dir)

message("Wrote six primary process-outcome analysis caches to: ", out_dir)
