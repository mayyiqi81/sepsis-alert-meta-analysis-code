# Small-study-effect diagnostics for the Figure 3 evidence-map overview.
#
# This script writes Egger/Begg diagnostics used by
# Figure_3_meta_analysis_summary.R. It does not assemble figures
# and does not calculate pooled meta-analytic summaries.

library(dplyr)

script_arg <- c(grep("^--file=", commandArgs(FALSE), value = TRUE), grep("\\.R$", commandArgs(FALSE), value = TRUE, ignore.case = TRUE))
script_dir <- if (length(script_arg) > 0) dirname(normalizePath(sub("^--file=", "", script_arg[1]), mustWork = FALSE)) else getwd()
project_root <- normalizePath(file.path(script_dir, "..", ".."), winslash = "/", mustWork = TRUE)
source(file.path(project_root, "analysis_scripts", "common", "project_paths.R"))
source(file.path(project_root, "analysis_scripts", "common", "small_study_effect_functions.R"))
options(sepsis_project_root = set_project_root_from_script(project_root))
script_dir <- project_root

data_out_dir <- Sys.getenv(
  "BAYES_EVIDENCE_MAP_DATA_DIR",
  unset = "computed_data/evidence_map"
)
if (!grepl("^([A-Za-z]:[\\\\/]|/)", data_out_dir)) data_out_dir <- file.path(script_dir, data_out_dir)
if (!dir.exists(data_out_dir)) dir.create(data_out_dir, recursive = TRUE)

plot_prefix <- "evidence_map_forest_funnel_overview"
diagnostics_file_name <- "evidence_map_forest_funnel_diagnostics.csv"

clean_value <- function(x) {
  x <- trimws(as.character(x))
  x[x %in% c("", "/", "NA", "Na", "na", "N/A")] <- NA_character_
  x
}

safe_num <- function(x) {
  suppressWarnings(as.numeric(trimws(as.character(x))))
}

clean_names <- function(x) {
  trimws(sub("^\ufeff", "", x))
}

read_mortality_es <- function() {
  path <- file.path(
    script_dir,
    "computed_data",
    "mortality_main",
    "canonical",
    "canonical_mortality_effect_sizes.csv"
  )
  if (!file.exists(path)) stop("Missing canonical mortality effect-size cache: ", path)
  raw <- read.csv(path, stringsAsFactors = FALSE, check.names = FALSE, encoding = "UTF-8")
  names(raw) <- clean_names(names(raw))
  required <- c("study", "design", "yi", "vi_original", "vi_bias_adjusted", "n_e", "n_c")
  missing <- setdiff(required, names(raw))
  if (length(missing) > 0) {
    stop("Canonical mortality effect-size cache lacks: ", paste(missing, collapse = ", "))
  }

  raw %>%
    transmute(
      study = clean_value(study),
      design = toupper(clean_value(design)),
      yi = safe_num(yi),
      vi_original = safe_num(vi_original),
      vi_bias_adjusted = safe_num(vi_bias_adjusted),
      vi = vi_bias_adjusted,
      vi_funnel = vi_original,
      n_total = safe_num(n_e) + safe_num(n_c)
    ) %>%
    filter(
      !is.na(study), design %in% c("RCT", "NRS"),
      is.finite(yi), is.finite(vi), vi > 0,
      is.finite(vi_funnel), vi_funnel > 0
    )
}

read_formal_continuous_effect_cache <- function(path, outcome_label) {
  if (!file.exists(path)) stop("Missing formal ", outcome_label, " effect-size cache: ", path)
  raw <- read.csv(path, stringsAsFactors = FALSE, check.names = FALSE, encoding = "UTF-8")
  names(raw) <- clean_names(names(raw))
  required <- c("study", "design", "yi", "vi_original", "n_e", "n_c")
  missing <- setdiff(required, names(raw))
  if (length(missing) > 0) {
    stop("Formal ", outcome_label, " effect-size cache lacks: ", paste(missing, collapse = ", "))
  }
  raw %>%
    transmute(
      study = study,
      design = design,
      yi = yi,
      vi = safe_num(vi_original),
      vi_funnel = vi,
      n_e = safe_num(n_e),
      n_c = safe_num(n_c),
      n_total = n_e + n_c
    ) %>%
    filter(!is.na(study), design %in% c("RCT", "NRS"), is.finite(yi), is.finite(vi), vi > 0)
}

read_effect_size_file <- function(path) {
  if (!file.exists(path)) {
    stop("Missing effect-size file: ", path)
  }

  raw <- read.csv(path, stringsAsFactors = FALSE, check.names = FALSE, encoding = "UTF-8")
  names(raw) <- clean_names(names(raw))
  required <- c("study", "design", "yi")
  missing <- setdiff(required, names(raw))
  if (length(missing) > 0 || (!"vi" %in% names(raw) && !"se" %in% names(raw))) {
    stop(
      "Effect-size cache lacks required study/design/yi and vi-or-se columns: ",
      path
    )
  }
  dat <- raw %>%
    transmute(
      study = clean_value(study),
      design = if ("design" %in% names(raw)) toupper(clean_value(design)) else NA_character_,
      yi = safe_num(yi),
      vi = if ("vi" %in% names(raw)) safe_num(vi) else safe_num(se)^2,
      vi_funnel = vi,
      n_e = if ("n_e" %in% names(raw)) safe_num(n_e) else NA_real_,
      n_c = if ("n_c" %in% names(raw)) safe_num(n_c) else NA_real_,
      n_total = n_e + n_c
    ) %>%
    filter(!is.na(study), !is.na(yi), !is.na(vi), vi > 0)

  if (any(!dat$design %in% c("RCT", "NRS"))) {
    stop("Every cached study must have Design equal to RCT or NRS: ", path)
  }
  dat
}

read_3_1_process_effect_size <- function(file_name) {
  read_effect_size_file(file.path(script_dir, "computed_data", "process_main", file_name))
}

diagnostic_row <- function(label, es) {
  es <- es %>%
    mutate(vi_funnel = ifelse(is.finite(vi_funnel) & vi_funnel > 0, vi_funnel, vi))
  cbind(
    data.frame(label = label, k = nrow(es), stringsAsFactors = FALSE),
    publication_bias_tests(es, vi_col = "vi_funnel")
  )
}

mortality_es <- read_mortality_es()

diagnostics <- bind_rows(
  diagnostic_row("Mortality: RCT", mortality_es %>% filter(design == "RCT")),
  diagnostic_row("Mortality: NRS", mortality_es %>% filter(design == "NRS")),
  diagnostic_row("Mortality: RCT + NRS", mortality_es),
  diagnostic_row("Hospital LOS", read_formal_continuous_effect_cache(file.path(script_dir, "computed_data", "los_hospital_main", "los_effect_sizes.csv"), "hospital LOS")),
  diagnostic_row("ICU LOS", read_formal_continuous_effect_cache(file.path(script_dir, "computed_data", "los_icu_main", "icu_los_effect_sizes.csv"), "ICU LOS")),
  diagnostic_row("Lactate measurement completion proportion", read_3_1_process_effect_size("3.1_lactate_measurement_completion_proportion_meta_analysis_effect_sizes.csv")),
  diagnostic_row("Blood culture collection completion proportion", read_3_1_process_effect_size("3.2_blood_culture_collection_completion_proportion_meta_analysis_effect_sizes.csv")),
  diagnostic_row("Timely antibiotic administration proportion", read_3_1_process_effect_size("3.3_timely_antibiotic_administration_proportion_meta_analysis_effect_sizes.csv")),
  diagnostic_row("Fluid resuscitation completion proportion", read_3_1_process_effect_size("3.4_fluid_resuscitation_completion_proportion_meta_analysis_effect_sizes.csv")),
  diagnostic_row("Overall sepsis bundle completion proportion", read_3_1_process_effect_size("3.5_overall_sepsis_bundle_completion_proportion_meta_analysis_effect_sizes.csv")),
  diagnostic_row("Time to antibiotic administration", read_3_1_process_effect_size("3.6_time_to_antibiotic_administration_meta_analysis_effect_sizes.csv"))
)

write.csv(
  diagnostics,
  file.path(data_out_dir, diagnostics_file_name),
  row.names = FALSE,
  na = ""
)

message("Wrote Figure 3 small-study diagnostics to: ", data_out_dir)
