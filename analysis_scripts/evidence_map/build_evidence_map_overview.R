# Build the Figure 3 evidence-map summary and study-level plotting caches.
#
# Pooled estimates are read from their canonical upstream analysis outputs. This
# script does not refit any meta-analysis model. Small-study-effect diagnostics
# are produced separately by build_small_study_effect_diagnostics.R.

library(dplyr)

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

data_out_dir <- Sys.getenv(
  "BAYES_EVIDENCE_MAP_DATA_DIR",
  unset = file.path(project_root, "computed_data", "evidence_map")
)
if (!grepl("^([A-Za-z]:[\\\\/]|/)", data_out_dir)) {
  data_out_dir <- file.path(project_root, data_out_dir)
}
if (!dir.exists(data_out_dir)) dir.create(data_out_dir, recursive = TRUE)

clean_value <- function(x) {
  x <- trimws(as.character(x))
  x <- gsub("<U\\+2020>", "", x)
  x <- gsub("\u2020", "", x, fixed = TRUE)
  x[x %in% c("", "/", "NA", "Na", "na", "N/A")] <- NA_character_
  x
}

safe_num <- function(x) suppressWarnings(as.numeric(trimws(as.character(x))))
clean_names <- function(x) trimws(sub("^\ufeff", "", x))

fmt_num <- function(x, digits = 3) {
  ifelse(is.na(x), "", sprintf(paste0("%.", digits, "f"), x))
}

fmt_effect <- function(est, lower, upper, digits = 2) {
  paste0(fmt_num(est, digits), "\n[", fmt_num(lower, digits), ", ", fmt_num(upper, digits), "]")
}

read_csv_checked <- function(path, required, description) {
  if (!file.exists(path)) stop("Missing ", description, ": ", path)
  out <- read.csv(path, stringsAsFactors = FALSE, check.names = FALSE, encoding = "UTF-8")
  names(out) <- clean_names(names(out))
  missing <- setdiff(required, names(out))
  if (length(missing) > 0) {
    stop(description, " lacks: ", paste(missing, collapse = ", "))
  }
  out
}

read_existing_summary_row <- function(path, label, effect_type, order) {
  raw <- read_csv_checked(path, "k", paste0("meta-analysis summary for ", label))
  row <- if ("label" %in% names(raw)) {
    raw %>% filter(clean_value(.data$label) == .env$label) %>% slice(1)
  } else if ("group" %in% names(raw)) {
    raw %>% filter(tolower(clean_value(.data$group)) == "overall") %>% slice(1)
  } else {
    raw %>% slice(1)
  }
  if (nrow(row) != 1) stop("Could not identify one summary row for: ", label)

  estimate <- if (effect_type == "MD" && "MD" %in% names(row)) {
    safe_num(row$MD)
  } else if ("estimate" %in% names(row)) {
    safe_num(row$estimate)
  } else {
    stop("Summary lacks estimate column for: ", label)
  }
  lower <- if ("Lower_95" %in% names(row)) safe_num(row$Lower_95) else safe_num(row$lower_95)
  upper <- if ("Upper_95" %in% names(row)) safe_num(row$Upper_95) else safe_num(row$upper_95)
  q_value <- if ("Q" %in% names(row)) safe_num(row$Q) else NA_real_
  q_p_value <- if ("Q_p_value" %in% names(row)) {
    safe_num(row$Q_p_value)
  } else if ("p_heterogeneity" %in% names(row)) {
    safe_num(row$p_heterogeneity)
  } else {
    NA_real_
  }

  data.frame(
    order = order,
    label = label,
    k = safe_num(row$k),
    effect_type = effect_type,
    estimate = estimate,
    lower_95 = lower,
    upper_95 = upper,
    tau2 = safe_num(row$tau2),
    I2_percent = safe_num(row$I2_percent),
    Q = q_value,
    Q_p_value = q_p_value,
    method = clean_value(row$method),
    stringsAsFactors = FALSE
  )
}

canonical_dir <- file.path(project_root, "computed_data", "mortality_main", "canonical")
mortality_summary <- read_csv_checked(
  file.path(canonical_dir, "canonical_mortality_posterior_summary.csv"),
  c("group_id", "k", "RR", "Lower_95_CrI", "Upper_95_CrI", "tau2"),
  "canonical mortality posterior summary"
)

make_mortality_summary_row <- function(label, group_id, order) {
  row <- mortality_summary %>% filter(.data$group_id == .env$group_id)
  if (nrow(row) != 1) stop("Expected one canonical mortality row for: ", group_id)
  data.frame(
    order = order,
    label = label,
    k = safe_num(row$k),
    effect_type = "RR",
    estimate = safe_num(row$RR),
    lower_95 = safe_num(row$Lower_95_CrI),
    upper_95 = safe_num(row$Upper_95_CrI),
    tau2 = safe_num(row$tau2),
    I2_percent = NA_real_,
    Q = NA_real_,
    Q_p_value = NA_real_,
    method = "Bayesian random-effects model",
    stringsAsFactors = FALSE
  )
}

process_summary_path <- function(file_name) {
  file.path(project_root, "computed_data", "process_main", file_name)
}

summary_df <- bind_rows(
  make_mortality_summary_row("Mortality: RCT", "rct_only", 1),
  make_mortality_summary_row("Mortality: NRS", "nrs_only", 2),
  make_mortality_summary_row("Mortality: RCT + NRS", "overall", 3),
  read_existing_summary_row(
    file.path(project_root, "computed_data", "los_hospital_main", "los_meta_summary.csv"),
    "Hospital LOS", "MD", 4
  ),
  read_existing_summary_row(
    file.path(project_root, "computed_data", "los_icu_main", "icu_los_meta_summary.csv"),
    "ICU LOS", "MD", 5
  ),
  read_existing_summary_row(
    process_summary_path("3.1_lactate_measurement_completion_proportion_meta_analysis_summary.csv"),
    "Lactate measurement completion proportion", "RD", 6
  ),
  read_existing_summary_row(
    process_summary_path("3.2_blood_culture_collection_completion_proportion_meta_analysis_summary.csv"),
    "Blood culture collection completion proportion", "RD", 7
  ),
  read_existing_summary_row(
    process_summary_path("3.3_timely_antibiotic_administration_proportion_meta_analysis_summary.csv"),
    "Timely antibiotic administration proportion", "RD", 8
  ),
  read_existing_summary_row(
    process_summary_path("3.4_fluid_resuscitation_completion_proportion_meta_analysis_summary.csv"),
    "Fluid resuscitation completion proportion", "RD", 9
  ),
  read_existing_summary_row(
    process_summary_path("3.5_overall_sepsis_bundle_completion_proportion_meta_analysis_summary.csv"),
    "Overall sepsis bundle completion proportion", "RD", 10
  ),
  read_existing_summary_row(
    process_summary_path("3.6_time_to_antibiotic_administration_meta_analysis_summary.csv"),
    "Time to antibiotic administration", "MD", 11
  )
) %>% arrange(order)

diagnostics <- read_csv_checked(
  file.path(data_out_dir, "evidence_map_forest_funnel_diagnostics.csv"),
  c("label", "egger_intercept", "egger_p", "begg_tau", "begg_p"),
  "Figure 3 small-study-effect diagnostics"
)

summary_df <- summary_df %>%
  left_join(
    diagnostics %>%
      transmute(
        label = clean_value(label),
        egger_intercept = safe_num(egger_intercept),
        egger_p = safe_num(egger_p),
        begg_tau = safe_num(begg_tau),
        begg_p = safe_num(begg_p)
      ),
    by = "label"
  ) %>%
  select(
    order, label, k, effect_type, estimate, lower_95, upper_95,
    tau2, I2_percent, Q, Q_p_value,
    egger_intercept, egger_p, begg_tau, begg_p, method
  )

read_effect_cache <- function(path, label, order, effect_type, mortality = FALSE, design_filter = NULL) {
  required <- if (mortality) {
    c("study", "design", "yi", "vi_original", "vi_bias_adjusted", "n_e", "n_c")
  } else {
    c("study", "design", "yi", "n_e", "n_c")
  }
  raw <- read_csv_checked(path, required, paste0("effect-size cache for ", label))
  vi <- if (mortality) {
    safe_num(raw$vi_bias_adjusted)
  } else if ("vi_original" %in% names(raw)) {
    safe_num(raw$vi_original)
  } else if ("vi" %in% names(raw)) {
    safe_num(raw$vi)
  } else if ("se" %in% names(raw)) {
    safe_num(raw$se)^2
  } else {
    stop("Effect-size cache lacks a variance column for: ", label)
  }
  vi_funnel <- if (mortality) safe_num(raw$vi_original) else vi
  out <- data.frame(
    order = order,
    label = label,
    effect_type = effect_type,
    study = clean_value(raw$study),
    design = toupper(clean_value(raw$design)),
    yi = safe_num(raw$yi),
    vi = vi,
    vi_funnel = vi_funnel,
    n_total = safe_num(raw$n_e) + safe_num(raw$n_c),
    stringsAsFactors = FALSE
  ) %>%
    filter(!is.na(study), design %in% c("RCT", "NRS"), is.finite(yi), is.finite(vi), vi > 0)
  if (!is.null(design_filter)) out <- out %>% filter(design == design_filter)

  summary_row <- summary_df %>% filter(.data$order == .env$order)
  if (nrow(summary_row) != 1) stop("Missing summary row for: ", label)
  if (nrow(out) != as.integer(summary_row$k)) {
    stop(
      "Effect-size rows do not match summary k for ", label,
      ": effect rows = ", nrow(out), ", summary k = ", summary_row$k
    )
  }
  tau2 <- safe_num(summary_row$tau2)
  if (!is.finite(tau2)) tau2 <- 0
  out %>% mutate(random_weight = 100 * (1 / (vi + tau2)) / sum(1 / (vi + tau2)))
}

mortality_effect_file <- file.path(canonical_dir, "canonical_mortality_effect_sizes.csv")
study_df <- bind_rows(
  read_effect_cache(mortality_effect_file, "Mortality: RCT", 1, "RR", mortality = TRUE, design_filter = "RCT"),
  read_effect_cache(mortality_effect_file, "Mortality: NRS", 2, "RR", mortality = TRUE, design_filter = "NRS"),
  read_effect_cache(mortality_effect_file, "Mortality: RCT + NRS", 3, "RR", mortality = TRUE),
  read_effect_cache(
    file.path(project_root, "computed_data", "los_hospital_main", "los_effect_sizes.csv"),
    "Hospital LOS", 4, "MD"
  ),
  read_effect_cache(
    file.path(project_root, "computed_data", "los_icu_main", "icu_los_effect_sizes.csv"),
    "ICU LOS", 5, "MD"
  ),
  read_effect_cache(
    process_summary_path("3.1_lactate_measurement_completion_proportion_meta_analysis_effect_sizes.csv"),
    "Lactate measurement completion proportion", 6, "RD"
  ),
  read_effect_cache(
    process_summary_path("3.2_blood_culture_collection_completion_proportion_meta_analysis_effect_sizes.csv"),
    "Blood culture collection completion proportion", 7, "RD"
  ),
  read_effect_cache(
    process_summary_path("3.3_timely_antibiotic_administration_proportion_meta_analysis_effect_sizes.csv"),
    "Timely antibiotic administration proportion", 8, "RD"
  ),
  read_effect_cache(
    process_summary_path("3.4_fluid_resuscitation_completion_proportion_meta_analysis_effect_sizes.csv"),
    "Fluid resuscitation completion proportion", 9, "RD"
  ),
  read_effect_cache(
    process_summary_path("3.5_overall_sepsis_bundle_completion_proportion_meta_analysis_effect_sizes.csv"),
    "Overall sepsis bundle completion proportion", 10, "RD"
  ),
  read_effect_cache(
    process_summary_path("3.6_time_to_antibiotic_administration_meta_analysis_effect_sizes.csv"),
    "Time to antibiotic administration", 11, "MD"
  )
) %>% arrange(order, study)

grade <- read_csv_checked(
  file.path(project_root, "source_csv", "grade_certainty.csv"),
  c("outcome", "final_certainty", "within_grade_differentiation", "relative_evidentiary_confidence"),
  "GRADE certainty source"
) %>%
  transmute(
    grade_key = tolower(clean_value(outcome)),
    final_certainty = clean_value(final_certainty),
    within_grade_differentiation = clean_value(within_grade_differentiation),
    relative_evidentiary_confidence = clean_value(relative_evidentiary_confidence)
  )

plot_summary <- summary_df %>%
  mutate(
    display_label = case_when(
      label == "Lactate measurement completion proportion" ~ "Lactate measurement",
      label == "Blood culture collection completion proportion" ~ "Blood culture",
      label == "Time to antibiotic administration" ~ "Antibiotic time",
      label == "Timely antibiotic administration proportion" ~ "Antibiotic proportion",
      label == "Fluid resuscitation completion proportion" ~ "Fluid resuscitation",
      label == "Overall sepsis bundle completion proportion" ~ "Sepsis bundle",
      TRUE ~ label
    ),
    grade_outcome = case_when(
      order <= 3 ~ label,
      display_label == "Hospital LOS" ~ "Hospital length of stay",
      display_label == "ICU LOS" ~ "ICU length of stay",
      TRUE ~ display_label
    ),
    grade_key = tolower(clean_value(grade_outcome))
  ) %>%
  left_join(grade, by = "grade_key") %>%
  select(-grade_key, -grade_outcome)

write.csv(
  summary_df,
  file.path(data_out_dir, "evidence_map_forest_funnel_overview_summary.csv"),
  row.names = FALSE,
  na = ""
)
write.csv(
  plot_summary,
  file.path(data_out_dir, "evidence_map_forest_funnel_plot_summary.csv"),
  row.names = FALSE,
  na = ""
)
write.csv(
  study_df,
  file.path(data_out_dir, "evidence_map_forest_funnel_plot_studies.csv"),
  row.names = FALSE,
  na = ""
)

message("Wrote Figure 3 evidence-map plotting caches to: ", data_out_dir)
