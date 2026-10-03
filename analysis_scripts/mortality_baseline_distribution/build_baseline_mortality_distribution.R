# Figure S6F descriptive-data builder.

library(dplyr)

script_arg <- c(grep("^--file=", commandArgs(FALSE), value = TRUE), grep("\\.R$", commandArgs(FALSE), value = TRUE, ignore.case = TRUE))
script_dir <- if (length(script_arg) > 0) dirname(normalizePath(sub("^--file=", "", script_arg[1]), mustWork = FALSE)) else getwd()
project_root <- normalizePath(file.path(script_dir, "..", ".."), winslash = "/", mustWork = TRUE)
workspace_dir <- project_root
source(file.path(workspace_dir, "analysis_scripts", "common", "bayes_mortality_population_helpers.R"))

baseline_mortality_to_percent_display <- function(x) x * 100
out_dir <- file.path(workspace_dir, "computed_data", "mortality_baseline_distribution")
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
output_prefix <- "Figure_S6F_baseline_mortality_by_analytic_population"
input_file <- file.path(workspace_dir, "source_csv", "mortality.csv")

# 2. Population mapping and audit data -------------------------------------

safe_num <- function(x) {
  suppressWarnings(as.numeric(gsub("%", "", trimws(as.character(x)), fixed = TRUE)))
}

population_levels <- c(
  "all_screened",
  "infection_targeted",
  "organ_dysfunction_severe_sepsis",
  "organ_dysfunction_sepsis3"
)
population_labels <- c(
  all_screened = "All screened",
  infection_targeted = "Infection-targeted",
  organ_dysfunction_severe_sepsis = "Organ dysfunction-targeted\n(severe sepsis)",
  organ_dysfunction_sepsis3 = "Organ dysfunction-targeted\n(Sepsis-3)"
)

derive_population <- function(scope) {
  dplyr::case_when(
    is_all_screened(scope) ~ "all_screened",
    is_infection_targeted(scope) ~ "infection_targeted",
    is_severe_sepsis_targeted(scope) ~ "organ_dysfunction_severe_sepsis",
    is_sepsis3_targeted(scope) ~ "organ_dysfunction_sepsis3",
    TRUE ~ NA_character_
  )
}

raw <- read.csv(input_file, stringsAsFactors = FALSE, check.names = FALSE, encoding = "UTF-8")
names(raw) <- trimws(sub("^\\ufeff", "", names(raw)))
required_cols <- c(
  "Study",
  "Design",
  "events_c",
  "n_c",
  "Baseline mortality rate in control group",
  "Clinical population scope"
)
missing_cols <- setdiff(required_cols, names(raw))
if (length(missing_cols) > 0) stop("Missing required mortality columns: ", paste(missing_cols, collapse = ", "))

plot_data <- raw %>%
  transmute(
    study = trimws(Study),
    design = toupper(trimws(Design)),
    events_c = safe_num(events_c),
    n_c = safe_num(n_c),
    baseline_mortality = safe_num(.data[["Baseline mortality rate in control group"]]),
    baseline_mortality_from_counts = events_c / n_c,
    clinical_population_scope = trimws(.data[["Clinical population scope"]]),
    population_id = derive_population(.data[["Clinical population scope"]])
  ) %>%
  mutate(
    baseline_mortality_percent = baseline_mortality_to_percent_display(baseline_mortality)
  )

invalid_rows <- plot_data %>%
  filter(
    is.na(study) | study == "" | !design %in% c("RCT", "NRS") |
      is.na(events_c) | !is.finite(events_c) | events_c < 0 |
      is.na(n_c) | !is.finite(n_c) | n_c <= 0 | events_c > n_c |
      is.na(baseline_mortality) | !is.finite(baseline_mortality) |
      baseline_mortality < 0 | baseline_mortality > 1 |
      is.na(baseline_mortality_from_counts) |
      abs(baseline_mortality - baseline_mortality_from_counts) > 1e-8 |
      is.na(population_id)
  )
if (nrow(invalid_rows) > 0) {
  stop(
    "S6F found invalid, inconsistent, or unclassified mortality rows; ",
    "the formal baseline mortality field must agree with events_c / n_c within 1e-8. Count: ",
    nrow(invalid_rows)
  )
}
if (nrow(plot_data) != 63) {
  stop("Expected 63 mortality rows for S6F, found ", nrow(plot_data), ".")
}
if (anyDuplicated(plot_data$study) > 0) {
  stop("S6F requires one current mortality row per study; duplicate study rows were found.")
}

population_counts <- table(factor(plot_data$population_id, levels = population_levels))
if (any(population_counts == 0)) {
  stop("S6F population layer is empty: ", paste(names(population_counts)[population_counts == 0], collapse = ", "))
}
if (sum(population_counts) != nrow(plot_data)) {
  stop("S6F population layers are not exhaustive.")
}

plot_data <- plot_data %>%
  mutate(
    population_id = factor(population_id, levels = rev(population_levels)),
    population_label = factor(
      unname(population_labels[as.character(population_id)]),
      levels = rev(unname(population_labels))
    )
  )

summary_table <- plot_data %>%
  group_by(population_id, population_label) %>%
  summarise(
    k = n(),
    n_rct = sum(design == "RCT"),
    n_nrs = sum(design == "NRS"),
    min = min(baseline_mortality),
    Q1 = as.numeric(quantile(baseline_mortality, 0.25, names = FALSE)),
    median = median(baseline_mortality),
    Q3 = as.numeric(quantile(baseline_mortality, 0.75, names = FALSE)),
    max = max(baseline_mortality),
    .groups = "drop"
  ) %>%
  mutate(
    group = as.character(population_label),
    population_id = as.character(population_id),
    median_label = paste0(
      "k=", k, "\nmedian ",
      sprintf("%.1f%%", baseline_mortality_to_percent_display(median))
    ),
    .before = 1
  ) %>%
  arrange(match(population_id, population_levels))

metadata_table <- data.frame(
  key = c(
    "input_file",
    "population_source",
    "baseline_mortality_source",
    "baseline_mortality_qa",
    "population_rule",
    "n_rows",
    "n_rows_in_layers"
  ),
  value = c(
    normalizePath(input_file, mustWork = FALSE),
    "Clinical population scope",
    "Baseline mortality rate in control group",
    "events_c / n_c agrees with the formal baseline field within 1e-8; the formal field is used for display",
    "All screened; Infection-targeted; Severe sepsis-targeted; Sepsis-3-targeted; the latter two are displayed as underlying definitions within organ dysfunction-targeted.",
    nrow(plot_data),
    sum(summary_table$k)
  ),
  stringsAsFactors = FALSE
)

utils::write.csv(summary_table, file.path(out_dir, paste0(output_prefix, "_summary.csv")), row.names = FALSE, fileEncoding = "UTF-8")
utils::write.csv(plot_data, file.path(out_dir, paste0(output_prefix, "_data.csv")), row.names = FALSE, fileEncoding = "UTF-8")
utils::write.csv(metadata_table, file.path(out_dir, paste0(output_prefix, "_metadata.csv")), row.names = FALSE, fileEncoding = "UTF-8")
message("Wrote Figure S6F descriptive cache to: ", out_dir)
