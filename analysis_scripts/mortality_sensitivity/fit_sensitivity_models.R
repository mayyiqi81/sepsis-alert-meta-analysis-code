# Mortality sensitivity analyses: final pooled results only.
# 1.Primary analysis(NRS bias SD = 0.20)
# No down-weighting reference
# NRS Bias SD sensitivity
# Power-prior sensitivity
# Exclusion analyses
# Effect-estimate availability
# Population sensitivity

library(metafor)
library(dplyr)
library(forestplot)
library(grid)

script_arg <- c(grep("^--file=", commandArgs(FALSE), value = TRUE), grep("\\.R$", commandArgs(FALSE), value = TRUE, ignore.case = TRUE))
script_dir <- if (length(script_arg) > 0) dirname(normalizePath(sub("^--file=", "", script_arg[1]), mustWork = FALSE)) else getwd()
project_root <- normalizePath(file.path(script_dir, "..", ".."), winslash = "/", mustWork = TRUE)
source(file.path(project_root, "analysis_scripts", "common", "project_paths.R"))
options(sepsis_project_root = set_project_root_from_script(project_root))
script_dir <- project_root
workspace_dir <- script_dir
source(file.path(workspace_dir, "analysis_scripts", "common", "bayes_prior_scale_spec.R"))
source(file.path(workspace_dir, "analysis_scripts", "common", "bayes_mortality_population_helpers.R"))
source(file.path(workspace_dir, "analysis_scripts", "common", "bayes_mortality_functions.R"))
input_file <- file.path(workspace_dir, "source_csv", "mortality.csv")
out_dir <- Sys.getenv(
  "MORTALITY_SENSITIVITY_OUT_DIR",
  unset = "computed_data/mortality_sensitivity"
)
if (!grepl("^([A-Za-z]:[\\\\/]|/)", out_dir)) out_dir <- file.path(workspace_dir, out_dir)
canonical_out_dir <- Sys.getenv(
  "BAYES_CANONICAL_OUT_DIR",
  unset = "computed_data/mortality_main/canonical"
)
if (!grepl("^([A-Za-z]:[\\\\/]|/)", canonical_out_dir)) canonical_out_dir <- file.path(workspace_dir, canonical_out_dir)
plot_font <- "Segoe UI"
fit_control <- list(maxiter = 10000, stepadj = 0.1)
bayes_settings <- bayes_meta_settings_from_env()
primary_bias_sd <- 0.20
primary_bias_sd_label <- sprintf("%.2f", primary_bias_sd)
color_text <- "#0A265E"
color_axis <- color_text
color_main <- "#5A5368"
color_model_sensitivity <- "#8A748F"
color_population_sensitivity <- "#CFC2D2"
color_section_title <- "#5A5368"
color_model_title <- "#8A748F"
color_population_title <- "#5A5368"

if (!dir.exists(out_dir)) {
  dir.create(out_dir, recursive = TRUE)
}
fits_dir <- file.path(out_dir, "fits")
dir.create(fits_dir, recursive = TRUE, showWarnings = FALSE)

clean_value <- function(x) {
  x <- trimws(as.character(x))
  x[x %in% c("", "/", "NA", "Na", "na", "N/A")] <- NA_character_
  x
}

safe_num <- function(x) {
  suppressWarnings(as.numeric(trimws(as.character(x))))
}

fmt_rr <- function(rr, lower, upper) {
  sprintf("%.3f [%.3f, %.3f]", rr, lower, upper)
}

safe_model_value <- function(x) {
  ifelse(is.null(x) || length(x) == 0, NA_real_, as.numeric(x)[1])
}

nrs_bias_variance <- function(dat, sigma_bias) {
  ifelse(dat$is_nrs, sigma_bias^2, 0)
}

nrs_bias_note <- function(prefix, retained = FALSE) {
  scope <- ifelse(retained, "retained NRS", "NRS")
  paste0(prefix, "; all ", scope, " variance inflated by bias SD = ", primary_bias_sd_label)
}

assert_expected_counts <- function(dat, label, expected_rct, expected_nrs, expected_total) {
  observed_rct <- sum(dat$design == "RCT", na.rm = TRUE)
  observed_nrs <- sum(dat$design == "NRS", na.rm = TRUE)
  observed_total <- nrow(dat)
  if (observed_rct != expected_rct ||
      observed_nrs != expected_nrs ||
      observed_total != expected_total) {
    stop(
      label, " count check failed: observed RCT=", observed_rct,
      ", NRS=", observed_nrs,
      ", total=", observed_total,
      "; expected RCT=", expected_rct,
      ", NRS=", expected_nrs,
      ", total=", expected_total,
      call. = FALSE
    )
  }
  invisible(TRUE)
}

reml_reference_results <- list()
canonical_bayes <- NULL
sensitivity_fit_inventory <- list()

sensitivity_fit_slug <- function(x) {
  slug <- tolower(trimws(as.character(x)))
  slug <- gsub("[^a-z0-9]+", "_", slug)
  gsub("(^_|_$)", "", slug)
}

canonical_group_for_sensitivity <- function(analysis, scenario) {
  if (!identical(analysis, "1.Primary analysis(NRS bias SD = 0.20)")) {
    return(NA_character_)
  }
  dplyr::case_when(
    identical(scenario, "RCT only") ~ "rct_only",
    identical(scenario, "NRS only: bias SD = 0.20") ~ "nrs_only",
    identical(scenario, "Overall: bias SD = 0.20") ~ "overall",
    TRUE ~ NA_character_
  )
}

fit_adjusted <- function(dat, vi_adjusted, analysis, scenario, note) {
  model_dat <- dat
  model_dat$vi_adjusted_current <- vi_adjusted
  canonical_group_id <- canonical_group_for_sensitivity(analysis, scenario)
  if (!is.na(canonical_group_id)) {
    if (is.null(canonical_bayes)) {
      canonical_bayes <<- load_canonical_bayes_mortality(
        canonical_out_dir,
        require_fits = FALSE,
        input_file = input_file,
        expected_settings = bayes_settings,
        expected_nrs_bias_sd = primary_bias_sd
      )
      assert_canonical_effect_sizes_match(es, canonical_bayes$effect_sizes)
    }
    bayes_summary <- canonical_summary_row(canonical_bayes, canonical_group_id)
    fit_path <- file.path(canonical_out_dir, "fits", paste0(canonical_group_id, ".brmsfit.rds"))
    fit_source <- "canonical_main"
  } else {
    fit_id <- paste(
      sensitivity_fit_slug(analysis),
      sensitivity_fit_slug(scenario),
      sep = "__"
    )
    fit_path <- file.path(fits_dir, paste0(fit_id, ".brmsfit.rds"))
    expected_metadata <- make_bayes_model_cache_metadata(
      data = model_dat,
      vi_col = "vi_adjusted_current",
      response_scale = "log_rr",
      prior_spec_id = "log_rr_shared_current",
      uncertainty_spec_id = "sensitivity_defined_variance_log_rr",
      settings = bayes_settings,
      study_col = "study",
      yi_col = "yi",
      formula_text = "yi | se(.bayes_se, sigma = FALSE) ~ 1 + (1 | .bayes_study)"
    )
    force_refit <- identical(Sys.getenv("BAYES_FORCE_REFIT", unset = ""), "1")
    fit <- NULL
    if (file.exists(fit_path) && !force_refit) {
      cached_fit <- readRDS(fit_path)
      if (isTRUE(assert_cached_bayes_model_matches(cached_fit, expected_metadata))) fit <- cached_fit
    }
    if (is.null(fit)) {
      fit <- fit_bayes_mortality_meta_with_quality_gate(
        model_dat,
        "vi_adjusted_current",
        response_scale = "log_rr",
        prior_spec_id = "log_rr_shared_current",
        uncertainty_spec_id = "sensitivity_defined_variance_log_rr",
        settings = bayes_settings
      )
      saveRDS(fit, fit_path)
    }
    bayes_summary <- summarise_bayes_mortality_fit(
      fit,
      group = paste(analysis, scenario, sep = " | "),
      k_rct = sum(dat$design == "RCT"),
      k_nrs = sum(dat$design == "NRS"),
      nrs_bias_sd = primary_bias_sd
    )
    fit_source <- "mortality_sensitivity_cache"
  }
  sensitivity_fit_inventory[[length(sensitivity_fit_inventory) + 1L]] <<- data.frame(
    analysis = analysis,
    scenario = scenario,
    k = nrow(dat),
    fit_source = fit_source,
    fit_file = gsub("\\\\", "/", substring(normalizePath(fit_path, winslash = "/", mustWork = FALSE), nchar(workspace_dir) + 2L)),
    stringsAsFactors = FALSE
  )
  fit_reml <- rma(
    yi = dat$yi,
    vi = vi_adjusted,
    method = "REML",
    control = fit_control
  )
  reml_reference_results[[length(reml_reference_results) + 1]] <<- summarise_reml_mortality_reference(
    fit_reml,
    paste(analysis, scenario, sep = " | "),
    sum(dat$design == "RCT"),
    sum(dat$design == "NRS"),
    primary_bias_sd
  )
  random_weights <- 1 / (vi_adjusted + bayes_summary$tau2)
  rct_weight_percent <- ifelse(
    sum(random_weights, na.rm = TRUE) > 0,
    100 * sum(random_weights[dat$design == "RCT"], na.rm = TRUE) / sum(random_weights, na.rm = TRUE),
    NA_real_
  )

  data.frame(
    analysis = analysis,
    scenario = scenario,
    k = bayes_summary$k,
    k_rct = sum(dat$design == "RCT"),
    k_nrs = sum(dat$design == "NRS"),
    RR = bayes_summary$RR,
    Lower_95_CrI = bayes_summary$Lower_95_CrI,
    Upper_95_CrI = bayes_summary$Upper_95_CrI,
    RR_label = format_bayes_rr(bayes_summary$RR, bayes_summary$Lower_95_CrI, bayes_summary$Upper_95_CrI),
    tau = bayes_summary$tau,
    tau2 = bayes_summary$tau2,
    rhat_max = bayes_summary$rhat_max,
    ess_bulk_min = bayes_summary$ess_bulk_min,
    ess_tail_min = bayes_summary$ess_tail_min,
    divergent_transitions = bayes_summary$divergent_transitions,
    max_treedepth = bayes_summary$max_treedepth,
    max_treedepth_hits = bayes_summary$max_treedepth_hits,
    method = "Bayesian random-effects model",
    note = note,
    rct_weight_percent = rct_weight_percent,
    rct_weight_label = ifelse(is.na(rct_weight_percent), "", sprintf("%.1f%%", rct_weight_percent)),
    tau_label = ifelse(is.na(bayes_summary$tau), "", sprintf("%.3f", bayes_summary$tau)),
    stringsAsFactors = FALSE
  )
}

add_result <- function(results, dat, vi_adjusted, analysis, scenario, note) {
  results[[length(results) + 1]] <- fit_adjusted(
    dat = dat,
    vi_adjusted = vi_adjusted,
    analysis = analysis,
    scenario = scenario,
    note = note
  )
  results
}

add_nrs_bias_result <- function(results, sigma_bias, analysis, scenario) {
  add_result(
    results,
    es,
    es$vi_original + nrs_bias_variance(es, sigma_bias),
    analysis,
    scenario,
    paste0("RCT unchanged; all NRS variance inflated by bias SD = ", sprintf("%.2f", sigma_bias))
  )
}

add_nrs_bias_dataset_result <- function(results, dat, sigma_bias, analysis, scenario, note) {
  add_result(
    results,
    dat,
    dat$vi_original + nrs_bias_variance(dat, sigma_bias),
    analysis,
    scenario,
    note
  )
}

add_power_prior_result <- function(results, alpha, scenario) {
  add_result(
    results,
    es,
    ifelse(es$is_nrs, es$vi_original / alpha, es$vi_original),
    "Power-prior sensitivity",
    scenario,
    paste0("NRS variance multiplied by 1/alpha; alpha = ", alpha)
  )
}

read_effect_data <- function(file, use_adjusted = FALSE) {
  raw <- read.csv(file, stringsAsFactors = FALSE, check.names = FALSE)
  names(raw) <- trimws(names(raw))
  raw <- raw[, names(raw) != "", drop = FALSE]

  required <- c("Study", "events_e", "n_e", "events_c", "n_c", "Design")
  missing_cols <- setdiff(required, names(raw))
  if (length(missing_cols) > 0) {
    stop("Missing required columns in ", file, ": ", paste(missing_cols, collapse = ", "))
  }

  validation_data <- data.frame(
    study = raw$Study,
    design = raw$Design,
    events_e = safe_num(raw$events_e),
    n_e = safe_num(raw$n_e),
    events_c = safe_num(raw$events_c),
    n_c = safe_num(raw$n_c),
    aRR = if ("aRR" %in% names(raw)) safe_num(raw$aRR) else NA_real_,
    aRR_lower_95 = if ("aRR_lower_95" %in% names(raw)) safe_num(raw$aRR_lower_95) else NA_real_,
    aRR_upper_95 = if ("aRR_upper_95" %in% names(raw)) safe_num(raw$aRR_upper_95) else NA_real_,
    stringsAsFactors = FALSE
  )
  validate_mortality_2x2(validation_data, context = basename(file))
  validate_adjusted_rr_fields(validation_data, context = basename(file))

  dat <- raw %>%
    mutate(
      study = clean_value(Study),
      design = toupper(clean_value(Design)),
      events_e = safe_num(events_e),
      n_e = safe_num(n_e),
      events_c = safe_num(events_c),
      n_c = safe_num(n_c)
    ) %>%
    filter(
      !is.na(study),
      !is.na(events_e), !is.na(n_e),
      !is.na(events_c), !is.na(n_c),
      design %in% c("RCT", "NRS")
    )

  if (use_adjusted && "aRR" %in% names(raw)) {
    dat$aRR <- safe_num(raw$aRR[match(dat$study, raw$Study)])
  } else {
    dat$aRR <- rep(NA_real_, nrow(dat))
  }

  dat
}

make_effect_data <- function(dat, use_adjusted = FALSE, adjusted_only = FALSE) {
  if (!"aRR" %in% names(dat)) {
    dat$aRR <- rep(NA_real_, nrow(dat))
  }
  if (adjusted_only) {
    dat <- dat %>% filter(!is.na(aRR))
  }

  escalc(
    measure = "RR",
    ai = events_e,
    bi = n_e - events_e,
    ci = events_c,
    di = n_c - events_c,
    data = dat,
    add = 0.5,
    to = "all"
  ) %>%
    mutate(
      yi_event = yi,
      yi = ifelse(use_adjusted & !is.na(aRR), log(aRR), yi_event),
      vi_original = vi,
      is_nrs = design == "NRS"
    )
}

df_raw <- read.csv(input_file, stringsAsFactors = FALSE, check.names = FALSE, fileEncoding = "UTF-8-BOM")
names(df_raw) <- tolower(trimws(gsub("^\ufeff", "", names(df_raw))))
df_raw <- df_raw[, names(df_raw) != "", drop = FALSE]

required_cols <- c(
  "study", "events_e", "n_e", "events_c", "n_c", "design",
  "publication type", "clinical population scope", "intervention package type",
  "comparator condition", "additional cohort restriction",
  "adjusted mortality estimate not used in primary synthesis"
)
missing_cols <- setdiff(required_cols, names(df_raw))
if (length(missing_cols) > 0) {
  stop("Missing required columns: ", paste(missing_cols, collapse = ", "))
}

availability_flag <- tolower(trimws(df_raw[["adjusted mortality estimate not used in primary synthesis"]]))
if (anyNA(availability_flag) || !all(availability_flag %in% c("yes", "no"))) {
  stop("Adjusted mortality estimate availability must be explicitly marked Yes or No.")
}

optional_cols <- intersect(c("arr", "arr_lower_95", "arr_upper_95"), names(df_raw))

validation_data <- data.frame(
  study = df_raw$study,
  design = df_raw$design,
  events_e = safe_num(df_raw$events_e),
  n_e = safe_num(df_raw$n_e),
  events_c = safe_num(df_raw$events_c),
  n_c = safe_num(df_raw$n_c),
  aRR = if ("arr" %in% names(df_raw)) safe_num(df_raw$arr) else NA_real_,
  aRR_lower_95 = if ("arr_lower_95" %in% names(df_raw)) safe_num(df_raw$arr_lower_95) else NA_real_,
  aRR_upper_95 = if ("arr_upper_95" %in% names(df_raw)) safe_num(df_raw$arr_upper_95) else NA_real_,
  stringsAsFactors = FALSE
)
validate_mortality_2x2(validation_data, context = basename(input_file))
validate_adjusted_rr_fields(validation_data, context = basename(input_file))

df <- df_raw[, c(required_cols, optional_cols)] %>%
  mutate(
    study = clean_value(study),
    design = toupper(clean_value(design)),
    publication_type = tolower(clean_value(`publication type`)),
    clinical_population_scope = clean_value(`clinical population scope`),
    intervention_package_type = clean_value(`intervention package type`),
    comparator_condition = clean_value(`comparator condition`),
    additional_cohort_restriction = tolower(clean_value(`additional cohort restriction`)),
    adjusted_estimate_not_used = tolower(trimws(`adjusted mortality estimate not used in primary synthesis`)) == "yes",
    events_e = safe_num(events_e),
    n_e = safe_num(n_e),
    events_c = safe_num(events_c),
    n_c = safe_num(n_c),
    aRR = if ("arr" %in% names(.)) safe_num(arr) else NA_real_,
    aRR_lower_95 = if ("arr_lower_95" %in% names(.)) safe_num(arr_lower_95) else NA_real_,
    aRR_upper_95 = if ("arr_upper_95" %in% names(.)) safe_num(arr_upper_95) else NA_real_
  ) %>%
  filter(
    !is.na(study),
    !is.na(events_e), !is.na(n_e),
    !is.na(events_c), !is.na(n_c),
    design %in% c("RCT", "NRS")
  )

if (nrow(df) != 63) {
  stop("Expected 63 mortality rows after filtering; found ", nrow(df), call. = FALSE)
}

es <- escalc(
  measure = "RR",
  ai = events_e,
  bi = n_e - events_e,
  ci = events_c,
  di = n_c - events_c,
  data = df,
  add = 0.5,
  to = "all"
) %>%
  mutate(
    yi_event = yi,
    uses_adjusted_rr = !is.na(aRR) & !is.na(aRR_lower_95) & !is.na(aRR_upper_95),
    yi = ifelse(uses_adjusted_rr, log(aRR), yi_event),
    vi_original = ifelse(uses_adjusted_rr, ((log(aRR_upper_95) - log(aRR_lower_95)) / (2 * 1.96))^2, vi),
    vi_bias_adjusted = ifelse(design == "NRS", vi_original + primary_bias_sd^2, vi_original),
    is_nrs = design == "NRS",
    is_conference_nrs = is_nrs & grepl("conference", publication_type, fixed = TRUE)
  )

rct_only <- es %>% filter(design == "RCT")
nrs_only <- es %>% filter(design == "NRS")
assert_expected_counts(es, "Primary RR dataset", 10, 53, 63)
results <- list()

results <- add_result(
  results,
  rct_only,
  rct_only$vi_original,
  "1.Primary analysis(NRS bias SD = 0.20)",
  "RCT only",
  "RCT only; no NRS included"
)

results <- add_result(
  results,
  nrs_only,
  nrs_only$vi_original + primary_bias_sd^2,
  "1.Primary analysis(NRS bias SD = 0.20)",
  "NRS only: bias SD = 0.20",
  nrs_bias_note("NRS only")
)

results <- add_nrs_bias_result(
  results,
  primary_bias_sd,
  "1.Primary analysis(NRS bias SD = 0.20)",
  "Overall: bias SD = 0.20"
)

results <- add_result(
  results,
  nrs_only,
  nrs_only$vi_original,
  "No down-weighting reference",
  "NRS only",
  "NRS only; original variances"
)

results <- add_result(
  results,
  es,
  es$vi_original,
  "No down-weighting reference",
  "Total",
  "RCT + NRS; original variances"
)

results <- add_nrs_bias_result(
  results,
  0.30,
  "NRS Bias SD sensitivity",
  "bias SD= 0.30"
)

results <- add_nrs_bias_result(
  results,
  0.40,
  "NRS Bias SD sensitivity",
  "bias SD= 0.40"
)

# Power-prior sensitivity.
results <- add_power_prior_result(results, 0.25, "alpha = 0.25")
results <- add_power_prior_result(results, 0.1, "alpha = 0.1")

results <- add_nrs_bias_dataset_result(
  results,
  es %>% filter(!is_conference_nrs),
  primary_bias_sd,
  "Exclusion analyses",
  "Exclude Conference-abstract NRS",
  nrs_bias_note("Exclude NRS conference abstracts", retained = TRUE)
)

exclude_co_intervention_es <- es %>% filter(intervention_package_type != "Co-intervention")
if (nrow(exclude_co_intervention_es) < nrow(es)) {
  results <- add_nrs_bias_dataset_result(
    results,
    exclude_co_intervention_es,
    primary_bias_sd,
    "Exclusion analyses",
    "Exclude Co-intervention",
    nrs_bias_note("Exclude studies with Co-intervention intervention package type", retained = TRUE)
  )
}

inactive_comparator_es <- es %>%
  filter(comparator_condition %in% c(
    "Usual care without active electronic alerting",
    "Usual care with silent electronic monitoring"
  ))
assert_expected_counts(inactive_comparator_es, "Inactive comparator only", 6, 42, 48)
results <- add_nrs_bias_dataset_result(
  results,
  inactive_comparator_es,
  primary_bias_sd,
  "Exclusion analyses",
  "Inactive comparator only",
  nrs_bias_note("Include only inactive comparator studies", retained = TRUE)
)

# This retains A or U-section-sign studies; F studies are not part of this subset.
adjusted_unused_es <- es %>% filter(adjusted_estimate_not_used)
assert_expected_counts(adjusted_unused_es, "Adjusted estimate not used in primary synthesis", 1, 21, 22)
if (any(adjusted_unused_es$uses_adjusted_rr) ||
    !isTRUE(all.equal(as.numeric(adjusted_unused_es$yi), as.numeric(adjusted_unused_es$yi_event))) ||
    !isTRUE(all.equal(as.numeric(adjusted_unused_es$vi_original), as.numeric(adjusted_unused_es$vi)))) {
  stop("The adjusted-unused sensitivity subset must retain primary group-level RR estimates and sampling variances.")
}
results <- add_nrs_bias_dataset_result(
  results,
  adjusted_unused_es,
  primary_bias_sd,
  "Effect-estimate availability",
  "Adjusted estimate not used in primary synthesis",
  nrs_bias_note("Retain studies with A or U-section-sign mortality evidence: adjusted/design-adjusted evidence was reported but not used as the primary quantitative effect; use the primary group-level RR", retained = TRUE)
)
write.csv(
  as.data.frame(adjusted_unused_es[, c(
    "study", "design", "events_e", "n_e", "events_c", "n_c",
    "adjusted_estimate_not_used", "uses_adjusted_rr", "yi", "vi_original", "vi_bias_adjusted"
  )]),
  file.path(out_dir, "mortality_adjusted_unused_effect_sizes.csv"),
  row.names = FALSE
)

no_additional_restriction_es <- es %>% filter(additional_cohort_restriction == "no")
assert_expected_counts(no_additional_restriction_es, "No additional cohort restriction", 10, 49, 59)
results <- add_nrs_bias_dataset_result(
  results,
  no_additional_restriction_es,
  primary_bias_sd,
  "Population sensitivity",
  "No additional cohort restriction",
  nrs_bias_note("Include only studies with no additional cohort restriction", retained = TRUE)
)

exclude_all_screened_es <- es %>%
  filter(!is_all_screened(clinical_population_scope))
assert_expected_counts(exclude_all_screened_es, "Exclude all-screened", 7, 53, 60)
if (nrow(exclude_all_screened_es) < nrow(es)) {
  results <- add_nrs_bias_dataset_result(
    results,
    exclude_all_screened_es,
    primary_bias_sd,
    "Population sensitivity",
    "Exclude all-screened population",
    nrs_bias_note("Exclude all-screened/all-population studies", retained = TRUE)
  )
}

restricted_severe_or_sepsis3_es <- es %>%
  filter(is_organ_dysfunction_targeted(clinical_population_scope))
assert_expected_counts(restricted_severe_or_sepsis3_es, "Severe sepsis/Sepsis-3", 4, 30, 34)
results <- add_nrs_bias_dataset_result(
  results,
  restricted_severe_or_sepsis3_es,
  primary_bias_sd,
  "Population sensitivity",
  "Restricted severe sepsis / Sepsis-3 population",
  nrs_bias_note("Include only severe sepsis or Sepsis-3 patients")
)

restricted_sepsis3_es <- es %>%
  filter(is_sepsis3_targeted(clinical_population_scope))
assert_expected_counts(restricted_sepsis3_es, "Sepsis-3 only", 3, 12, 15)
results <- add_nrs_bias_dataset_result(
  results,
  restricted_sepsis3_es,
  primary_bias_sd,
  "Population sensitivity",
  "Restricted Sepsis-3 patients",
  nrs_bias_note("Include only Sepsis-3 patients")
)

summary_table <- bind_rows(results)
summary_table <- summary_table %>%
  mutate(
    display_label = case_when(
      scenario == "RCT only" ~ "RCT only",
      scenario == "NRS only" ~ "NRS only",
      scenario == "Total" ~ "Overall",
      scenario == "NRS only: bias SD = 0.20" ~ "NRS only",
      scenario == "Overall: bias SD = 0.20" ~ "Overall",
      scenario == "bias SD= 0.30" ~ "bias SD= 0.30",
      scenario == "bias SD= 0.40" ~ "bias SD= 0.40",
      scenario == "alpha = 0.25" ~ "alpha = 0.25",
      scenario == "alpha = 0.1" ~ "alpha = 0.10",
      scenario == "Exclude Conference-abstract NRS" ~ "Exclude Conference NRS",
      scenario == "Exclude Co-intervention" ~ "Exclude Co-intervention",
      scenario == "Inactive comparator only" ~ "Inactive comparator only",
      scenario == "No additional cohort restriction" ~ "No additional cohort restriction",
      scenario == "Exclude all-screened population" ~ "Exclude all-screened",
      scenario == "Restricted severe sepsis / Sepsis-3 population" ~ "Severe sepsis/Sepsis-3",
      scenario == "Restricted Sepsis-3 patients" ~ "Sepsis-3 only",
      TRUE ~ scenario
    )
  )

write.csv(
  summary_table,
  file.path(out_dir, "mortality_sensitivity_final_results_summary.csv"),
  row.names = FALSE
)
write.csv(
  bind_rows(reml_reference_results),
  file.path(out_dir, "mortality_sensitivity_final_results_reml_reference.csv"),
  row.names = FALSE
)
fit_inventory <- bind_rows(sensitivity_fit_inventory)
write.csv(
  fit_inventory,
  file.path(out_dir, "mortality_sensitivity_fit_inventory.csv"),
  row.names = FALSE
)
expected_local_fits <- basename(fit_inventory$fit_file[fit_inventory$fit_source == "mortality_sensitivity_cache"])
stale_local_fits <- setdiff(list.files(fits_dir, pattern = "\\.brmsfit\\.rds$"), expected_local_fits)
if (length(stale_local_fits) > 0L) unlink(file.path(fits_dir, stale_local_fits))
message("Saved mortality sensitivity Bayesian cache to: ", out_dir)
