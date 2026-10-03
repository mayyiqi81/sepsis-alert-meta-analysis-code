# Bayesian sensitivity analyses for mortality-process meta-regression.

library(metafor)
library(dplyr)

script_arg <- c(grep("^--file=", commandArgs(FALSE), value = TRUE), grep("\\.R$", commandArgs(FALSE), value = TRUE, ignore.case = TRUE))
script_dir <- if (length(script_arg) > 0) dirname(normalizePath(sub("^--file=", "", script_arg[1]), mustWork = FALSE)) else getwd()
project_root <- normalizePath(file.path(script_dir, "..", ".."), winslash = "/", mustWork = TRUE)
source(file.path(project_root, "analysis_scripts", "common", "project_paths.R"))
options(sepsis_project_root = set_project_root_from_script(project_root))
script_dir <- project_root
workspace_dir <- script_dir
source(file.path(workspace_dir, "analysis_scripts", "common", "bayes_prior_scale_spec.R"))
source(file.path(workspace_dir, "analysis_scripts", "common", "antibiotic_time_direction_contract.R"))
source(file.path(workspace_dir, "analysis_scripts", "common", "bayes_mortality_population_helpers.R"))
source(file.path(workspace_dir, "analysis_scripts", "common", "bayes_mortality_functions.R"))
source(file.path(workspace_dir, "analysis_scripts", "common", "mortality_rd_functions.R"))

input_file <- file.path(workspace_dir, "source_csv", "mortality.csv")
rd_inputs_file <- file.path(workspace_dir, "source_csv", "mortality_rd_inputs.csv")
nrs_bias_sd <- 0.20
bayes_settings <- bayes_meta_settings_from_env()
prior_variant <- Sys.getenv("BAYES_PRIOR_VARIANT", unset = "canonical")
out_dir <- Sys.getenv(
  "MORTALITY_PROCESS_METAREG_DATA_DIR",
  unset = "computed_data/mortality_process_metareg"
)
if (!grepl("^([A-Za-z]:[\\\\/]|/)", out_dir)) out_dir <- file.path(workspace_dir, out_dir)
canonical_out_dir <- Sys.getenv(
  "BAYES_CANONICAL_OUT_DIR",
  unset = "computed_data/mortality_main/canonical"
)
if (!grepl("^([A-Za-z]:[\\\\/]|/)", canonical_out_dir)) canonical_out_dir <- file.path(workspace_dir, canonical_out_dir)

fits_dir <- file.path(out_dir, "fits")
output_prefix <- "mortality_process_meta_regression_sensitivity"
legacy_prefix <- "mortality_univariable_process_meta_regression_nrs_bias_sd"

MIN_POINT_SIZE <- 0.35
POINT_WEIGHT_SCALE <- 2

if (!dir.exists(out_dir)) {
  dir.create(out_dir, recursive = TRUE)
}
if (!dir.exists(fits_dir)) {
  dir.create(fits_dir, recursive = TRUE)
}

global_palette <- list(
  blue = c(light = "#DCECF9", main = "#AED7F3", deep = "#5A9CCB"),
  pink = c(light = "#F7D5DC", main = "#E8A3B1", deep = "#D76A84")
)

colors <- list(
  rct = global_palette$blue["deep"],
  nrs = global_palette$pink["deep"]
)

covariates <- data.frame(
  covariate = c(
    "lactate_measurement",
    "blood_culture",
    "antibiotic_time",
    "antibiotic_proportion",
    "fluid_resuscitation",
    "sepsis_bundle"
  ),
  column = c(
    "Lactate measurement RD, percentage points",
    "Blood culture RD, percentage points",
    "Antibiotic time MD, hours",
    "Antibiotic proportion RD, percentage points",
    "Fluid resuscitation RD, percentage points",
    "Sepsis bundle RD, percentage points"
  ),
  label = c(
    "Lactate measurement",
    "Blood culture",
    "Antibiotic time",
    "Antibiotic proportion",
    "Fluid resuscitation",
    "Sepsis bundle"
  ),
  x_unit = c(
    "Process-compliance RD, percentage points",
    "Process-compliance RD, percentage points",
    "Change in antibiotic time, hours",
    "Process-compliance RD, percentage points",
    "Process-compliance RD, percentage points",
    "Process-compliance RD, percentage points"
  ),
  strict_flag = c(
    "Strict_same_cohort_regression__Lactate_measurement",
    "Strict_same_cohort_regression__Blood_culture",
    "Strict_same_cohort_regression__Antibiotic_time",
    "Strict_same_cohort_regression__Antibiotic_proportion",
    "Strict_same_cohort_regression__Fluid_resuscitation",
    "Strict_same_cohort_regression__Sepsis_bundle"
  ),
  stringsAsFactors = FALSE
)

baseline_covariate <- data.frame(
  covariate = c(
    "baseline_mortality_rate_percent_all",
    "baseline_mortality_rate_percent_severe"
  ),
  column = c(
    "Baseline mortality rate in control group",
    "Baseline mortality rate in control group"
  ),
  label = c(
    "All studies",
    "Severe sepsis/Sepsis-3"
  ),
  x_unit = c(
    "Baseline mortality (%)",
    "Baseline mortality (%)"
  ),
  strict_flag = c(NA_character_, NA_character_),
  stringsAsFactors = FALSE
)

all_covariates <- bind_rows(covariates, baseline_covariate)

analysis_covariates <- function(spec) {
  if (identical(spec$analysis, "Bias-unadjusted log RR")) {
    bind_rows(covariates, baseline_covariate)
  } else {
    covariates
  }
}

analysis_specs <- list(
  list(
    analysis = "Bias-unadjusted log RR",
    analysis_id = "bias_unadjusted_logrr",
    response_yi = "yi_logrr",
    response_vi = "vi_original",
    strict = FALSE,
    y_mode = "rr",
    include_baseline = TRUE,
    subtitle = "Sensitivity analysis: without NRS bias adjustment",
    caption = "NRS variance was not inflated by the prespecified bias SD.",
    pdf = "1.5A.mortality_process_meta_regression_bias_unadjusted_logrr.pdf",
    note = "Mortality log RR retained the 1.6 point-estimate selection; NRS variance was not inflated by the prespecified bias SD."
  ),
  list(
    analysis = "Mortality RD",
    analysis_id = "mortality_rd",
    response_yi = "yi_rd",
    response_vi = "vi_rd",
    strict = FALSE,
    y_mode = "rd",
    include_baseline = FALSE,
    subtitle = "Sensitivity analysis: mortality risk-difference scale",
    caption = "Mortality RD used event-count estimates except for the adjusted absolute-effect inputs specified for Adams 2022 and Arabi 2025.",
    pdf = "1.5B.mortality_process_meta_regression_mortality_rd.pdf",
    note = "Mortality RD used event-count estimates except for Adams 2022 (reported adjusted RD) and Arabi 2025 (adjusted RD derived from aRR and control risk, with delta-method variance incorporating adjusted-log-RR and control-risk sampling uncertainty under an independence assumption); no NRS bias-SD inflation was applied on the RD scale."
  ),
  list(
    analysis = "Strict same-cohort",
    analysis_id = "strict_same_cohort",
    response_yi = "yi_logrr",
    response_vi = "vi_adjusted",
    strict = TRUE,
    y_mode = "rr",
    include_baseline = FALSE,
    subtitle = "Sensitivity analysis: strict mortality-process cohort alignment",
    caption = "Only studies with outcome-specific strict mortality-process cohort alignment were included.",
    pdf = "1.5C.mortality_process_meta_regression_strict_same_cohort.pdf",
    note = paste0("Mortality log RR retained the 1.6 point-estimate selection; NRS variance inflated by sigma_bias^2 = ", nrs_bias_sd^2, "; outcome-specific strict same-cohort flag required.")
  )
)

fmt_p <- function(p) {
  ifelse(is.na(p), NA_character_, ifelse(p < 0.001, "<0.001", sprintf("%.3f", p)))
}

clean_value <- function(x) {
  clean_mortality_value(x)
}

safe_num <- function(x) {
  safe_mortality_num(x)
}

strict_flag_column <- function(covariate_name) {
  paste0("strict_", covariate_name)
}

# The current source file uses the `*-targeted` labels, while the historical
# Figure 6 cache used equivalent parenthetical labels. Both encodings map to
# the same 34-study organ-dysfunction set; this compatibility pattern keeps
# cache regeneration aligned with the established membership.
read_mortality_data <- function() {
  if (!file.exists(input_file)) {
    stop("Missing input file: ", input_file)
  }

  df_raw <- read_mortality_csv(input_file, lower_names = FALSE)
assert_antibiotic_time_source_direction(df_raw)
  required_cols <- c(
    "Study", "Design", "events_e", "n_e", "events_c", "n_c",
    "Baseline mortality rate in control group",
    "Clinical population scope",
    covariates$column,
    covariates$strict_flag
  )
  missing_cols <- setdiff(required_cols, names(df_raw))
  if (length(missing_cols) > 0) {
    stop("Missing required columns: ", paste(missing_cols, collapse = ", "))
  }
  for (col in c("aRR", "aRR_lower_95", "aRR_upper_95")) {
    if (!col %in% names(df_raw)) {
      df_raw[[col]] <- NA_real_
    }
  }

  df <- df_raw %>%
    transmute(
      study = clean_value(Study),
      design = toupper(clean_value(Design)),
      events_e = safe_num(events_e),
      n_e = safe_num(n_e),
      events_c = safe_num(events_c),
      n_c = safe_num(n_c),
      aRR = safe_num(aRR),
      aRR_lower_95 = safe_num(aRR_lower_95),
      aRR_upper_95 = safe_num(aRR_upper_95),
      clinical_population_scope = clean_value(.data[["Clinical population scope"]]),
      # Keep the S6A baseline moderator on the 0-1 proportion scale used by
      # Figure 6. The plotting layer converts only displayed x coordinates to
      # percentages; the model and saved prediction scale stay unchanged.
      baseline_mortality_rate_percent_all = safe_num(.data[["Baseline mortality rate in control group"]]),
      baseline_mortality_rate_percent_severe = safe_num(.data[["Baseline mortality rate in control group"]]),
      blood_culture = safe_num(.data[["Blood culture RD, percentage points"]]),
      lactate_measurement = safe_num(.data[["Lactate measurement RD, percentage points"]]),
      antibiotic_proportion = safe_num(.data[["Antibiotic proportion RD, percentage points"]]),
      fluid_resuscitation = safe_num(.data[["Fluid resuscitation RD, percentage points"]]),
      sepsis_bundle = safe_num(.data[["Sepsis bundle RD, percentage points"]]),
      antibiotic_time = safe_num(.data[["Antibiotic time MD, hours"]]),
      strict_blood_culture = clean_value(.data[["Strict_same_cohort_regression__Blood_culture"]]),
      strict_lactate_measurement = clean_value(.data[["Strict_same_cohort_regression__Lactate_measurement"]]),
      strict_antibiotic_proportion = clean_value(.data[["Strict_same_cohort_regression__Antibiotic_proportion"]]),
      strict_fluid_resuscitation = clean_value(.data[["Strict_same_cohort_regression__Fluid_resuscitation"]]),
      strict_sepsis_bundle = clean_value(.data[["Strict_same_cohort_regression__Sepsis_bundle"]]),
      strict_antibiotic_time = clean_value(.data[["Strict_same_cohort_regression__Antibiotic_time"]])
    ) %>%
    filter(
      !is.na(study), study != "",
      design %in% c("RCT", "NRS"),
      !is.na(events_e), !is.na(n_e),
      !is.na(events_c), !is.na(n_c),
      !is.na(baseline_mortality_rate_percent_all),
      n_e > 0, n_c > 0
    )

  validate_mortality_2x2(df, context = "fit_process_metareg mortality data")
  validate_adjusted_rr_fields(df, context = "fit_process_metareg mortality data")

  if (nrow(df) != 63) {
    stop("Expected 63 mortality rows after filtering, found ", nrow(df), ".")
  }

  df
}

construct_effects <- function(df) {
  es_rr <- metafor::escalc(
    measure = "RR",
    ai = events_e,
    bi = n_e - events_e,
    ci = events_c,
    di = n_c - events_c,
    data = df,
    add = 0.5,
    to = "all"
  )

  out <- es_rr %>%
    mutate(
      yi_crude_logrr = yi,
      vi_crude_logrr = vi,
      uses_adjusted_rr = !is.na(aRR) & !is.na(aRR_lower_95) & !is.na(aRR_upper_95),
      p_c = events_c / n_c,
      yi_logrr = ifelse(uses_adjusted_rr, log(aRR), yi_crude_logrr),
      vi_original = ifelse(
        uses_adjusted_rr,
        ((log(aRR_upper_95) - log(aRR_lower_95)) / (2 * 1.96))^2,
        vi_crude_logrr
      ),
      vi_adjusted = ifelse(design == "NRS", vi_original + nrs_bias_sd^2, vi_original),
      vi_bias_adjusted = vi_adjusted,
      se_original = sqrt(vi_original),
      se_adjusted = sqrt(vi_adjusted),
      yi = yi_logrr,
      RR = exp(yi_logrr),
      lower_RR = exp(yi_logrr - 1.96 * se_adjusted),
      upper_RR = exp(yi_logrr + 1.96 * se_adjusted)
    )
  out <- add_mortality_rd_effects(out, rd_inputs_file)
  out$yi <- out$yi_logrr
  out
}

make_dataset <- function(es_all, spec, covariate_name) {
  dat <- es_all %>%
    mutate(covariate_value = .data[[covariate_name]]) %>%
    filter(!is.na(covariate_value), !is.na(.data[[spec$response_yi]]), !is.na(.data[[spec$response_vi]]))

  if (isTRUE(spec$strict)) {
    flag_col <- strict_flag_column(covariate_name)
    dat <- dat %>%
      filter(tolower(trimws(.data[[flag_col]])) == "yes")
  }

  if (identical(covariate_name, "baseline_mortality_rate_percent_severe")) {
    dat <- dat %>%
      filter(is_organ_dysfunction_targeted(clinical_population_scope))
  }

  dat
}

skip_result <- function(spec, covariate_name, dat, reason) {
  # Display-only conversions: RD response values are written as percentage
  # points and log-RR values as RRs; model fitting remains on the analytic
  # response scale.
  dat <- dat %>%
    mutate(
      random_weight = NA_real_,
      y_plot = if (spec$y_mode == "rd") .data[[spec$response_yi]] * 100 else exp(.data[[spec$response_yi]]),
      point_size = NA_real_,
      shape_value = ifelse(design == "RCT", 17, 16),
      point_color = ifelse(design == "RCT", colors$rct, colors$nrs)
    )
  cov_row <- all_covariates[all_covariates$covariate == covariate_name, ]
  list(
    summary = data.frame(
      analysis = spec$analysis,
      analysis_id = spec$analysis_id,
      covariate = covariate_name,
      label = cov_row$label,
      column = cov_row$column,
      x_unit = cov_row$x_unit,
      k = nrow(dat),
      n_rct = sum(dat$design == "RCT"),
      n_nrs = sum(dat$design == "NRS"),
      min_value = ifelse(nrow(dat) == 0, NA_real_, min(dat$covariate_value, na.rm = TRUE)),
      max_value = ifelse(nrow(dat) == 0, NA_real_, max(dat$covariate_value, na.rm = TRUE)),
      p_value = NA_real_,
      p_value_formatted = NA_character_,
      prob_direction = NA_real_,
      beta = NA_real_,
      beta_lower = NA_real_,
      beta_upper = NA_real_,
      tau = NA_real_,
      tau2 = NA_real_,
      rhat_max = NA_real_,
      ess_bulk_min = NA_real_,
      ess_tail_min = NA_real_,
      divergent_transitions = NA_real_,
      QM = NA_real_,
      df = NA_real_,
      I2_percent = NA_real_,
      AIC = NA_real_,
      BIC = NA_real_,
      method = "Skipped",
      note = reason,
      y_mode = spec$y_mode,
      response_yi = spec$response_yi,
      response_vi = spec$response_vi,
      fit_file = NA_character_,
      stringsAsFactors = FALSE
    ),
    predictions = data.frame(),
    data = dat
  )
}

moderator_id_for_covariate <- function(covariate_name) {
  if (covariate_name %in% baseline_covariate$covariate) {
    "baseline_mortality"
  } else if (identical(covariate_name, "antibiotic_time")) {
    "antibiotic_time"
  } else {
    "process_completion"
  }
}

fit_meta_regression <- function(dat, spec, covariate_name) {
  model_dat <- dat %>%
    mutate(
      yi = .data[[spec$response_yi]],
      vi_model = .data[[spec$response_vi]]
    )
  moderator_id <- moderator_id_for_covariate(covariate_name)
  response_scale <- if (identical(spec$y_mode, "rd")) "rd" else "log_rr"
  prior_spec_id <- if (identical(response_scale, "rd")) {
    "rd_metareg_canonical"
  } else {
    "log_rr_metareg_canonical"
  }
  uncertainty_spec_id <- switch(
    spec$analysis_id,
    bias_unadjusted_logrr = "bias_unadjusted_no_nrs_additional_uncertainty",
    mortality_rd = "mortality_rd_no_nrs_additional_uncertainty",
    strict_same_cohort = "strict_same_cohort_primary_nrs_bias_sd_0.20",
    stop("No uncertainty specification registered for analysis_id: ", spec$analysis_id)
  )
  fit_path <- file.path(fits_dir, paste0(output_prefix, "_", spec$analysis_id, "_", covariate_name, ".brmsfit.rds"))
  force_refit <- identical(Sys.getenv("BAYES_FORCE_REFIT", unset = ""), "1")
  expected_metadata <- make_bayes_metareg_cache_metadata(
    data = model_dat,
    moderator_col = "covariate_value",
    vi_col = "vi_model",
    moderator_id = moderator_id,
    response_scale = response_scale,
    prior_spec_id = prior_spec_id,
    uncertainty_spec_id = uncertainty_spec_id,
    settings = bayes_settings,
    prior_variant = prior_variant
  )
  reuse_cached_fit <- file.exists(fit_path) && !force_refit
  if (reuse_cached_fit) {
    cached_fit <- readRDS(fit_path)
    reuse_cached_fit <- isTRUE(assert_cached_metareg_fit_matches(cached_fit, expected_metadata))
    if (!reuse_cached_fit) {
      upgraded_fit <- upgrade_legacy_cached_metareg_fit(cached_fit, expected_metadata)
      if (!is.null(upgraded_fit)) {
        cached_fit <- upgraded_fit
        saveRDS(cached_fit, fit_path)
        reuse_cached_fit <- TRUE
      }
    }
  }
  if (reuse_cached_fit) {
    fit <- if (exists("cached_fit", inherits = FALSE)) cached_fit else readRDS(fit_path)
  } else {
    fit <- fit_bayes_mortality_metareg_with_quality_gate(
      model_dat,
      moderator_col = "covariate_value",
      moderator_id = moderator_id,
      response_scale = response_scale,
      prior_spec_id = prior_spec_id,
      uncertainty_spec_id = uncertainty_spec_id,
      vi_col = "vi_model",
      settings = bayes_settings,
      prior_variant = prior_variant
    )
    saveRDS(fit, fit_path)
  }
  list(fit = fit, fit_path = fit_path)
}

run_univariable_sensitivity <- function(es_all, spec, covariate_name) {
  dat <- make_dataset(es_all, spec, covariate_name)
  moderator_id <- moderator_id_for_covariate(covariate_name)

  if (nrow(dat) < 3) {
    return(skip_result(spec, covariate_name, dat, "Skipped: fewer than three studies"))
  }
  if (length(unique(dat$covariate_value)) < 2) {
    return(skip_result(spec, covariate_name, dat, "Skipped: fewer than two unique covariate values"))
  }

  fit_result <- fit_meta_regression(dat, spec, covariate_name)
  fit <- fit_result$fit
  fit_summary <- summarise_bayes_mortality_metareg(fit, covariate_name)
  response_scale <- if (identical(spec$y_mode, "rd")) "rd" else "log_rr"
  prior_spec_id <- if (identical(response_scale, "rd")) {
    "rd_metareg_canonical"
  } else {
    "log_rr_metareg_canonical"
  }

  dat <- dat %>%
    mutate(
      random_weight = 100 * (1 / (.data[[spec$response_vi]] + fit_summary$tau2)) /
        sum(1 / (.data[[spec$response_vi]] + fit_summary$tau2), na.rm = TRUE),
      y_plot = if (spec$y_mode == "rd") .data[[spec$response_yi]] * 100 else exp(.data[[spec$response_yi]]),
      point_size = pmax(MIN_POINT_SIZE, POINT_WEIGHT_SCALE * sqrt(random_weight)),
      shape_value = ifelse(design == "RCT", 17, 16),
      point_color = ifelse(design == "RCT", colors$rct, colors$nrs)
    )

  x_grid_raw <- seq(min(dat$covariate_value), max(dat$covariate_value), length.out = 100)
  pred_df <- predict_bayes_mortality_metareg_raw(fit, x_grid_raw, moderator_id) %>%
    mutate(
      x = x_grid_raw,
      analysis = spec$analysis,
      analysis_id = spec$analysis_id,
      covariate = covariate_name,
      y = if (spec$y_mode == "rd") pred * 100 else RR,
      y_lower = if (spec$y_mode == "rd") ci_lb * 100 else lower_RR,
      y_upper = if (spec$y_mode == "rd") ci_ub * 100 else upper_RR,
      RR = if (spec$y_mode == "rd") NA_real_ else RR,
      lower_RR = if (spec$y_mode == "rd") NA_real_ else lower_RR,
      upper_RR = if (spec$y_mode == "rd") NA_real_ else upper_RR,
      .before = 1
    )

  cov_row <- all_covariates[all_covariates$covariate == covariate_name, ]
  list(
    summary = data.frame(
      analysis = spec$analysis,
      analysis_id = spec$analysis_id,
      covariate = covariate_name,
      label = cov_row$label,
      column = cov_row$column,
      x_unit = cov_row$x_unit,
      k = nrow(dat),
      n_rct = sum(dat$design == "RCT"),
      n_nrs = sum(dat$design == "NRS"),
      min_value = min(dat$covariate_value, na.rm = TRUE),
      max_value = max(dat$covariate_value, na.rm = TRUE),
      p_value = fit_summary$p_value,
      p_value_formatted = fmt_p(fit_summary$p_value),
      prob_direction = fit_summary$prob_direction,
       beta = fit_summary$beta_model,
       beta_lower = fit_summary$beta_lower_model,
       beta_upper = fit_summary$beta_upper_model,
       beta_model = fit_summary$beta_model,
       beta_lower_model = fit_summary$beta_lower_model,
       beta_upper_model = fit_summary$beta_upper_model,
       beta_reported = fit_summary$beta_reported,
       beta_lower_reported = fit_summary$beta_lower_reported,
       beta_upper_reported = fit_summary$beta_upper_reported,
       reported_increment_label = fit_summary$reported_increment_label,
       moderator_id = fit_summary$moderator_id,
       moderator_centre_raw = fit_summary$moderator_centre_raw,
       moderator_increment_raw = fit_summary$moderator_increment_raw,
       prior_spec_id = prior_spec_id,
       prior_variant = prior_variant,
       tau = fit_summary$tau,
      tau2 = fit_summary$tau2,
      rhat_max = fit_summary$rhat_max,
      ess_bulk_min = fit_summary$ess_bulk_min,
      ess_tail_min = fit_summary$ess_tail_min,
      divergent_transitions = fit_summary$divergent_transitions,
      QM = NA_real_,
      df = NA_real_,
      I2_percent = NA_real_,
      AIC = NA_real_,
      BIC = NA_real_,
      method = fit_summary$method,
      note = spec$note,
      y_mode = spec$y_mode,
      response_yi = spec$response_yi,
      response_vi = spec$response_vi,
      fit_file = fit_result$fit_path,
      stringsAsFactors = FALSE
    ),
    predictions = pred_df,
    data = dat
  )
}

run_analysis <- function(es_all, spec) {
  plot_covariates <- analysis_covariates(spec)
  results <- lapply(plot_covariates$covariate, function(x) run_univariable_sensitivity(es_all, spec, x))

  summary_table <- bind_rows(lapply(results, `[[`, "summary")) %>%
    arrange(match(covariate, plot_covariates$covariate))

  predictions <- bind_rows(lapply(results, `[[`, "predictions"))
  if (nrow(predictions) > 0) {
    predictions <- predictions %>%
      left_join(all_covariates[, c("covariate", "label", "column", "x_unit")], by = "covariate")
  }

  plot_data <- bind_rows(lapply(seq_along(results), function(i) {
    covariate_name <- plot_covariates$covariate[i]
    results[[i]]$data %>%
      mutate(
        analysis = spec$analysis,
        analysis_id = spec$analysis_id,
        covariate = covariate_name,
        covariate_value = .data[[covariate_name]]
      ) %>%
      select(
        analysis, analysis_id, covariate, study, design, covariate_value,
        events_e, n_e, events_c, n_c,
        yi_logrr, vi_original, vi_adjusted, vi_bias_adjusted,
        yi_rd, vi_rd, yi_rd_event_counts, vi_rd_event_counts, rd_lower_95, rd_upper_95, rd_source, uses_adjusted_rr,
        random_weight, point_size, shape_value, point_color, y_plot
      )
  })) %>%
    left_join(all_covariates[, c("covariate", "label", "column", "x_unit")], by = "covariate")

  list(
    spec = spec,
    summary = summary_table,
    predictions = predictions,
    plot_data = plot_data
  )
}

run_qa_checks <- function(es_all, analysis_results) {
  process_missing <- setdiff(covariates$covariate, names(es_all))
  if (length(process_missing) > 0) {
    stop("Process columns missing after import: ", paste(process_missing, collapse = ", "))
  }

  strict_missing <- setdiff(paste0("strict_", covariates$covariate), names(es_all))
  if (length(strict_missing) > 0) {
    stop("Strict same-cohort flag columns missing after import: ", paste(strict_missing, collapse = ", "))
  }

  if (any(abs(es_all$vi_adjusted[es_all$design == "NRS"] - (es_all$vi_original[es_all$design == "NRS"] + nrs_bias_sd^2)) > 1e-12)) {
    stop("NRS log-RR bias variance inflation check failed.")
  }
  if (any(abs(es_all$vi_adjusted[es_all$design == "RCT"] - es_all$vi_original[es_all$design == "RCT"]) > 1e-12)) {
    stop("RCT log-RR variance should not be bias-inflated.")
  }

  rd_check <- (es_all$events_e / es_all$n_e) - (es_all$events_c / es_all$n_c)
  event_count_rows <- es_all$rd_source == "event_counts"
  if (max(abs(es_all$yi_rd[event_count_rows] - rd_check[event_count_rows]), na.rm = TRUE) > 1e-12) {
    stop("Mortality event-count RD check failed.")
  }
  adams_row <- es_all[es_all$study == "Adams 2022", , drop = FALSE]
  arabi_row <- es_all[es_all$study == "Arabi 2025", , drop = FALSE]
  rd_inputs <- read_mortality_rd_inputs(rd_inputs_file)
  adams_input <- rd_inputs[rd_inputs$Study == "Adams 2022", , drop = FALSE]
  if (nrow(adams_row) != 1 || nrow(adams_input) != 1 || adams_row$rd_source != "adjusted_rd_reported" ||
      abs(adams_row$yi_rd - adams_input$rd_percentage_points / 100) > 1e-12) {
    stop("Adams 2022 adjusted RD override check failed.")
  }
  arabi_expected_rd <- arabi_row$p_c * (arabi_row$aRR - 1)
  if (nrow(arabi_row) != 1 || arabi_row$rd_source != "adjusted_rd_from_arr" || abs(arabi_row$yi_rd - arabi_expected_rd) > 1e-12) {
    stop("Arabi 2025 aRR-derived RD override check failed.")
  }

  main_sets <- lapply(covariates$covariate, function(covariate_name) {
    es_all %>%
      filter(!is.na(.data[[covariate_name]])) %>%
      arrange(study) %>%
      pull(study)
  })
  names(main_sets) <- covariates$covariate

  sens1 <- analysis_results[["Bias-unadjusted log RR"]]
  for (covariate_name in covariates$covariate) {
    sens1_set <- sens1$plot_data %>%
      filter(covariate == covariate_name) %>%
      arrange(study) %>%
      pull(study)
    if (!identical(sens1_set, main_sets[[covariate_name]])) {
      stop("Sensitivity 1 study set differs from 1.6 complete-case set for ", covariate_name)
    }

    sens1_yi <- sens1$plot_data %>%
      filter(covariate == covariate_name) %>%
      arrange(study) %>%
      select(study, yi_logrr)
    main_yi <- es_all %>%
      filter(study %in% sens1_set) %>%
      arrange(study) %>%
      select(study, yi_logrr)
    if (!isTRUE(all.equal(sens1_yi, main_yi, tolerance = 1e-12, check.attributes = FALSE))) {
      stop("Sensitivity 1 mortality point estimates differ from 1.6 logic for ", covariate_name)
    }
  }

  strict_res <- analysis_results[["Strict same-cohort"]]
  for (covariate_name in covariates$covariate) {
    flag_col <- strict_flag_column(covariate_name)
    strict_set <- strict_res$plot_data %>%
      filter(covariate == covariate_name) %>%
      arrange(study) %>%
      pull(study)
    if (!all(strict_set %in% main_sets[[covariate_name]])) {
      stop("Strict same-cohort set includes rows outside the 1.6 complete-case set for ", covariate_name)
    }
    flag_values <- es_all %>%
      filter(study %in% strict_set) %>%
      pull(.data[[flag_col]])
    if (length(flag_values) > 0 && any(tolower(trimws(flag_values)) != "yes")) {
      stop("Strict same-cohort set includes non-Yes flag values for ", covariate_name)
    }
  }

  invisible(TRUE)
}

df <- read_mortality_data()
es_all <- construct_effects(df)

canonical_bayes <- load_canonical_bayes_mortality(
  canonical_out_dir,
  require_fits = FALSE,
  input_file = input_file,
  expected_settings = bayes_settings,
  expected_nrs_bias_sd = nrs_bias_sd
)
assert_canonical_effect_sizes_match(
  current_es = es_all %>% select(study, design, yi, vi_original, vi_bias_adjusted),
  canonical_es = canonical_bayes$effect_sizes
)
canonical_overall <- canonical_summary_row(canonical_bayes, "overall")
main_tau2_for_weights <- as.numeric(canonical_overall$tau2)
if (!is.finite(main_tau2_for_weights)) {
  stop("Canonical Bayesian overall tau2 is missing or non-finite: ", canonical_out_dir)
}
es_all <- es_all %>%
  mutate(
    main_random_weight = 100 * (1 / (vi_adjusted + main_tau2_for_weights)) /
      sum(1 / (vi_adjusted + main_tau2_for_weights), na.rm = TRUE)
  )

analysis_results <- lapply(analysis_specs, function(spec) run_analysis(es_all, spec))
names(analysis_results) <- vapply(analysis_specs, `[[`, character(1), "analysis")

run_qa_checks(es_all, analysis_results)

summary_table <- bind_rows(lapply(analysis_results, `[[`, "summary")) %>%
  select(
    analysis, analysis_id, covariate, label, k, n_rct, n_nrs,
    beta, beta_lower, beta_upper, p_value, p_value_formatted, prob_direction,
    tau, tau2, rhat_max, ess_bulk_min, ess_tail_min, divergent_transitions,
    QM, df, I2_percent, AIC, BIC, method, note, everything()
  )

prediction_table <- bind_rows(lapply(analysis_results, `[[`, "predictions"))
analysis_data <- bind_rows(lapply(analysis_results, `[[`, "plot_data"))

baseline_panel_ids <- baseline_covariate$covariate
baseline_analysis <- "Bias-unadjusted log RR"
baseline_data <- analysis_data %>%
  filter(analysis == baseline_analysis, covariate %in% baseline_panel_ids)
baseline_common_x <- seq(
  min(baseline_data$covariate_value, na.rm = TRUE),
  max(baseline_data$covariate_value, na.rm = TRUE),
  length.out = 100
)
baseline_common_predictions <- bind_rows(lapply(baseline_panel_ids, function(panel_id) {
  summary_row <- summary_table %>%
    filter(analysis == baseline_analysis, covariate == panel_id)
  if (nrow(summary_row) != 1) stop("Expected one baseline sensitivity summary row for: ", panel_id)
  fit_name <- basename(gsub("\\\\", "/", summary_row$fit_file[[1]]))
  fit_file <- file.path(fits_dir, fit_name)
  if (!file.exists(fit_file)) stop("Missing accepted baseline sensitivity fit: ", fit_file)
  fit <- readRDS(fit_file)
  predict_bayes_mortality_metareg_raw(fit, baseline_common_x, "baseline_mortality") %>%
    transmute(
      analysis = baseline_analysis,
      analysis_id = "bias_unadjusted_logrr",
      covariate = panel_id,
      y = RR,
      y_lower = lower_RR,
      y_upper = upper_RR,
      x = 100 * x,
      pred,
      ci_lb,
      ci_ub,
      RR,
      lower_RR,
      upper_RR
    )
}))

write.csv(
  summary_table,
  file.path(out_dir, paste0(output_prefix, "_summary.csv")),
  row.names = FALSE,
  fileEncoding = "UTF-8"
)
write.csv(
  prediction_table,
  file.path(out_dir, paste0(output_prefix, "_predictions.csv")),
  row.names = FALSE,
  fileEncoding = "UTF-8"
)
write.csv(
  analysis_data,
  file.path(out_dir, paste0(output_prefix, "_data.csv")),
  row.names = FALSE,
  fileEncoding = "UTF-8"
)
write.csv(
  baseline_common_predictions,
  file.path(out_dir, paste0(output_prefix, "_baseline_common_predictions.csv")),
  row.names = FALSE,
  fileEncoding = "UTF-8"
)
write.csv(
  data.frame(
    key = c(
      "created_at", "iter", "warmup", "chains", "seed", "adapt_delta", "max_treedepth",
      "canonical_out_dir", "canonical_overall_tau2", "output_prefix"
    ),
    value = c(
      format(Sys.time(), "%Y-%m-%d %H:%M:%S %z"),
      bayes_settings$iter,
      bayes_settings$warmup,
      bayes_settings$chains,
      bayes_settings$seed,
      bayes_settings$adapt_delta,
      bayes_settings$max_treedepth,
      canonical_out_dir,
      main_tau2_for_weights,
      output_prefix
    )
  ),
  file.path(out_dir, paste0(output_prefix, "_metadata.csv")),
  row.names = FALSE,
  fileEncoding = "UTF-8"
)

legacy_files <- c(
  list.files(out_dir, pattern = paste0("^", legacy_prefix), full.names = TRUE),
  list.files(fits_dir, pattern = paste0("^", legacy_prefix), full.names = TRUE)
)
if (length(legacy_files) > 0) {
  unlink(legacy_files, force = TRUE)
}

message("Saved 1.5 Bayesian mortality-process sensitivity cache to: ", out_dir)
