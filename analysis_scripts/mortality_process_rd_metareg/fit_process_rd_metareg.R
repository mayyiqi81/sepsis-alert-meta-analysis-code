# Univariable mortality meta-regression for process completion RD/time covariates.

library(metafor)
library(dplyr)
library(ggplot2)
library(ggtext)
library(cowplot)
library(gridExtra)
library(grid)

script_arg <- c(grep("^--file=", commandArgs(FALSE), value = TRUE), grep("\\.R$", commandArgs(FALSE), value = TRUE, ignore.case = TRUE))
script_dir <- if (length(script_arg) > 0) dirname(normalizePath(sub("^--file=", "", script_arg[1]), mustWork = FALSE)) else getwd()
project_root <- normalizePath(file.path(script_dir, "..", ".."), winslash = "/", mustWork = TRUE)
source(file.path(project_root, "analysis_scripts", "common", "project_paths.R"))
options(sepsis_project_root = set_project_root_from_script(project_root))
script_dir <- project_root
workspace_dir <- script_dir
source(file.path(workspace_dir, "analysis_scripts", "common", "bayes_prior_scale_spec.R"))
source(file.path(workspace_dir, "analysis_scripts", "common", "antibiotic_time_direction_contract.R"))
source(file.path(workspace_dir, "analysis_scripts", "common", "bayes_mortality_functions.R"))
source(file.path(workspace_dir, "analysis_scripts", "common", "bayes_mortality_population_helpers.R"))
input_file <- file.path(workspace_dir, "source_csv", "mortality.csv")
nrs_bias_sd <- 0.20
bayes_settings <- bayes_meta_settings_from_env()
prior_variant <- Sys.getenv("BAYES_PRIOR_VARIANT", unset = "canonical")
analysis_label <- sprintf("NRS bias SD = %.2f", nrs_bias_sd)
weight_label <- "adjusted inverse-variance weight"
out_dir <- Sys.getenv(
  "MORTALITY_PROCESS_RD_METAREG_DATA_DIR",
  unset = "computed_data/mortality_process_rd_metareg"
)
if (!grepl("^([A-Za-z]:[\\\\/]|/)", out_dir)) out_dir <- file.path(workspace_dir, out_dir)
canonical_out_dir <- Sys.getenv(
  "BAYES_CANONICAL_OUT_DIR",
  unset = "computed_data/mortality_main/canonical"
)
if (!grepl("^([A-Za-z]:[\\\\/]|/)", canonical_out_dir)) canonical_out_dir <- file.path(workspace_dir, canonical_out_dir)
fits_dir <- file.path(out_dir, "fits")
output_prefix <- "mortality_univariable_process_rd_meta_regression_nrs_bias_sd"
MIN_POINT_SIZE <- 0.35
POINT_WEIGHT_SCALE <- 2
FIG_WIDTH_PX <- 1281
ROW_HEIGHT_PX <- 380
ROW_SPACER_PX <- 22
BASELINE_ROW_HEIGHT_PX <- 462
BASELINE_ROW_WIDTH_PX <- round(FIG_WIDTH_PX * 0.80)
FIG_HEIGHT_PX <- ROW_HEIGHT_PX * 2 + ROW_SPACER_PX * 2 + BASELINE_ROW_HEIGHT_PX
FIG_WIDTH <- 12.81
FIG_HEIGHT <- 12.66
plot_font <- "Segoe UI"

if (.Platform$OS.type == "windows") {
  windowsFonts("Segoe UI" = windowsFont("Segoe UI"))
}

global_palette <- list(
  blue = c(light = "#DCECF9", main = "#AED7F3", deep = "#5A9CCB"),
  pink = c(light = "#F7D5DC", main = "#E8A3B1", deep = "#D76A84"),
  transition = c(light = "#E8E0F4", main = "#C7C1E6", deep = "#9E8BC8"),
  all = c(light = "#CFC2D2", main = "#8A748F", deep = "#5A5368"),
  blue_gray = c(light = "#E2E7EE", main = "#CFD6E3", deep = "#7A879C")
)

colors <- list(
  text = "#0A265E",
  muted = global_palette$blue_gray["deep"],
  subtext = "#5F6E84",
  faint = "#F4F7FA",
  rule = global_palette$blue_gray["main"],
  pooled = global_palette$all["deep"],
  rct = global_palette$blue["deep"],
  nrs = global_palette$pink["deep"]
)

if (!dir.exists(out_dir)) {
  dir.create(out_dir, recursive = TRUE)
}
if (!dir.exists(fits_dir)) {
  dir.create(fits_dir, recursive = TRUE)
}

redundant_outputs <- c(
  paste0(output_prefix, "_coefficients.csv"),
  paste0(output_prefix, "_predictions.csv"),
  paste0(output_prefix, "_data.csv"),
  paste0("mortality_meta_regression_nrs_bias_sd_", c(
    "blood_culture",
    "lactate_measurement",
    "antibiotic_proportion",
    "fluid_resuscitation",
    "sepsis_bundle",
    "antibiotic_time"
  ), rep(c(".pdf", ".png"), each = 6))
)
unlink(file.path(out_dir, redundant_outputs), force = TRUE)

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
  stringsAsFactors = FALSE
)

covariate_order <- covariates$label

fmt_p <- function(p) {
  ifelse(is.na(p), NA_character_, ifelse(p < 0.001, "<0.001", sprintf("%.3f", p)))
}

fmt_number <- function(x) {
  ifelse(is.na(x), "", sprintf("%.3f", x))
}

fmt_number_display <- function(x) {
  gsub("-", "\u2212", fmt_number(x), fixed = TRUE)
}

safe_num <- function(x) {
  x <- gsub("%", "", trimws(as.character(x)), fixed = TRUE)
  suppressWarnings(as.numeric(x))
}

safe_model_value <- function(x) {
  ifelse(is.null(x) || length(x) == 0, NA_real_, as.numeric(x)[1])
}

fit_meta_regression <- function(dat, covariate_name, fit_path) {
  moderator_id <- if (identical(covariate_name, "baseline_mortality")) {
    "baseline_mortality"
  } else if (identical(covariate_name, "antibiotic_time")) {
    "antibiotic_time"
  } else {
    "process_completion"
  }
  expected_metadata <- make_bayes_metareg_cache_metadata(
    data = dat,
    moderator_col = covariate_name,
    vi_col = "vi_adjusted",
    moderator_id = moderator_id,
    response_scale = "log_rr",
    prior_spec_id = "log_rr_metareg_canonical",
    uncertainty_spec_id = "primary_nrs_log_rr_bias_sd_0.20",
    settings = bayes_settings,
    prior_variant = prior_variant
  )
  force_refit <- identical(Sys.getenv("BAYES_FORCE_REFIT", unset = ""), "1")
  fit <- NULL
  if (file.exists(fit_path) && !force_refit) {
    cached_fit <- readRDS(fit_path)
    if (isTRUE(assert_cached_metareg_fit_matches(cached_fit, expected_metadata))) {
      fit <- cached_fit
    } else {
      fit <- upgrade_legacy_cached_metareg_fit(cached_fit, expected_metadata)
      if (!is.null(fit)) saveRDS(fit, fit_path)
    }
  }
  if (is.null(fit)) {
    fit <- fit_bayes_mortality_metareg_with_quality_gate(
      dat,
      moderator_col = covariate_name,
      moderator_id = moderator_id,
      response_scale = "log_rr",
      prior_spec_id = "log_rr_metareg_canonical",
      uncertainty_spec_id = "primary_nrs_log_rr_bias_sd_0.20",
      vi_col = "vi_adjusted",
      settings = bayes_settings,
      prior_variant = prior_variant
    )
    saveRDS(fit, fit_path)
  }
  fit
}

df_raw <- read.csv(
  input_file,
  stringsAsFactors = FALSE,
  check.names = FALSE,
  encoding = "UTF-8"
)
names(df_raw) <- trimws(sub("^\ufeff", "", names(df_raw)))
df_raw <- df_raw[, names(df_raw) != "", drop = FALSE]
assert_antibiotic_time_source_direction(df_raw)

required_cols <- c(
  "Study", "Design", "events_e", "n_e", "events_c", "n_c",
  "Baseline mortality rate in control group", "Clinical population scope",
  covariates$column
)
missing_cols <- setdiff(required_cols, names(df_raw))
if (length(missing_cols) > 0) {
  stop("Missing required columns: ", paste(missing_cols, collapse = ", "))
}

df <- df_raw %>%
  transmute(
    study = trimws(Study),
    design = toupper(trimws(Design)),
    events_e = safe_num(events_e),
    n_e = safe_num(n_e),
    events_c = safe_num(events_c),
    n_c = safe_num(n_c),
    aRR = if ("aRR" %in% names(df_raw)) safe_num(aRR) else NA_real_,
    aRR_lower_95 = if ("aRR_lower_95" %in% names(df_raw)) safe_num(aRR_lower_95) else NA_real_,
    aRR_upper_95 = if ("aRR_upper_95" %in% names(df_raw)) safe_num(aRR_upper_95) else NA_real_,
    clinical_population_scope = trimws(.data[["Clinical population scope"]]),
    baseline_mortality_rate_in_control_group = safe_num(.data[["Baseline mortality rate in control group"]]),
    baseline_mortality = baseline_mortality_rate_in_control_group,
    blood_culture = safe_num(.data[["Blood culture RD, percentage points"]]),
    lactate_measurement = safe_num(.data[["Lactate measurement RD, percentage points"]]),
    antibiotic_proportion = safe_num(.data[["Antibiotic proportion RD, percentage points"]]),
    fluid_resuscitation = safe_num(.data[["Fluid resuscitation RD, percentage points"]]),
    sepsis_bundle = safe_num(.data[["Sepsis bundle RD, percentage points"]]),
    antibiotic_time = safe_num(.data[["Antibiotic time MD, hours"]])
  ) %>%
  filter(
    !is.na(study), study != "",
    design %in% c("RCT", "NRS"),
    !is.na(events_e), !is.na(n_e),
    !is.na(events_c), !is.na(n_c),
    !is.na(baseline_mortality),
    n_c > 0
  )

validate_mortality_2x2(df, context = "fit_process_rd_metareg mortality data")
validate_adjusted_rr_fields(df, context = "fit_process_rd_metareg mortality data")

es_all <- escalc(
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
    yi_crude = yi,
    vi_crude = vi,
    uses_adjusted_rr = !is.na(aRR) & !is.na(aRR_lower_95) & !is.na(aRR_upper_95),
    yi = ifelse(uses_adjusted_rr, log(aRR), yi_crude),
    vi_original = ifelse(uses_adjusted_rr, ((log(aRR_upper_95) - log(aRR_lower_95)) / (2 * 1.96))^2, vi_crude),
    vi_adjusted = ifelse(design == "NRS", vi_original + nrs_bias_sd^2, vi_original),
    vi_bias_adjusted = vi_adjusted,
    se_adjusted = sqrt(vi_adjusted),
    RR = exp(yi),
    lower_RR = exp(yi - 1.96 * se_adjusted),
    upper_RR = exp(yi + 1.96 * se_adjusted),
    inverse_variance_weight = 1 / vi_adjusted
  )

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

run_univariable <- function(covariate_name) {
  dat <- es_all %>%
    filter(!is.na(.data[[covariate_name]]))
  moderator_id <- if (identical(covariate_name, "antibiotic_time")) {
    "antibiotic_time"
  } else {
    "process_completion"
  }

  if (nrow(dat) < 3 || length(unique(dat[[covariate_name]])) < 2) {
    dat <- dat %>% mutate(random_weight = NA_real_)
    return(list(
      summary = data.frame(
        covariate = covariate_name,
        k = nrow(dat),
        n_rct = sum(dat$design == "RCT"),
        n_nrs = sum(dat$design == "NRS"),
        min_value = ifelse(nrow(dat) == 0, NA_real_, min(dat[[covariate_name]], na.rm = TRUE)),
        max_value = ifelse(nrow(dat) == 0, NA_real_, max(dat[[covariate_name]], na.rm = TRUE)),
        QM = NA_real_,
        df = NA_real_,
        p_value = NA_real_,
        p_value_formatted = NA_character_,
        tau2 = NA_real_,
        I2_percent = NA_real_,
        AIC = NA_real_,
        BIC = NA_real_,
        method = NA_character_,
        note = "Skipped: fewer than three studies or fewer than two unique covariate values",
        stringsAsFactors = FALSE
      ),
      predictions = data.frame(),
      data = dat
    ))
  }

  fit_path <- file.path(fits_dir, paste0(output_prefix, "_", covariate_name, ".brmsfit.rds"))
  fit <- fit_meta_regression(dat, covariate_name, fit_path)
  fit_summary <- summarise_bayes_mortality_metareg(fit, covariate_name)
  dat <- dat %>%
    mutate(
      random_weight = 100 * (1 / (vi_adjusted + fit_summary$tau2)) /
        sum(1 / (vi_adjusted + fit_summary$tau2), na.rm = TRUE)
    )

  x_grid_raw <- seq(min(dat[[covariate_name]]), max(dat[[covariate_name]]), length.out = 100)
  pred_df <- predict_bayes_mortality_metareg_raw(fit, x_grid_raw, moderator_id) %>%
    mutate(x = x_grid_raw, covariate = covariate_name, .before = 1)

  list(
    summary = data.frame(
      covariate = covariate_name,
      k = nrow(dat),
      n_rct = sum(dat$design == "RCT"),
      n_nrs = sum(dat$design == "NRS"),
      min_value = min(dat[[covariate_name]], na.rm = TRUE),
      max_value = max(dat[[covariate_name]], na.rm = TRUE),
      QM = NA_real_,
      df = NA_real_,
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
       prior_spec_id = "log_rr_metareg_canonical",
       prior_variant = prior_variant,
      tau = fit_summary$tau,
      tau2 = fit_summary$tau2,
      rhat_max = fit_summary$rhat_max,
      ess_bulk_min = fit_summary$ess_bulk_min,
      ess_tail_min = fit_summary$ess_tail_min,
      divergent_transitions = fit_summary$divergent_transitions,
      AIC = NA_real_,
      BIC = NA_real_,
      method = fit_summary$method,
      note = paste0("NRS variance inflated by sigma_bias^2 = ", nrs_bias_sd^2),
      stringsAsFactors = FALSE
    ),
    predictions = pred_df,
    data = dat
  )
}

results <- lapply(covariates$covariate, run_univariable)
summary_table <- bind_rows(lapply(results, `[[`, "summary")) %>%
  left_join(covariates, by = "covariate") %>%
  select(covariate, label, column, everything()) %>%
  arrange(match(covariate, covariates$covariate))
prediction_table <- bind_rows(lapply(results, `[[`, "predictions")) %>%
  left_join(covariates, by = "covariate") %>%
  select(covariate, label, column, everything())

analysis_data <- bind_rows(lapply(seq_along(results), function(i) {
  covariate_name <- covariates$covariate[i]
  results[[i]]$data %>%
    mutate(
      covariate = covariate_name,
      covariate_value = .data[[covariate_name]]
    ) %>%
    select(
      covariate, study, design, clinical_population_scope,
      baseline_mortality_rate_in_control_group, baseline_mortality, covariate_value,
      events_e, n_e, events_c, n_c,
      yi, vi_original, vi_adjusted, se_adjusted,
      RR, lower_RR, upper_RR, inverse_variance_weight, random_weight, main_random_weight
    )
})) %>%
  left_join(covariates, by = "covariate") %>%
  select(covariate, label, column, everything())

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
baseline_selected <- c("All studies", "Severe sepsis/Sepsis-3")
baseline_subgroups <- list(
  list(
    subgroup = "All studies",
    data = es_all
  ),
  list(
    subgroup = "Severe sepsis/Sepsis-3",
    data = es_all %>%
      filter(is_organ_dysfunction_targeted(clinical_population_scope))
  )
)
baseline_common_grid <- seq(
  0,
  max(unlist(lapply(baseline_subgroups, function(x) x$data$baseline_mortality)), na.rm = TRUE),
  length.out = 100
)

run_baseline_regression <- function(dat, subgroup_label, x_grid) {
  if (nrow(dat) < 3 || length(unique(dat$baseline_mortality)) < 2) {
    stop("Baseline mortality meta-regression requires at least three studies and two unique baseline mortality values: ", subgroup_label)
  }

  safe_subgroup_id <- gsub("[^A-Za-z0-9]+", "_", tolower(subgroup_label))
  safe_subgroup_id <- gsub("(^_|_$)", "", safe_subgroup_id)
  fit_path <- file.path(fits_dir, paste0(output_prefix, "_baseline_mortality_", safe_subgroup_id, ".brmsfit.rds"))
  fit <- fit_meta_regression(dat, "baseline_mortality", fit_path)
  fit_summary <- summarise_bayes_mortality_metareg(fit, "baseline_mortality")

  dat <- dat %>%
    mutate(
      subgroup = subgroup_label,
      random_weight = 100 * (1 / (vi_adjusted + fit_summary$tau2)) /
        sum(1 / (vi_adjusted + fit_summary$tau2), na.rm = TRUE)
    )

  predictions <- predict_bayes_mortality_metareg_raw(fit, x_grid, "baseline_mortality") %>%
    transmute(
      subgroup = subgroup_label,
      baseline_mortality = x_grid,
      pred,
      ci_lb,
      ci_ub,
      RR,
      lower_RR,
      upper_RR
    )

  summary <- data.frame(
    subgroup = subgroup_label,
    covariate = "baseline_mortality",
    label = "Baseline mortality",
    k = nrow(dat),
    n_rct = sum(dat$design == "RCT"),
    n_nrs = sum(dat$design == "NRS"),
    min_value = min(dat$baseline_mortality, na.rm = TRUE),
    max_value = max(dat$baseline_mortality, na.rm = TRUE),
    QM = NA_real_,
    df = NA_real_,
    p_value = fit_summary$p_value,
    p_value_formatted = fmt_p(fit_summary$p_value),
    prob_direction = fit_summary$prob_direction,
    beta = fit_summary$beta,
    beta_lower = fit_summary$beta_lower,
    beta_upper = fit_summary$beta_upper,
    tau = fit_summary$tau,
    tau2 = fit_summary$tau2,
    rhat_max = fit_summary$rhat_max,
    ess_bulk_min = fit_summary$ess_bulk_min,
    ess_tail_min = fit_summary$ess_tail_min,
    divergent_transitions = fit_summary$divergent_transitions,
    AIC = NA_real_,
    BIC = NA_real_,
    method = fit_summary$method,
    note = paste0("NRS variance inflated by sigma_bias^2 = ", nrs_bias_sd^2),
    stringsAsFactors = FALSE
  )

  list(summary = summary, predictions = predictions, data = dat)
}

baseline_results <- lapply(
  baseline_subgroups,
  function(x) run_baseline_regression(x$data, x$subgroup, baseline_common_grid)
)
baseline_summary <- bind_rows(lapply(baseline_results, `[[`, "summary")) %>%
  mutate(subgroup = factor(subgroup, levels = baseline_selected)) %>%
  arrange(subgroup)
baseline_predictions <- bind_rows(lapply(baseline_results, `[[`, "predictions")) %>%
  mutate(subgroup = factor(subgroup, levels = baseline_selected)) %>%
  arrange(subgroup, baseline_mortality)
baseline_data <- bind_rows(lapply(baseline_results, `[[`, "data")) %>%
  mutate(subgroup = factor(subgroup, levels = baseline_selected)) %>%
  arrange(subgroup)

write.csv(
  baseline_summary,
  file.path(out_dir, paste0(output_prefix, "_baseline_mortality_summary.csv")),
  row.names = FALSE,
  fileEncoding = "UTF-8"
)
write.csv(
  baseline_predictions,
  file.path(out_dir, paste0(output_prefix, "_baseline_mortality_predictions.csv")),
  row.names = FALSE,
  fileEncoding = "UTF-8"
)
write.csv(
  baseline_data,
  file.path(out_dir, paste0(output_prefix, "_baseline_mortality_data.csv")),
  row.names = FALSE,
  fileEncoding = "UTF-8"
)

write.csv(
  data.frame(
    key = c("created_at", "iter", "warmup", "chains", "seed", "adapt_delta", "max_treedepth", "canonical_out_dir", "canonical_overall_tau2"),
    value = c(
      format(Sys.time(), "%Y-%m-%d %H:%M:%S %z"),
      bayes_settings$iter,
      bayes_settings$warmup,
      bayes_settings$chains,
      bayes_settings$seed,
      bayes_settings$adapt_delta,
      bayes_settings$max_treedepth,
      canonical_out_dir,
      main_tau2_for_weights
    )
  ),
  file.path(out_dir, paste0(output_prefix, "_metadata.csv")),
  row.names = FALSE,
  fileEncoding = "UTF-8"
)

message("Saved mortality_univariable_process_rd_meta_regression_nrs_bias_sd Bayesian cache to: ", out_dir)
