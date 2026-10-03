# Mortality subgroup analyses with subgroup-difference p-values.
# NRS studies use bias variance sigma_bias = 0.20.

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
  "MORTALITY_SUBGROUP_DATA_DIR",
  unset = "computed_data/mortality_subgroup"
)
if (!grepl("^([A-Za-z]:[\\\\/]|/)", out_dir)) out_dir <- file.path(workspace_dir, out_dir)
fits_dir <- file.path(out_dir, "fits")
nrs_bias_sd <- 0.20
plot_font <- "Segoe UI"
fit_control <- list(maxiter = 10000, stepadj = 0.1)
bayes_settings <- bayes_meta_settings_from_env()
color_text <- "#0A265E"
color_axis <- color_text
color_estimate <- "#5A5368"
color_section_title <- "#5A5368"
color_population_title <- "#8A748F"
color_alert_title <- "#8A748F"
color_subgroup_p <- "#8A748F"

if (!dir.exists(out_dir)) {
  dir.create(out_dir, recursive = TRUE)
}
if (!dir.exists(fits_dir)) {
  dir.create(fits_dir, recursive = TRUE)
}

clean_value <- function(x) {
  x <- trimws(as.character(x))
  x[x %in% c("", "/", "NA", "Na", "na", "N/A")] <- NA_character_
  x
}

subgroup_level_mappings <- list(
  population_ascertainment = c(
    "alert-triggered" = "Alert-triggered",
    "alert-triggered + diagnosis-defined" = "Alert-triggered + diagnosis-defined",
    "all screened" = "All screened",
    "diagnosis-defined" = "Diagnosis-defined"
  ),
  clinical_population_scope = c(
    "infection-targeted" = "Infection-targeted",
    "all screened" = "All screened",
    "severe sepsis-targeted" = "Severe sepsis-targeted",
    "(severe sepsis)" = "Severe sepsis-targeted",
    "sepsis-3-targeted" = "Sepsis-3-targeted",
    "(sepsis-3)" = "Sepsis-3-targeted"
  ),
  baseline_mortality_in_control_group = c(
    "<10%" = "<10%",
    "10% to <20%" = "10% to <20%",
    ">=20%" = ">=20%"
  ),
  sample_size_1000 = c(
    "small sample (<1000)" = "Small sample (<1000)",
    "large sample (>=1000)" = "Large sample (>=1000)"
  ),
  centers_sites = c("single" = "Single", "multiple" = "Multiple"),
  deployment_setting = c(
    "icu" = "ICU", "wards" = "Wards", "ed" = "ED", "multiple settings" = "Multiple settings", "nr" = "NR"
  ),
  intervention_package_type = c(
    "alert only" = "Alert only", "actionable alert" = "Actionable alert", "co-intervention" = "Co-intervention"
  ),
  alert_targeting_approach = c(
    "suspected infection" = "Suspected infection",
    "sepsis with organ dysfunction" = "Sepsis with organ dysfunction",
    "tiered sepsis alerts" = "Tiered sepsis alerts"
  ),
  algorithm_class = c("rule" = "Rule", "ml" = "ML"),
  recipient = c(
    "physician" = "Physician", "nurse" = "Nurse", "clinical care team" = "Clinical care team",
    "dedicated response service" = "Dedicated response service", "nr" = "NR"
  ),
  delivery_mode = c("ehr display" = "EHR display", "page/message" = "Page/message", "nr" = "NR")
)

canonicalize_subgroup_level <- function(var_name, x) {
  if (!var_name %in% names(subgroup_level_mappings)) {
    stop("No explicit subgroup mapping registered for variable: ", var_name)
  }
  raw <- trimws(as.character(x))
  keys <- tolower(raw)
  mapping <- subgroup_level_mappings[[var_name]]
  canonical <- unname(mapping[keys])
  canonical[is.na(raw) | raw == ""] <- NA_character_
  unmatched <- !is.na(raw) & nzchar(raw) & is.na(canonical)
  if (any(unmatched)) {
    stop(
      "Unmapped subgroup level(s) for ", var_name, ": ",
      paste(unique(raw[unmatched]), collapse = ", ")
    )
  }
  canonical
}

fmt_rr <- function(rr, lower, upper) {
  sprintf("%.3f [%.3f, %.3f]", rr, lower, upper)
}

fmt_p <- function(p) {
  ifelse(is.na(p), "", ifelse(p < 0.001, "<0.001", sprintf("%.3f", p)))
}

fmt_p_label <- function(p) {
  p_text <- fmt_p(p)
  ifelse(p_text == "", "", sprintf("%-6s%s", p_text, ifelse(p < 0.05, "*", " ")))
}

fmt_p_subgroup_label <- function(p_subgroup) {
  p_subgroup_text <- fmt_p(p_subgroup)
  if (p_subgroup_text == "") {
    return("")
  }
  paste0("P for subgroup difference = ", p_subgroup_text)
}

subgroup_level_order <- function(var_name, x) {
  observed <- unique(as.character(x))
  preferred <- switch(
    var_name,
    population_ascertainment = c("Alert-triggered", "Alert-triggered + diagnosis-defined", "Diagnosis-defined"),
    clinical_population_scope = c("Infection-targeted", "Severe sepsis-targeted", "Sepsis-3-targeted"),
    baseline_mortality_in_control_group = c("<10%", "10% to <20%", ">=20%"),
    sample_size_1000 = c("Small sample (<1000)", "Large sample (>=1000)"),
    centers_sites = c("Single", "Multiple"),
    deployment_setting = c("ICU", "Wards", "ED", "Multiple settings"),
    intervention_package_type = c("Alert only", "Actionable alert", "Co-intervention"),
    alert_targeting_approach = c("Suspected infection", "Sepsis with organ dysfunction", "Tiered sepsis alerts"),
    algorithm_class = c("Rule", "ML"),
    recipient = c("Physician", "Nurse", "Clinical care team", "Dedicated response service"),
    delivery_mode = c("EHR display", "Page/message"),
    character()
  )
  c(preferred[preferred %in% observed], setdiff(observed, preferred))
}

safe_model_value <- function(x) {
  ifelse(is.null(x) || length(x) == 0, NA_real_, as.numeric(x)[1])
}

subgroup_fit_slug <- function(x) {
  slug <- tolower(trimws(as.character(x)))
  slug <- gsub(">=", " ge ", slug, fixed = TRUE)
  slug <- gsub("<=", " le ", slug, fixed = TRUE)
  slug <- gsub(">", " gt ", slug, fixed = TRUE)
  slug <- gsub("<", " lt ", slug, fixed = TRUE)
  slug <- gsub("%", " pct ", slug, fixed = TRUE)
  slug <- gsub("[^a-z0-9]+", "_", slug)
  gsub("(^_|_$)", "", slug)
}

legacy_subgroup_fit_ids <- function(var_name, level) {
  if (identical(var_name, "baseline_mortality_in_control_group")) {
    return(paste0(
      "mortality_subgroup_", var_name, "_",
      switch(level, "<10%" = "_10", "10% to <20%" = "10_to_20", ">=20%" = "_20", "")
    ))
  }
  if (identical(var_name, "sample_size_1000")) {
    return(paste0(
      "mortality_subgroup_", var_name, "_",
      switch(level, "Small sample (<1000)" = "small_sample_1000", "Large sample (>=1000)" = "large_sample_1000", "")
    ))
  }
  character()
}

fit_pool <- function(dat, fit_id = NULL, legacy_fit_ids = character()) {
  random_weight_summary <- function(weights) {
    ifelse(
      sum(weights, na.rm = TRUE) > 0,
      100 * sum(weights[dat$design == "RCT"], na.rm = TRUE) / sum(weights, na.rm = TRUE),
      NA_real_
    )
  }

  if (nrow(dat) == 1) {
    yi <- dat$yi[1]
    se <- sqrt(dat$vi_adjusted[1])
    rct_weight_percent <- ifelse(dat$design[1] == "RCT", 100, 0)
    return(data.frame(
      k = 1,
      k_rct = sum(dat$design == "RCT"),
      k_nrs = sum(dat$design == "NRS"),
      events_e = sum(dat$events_e),
      n_e = sum(dat$n_e),
      events_c = sum(dat$events_c),
      n_c = sum(dat$n_c),
      RR = exp(yi),
      Lower_95_CrI = exp(yi - 1.96 * se),
      Upper_95_CrI = exp(yi + 1.96 * se),
      tau = NA_real_,
      tau2 = NA_real_,
      rhat_max = NA_real_,
      ess_bulk_min = NA_real_,
      ess_tail_min = NA_real_,
      divergent_transitions = NA_real_,
      max_treedepth = NA_real_,
      max_treedepth_hits = NA_real_,
      method = "single study normal approximation",
      rct_weight_percent = rct_weight_percent,
      rct_weight_label = sprintf("%.1f%%", rct_weight_percent),
      tau_label = "",
      stringsAsFactors = FALSE
    ))
  }

  fit_path <- if (!is.null(fit_id) && nzchar(fit_id)) file.path(fits_dir, paste0(fit_id, ".brmsfit.rds")) else NA_character_
  expected_metadata <- make_bayes_model_cache_metadata(
    data = dat,
    vi_col = "vi_adjusted",
    response_scale = "log_rr",
    prior_spec_id = "log_rr_shared_current",
    uncertainty_spec_id = "primary_nrs_log_rr_bias_sd_0.20",
    settings = bayes_settings,
    study_col = "study",
    yi_col = "yi",
    formula_text = "yi | se(.bayes_se, sigma = FALSE) ~ 1 + (1 | .bayes_study)"
  )
  force_refit <- identical(Sys.getenv("BAYES_FORCE_REFIT", unset = ""), "1")
  candidate_paths <- unique(c(
    fit_path,
    file.path(fits_dir, paste0(legacy_fit_ids[nzchar(legacy_fit_ids)], ".brmsfit.rds"))
  ))
  candidate_paths <- candidate_paths[!is.na(candidate_paths) & file.exists(candidate_paths)]
  reuse <- FALSE
  cached_fit <- NULL
  if (!force_refit && length(candidate_paths) > 0L) {
    for (candidate_path in candidate_paths) {
      candidate_fit <- readRDS(candidate_path)
      if (isTRUE(assert_cached_bayes_model_matches(candidate_fit, expected_metadata))) {
        cached_fit <- candidate_fit
        reuse <- TRUE
        break
      }
      upgraded_fit <- upgrade_unattributed_cached_intercept_fit(
        candidate_fit,
        expected_data = dat,
        expected_metadata = expected_metadata
      )
      if (!is.null(upgraded_fit)) {
        cached_fit <- upgraded_fit
        reuse <- TRUE
        if (!is.na(fit_path)) saveRDS(cached_fit, fit_path)
        break
      }
    }
  }
  if (reuse) {
    fit <- cached_fit
  } else {
    fit <- fit_bayes_mortality_meta_with_quality_gate(
      dat,
      "vi_adjusted",
      response_scale = "log_rr",
      prior_spec_id = "log_rr_shared_current",
      uncertainty_spec_id = "primary_nrs_log_rr_bias_sd_0.20",
      settings = bayes_settings
    )
    if (!is.na(fit_path)) saveRDS(fit, fit_path)
  }
  bayes_summary <- summarise_bayes_mortality_fit(
    fit,
    group = "subgroup",
    k_rct = sum(dat$design == "RCT"),
    k_nrs = sum(dat$design == "NRS"),
    nrs_bias_sd = nrs_bias_sd
  )
  random_weights <- 1 / (dat$vi_adjusted + bayes_summary$tau2)
  rct_weight_percent <- random_weight_summary(random_weights)

  data.frame(
    k = bayes_summary$k,
    k_rct = sum(dat$design == "RCT"),
    k_nrs = sum(dat$design == "NRS"),
    events_e = sum(dat$events_e),
    n_e = sum(dat$n_e),
    events_c = sum(dat$events_c),
    n_c = sum(dat$n_c),
    RR = bayes_summary$RR,
    Lower_95_CrI = bayes_summary$Lower_95_CrI,
    Upper_95_CrI = bayes_summary$Upper_95_CrI,
    tau = bayes_summary$tau,
    tau2 = bayes_summary$tau2,
    rhat_max = bayes_summary$rhat_max,
    ess_bulk_min = bayes_summary$ess_bulk_min,
    ess_tail_min = bayes_summary$ess_tail_min,
    divergent_transitions = bayes_summary$divergent_transitions,
    max_treedepth = bayes_summary$max_treedepth,
    max_treedepth_hits = bayes_summary$max_treedepth_hits,
    method = "Bayesian random-effects model",
    rct_weight_percent = rct_weight_percent,
    rct_weight_label = ifelse(is.na(rct_weight_percent), "", sprintf("%.1f%%", rct_weight_percent)),
    tau_label = ifelse(is.na(bayes_summary$tau), "", sprintf("%.3f", bayes_summary$tau)),
    stringsAsFactors = FALSE
  )
}

fit_subgroup_p <- function(dat, var_name) {
  dat <- dat %>% filter(!is.na(.data[[var_name]]))
  if (nrow(dat) < 3 || length(unique(dat[[var_name]])) < 2) {
    return(list(p = NA_real_, method = "skipped"))
  }

  dat[[var_name]] <- factor(dat[[var_name]], levels = unique(dat[[var_name]]))
  fit <- rma(
    yi = yi,
    vi = vi_adjusted,
    mods = as.formula(paste0("~ ", var_name)),
    data = dat,
    method = "REML",
    control = fit_control
  )

  list(p = safe_model_value(fit$QMp), method = "REML reference")
}

df_raw <- read.csv(input_file, stringsAsFactors = FALSE, check.names = FALSE)
names(df_raw) <- trimws(names(df_raw))
if (length(names(df_raw)) > 0 && grepl("study$", names(df_raw)[1], ignore.case = TRUE)) names(df_raw)[1] <- "Study"
df_raw <- df_raw[, names(df_raw) != "", drop = FALSE]

control_risk_col <- if ("Control-group mortality risk" %in% names(df_raw)) {
  "Control-group mortality risk"
} else {
  "Baseline mortality in control group"
}
first_existing_col <- function(candidates) {
  hit <- candidates[candidates %in% names(df_raw)]
  if (length(hit) == 0) NA_character_ else hit[1]
}
population_ascertainment_col <- first_existing_col(c("Population ascertainment"))
clinical_population_scope_col <- first_existing_col(c("Clinical population scope"))
deployment_setting_col <- first_existing_col(c("Deployment setting"))
algorithm_class_col <- first_existing_col(c("Algorithm class"))
alert_targeting_approach_col <- first_existing_col(c("Alert targeting approach"))

required_cols <- c(
  "Study", "events_e", "n_e", "events_c", "n_c", "Design",
  "Intervention package type", "Recipient", "Delivery mode", "Centers / sites", control_risk_col,
  population_ascertainment_col, clinical_population_scope_col, deployment_setting_col,
  alert_targeting_approach_col, algorithm_class_col
)
required_cols <- required_cols[!is.na(required_cols)]
missing_cols <- setdiff(required_cols, names(df_raw))
if (length(missing_cols) > 0) {
  stop("Missing required columns: ", paste(missing_cols, collapse = ", "))
}

df <- df_raw %>%
  transmute(
    study = clean_value(Study),
    events_e = safe_mortality_num(events_e),
    n_e = safe_mortality_num(n_e),
    events_c = safe_mortality_num(events_c),
    n_c = safe_mortality_num(n_c),
    aRR = if ("aRR" %in% names(df_raw)) safe_mortality_num(aRR) else NA_real_,
    aRR_lower_95 = if ("aRR_lower_95" %in% names(df_raw)) safe_mortality_num(aRR_lower_95) else NA_real_,
    aRR_upper_95 = if ("aRR_upper_95" %in% names(df_raw)) safe_mortality_num(aRR_upper_95) else NA_real_,
    design = toupper(clean_value(Design)),
    total_n = safe_mortality_num(n_e) + safe_mortality_num(n_c),
    population_ascertainment = if (!is.na(population_ascertainment_col)) clean_value(.data[[population_ascertainment_col]]) else NA_character_,
    clinical_population_scope = if (!is.na(clinical_population_scope_col)) clean_value(.data[[clinical_population_scope_col]]) else NA_character_,
    deployment_setting = if (!is.na(deployment_setting_col)) clean_value(.data[[deployment_setting_col]]) else NA_character_,
    intervention_package_type = clean_value(.data[["Intervention package type"]]),
    alert_targeting_approach = if (!is.na(alert_targeting_approach_col)) clean_value(.data[[alert_targeting_approach_col]]) else NA_character_,
    algorithm_class = if (!is.na(algorithm_class_col)) clean_value(.data[[algorithm_class_col]]) else NA_character_,
    recipient = clean_value(Recipient),
    delivery_mode = clean_value(.data[["Delivery mode"]]),
    centers_sites = clean_value(.data[["Centers / sites"]]),
    baseline_mortality_in_control_group = clean_value(.data[[control_risk_col]])
  ) %>%
  filter(
    !is.na(study),
    !is.na(events_e), !is.na(n_e),
    !is.na(events_c), !is.na(n_c),
    design %in% c("RCT", "NRS")
  )

validate_mortality_2x2(df, context = "fit_subgroup_analysis mortality data")
validate_adjusted_rr_fields(df, context = "fit_subgroup_analysis mortality data")

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
    yi_crude = yi,
    vi_crude = vi,
    uses_adjusted_rr = !is.na(aRR) & !is.na(aRR_lower_95) & !is.na(aRR_upper_95),
    yi = ifelse(uses_adjusted_rr, log(aRR), yi_crude),
    vi_original = ifelse(uses_adjusted_rr, ((log(aRR_upper_95) - log(aRR_lower_95)) / (2 * 1.96))^2, vi_crude),
    vi_adjusted = ifelse(design == "NRS", vi_original + nrs_bias_sd^2, vi_original),
    sample_size_1000 = factor(
      ifelse(total_n < 1000, "Small sample (<1000)", "Large sample (>=1000)"),
      levels = c("Small sample (<1000)", "Large sample (>=1000)")
    ),
    centers_sites = factor(
      centers_sites,
      levels = c("Single", "Multiple")
    ),
    baseline_mortality_in_control_group = factor(
      baseline_mortality_in_control_group,
      levels = c("<10%", "10% to <20%", ">=20%")
    )
  )

mapping_audit_rows <- list()
for (var_name in names(subgroup_level_mappings)) {
  if (!var_name %in% names(es)) next
  raw_values <- as.character(es[[var_name]])
  canonical_values <- canonicalize_subgroup_level(var_name, raw_values)
  keep <- !is.na(raw_values) & nzchar(trimws(raw_values))
  if (any(keep)) {
    mapping_audit_rows[[length(mapping_audit_rows) + 1]] <- data.frame(
      variable = var_name,
      raw_value = raw_values[keep],
      canonical_value = canonical_values[keep],
      n_rows = 1L,
      mapping_reason = ifelse(
        tolower(trimws(raw_values[keep])) == tolower(trimws(canonical_values[keep])),
        "canonical level or case/format normalization",
        "explicit documented alias"
      ),
      stringsAsFactors = FALSE
    ) %>%
      group_by(variable, raw_value, canonical_value, mapping_reason) %>%
      summarise(n_rows = sum(n_rows), .groups = "drop")
  }
  es[[var_name]] <- canonical_values
}
write.csv(
  bind_rows(mapping_audit_rows),
  file.path(out_dir, "mortality_subgroup_level_mapping_audit.csv"),
  row.names = FALSE
)

subgroups <- data.frame(
  analysis = c(
    rep("A. Population- and study-level subgroup analyses", 5),
    rep("B. Alert-related subgroup analyses", 6)
  ),
  variable = c(
    "population_ascertainment",
    "clinical_population_scope",
    "baseline_mortality_in_control_group",
    "sample_size_1000",
    "centers_sites",
    "deployment_setting",
    "intervention_package_type",
    "alert_targeting_approach",
    "algorithm_class",
    "recipient",
    "delivery_mode"
  ),
  label = c(
    "Population ascertainment",
    "Clinical population scope",
    "Baseline mortality in control group",
    "Total analytic sample size",
    "Centres/sites",
    "Deployment setting",
    "Intervention package",
    "Alert targeting approach",
    "Algorithm type",
    "Recipient",
    "Delivery mode"
  ),
  stringsAsFactors = FALSE
)

summary_rows <- list()
p_rows <- list()
fit_inventory_rows <- list()

for (i in seq_len(nrow(subgroups))) {
  var_name <- subgroups$variable[i]
  label <- subgroups$label[i]
  analysis <- subgroups$analysis[i]
  dat <- es %>% filter(!is.na(.data[[var_name]]))

  dat <- dat %>% filter(.data[[var_name]] != "NR")
  dat <- dat %>% filter(tolower(as.character(.data[[var_name]])) != "all")
  if (var_name %in% c("population_ascertainment", "clinical_population_scope")) {
    dat <- dat %>% filter(!is_all_screened(.data[[var_name]]))
  }
  if (var_name == "recipient") {
    dat <- dat %>% filter(.data[[var_name]] != "Specialized response team")
  }

  p_res <- fit_subgroup_p(dat, var_name)
  ordered_levels <- subgroup_level_order(var_name, dat[[var_name]])
  p_rows[[length(p_rows) + 1]] <- data.frame(
    analysis = analysis,
    variable = var_name,
    label = label,
    k = nrow(dat),
    levels = paste(ordered_levels, collapse = "; "),
    p_subgroup = p_res$p,
    p_subgroup_formatted = fmt_p(p_res$p),
    p_method = p_res$method,
    nrs_bias_sd = nrs_bias_sd,
    stringsAsFactors = FALSE
  )

  for (level in ordered_levels) {
    level_dat <- dat %>% filter(.data[[var_name]] == level)
    fit_id <- paste0("mortality_subgroup_", var_name, "_", subgroup_fit_slug(level))
    pooled <- fit_pool(
      level_dat,
      fit_id = fit_id,
      legacy_fit_ids = legacy_subgroup_fit_ids(var_name, level)
    )
    fit_inventory_rows[[length(fit_inventory_rows) + 1]] <- data.frame(
      variable = var_name,
      subgroup = level,
      k = nrow(level_dat),
      fit_required = nrow(level_dat) >= 2,
      fit_file = if (nrow(level_dat) >= 2) paste0(fit_id, ".brmsfit.rds") else NA_character_,
      stringsAsFactors = FALSE
    )
    summary_rows[[length(summary_rows) + 1]] <- cbind(
      data.frame(
        analysis = analysis,
        variable = var_name,
        label = label,
        subgroup = level,
        p_subgroup = p_res$p,
        p_subgroup_formatted = fmt_p(p_res$p),
        stringsAsFactors = FALSE
      ),
      pooled
    )
  }
}

summary_table <- bind_rows(summary_rows)
p_table <- bind_rows(p_rows)
fit_inventory <- bind_rows(fit_inventory_rows)

expected_fit_files <- na.omit(fit_inventory$fit_file)
existing_fit_files <- list.files(fits_dir, pattern = "^mortality_subgroup_.*\\.brmsfit\\.rds$", full.names = FALSE)
stale_fit_files <- setdiff(existing_fit_files, expected_fit_files)
if (length(stale_fit_files) > 0L) {
  unlink(file.path(fits_dir, stale_fit_files))
}

write.csv(
  summary_table,
  file.path(out_dir, "mortality_subgroup_summary.csv"),
  row.names = FALSE
)
write.csv(
  p_table,
  file.path(out_dir, "mortality_subgroup_pvalues.csv"),
  row.names = FALSE
)
write.csv(
  fit_inventory,
  file.path(out_dir, "mortality_subgroup_fit_inventory.csv"),
  row.names = FALSE,
  na = ""
)
write.csv(
  data.frame(
    key = c("created_at", "iter", "warmup", "chains", "seed", "adapt_delta", "max_treedepth", "nrs_bias_sd"),
    value = c(
      format(Sys.time(), "%Y-%m-%d %H:%M:%S %z"),
      bayes_settings$iter,
      bayes_settings$warmup,
      bayes_settings$chains,
      bayes_settings$seed,
      bayes_settings$adapt_delta,
      bayes_settings$max_treedepth,
      nrs_bias_sd
    )
  ),
  file.path(out_dir, "mortality_subgroup_metadata.csv"),
  row.names = FALSE
)

message("Saved mortality subgroup Bayesian cache to: ", out_dir)
