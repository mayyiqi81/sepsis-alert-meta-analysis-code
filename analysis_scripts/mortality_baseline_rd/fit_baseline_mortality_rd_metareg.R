# Figure S6E data and model generation.
#
# This script owns the S6E analysis cache: mortality RD construction, Bayesian
# random-effects meta-regression fits, predictions, audit tables, and fit RDS
# files. The companion Figure_S6E_baseline_mortality_and_mortality_rd.R script is display-only.

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
source(file.path(project_root, "analysis_scripts", "common", "project_paths.R"))
options(sepsis_project_root = set_project_root_from_script(project_root))
script_dir <- project_root
workspace_dir <- script_dir

source(file.path(workspace_dir, "analysis_scripts", "common", "bayes_prior_scale_spec.R"))
source(file.path(workspace_dir, "analysis_scripts", "common", "bayes_mortality_population_helpers.R"))
source(file.path(workspace_dir, "analysis_scripts", "common", "bayes_mortality_functions.R"))
source(file.path(workspace_dir, "analysis_scripts", "common", "mortality_rd_functions.R"))

out_dir <- Sys.getenv(
  "S6E_DATA_DIR",
  unset = file.path("computed_data", "mortality_baseline_rd")
)
if (!grepl("^([A-Za-z]:[\\\\/]|/)", out_dir)) out_dir <- file.path(workspace_dir, out_dir)
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

output_prefix <- "Figure_S6E_baseline_mortality_and_mortality_rd"
fits_dir <- file.path(out_dir, "fits")
if (!dir.exists(fits_dir)) dir.create(fits_dir, recursive = TRUE)

input_file <- file.path(workspace_dir, "source_csv", "mortality.csv")
rd_inputs_file <- file.path(workspace_dir, "source_csv", "mortality_rd_inputs.csv")
nrs_bias_sd <- 0.20
bayes_settings <- bayes_meta_settings_from_env()
prior_variant <- Sys.getenv("BAYES_PRIOR_VARIANT", unset = "canonical")

study_point_size <- function(
  n_total,
  denominator = 29852,
  min_size = 0.95,
  max_size = 5.00
) {
  n_for_size <- ifelse(
    is.finite(n_total) & n_total > 0,
    pmin(n_total, denominator),
    NA_real_
  )
  ifelse(
    is.finite(n_for_size),
    max_size * (n_for_size / denominator)^(1 / 5),
    min_size
  )
}

safe_num <- function(x) {
  safe_mortality_num(x)
}

fmt_p_tail_value <- function(p) {
  ifelse(is.na(p), NA_character_, ifelse(p < 0.001, "<0.001", sprintf("%.3f", p)))
}

read_mortality_input <- function(path) {
  raw <- read_mortality_csv(path, lower_names = FALSE)
  required_cols <- c(
    "Study", "Design", "events_e", "n_e", "events_c", "n_c",
    "Baseline mortality rate in control group", "Clinical population scope"
  )
  missing_cols <- setdiff(required_cols, names(raw))
  if (length(missing_cols) > 0) {
    stop("Missing required mortality columns: ", paste(missing_cols, collapse = ", "))
  }
  for (col in c("aRR", "aRR_lower_95", "aRR_upper_95")) {
    if (!col %in% names(raw)) raw[[col]] <- NA_real_
  }
  validate_mortality_2x2(raw, context = basename(path))
  validate_adjusted_rr_fields(raw, context = basename(path))

  raw %>%
    transmute(
      study = clean_mortality_value(Study),
      design = toupper(clean_mortality_value(Design)),
      events_e = safe_num(events_e),
      n_e = safe_num(n_e),
      events_c = safe_num(events_c),
      n_c = safe_num(n_c),
      aRR = safe_num(aRR),
      aRR_lower_95 = safe_num(aRR_lower_95),
      aRR_upper_95 = safe_num(aRR_upper_95),
      clinical_population_scope = clean_mortality_value(.data[["Clinical population scope"]]),
      baseline_mortality = safe_num(.data[["Baseline mortality rate in control group"]])
    ) %>%
    filter(
      !is.na(study), study != "",
      design %in% c("RCT", "NRS"),
      !is.na(events_e), !is.na(n_e), !is.na(events_c), !is.na(n_c),
      !is.na(baseline_mortality), n_e > 0, n_c > 0,
      events_e >= 0, events_e <= n_e,
      events_c >= 0, events_c <= n_c
    )
}

construct_mortality_rd_effects <- function(df) {
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
      n_total = n_e + n_c
    )
  add_mortality_rd_effects(out, rd_inputs_file)
}

es_all <- construct_mortality_rd_effects(read_mortality_input(input_file)) %>%
  filter(is.finite(baseline_mortality), is.finite(yi_rd), is.finite(vi_rd), vi_rd > 0)

if (nrow(es_all) != 63) {
  stop("Expected 63 mortality rows for S6E, found ", nrow(es_all), ".")
}

rd_check <- es_all$events_e / es_all$n_e - es_all$events_c / es_all$n_c
event_rows <- es_all$rd_source == "event_counts"
if (max(abs(es_all$yi_rd[event_rows] - rd_check[event_rows]), na.rm = TRUE) > 1e-12) {
  stop("S6E event-count RD construction check failed.")
}
adams <- es_all[es_all$study == "Adams 2022", , drop = FALSE]
arabi <- es_all[es_all$study == "Arabi 2025", , drop = FALSE]
rd_inputs <- read_mortality_rd_inputs(rd_inputs_file)
adams_input <- rd_inputs[rd_inputs$Study == "Adams 2022", , drop = FALSE]
if (nrow(adams) != 1 || nrow(adams_input) != 1 || adams$rd_source != "adjusted_rd_reported" ||
    abs(adams$yi_rd - adams_input$rd_percentage_points / 100) > 1e-12) {
  stop("S6E Adams 2022 RD override check failed.")
}
if (nrow(arabi) != 1 || arabi$rd_source != "adjusted_rd_from_arr" || abs(arabi$yi_rd - arabi$p_c * (arabi$aRR - 1)) > 1e-12) {
  stop("S6E Arabi 2025 RD override check failed.")
}

# Derive the population directly from the same mortality input and shared
# population helper used by Figure 6. S6E therefore has no dependency on a
# display cache, and cannot silently inherit a stale or differently filtered
# subgroup membership.
organ_dysfunction_studies <- unique(
  es_all$study[is_organ_dysfunction_targeted(es_all$clinical_population_scope)]
)
if (length(organ_dysfunction_studies) != 34 || anyDuplicated(organ_dysfunction_studies)) {
  stop("S6E organ dysfunction-targeted membership must contain 34 unique mortality studies.")
}

group_specs <- list(
  all_studies = list(
    label = "All studies",
    data = es_all
  ),
  organ_dysfunction_targeted = list(
    label = "Organ dysfunction-targeted",
    data = es_all %>% filter(study %in% organ_dysfunction_studies)
  )
)
if (nrow(group_specs$organ_dysfunction_targeted$data) != 34) {
  stop("Expected 34 organ dysfunction-targeted mortality rows for S6E.")
}

run_group <- function(group_id, group_spec) {
  dat <- group_spec$data %>%
    filter(is.finite(baseline_mortality), is.finite(yi_rd), is.finite(vi_rd), vi_rd > 0)
  model_dat <- dat %>% mutate(yi = yi_rd)
  fit_path <- file.path(fits_dir, paste0(output_prefix, "_", group_id, ".brmsfit.rds"))
  force_refit <- identical(Sys.getenv("BAYES_FORCE_REFIT", unset = ""), "1")
  expected_metadata <- make_bayes_metareg_cache_metadata(
    data = model_dat,
    moderator_col = "baseline_mortality",
    vi_col = "vi_rd",
    moderator_id = "baseline_mortality",
    response_scale = "rd",
    prior_spec_id = "rd_metareg_canonical",
    uncertainty_spec_id = "mortality_rd_no_nrs_additional_uncertainty",
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
    fit <- cached_fit
  } else {
    fit <- fit_bayes_mortality_metareg_with_quality_gate(
      model_dat,
      moderator_col = "baseline_mortality",
      moderator_id = "baseline_mortality",
      response_scale = "rd",
      prior_spec_id = "rd_metareg_canonical",
      uncertainty_spec_id = "mortality_rd_no_nrs_additional_uncertainty",
      vi_col = "vi_rd",
      settings = bayes_settings,
      prior_variant = prior_variant
    )
    saveRDS(fit, fit_path)
  }

  fit_summary <- summarise_bayes_mortality_metareg(fit, "baseline_mortality")
  # Display-only conversions: the fitted RD and baseline moderator remain on
  # their analytic scales; these fields are written in percentage units for
  # downstream display tables.
  dat <- dat %>%
    mutate(
      group_id = group_id,
      group = group_spec$label,
      covariate_value = baseline_mortality * 100,
      y_plot = yi_rd * 100
    )

  x_grid_raw <- seq(min(dat$baseline_mortality), max(dat$baseline_mortality), length.out = 100)
  predictions <- predict_bayes_mortality_metareg_raw(fit, x_grid_raw, "baseline_mortality") %>%
    transmute(
      group_id = group_id,
      group = group_spec$label,
      baseline_mortality = x_grid_raw,
      pred = pred,
      ci_lb = ci_lb,
      ci_ub = ci_ub,
      # Display-only conversion from RD proportions to percentage points.
      y = pred * 100,
      y_lower = ci_lb * 100,
      y_upper = ci_ub * 100
    )

  summary <- data.frame(
    group_id = group_id,
    group = group_spec$label,
    covariate = "baseline_mortality",
    k = nrow(dat),
    n_rct = sum(dat$design == "RCT"),
    n_nrs = sum(dat$design == "NRS"),
    min_value = min(dat$baseline_mortality),
    max_value = max(dat$baseline_mortality),
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
       prior_spec_id = "rd_metareg_canonical",
       prior_variant = prior_variant,
    p_value = fit_summary$p_value,
    p_value_formatted = fmt_p_tail_value(fit_summary$p_value),
    tau = fit_summary$tau,
    tau2 = fit_summary$tau2,
    rhat_max = fit_summary$rhat_max,
    ess_bulk_min = fit_summary$ess_bulk_min,
    ess_tail_min = fit_summary$ess_tail_min,
    divergent_transitions = fit_summary$divergent_transitions,
    method = fit_summary$method,
    response_scale = "Mortality RD proportion in model; displayed in percentage points",
    moderator_scale = "Baseline mortality proportion in model; displayed in percentage points",
    nrs_uncertainty = "None on the RD scale",
    rd_construction = "Event-count RD except explicit inputs in source_csv/mortality_rd_inputs.csv; Arabi 2025 delta-method variance includes adjusted-log-RR and control-risk sampling uncertainty under an independence assumption",
    fit_file = fit_path,
    stringsAsFactors = FALSE
  )

  list(
    summary = summary,
    predictions = predictions,
    data = dat,
    diagnostics = extract_bayes_diagnostics(fit),
    fit = fit
  )
}

results <- lapply(names(group_specs), function(group_id) run_group(group_id, group_specs[[group_id]]))
names(results) <- names(group_specs)

summary_table <- bind_rows(lapply(results, `[[`, "summary"))
prediction_table <- bind_rows(lapply(results, `[[`, "predictions"))
analysis_data <- bind_rows(lapply(results, `[[`, "data")) %>%
  mutate(point_size = study_point_size(n_e + n_c))
diagnostics_table <- bind_rows(lapply(seq_along(results), function(i) {
  results[[i]]$diagnostics %>% mutate(group_id = names(results)[i], group = group_specs[[i]]$label, .before = 1)
}))
metadata_table <- data.frame(
  key = c(
    "input_file", "response", "moderator", "nrs_additional_uncertainty",
    "rd_construction", "organ_dysfunction_scope_rule", "bayes_settings_iter",
    "bayes_settings_warmup", "bayes_settings_chains", "bayes_settings_seed"
  ),
  value = c(
    normalizePath(input_file, mustWork = FALSE),
    "Mortality risk difference; modelled as proportion and displayed as percentage points",
    "Baseline mortality in the control group; modelled as proportion",
    "None on RD scale",
    "Shared mortality RD construction from source_csv/mortality_rd_inputs.csv; Arabi 2025 delta-method variance includes adjusted-log-RR and control-risk sampling uncertainty under an independence assumption",
    "Shared population helper: is_organ_dysfunction_targeted(clinical_population_scope)",
    bayes_settings$iter, bayes_settings$warmup, bayes_settings$chains, bayes_settings$seed
  ),
  stringsAsFactors = FALSE
)

# Preserve the original per-group prediction cache and add the common grid
# used by the display script, so no posterior prediction is done during plot
# assembly.
common_x_grid <- seq(
  min(es_all$baseline_mortality, na.rm = TRUE),
  max(es_all$baseline_mortality, na.rm = TRUE),
  length.out = 100
)
plot_prediction_table <- bind_rows(lapply(names(group_specs), function(group_id) {
  fit <- results[[group_id]]$fit
  group_id_value <- group_id
  group_label <- group_specs[[group_id]]$label
  predict_bayes_mortality_metareg_raw(fit, common_x_grid, "baseline_mortality") %>%
    transmute(
      group_id = group_id_value,
      group = group_label,
      baseline_mortality = common_x_grid,
      pred = pred,
      ci_lb = ci_lb,
      ci_ub = ci_ub,
      # Display-only conversion from RD proportions to percentage points.
      y = pred * 100,
      y_lower = ci_lb * 100,
      y_upper = ci_ub * 100
    )
}))

utils::write.csv(summary_table, file.path(out_dir, paste0(output_prefix, "_summary.csv")), row.names = FALSE, fileEncoding = "UTF-8")
utils::write.csv(prediction_table, file.path(out_dir, paste0(output_prefix, "_predictions.csv")), row.names = FALSE, fileEncoding = "UTF-8")
utils::write.csv(analysis_data, file.path(out_dir, paste0(output_prefix, "_data.csv")), row.names = FALSE, fileEncoding = "UTF-8")
utils::write.csv(diagnostics_table, file.path(out_dir, paste0(output_prefix, "_diagnostics.csv")), row.names = FALSE, fileEncoding = "UTF-8")
utils::write.csv(metadata_table, file.path(out_dir, paste0(output_prefix, "_metadata.csv")), row.names = FALSE, fileEncoding = "UTF-8")
utils::write.csv(
  plot_prediction_table,
  file.path(out_dir, paste0(output_prefix, "_plot_predictions.csv")),
  row.names = FALSE,
  fileEncoding = "UTF-8"
)

message("Saved S6E Bayesian data and model cache to: ", out_dir)
