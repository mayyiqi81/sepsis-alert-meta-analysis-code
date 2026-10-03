# Figure S6D analysis cache builder.

library(metafor)
library(dplyr)

script_arg <- c(grep("^--file=", commandArgs(FALSE), value = TRUE), grep("\\.R$", commandArgs(FALSE), value = TRUE, ignore.case = TRUE))
script_dir <- if (length(script_arg) > 0) dirname(normalizePath(sub("^--file=", "", script_arg[1]), mustWork = FALSE)) else getwd()
project_root <- normalizePath(file.path(script_dir, "..", ".."), winslash = "/", mustWork = TRUE)
workspace_dir <- project_root
source(file.path(workspace_dir, "analysis_scripts", "common", "project_paths.R"))
options(sepsis_project_root = set_project_root_from_script(workspace_dir))
source(file.path(workspace_dir, "analysis_scripts", "common", "bayes_prior_scale_spec.R"))

out_dir <- file.path(workspace_dir, "computed_data", "mortality_baseline_process")
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
output_prefix <- "Figure_S6D_baseline_process_level_and_process_improvement"
fits_dir <- file.path(out_dir, "fits")
if (!dir.exists(fits_dir)) dir.create(fits_dir, recursive = TRUE)

# 2. Display and analysis specifications -----------------------------------

process_covariates <- data.frame(
  outcome_id = c(
    "lactate_measurement",
    "blood_culture",
    "antibiotic_proportion",
    "fluid_resuscitation",
    "sepsis_bundle",
    "antibiotic_time"
  ),
  label = c(
    "Lactate measurement",
    "Blood culture",
    "Antibiotic proportion",
    "Fluid resuscitation",
    "Sepsis bundle",
    "Antibiotic time"
  ),
  file = c(
    "process_lactate_measurement.csv",
    "process_blood_culture.csv",
    "process_antibiotic_proportion.csv",
    "process_fluid_resuscitation.csv",
    "process_sepsis_bundle.csv",
    "process_antibiotic_time.csv"
  ),
  outcome_type = c("binary", "binary", "binary", "binary", "binary", "time"),
  stringsAsFactors = FALSE
)
process_covariates$file <- file.path(workspace_dir, "source_csv", process_covariates$file)
process_panel_ids <- matrix(process_covariates$outcome_id, nrow = 2, byrow = TRUE)

capitalise_reporting_label <- function(label) {
  paste0(toupper(substr(label, 1, 1)), substr(label, 2, nchar(label)))
}

figure_caption_title <- paste0(
  "Figure S6D. Supporting analysis: baseline process level and process improvement."
)
figure_caption_text <- paste0(
  "Separate univariable random-effects meta-regressions were fitted within each process outcome to examine the association between the baseline process level in the control group and the observed process effect. ",
  "Points show study-specific process effects, solid lines show fitted relationships, and shaded bands show 95% confidence intervals (CIs). ",
  "Blue triangles indicate randomised controlled trials (RCTs) and pink circles indicate non-randomised studies (NRSs). ",
  "Point size is scaled according to study sample size, with extreme sample sizes capped for visualisation; marker size does not represent meta-regression weight. ",
  "For completion outcomes, process effects are expressed as risk differences (RDs) in percentage points; antibiotic time is expressed as the intervention-minus-control time difference in hours, with negative values indicating shorter time with the intervention. Adjusted time differences were selected when available; otherwise, differences in harmonised group-level values were used. ",
  "These analyses are descriptive supporting analyses of baseline dependence and potential ceiling effects; because baseline level contributes directly to unadjusted process differences, the observed associations should not be interpreted as independent causal effect modification."
)


# 3. Data construction and frequentist model -------------------------------

safe_num <- function(x) {
  x <- gsub("%", "", trimws(as.character(x)), fixed = TRUE)
  suppressWarnings(as.numeric(x))
}

fmt_change <- function(x) {
  ifelse(is.na(x), "", sprintf("%.2f", x))
}

fmt_classic_p <- function(p) {
  ifelse(
    is.na(p),
    "",
    ifelse(p < 0.001, "<0.001", sprintf("%.3f", p))
  )
}

add_binary_process_rd_effect <- function(dat) {
  dat %>%
    mutate(
      crude_rd = 100 * (events_e / n_e - events_c / n_c),
      crude_se = 100 * sqrt(
        (events_e / n_e) * (1 - events_e / n_e) / n_e +
          (events_c / n_c) * (1 - events_c / n_c) / n_c
      ),
      adjusted_se = (adjusted_upper_95 - adjusted_lower_95) / (2 * 1.96),
      use_adjusted_rd = !is.na(adjusted_rd) &
        !is.na(adjusted_lower_95) &
        !is.na(adjusted_upper_95) &
        adjusted_upper_95 >= adjusted_lower_95 &
        adjusted_se > 0,
      selected_rd = ifelse(use_adjusted_rd, adjusted_rd, crude_rd),
      selected_se = ifelse(use_adjusted_rd, adjusted_se, crude_se),
      effect_basis = ifelse(
        use_adjusted_rd,
        "Adjusted RD from marked adjusted process row",
        "Crude RD from 2x2 counts"
      )
    )
}

read_binary_process <- function(spec) {
  raw <- read.csv(spec$file, stringsAsFactors = FALSE, check.names = FALSE)
  raw %>%
    transmute(
      outcome_id = spec$outcome_id,
      label = spec$label,
      outcome_type = spec$outcome_type,
      study = trimws(Study),
      design = toupper(trimws(Design)),
      events_e = safe_num(`Intervention events`),
      n_e = safe_num(`Intervention total`),
      events_c = safe_num(`Control events`),
      n_c = safe_num(`Control total`),
      adjusted_rd = safe_num(`Adjusted RD, percentage points`),
      adjusted_lower_95 = safe_num(`Adjusted lower 95, percentage points`),
      adjusted_upper_95 = safe_num(`Adjusted upper 95, percentage points`),
      baseline_level = 100 * events_c / n_c,
      n_total = n_e + n_c
    ) %>%
    filter(
      !is.na(study), study != "",
      design %in% c("RCT", "NRS"),
      !is.na(events_e), !is.na(n_e), !is.na(events_c), !is.na(n_c),
      n_e > 0, n_c > 0,
      events_e >= 0, events_e <= n_e,
      events_c >= 0, events_c <= n_c
    ) %>%
    add_binary_process_rd_effect() %>%
    mutate(
      yi = selected_rd,
      vi = selected_se^2,
      effect = yi,
      se = sqrt(vi),
      raw_effect_label = sprintf("RD %.2f percentage points", yi)
    )
}

read_time_process <- function(spec) {
  # Reuse the primary process analysis selection, including valid adjusted MDs.
  effect_path <- file.path(workspace_dir, "computed_data", "process_main",
    "3.6_time_to_antibiotic_administration_meta_analysis_effect_sizes.csv")
  if (!file.exists(effect_path)) {
    stop("Build process_main before the S6D baseline-process analysis.")
  }
  selected <- read.csv(effect_path, stringsAsFactors = FALSE, check.names = FALSE)
  selected %>%
    transmute(
      outcome_id = spec$outcome_id,
      label = spec$label,
      outcome_type = spec$outcome_type,
      study = study,
      design = design,
      baseline_level = mean_c,
      intervention_mean = mean_e,
      control_mean = mean_c,
      intervention_sd = sd_e,
      control_sd = sd_c,
      intervention_n = n_e,
      control_n = n_c,
      n_e = n_e,
      n_c = n_c,
      n_total = n_e + n_c,
      yi = yi,
      vi = vi,
      effect = yi,
      se = sqrt(vi),
      raw_effect_label = sprintf("Time difference %.3f h", yi),
      effect_basis = ifelse(use_adjusted_md,
        "Adjusted intervention-minus-control time difference from process_main",
        "Intervention-minus-control time difference from harmonised group values in process_main")
    ) %>%
    filter(
      !is.na(study), study != "",
      design %in% c("RCT", "NRS"),
      is.finite(baseline_level), is.finite(yi), is.finite(vi), vi > 0
    )
}

analysis_data <- bind_rows(lapply(seq_len(nrow(process_covariates)), function(i) {
  spec <- process_covariates[i, , drop = FALSE]
  if (identical(spec$outcome_type, "time")) read_time_process(spec) else read_binary_process(spec)
})) %>%
  filter(is.finite(baseline_level), is.finite(yi), is.finite(vi), vi > 0)

fit_meta_regression <- function(dat) {
  fit <- try(
    metafor::rma(
      yi = yi,
      vi = vi,
      mods = ~ baseline_level,
      data = dat,
      method = "REML",
      control = list(maxiter = 1000)
    ),
    silent = TRUE
  )
  if (inherits(fit, "try-error")) {
    fit <- try(
      metafor::rma(
        yi = yi,
        vi = vi,
        mods = ~ baseline_level,
        data = dat,
        method = "REML",
        control = list(maxiter = 10000)
      ),
      silent = TRUE
    )
    if (inherits(fit, "try-error")) {
      stop(
        "Figure S6D baseline-process meta-regression REML model failed to converge; ",
        "no alternative heterogeneity estimator was substituted. ",
        "outcome = ", paste(unique(dat$label), collapse = "; "),
        "; k = ", nrow(dat), "; retained studies = ",
        paste(dat$study, collapse = "; ")
      )
    }
  }
  attr(fit, "fit_method") <- "REML"
  fit
}

nice_y_config <- function(dat, pred) {
  y_values <- c(dat$y_plot, pred$y_lower, pred$y_upper)
  y_values <- y_values[is.finite(y_values)]
  reference_y <- 0
  y_values <- c(y_values, reference_y)

  if (unique(dat$outcome_type) == "binary") {
    lower <- floor(min(y_values) / 5) * 5
    upper <- ceiling(max(y_values) / 5) * 5
    if (lower == upper) {
      lower <- lower - 5
      upper <- upper + 5
    }
    breaks <- pretty(c(lower, upper), n = 5)
    return(list(
      y_limits = c(lower, upper),
      y_breaks = breaks[breaks >= lower & breaks <= upper],
      y_title = "Process improvement RD (percentage points)",
      reference_y = 0
    ))
  }

  breaks <- pretty(range(y_values), n = 5)
  list(
    y_limits = range(breaks),
    y_breaks = breaks,
    y_title = "Antibiotic time difference (h)",
    reference_y = 0
  )
}

run_outcome <- function(dat, outcome_id) {
  if (nrow(dat) < 3 || length(unique(dat$baseline_level)) < 2) {
    stop("S6D requires at least three studies and two baseline values for ", outcome_id)
  }

  fit <- fit_meta_regression(dat)
  fit_method <- attr(fit, "fit_method")
  fit_path <- file.path(fits_dir, paste0(output_prefix, "_", outcome_id, ".rds"))
  saveRDS(fit, fit_path)

  coefficient_table <- coef(summary(fit))
  slope <- coefficient_table["baseline_level", "estimate"]
  slope_lower <- coefficient_table["baseline_level", "ci.lb"]
  slope_upper <- coefficient_table["baseline_level", "ci.ub"]
  p_value <- as.numeric(fit$QMp)
  moderator_id <- if (identical(unique(dat$outcome_type), "binary")) {
    "process_completion"
  } else {
    "antibiotic_time"
  }
  moderator_spec <- get_bayes_metareg_moderator_spec(moderator_id)
  reporting_increment <- as.numeric(moderator_spec$reported_increment)
  reporting_label <- capitalise_reporting_label(
    sub("( increase| change)$", " higher baseline", moderator_spec$reported_increment_label)
  )

  dat <- dat %>%
    mutate(
      random_weight = 100 * (1 / (vi + fit$tau2)) /
        sum(1 / (vi + fit$tau2), na.rm = TRUE),
      y_plot = yi
    )

  x_grid <- seq(min(dat$baseline_level), max(dat$baseline_level), length.out = 100)
  pred <- predict(fit, newmods = x_grid)
  pred_df <- data.frame(
    outcome_id = outcome_id,
    label = unique(dat$label),
    outcome_type = unique(dat$outcome_type),
    baseline_level = x_grid,
    pred = pred$pred,
    ci_lb = pred$ci.lb,
    ci_ub = pred$ci.ub,
    y = pred$pred,
    y_lower = pred$ci.lb,
    y_upper = pred$ci.ub,
    stringsAsFactors = FALSE
  )

  y_config <- nice_y_config(dat, pred_df)
  annotation_label <- if (unique(dat$outcome_type) == "binary") {
    display_beta <- slope * reporting_increment
    display_lower <- slope_lower * reporting_increment
    display_upper <- slope_upper * reporting_increment
    paste0(
      "Change in process RD: ", fmt_change(display_beta),
      " (95% CI ", fmt_change(display_lower), " to ", fmt_change(display_upper), ")\n",
      reporting_label, "; P = ", fmt_classic_p(p_value)
    )
  } else {
    paste0(
      "Change in time MD: ", fmt_change(slope * reporting_increment),
      " (95% CI ", fmt_change(slope_lower * reporting_increment), " to ", fmt_change(slope_upper * reporting_increment), ")\n",
      reporting_label, if (p_value < 0.001) "; P < 0.001" else paste0("; P = ", fmt_classic_p(p_value))
    )
  }

  summary <- data.frame(
    outcome_id = outcome_id,
    label = unique(dat$label),
    outcome_type = unique(dat$outcome_type),
    k = nrow(dat),
    n_rct = sum(dat$design == "RCT"),
    n_nrs = sum(dat$design == "NRS"),
    min_baseline = min(dat$baseline_level),
    max_baseline = max(dat$baseline_level),
    beta = unname(slope),
    beta_lower = unname(slope_lower),
    beta_upper = unname(slope_upper),
    p_value = p_value,
    p_value_formatted = fmt_classic_p(p_value),
    tau2 = as.numeric(fit$tau2),
    I2_percent = as.numeric(fit$I2),
    AIC = AIC(fit),
    BIC = BIC(fit),
    method = fit_method,
    fit_file = fit_path,
    effect_source = if (unique(dat$outcome_type) == "binary") {
      "Adjusted RD when a valid adjusted estimate was available; otherwise crude RD from 2x2 counts"
    } else {
      "Selected intervention-minus-control time differences and variances from process_main"
    },
    stringsAsFactors = FALSE
  )

  list(
    summary = summary,
    predictions = pred_df,
    data = dat,
    y_config = y_config,
    annotation_label = annotation_label
  )
}

results <- lapply(split(analysis_data, analysis_data$outcome_id), function(dat) {
  run_outcome(dat, unique(dat$outcome_id))
})
names(results) <- vapply(results, function(x) x$summary$outcome_id, character(1))

summary_table <- bind_rows(lapply(results, `[[`, "summary")) %>%
  mutate(outcome_order = match(outcome_id, process_covariates$outcome_id)) %>%
  arrange(outcome_order) %>%
  select(-outcome_order)
prediction_table <- bind_rows(lapply(results, `[[`, "predictions"))
analysis_data_out <- bind_rows(lapply(results, `[[`, "data"))
# 4. Audit outputs -----------------------------------------------------------

summary_table$display_annotation <- vapply(
  summary_table$outcome_id,
  function(id) results[[id]]$annotation_label,
  character(1)
)
diagnostics_table <- summary_table %>%
  transmute(
    outcome_id,
    k_model = k,
    k_displayed = k,
    n_rct,
    n_nrs,
    method,
    prediction_min = min_baseline,
    prediction_max = max_baseline,
    hidden_after_fit = FALSE,
    note = "No post-fit display-only study exclusion; model and displayed study sets are identical."
  )
metadata_table <- data.frame(
  key = c("input_directory", "output_prefix", "antibiotic_time_display_rule", "legacy_filter_removed"),
  value = c(
    normalizePath(workspace_dir, mustWork = FALSE),
    output_prefix,
    "Intervention-minus-control time difference in hours; baseline is the harmonised control-group time in hours; adjusted estimates selected as in process_main.",
    "Jung 2018 and baseline_level > 20 were not excluded from fitting or display."
  ),
  stringsAsFactors = FALSE
)

utils::write.csv(summary_table, file.path(out_dir, paste0(output_prefix, "_summary.csv")), row.names = FALSE, fileEncoding = "UTF-8")
utils::write.csv(prediction_table, file.path(out_dir, paste0(output_prefix, "_predictions.csv")), row.names = FALSE, fileEncoding = "UTF-8")
utils::write.csv(analysis_data_out, file.path(out_dir, paste0(output_prefix, "_data.csv")), row.names = FALSE, fileEncoding = "UTF-8")
utils::write.csv(diagnostics_table, file.path(out_dir, paste0(output_prefix, "_diagnostics.csv")), row.names = FALSE, fileEncoding = "UTF-8")
utils::write.csv(metadata_table, file.path(out_dir, paste0(output_prefix, "_metadata.csv")), row.names = FALSE, fileEncoding = "UTF-8")
message("Wrote Figure S6D analysis cache to: ", out_dir)
