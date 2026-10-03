# Process-outcome sensitivity analyses (calculation only).

library(metafor)
library(dplyr)

script_arg <- c(grep("^--file=", commandArgs(FALSE), value = TRUE), grep("\\.R$", commandArgs(FALSE), value = TRUE, ignore.case = TRUE))
script_dir <- if (length(script_arg) > 0) dirname(normalizePath(sub("^--file=", "", script_arg[1]), mustWork = FALSE)) else getwd()
project_root <- normalizePath(file.path(script_dir, "..", ".."), winslash = "/", mustWork = TRUE)
script_dir <- project_root
out_dir <- file.path(project_root, "computed_data", "process_sensitivity")
fit_control <- list(maxiter = 10000, stepadj = 0.1)
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

clean_value <- function(x) {
  x <- trimws(as.character(x))
  x[x %in% c("", "/", "NA", "Na", "na", "N/A")] <- NA_character_
  x
}

safe_num <- function(x) suppressWarnings(as.numeric(trimws(as.character(x))))

is_include_flag <- function(x) {
  x <- tolower(trimws(as.character(x)))
  x == "yes" | grepl("^include", x)
}

is_yes_flag <- function(x) {
  tolower(trimws(as.character(x))) == "yes"
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
      effect_basis = ifelse(use_adjusted_rd, "Adjusted RD from 4.3 marked \u2021 row", "Crude RD from 2x2 counts")
    )
}

safe_model_value <- function(x) {
  ifelse(is.null(x) || length(x) == 0, NA_real_, as.numeric(x)[1])
}

fmt_effect <- function(x, lower, upper, digits = 3) {
  sprintf(paste0("%.", digits, "f [%.", digits, "f, %.", digits, "f]"), x, lower, upper)
}

fit_random <- function(dat, outcome_label = "unknown outcome", scenario_label = "unknown scenario") {
  if (nrow(dat) < 1) stop("No studies retained for scenario: ", outcome_label, " / ", scenario_label)
  if (nrow(dat) == 1) {
    yi <- dat$yi[1]
    se <- sqrt(dat$vi[1])
    return(list(
      k = 1,
      b = yi,
      ci.lb = yi - 1.96 * se,
      ci.ub = yi + 1.96 * se,
      tau2 = NA_real_,
      I2 = NA_real_,
      QE = NA_real_,
      QEp = NA_real_,
      p = 1,
      method_used = "single study"
    ))
  }
  fit <- try(rma(yi = yi, vi = vi, data = dat, method = "REML", control = fit_control), silent = TRUE)
  if (inherits(fit, "try-error")) {
    stop(
      "Process sensitivity REML model failed to converge; ",
      "no alternative heterogeneity estimator was substituted. ",
      "outcome = ", outcome_label, "; scenario = ", scenario_label,
      "; k = ", nrow(dat), "; retained studies = ",
      paste(dat$study, collapse = "; ")
    )
  }
  fit$method_used <- "REML"
  fit
}

read_binary_process <- function(file, outcome_label) {
  raw <- read.csv(file, check.names = FALSE, encoding = "UTF-8")
  names(raw) <- trimws(sub("^\ufeff", "", names(raw)))
  dat <- raw %>%
    transmute(
      outcome = outcome_label,
      study = clean_value(Study),
      design = toupper(clean_value(Design)),
      events_e = safe_num(`Intervention events`),
      n_e = safe_num(`Intervention total`),
      events_c = safe_num(`Control events`),
      n_c = safe_num(`Control total`),
      adjusted_rd = safe_num(`Adjusted RD, percentage points`),
      adjusted_lower_95 = safe_num(`Adjusted lower 95, percentage points`),
      adjusted_upper_95 = safe_num(`Adjusted upper 95, percentage points`),
      `Main analysis eligibility` = "Yes",
      `SA1 exclude substantial co-intervention`,
      `SA2 strict outcome definition`,
      sa3_potential_outliers_display = clean_value(`SA3 potential outliers display`),
      sa3_potential_outlier_flag = is_yes_flag(`SA3 potential outliers display`)
    ) %>%
    filter(
      !is.na(study),
      design %in% c("RCT", "NRS"),
      !is.na(events_e), !is.na(n_e), !is.na(events_c), !is.na(n_c),
      n_e > 0, n_c > 0,
      events_e >= 0, events_e <= n_e,
      events_c >= 0, events_c <= n_c
    ) %>%
    add_binary_process_rd_effect() %>%
    mutate(
      yi = selected_rd,
      vi = selected_se^2
    ) %>%
    filter(!is.na(yi), !is.na(vi), vi > 0)
  out <- dat %>% as.data.frame()
  for (nm in c(
    "Main analysis eligibility",
    "SA1 exclude substantial co-intervention",
    "SA2 strict outcome definition"
  )) {
    out[[nm]] <- dat[[nm]]
  }
  out
}

read_antibiotic_time <- function(file, outcome_label) {
  raw <- read.csv(file, check.names = FALSE, encoding = "UTF-8")
  names(raw) <- trimws(sub("^\ufeff", "", names(raw)))
  raw %>%
    transmute(
      outcome = outcome_label,
      study = clean_value(Study),
      design = toupper(clean_value(Design)),
      mean_e = safe_num(`Intervention value (h)`),
      sd_e = safe_num(`Intervention SD (h)`),
      n_e = safe_num(`Intervention n`),
      mean_c = safe_num(`Control value (h)`),
      sd_c = safe_num(`Control SD (h)`),
      n_c = safe_num(`Control n`),
      adjusted_md = safe_num(`Adjusted MD, hours`),
      adjusted_lower_95 = safe_num(`Adjusted lower 95, hours`),
      adjusted_upper_95 = safe_num(`Adjusted upper 95, hours`),
      `Main analysis eligibility` = "Yes",
      `SA1 exclude substantial co-intervention`,
      `SA2 strict outcome definition`,
      sa3_potential_outliers_display = clean_value(`SA3 potential outliers display`),
      sa3_potential_outlier_flag = is_yes_flag(`SA3 potential outliers display`)
    ) %>%
    mutate(
      crude_md = mean_e - mean_c,
      crude_se = sqrt(sd_e^2 / n_e + sd_c^2 / n_c),
      adjusted_se = (adjusted_upper_95 - adjusted_lower_95) / (2 * 1.96),
      use_adjusted_md = !is.na(adjusted_md) &
        !is.na(adjusted_lower_95) &
        !is.na(adjusted_upper_95) &
        adjusted_upper_95 >= adjusted_lower_95 &
        adjusted_se > 0,
      yi = ifelse(use_adjusted_md, adjusted_md, crude_md),
      vi = ifelse(use_adjusted_md, adjusted_se^2, crude_se^2)
    ) %>%
    filter(!is.na(study), design %in% c("RCT", "NRS"), !is.na(yi), !is.na(vi), vi > 0)
}

process_specs <- list(
  list(file = "process_blood_culture.csv", outcome = "Blood culture", type = "rd", prefix = "3.2_blood_culture_collection_completion_proportion_meta_analysis"),
  list(file = "process_lactate_measurement.csv", outcome = "Lactate measurement", type = "rd", prefix = "3.1_lactate_measurement_completion_proportion_meta_analysis"),
  list(file = "process_antibiotic_proportion.csv", outcome = "Antibiotic proportion", type = "rd", prefix = "3.3_timely_antibiotic_administration_proportion_meta_analysis"),
  list(file = "process_fluid_resuscitation.csv", outcome = "Fluid resuscitation", type = "rd", prefix = "3.4_fluid_resuscitation_completion_proportion_meta_analysis"),
  list(file = "process_sepsis_bundle.csv", outcome = "Sepsis bundle completion", type = "rd", prefix = "3.5_overall_sepsis_bundle_completion_proportion_meta_analysis"),
  list(file = "process_antibiotic_time.csv", outcome = "Antibiotic time", type = "md", prefix = "3.6_time_to_antibiotic_administration_meta_analysis")
)
process_specs <- lapply(process_specs, function(spec) {
  spec$file <- file.path(script_dir, "source_csv", spec$file)
  spec
})

scenario_specs <- list(
  list(analysis = "Main analysis", scenario = "Main analysis", column = "Main analysis eligibility"),
  list(analysis = "Sensitivity analyses", scenario = "SA1: exclude substantial co-intervention", column = "SA1 exclude substantial co-intervention"),
  list(analysis = "Sensitivity analyses", scenario = "SA2: strict outcome definition", column = "SA2 strict outcome definition"),
  list(analysis = "Sensitivity analyses", scenario = "SA3: outlier exclusion sensitivity analysis", column = NA_character_, filter = "sa3")
)

diagnose_outliers <- function(dat) {
  if (nrow(dat) < 3) {
    return(data.frame(
      study = dat$study,
      yi = dat$yi,
      rstudent = NA_real_,
      cook_d = NA_real_,
      dffits = NA_real_,
      hat = NA_real_,
      influential = FALSE,
      potential_outlier = FALSE,
      bonferroni_outlier = FALSE,
      eligible = FALSE,
      stringsAsFactors = FALSE
    ))
  }
  fit <- fit_random(dat, paste(unique(dat$outcome), collapse = "; "), "SA3 outlier diagnostics")
  rs <- try(rstudent(fit), silent = TRUE)
  infl <- try(influence(fit), silent = TRUE)
  if (inherits(rs, "try-error") || inherits(infl, "try-error")) {
    return(data.frame(
      study = dat$study,
      yi = dat$yi,
      rstudent = NA_real_,
      cook_d = NA_real_,
      dffits = NA_real_,
      hat = NA_real_,
      influential = FALSE,
      potential_outlier = FALSE,
      bonferroni_outlier = FALSE,
      eligible = FALSE,
      stringsAsFactors = FALSE
    ))
  }
  infl_df <- as.data.frame(infl$inf)
  z <- rs$z
  bonf_cutoff <- qnorm(1 - 0.05 / (2 * fit$k))
  influential_flag <- if ("inf" %in% names(infl$inf)) infl$inf$inf else if ("inf" %in% names(infl_df)) infl_df$inf else rep(FALSE, nrow(dat))
  if (is.factor(influential_flag)) influential_flag <- as.character(influential_flag)
  if (!is.logical(influential_flag)) influential_flag <- influential_flag %in% c(TRUE, "TRUE", "TRUE ", "true", "1", 1)
  out <- data.frame(
    study = dat$study,
    yi = dat$yi,
    rstudent = z,
    cook_d = infl_df$cook.d,
    dffits = infl_df$dffits,
    hat = infl_df$hat,
    influential = influential_flag,
    potential_outlier = abs(z) > 2,
    bonferroni_outlier = abs(z) > bonf_cutoff,
    bonferroni_cutoff = bonf_cutoff,
    stringsAsFactors = FALSE
  )
  out$eligible <- out$potential_outlier
  out$classification <- ifelse(
    out$potential_outlier & out$influential,
    "Outlier and influential",
    ifelse(
      out$potential_outlier,
      "Outlier only",
      ifelse(out$influential, "Influential only", "Neither")
    )
  )
  out[order(abs(out$rstudent), decreasing = TRUE), , drop = FALSE]
}

apply_sa3 <- function(dat) {
  dx <- diagnose_outliers(dat)
  dx$step <- 1L
  stored_removed <- dat %>%
    filter(sa3_potential_outlier_flag) %>%
    pull(study)
  computed_removed <- dx %>%
    filter(isTRUE(potential_outlier) | potential_outlier %in% TRUE) %>%
    pull(study)
  if (!setequal(stored_removed, computed_removed)) {
    discrepancy <- paste(
      c(
        paste0("stored_only=", paste(setdiff(stored_removed, computed_removed), collapse = ";")),
        paste0("computed_only=", paste(setdiff(computed_removed, stored_removed), collapse = ";"))
      ),
      collapse = " | "
    )
    stop("SA3 stored/display flags do not match computed externally studentized residuals for ",
         unique(dat$outcome), ": ", discrepancy)
  }
  current <- dat %>% filter(!study %in% computed_removed)
  attr(current, "removed") <- computed_removed
  attr(current, "diagnostics") <- dx
  current
}

canonical_process_main_row <- function(kept, spec) {
  cache_dir <- file.path(project_root, "computed_data", "process_main")
  summary_table <- read.csv(
    file.path(cache_dir, paste0(spec$prefix, "_plot_summary.csv")),
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  effects <- read.csv(
    file.path(cache_dir, paste0(spec$prefix, "_effect_sizes.csv")),
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  overall <- summary_table[summary_table$analysis_scope == "Overall", , drop = FALSE]
  if (nrow(overall) != 1L) stop("Canonical process summary must contain one Overall row: ", spec$outcome)

  current <- kept[, c("study", "design", "yi", "vi"), drop = FALSE]
  canonical <- effects[, c("study", "design", "yi", "vi"), drop = FALSE]
  rownames(current) <- NULL
  rownames(canonical) <- NULL
  if (!identical(current$study, canonical$study) || !identical(current$design, canonical$design) ||
      !isTRUE(all.equal(current$yi, canonical$yi, tolerance = 1e-12)) ||
      !isTRUE(all.equal(current$vi, canonical$vi, tolerance = 1e-12))) {
    stop("Process sensitivity main-analysis inputs do not match the canonical main analysis: ", spec$outcome)
  }

  weights <- 1 / (effects$vi + overall$tau2)
  rct_weight <- 100 * sum(weights[effects$design == "RCT"]) / sum(weights)
  effect_name <- if (identical(spec$type, "rd")) "RD percentage points" else "MD hours"
  data.frame(
    outcome = spec$outcome,
    analysis = "Main analysis",
    scenario = "Main analysis",
    k = overall$k,
    k_rct = sum(effects$design == "RCT"),
    k_nrs = sum(effects$design == "NRS"),
    effect = effect_name,
    estimate = overall$estimate,
    lower_95 = overall$lower_95,
    upper_95 = overall$upper_95,
    effect_label = fmt_effect(overall$estimate, overall$lower_95, overall$upper_95, 2),
    tau2 = overall$tau2,
    I2_percent = overall$I2_percent,
    p_heterogeneity = overall$Q_p_value,
    method = overall$method,
    note = "",
    rct_weight_percent = rct_weight,
    rct_weight_label = sprintf("%.1f%%", rct_weight),
    i2_label = paste0(sprintf("%.1f", overall$I2_percent), "%"),
    retained_studies = paste(effects$study, collapse = "; "),
    stringsAsFactors = FALSE
  )
}

fit_scenario <- function(dat, spec, scenario) {
  if (!is.na(scenario$column)) {
    kept <- dat %>% filter(is_include_flag(.data[[scenario$column]]))
  } else {
    kept <- dat %>% filter(is_include_flag(.data[["Main analysis eligibility"]]))
    kept <- apply_sa3(kept)
  }
  if (identical(scenario$scenario, "Main analysis")) {
    return(canonical_process_main_row(kept, spec))
  }
  fit <- fit_random(kept, spec$outcome, scenario$scenario)
  if (spec$type == "rd") {
    estimate <- safe_model_value(fit$b)
    lower <- safe_model_value(fit$ci.lb)
    upper <- safe_model_value(fit$ci.ub)
    effect <- "RD percentage points"
    effect_label <- fmt_effect(estimate, lower, upper, 2)
  } else {
    estimate <- safe_model_value(fit$b)
    lower <- safe_model_value(fit$ci.lb)
    upper <- safe_model_value(fit$ci.ub)
    effect <- "MD hours"
    effect_label <- fmt_effect(estimate, lower, upper, 2)
  }
  weights <- 1 / (kept$vi + safe_model_value(fit$tau2))
  rct_weight <- ifelse(sum(weights, na.rm = TRUE) > 0, 100 * sum(weights[kept$design == "RCT"], na.rm = TRUE) / sum(weights, na.rm = TRUE), NA_real_)
  removed <- attr(kept, "removed")
  note <- if (!is.null(removed) && length(removed) > 0) {
    paste0("Excluded by SA3 potential-outlier flag: ", paste(removed, collapse = "; "))
  } else if (!is.na(scenario$column)) {
    ""
  } else {
    "No study was flagged for SA3 potential-outlier exclusion."
  }
  data.frame(
    outcome = spec$outcome,
    analysis = scenario$analysis,
    scenario = scenario$scenario,
    k = nrow(kept),
    k_rct = sum(kept$design == "RCT"),
    k_nrs = sum(kept$design == "NRS"),
    effect = effect,
    estimate = estimate,
    lower_95 = lower,
    upper_95 = upper,
    effect_label = effect_label,
    tau2 = safe_model_value(fit$tau2),
    I2_percent = safe_model_value(fit$I2),
    p_heterogeneity = safe_model_value(fit$QEp),
    method = fit$method_used,
    note = note,
    rct_weight_percent = rct_weight,
    rct_weight_label = ifelse(is.na(rct_weight), "", sprintf("%.1f%%", rct_weight)),
    i2_label = ifelse(is.na(safe_model_value(fit$I2)), "", paste0(sprintf("%.1f", safe_model_value(fit$I2)), "%")),
    retained_studies = paste(kept$study, collapse = "; "),
    stringsAsFactors = FALSE
  )
}

all_results <- list()
diagnostics_results <- list()
for (spec in process_specs) {
  file <- spec$file
  dat <- if (spec$type == "rd") read_binary_process(file, spec$outcome) else read_antibiotic_time(file, spec$outcome)
  for (scenario in scenario_specs) {
    all_results[[length(all_results) + 1]] <- fit_scenario(dat, spec, scenario)
    if (!is.null(scenario$filter) && scenario$filter == "sa3") {
      main_dat <- dat %>% filter(is_include_flag(.data[["Main analysis eligibility"]]))
      sa3_dat <- apply_sa3(main_dat)
      dx <- attr(sa3_dat, "diagnostics")
      if (!is.null(dx) && nrow(dx) > 0) {
        dx$outcome <- spec$outcome
        dx$scenario <- scenario$scenario
        diagnostics_results[[length(diagnostics_results) + 1]] <- dx
      }
    }
  }
}

summary_table <- bind_rows(all_results)
write.csv(
  summary_table,
  file.path(out_dir, "process_outcomes_sensitivity_final_results_summary.csv"),
  row.names = FALSE
)
if (length(diagnostics_results) > 0) {
  write.csv(
    bind_rows(diagnostics_results) %>% select(outcome, scenario, step, everything()),
    file.path(out_dir, "process_outcomes_sensitivity_outlier_diagnostics.csv"),
    row.names = FALSE
  )
}

