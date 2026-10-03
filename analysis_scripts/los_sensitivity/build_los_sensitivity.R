# Build length-of-stay sensitivity-analysis caches.

suppressPackageStartupMessages({
  library(metafor)
  library(dplyr)
})

script_arg <- c(grep("^--file=", commandArgs(FALSE), value = TRUE), grep("\\.R$", commandArgs(FALSE), value = TRUE, ignore.case = TRUE))
script_dir <- if (length(script_arg) > 0) dirname(normalizePath(sub("^--file=", "", script_arg[1]), mustWork = FALSE)) else getwd()
project_root <- normalizePath(file.path(script_dir, "..", ".."), winslash = "/", mustWork = TRUE)
out_dir <- file.path(project_root, "computed_data", "los_sensitivity")
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
script_dir <- project_root
fit_control <- list(maxiter = 10000, stepadj = 0.1)

clean_value <- function(x) {
  x <- trimws(as.character(x))
  x[x %in% c("", "/", "NA", "Na", "na", "N/A")] <- NA_character_
  x
}

safe_num <- function(x) suppressWarnings(as.numeric(trimws(as.character(x))))

safe_model_value <- function(x) {
  ifelse(is.null(x) || length(x) == 0, NA_real_, as.numeric(x)[1])
}

is_yes <- function(x) tolower(clean_value(x)) == "yes"

fmt_effect <- function(x, lower, upper, digits = 2) {
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
      "LOS sensitivity REML model failed to converge; ",
      "no alternative heterogeneity estimator was substituted. ",
      "outcome = ", outcome_label, "; scenario = ", scenario_label,
      "; k = ", nrow(dat), "; retained studies = ",
      paste(dat$study, collapse = "; ")
    )
  }
  fit$method_used <- "REML"
  fit
}

read_los_file <- function(file, outcome_label) {
  raw <- read.csv(file, check.names = FALSE, stringsAsFactors = FALSE, encoding = "UTF-8")
  names(raw) <- trimws(sub("^\ufeff", "", names(raw)))
  raw <- raw[, names(raw) != "", drop = FALSE]
  required_cols <- c(
    "Study", "Design", "Intervention package type",
    "Control mean", "Control SD", "Control n",
    "Intervention mean", "Intervention SD", "Intervention n",
    "Adjusted effect estimate, days", "Adjusted lower 95, days", "Adjusted upper 95, days",
    "special LOS construction", "SA3 potential outliers display"
  )
  missing_cols <- setdiff(required_cols, names(raw))
  if (length(missing_cols) > 0) {
    stop("Missing required columns in ", basename(file), ": ", paste(missing_cols, collapse = ", "))
  }

  has_icu_denominator <- "ICU denominator population only" %in% names(raw)
  has_icu_setting <- "studies conducted in ICU only" %in% names(raw)
  has_source_note <- "Source note" %in% names(raw)
  has_original_dispersion_type <- "Original dispersion type" %in% names(raw)
  has_conversion_rule <- "Conversion rule" %in% names(raw)

  dat_all <- raw %>%
    transmute(
      outcome = outcome_label,
      row_order = row_number(),
      study = clean_value(Study),
      design_raw = clean_value(Design),
      design = toupper(design_raw),
      mean_c = safe_num(.data[["Control mean"]]),
      sd_c = safe_num(.data[["Control SD"]]),
      n_c = safe_num(.data[["Control n"]]),
    mean_e = safe_num(.data[["Intervention mean"]]),
    sd_e = safe_num(.data[["Intervention SD"]]),
    n_e = safe_num(.data[["Intervention n"]]),
    intervention_package_type_source = clean_value(.data[["Intervention package type"]]),
    adjusted_yi = safe_num(.data[["Adjusted effect estimate, days"]]),
    adjusted_lower = safe_num(.data[["Adjusted lower 95, days"]]),
    adjusted_upper = safe_num(.data[["Adjusted upper 95, days"]]),
    adjusted_se = (adjusted_upper - adjusted_lower) / (2 * 1.96),
    adjusted_effect_notes = if ("Adjusted effect notes" %in% names(raw)) clean_value(.data[["Adjusted effect notes"]]) else NA_character_,
      original_dispersion_type = if (has_original_dispersion_type) clean_value(.data[["Original dispersion type"]]) else NA_character_,
      conversion_rule = if (has_conversion_rule) clean_value(.data[["Conversion rule"]]) else NA_character_,
      source_note = if (has_source_note) clean_value(.data[["Source note"]]) else NA_character_,
      special_los_construction_raw = clean_value(.data[["special LOS construction"]]),
      special_los_construction = is_yes(.data[["special LOS construction"]]),
      sa3_potential_outliers_display = clean_value(.data[["SA3 potential outliers display"]]),
      sa3_potential_outlier_flag = is_yes(.data[["SA3 potential outliers display"]]),
      denominator_strictly_icu_population = if (has_icu_denominator) clean_value(.data[["ICU denominator population only"]]) else NA_character_,
      study_conducted_in_icu = if (has_icu_setting) clean_value(.data[["studies conducted in ICU only"]]) else NA_character_
    ) %>%
    mutate(
      intervention_package_type = intervention_package_type_source,
      qualitative_setting = NA_character_,
      obvious_multicomponent_cointervention = tolower(intervention_package_type) == "co-intervention",
      obvious_multicomponent_cointervention = ifelse(is.na(obvious_multicomponent_cointervention), FALSE, obvious_multicomponent_cointervention),
      eligible_adjusted = !is.na(adjusted_yi) & !is.na(adjusted_lower) & !is.na(adjusted_upper) &
        !is.na(adjusted_se) & adjusted_se > 0 & adjusted_upper > adjusted_lower & design %in% c("RCT", "NRS"),
      eligible_group =
        !is.na(study) & !is.na(mean_e) & !is.na(sd_e) & !is.na(n_e) &
        !is.na(mean_c) & !is.na(sd_c) & !is.na(n_c) &
        sd_e > 0 & sd_c > 0 & n_e > 1 & n_c > 1 &
        design %in% c("RCT", "NRS")
    )

  excluded <- dat_all %>%
    filter(!eligible_adjusted & !eligible_group)

  dat <- dat_all %>%
    filter(eligible_group, !eligible_adjusted)

  es_parts <- list()
  if (nrow(dat) > 0) {
    es_parts[[length(es_parts) + 1]] <- escalc(
      measure = "MD",
      m1i = mean_e,
      sd1i = sd_e,
      n1i = n_e,
      m2i = mean_c,
      sd2i = sd_c,
      n2i = n_c,
      data = dat
    ) %>%
      as.data.frame() %>%
      mutate(
        vi_original = vi,
        se_original = sqrt(vi_original),
        MD = yi,
        Lower_95 = yi - 1.96 * se_original,
        Upper_95 = yi + 1.96 * se_original,
        selected_effect_basis = "U"
      )
  }
  adjusted_es <- dat_all %>%
    filter(eligible_adjusted) %>%
    transmute(
      outcome, row_order, study, design_raw, design,
      mean_c, sd_c, n_c, mean_e, sd_e, n_e,
      adjusted_yi, adjusted_lower, adjusted_upper, adjusted_se, adjusted_effect_notes,
      original_dispersion_type, conversion_rule, source_note,
      special_los_construction_raw, special_los_construction,
      sa3_potential_outliers_display, sa3_potential_outlier_flag,
      denominator_strictly_icu_population, study_conducted_in_icu,
      intervention_package_type, qualitative_setting,
      obvious_multicomponent_cointervention,
      yi = adjusted_yi,
      vi = adjusted_se^2,
      vi_original = adjusted_se^2,
      se_original = adjusted_se,
      MD = adjusted_yi,
      Lower_95 = adjusted_lower,
      Upper_95 = adjusted_upper,
      selected_effect_basis = "F"
    )
  if (nrow(adjusted_es) > 0) es_parts[[length(es_parts) + 1]] <- adjusted_es
  es <- bind_rows(es_parts)

  sa3_studies <- dat_all %>%
    filter(sa3_potential_outlier_flag) %>%
    pull(study) %>%
    unique()
  list(es = es, excluded = excluded, sa3_studies = sa3_studies)
}

diagnose_outliers <- function(dat) {
  base_cols <- data.frame(
    study = dat$study,
    MD = dat$MD,
    yi = dat$yi,
    rstudent = NA_real_,
    abs_rstudent = NA_real_,
    cook_d = NA_real_,
    dffits = NA_real_,
    hat = NA_real_,
    Potential_outlier = FALSE,
    influential = FALSE,
    bonferroni_outlier = FALSE,
    bonferroni_cutoff = NA_real_,
    eligible = FALSE,
    classification = "Not assessed",
    stringsAsFactors = FALSE
  )
  if (nrow(dat) < 3) return(base_cols)
  fit <- fit_random(dat, outcome_label, group)
  rs <- try(rstudent(fit), silent = TRUE)
  infl <- try(influence(fit), silent = TRUE)
  if (inherits(rs, "try-error") || inherits(infl, "try-error")) return(base_cols)
  infl_df <- as.data.frame(infl$inf)
  z <- rs$z
  bonf_cutoff <- qnorm(1 - 0.05 / (2 * fit$k))
  influential_flag <- if ("inf" %in% names(infl$inf)) infl$inf$inf else if ("inf" %in% names(infl_df)) infl_df$inf else rep(FALSE, nrow(dat))
  if (is.factor(influential_flag)) influential_flag <- as.character(influential_flag)
  if (!is.logical(influential_flag)) influential_flag <- influential_flag %in% c(TRUE, "TRUE", "TRUE ", "true", "1", 1)
  out <- data.frame(
    study = dat$study,
    MD = dat$MD,
    yi = dat$yi,
    rstudent = z,
    abs_rstudent = abs(z),
    cook_d = infl_df$cook.d,
    dffits = infl_df$dffits,
    hat = infl_df$hat,
    Potential_outlier = abs(z) > 2,
    influential = influential_flag,
    bonferroni_outlier = abs(z) > bonf_cutoff,
    bonferroni_cutoff = bonf_cutoff,
    stringsAsFactors = FALSE
  )
  out$eligible <- out$Potential_outlier
  out$classification <- ifelse(
    out$Potential_outlier & out$influential,
    "Outlier and influential",
    ifelse(
      out$Potential_outlier,
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
    filter(isTRUE(Potential_outlier) | Potential_outlier %in% TRUE) %>%
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

scenario_specs <- list(
  list(analysis = "Main analysis", scenario = "Main analysis", filter = "main"),
  list(analysis = "Sensitivity analyses", scenario = "SA1: exclude survivor-only, free-day-derived, or non-LOS-specific denominator", filter = "special_los"),
  list(analysis = "Sensitivity analyses", scenario = "SA2: exclude obvious multicomponent cointervention", filter = "cointervention"),
  list(analysis = "Sensitivity analyses", scenario = "SA3: Outlier exclusion sensitivity analysis", filter = "sa3")
)

icu_scenario_specs <- list(
  list(analysis = "Sensitivity analyses", scenario = "SA4: restricted to ICU denominator population", filter = "icu_denominator"),
  list(analysis = "Sensitivity analyses", scenario = "SA5: restricted to studies conducted in ICU", filter = "icu_setting")
)

scenario_data <- function(dat, scenario) {
  if (scenario$filter == "main") {
    dat
  } else if (scenario$filter == "special_los") {
    dat %>% filter(!special_los_construction)
  } else if (scenario$filter == "cointervention") {
    dat %>% filter(!obvious_multicomponent_cointervention)
  } else if (scenario$filter == "sa3") {
    apply_sa3(dat)
  } else if (scenario$filter == "icu_denominator") {
    dat %>% filter(tolower(denominator_strictly_icu_population) == "yes")
  } else if (scenario$filter == "icu_setting") {
    dat %>% filter(tolower(study_conducted_in_icu) == "yes")
  } else {
    stop("Unknown scenario filter: ", scenario$filter)
  }
}

scenario_note <- function(original_dat, kept, scenario) {
  if (scenario$filter == "special_los") {
    excluded <- setdiff(original_dat$study, kept$study)
    ifelse(
      length(excluded) > 0,
      paste0("Excluded survivor-only/free-day-derived/non-LOS-specific denominator studies: ", paste(excluded, collapse = "; ")),
      "No study met the special LOS construction exclusion rule."
    )
  } else if (scenario$filter == "cointervention") {
    excluded <- setdiff(original_dat$study, kept$study)
    ifelse(
      length(excluded) > 0,
      paste0("Excluded obvious multicomponent cointervention studies: ", paste(excluded, collapse = "; ")),
      "No study met the obvious multicomponent cointervention exclusion rule."
    )
  } else if (scenario$filter == "sa3") {
    removed <- attr(kept, "removed")
    if (!is.null(removed) && length(removed) > 0) {
      paste0("Excluded by SA3 potential-outlier flag: ", paste(removed, collapse = "; "))
    } else {
      "No study was flagged for SA3 potential-outlier exclusion."
    }
  } else if (scenario$filter == "icu_denominator") {
    "Restricted to studies with denominator strictly ICU population."
  } else if (scenario$filter == "icu_setting") {
    "Restricted to studies conducted in ICU."
  } else {
    ""
  }
}

canonical_los_main_row <- function(kept, outcome_label) {
  spec <- switch(
    outcome_label,
    "Hospital LOS" = list(dir = "los_hospital_main", summary = "los_meta_summary.csv", effects = "los_effect_sizes.csv"),
    "ICU LOS" = list(dir = "los_icu_main", summary = "icu_los_meta_summary.csv", effects = "icu_los_effect_sizes.csv"),
    stop("No canonical LOS main-analysis cache registered for: ", outcome_label)
  )
  cache_dir <- file.path(project_root, "computed_data", spec$dir)
  summary_table <- read.csv(file.path(cache_dir, spec$summary), stringsAsFactors = FALSE, check.names = FALSE)
  effects <- read.csv(file.path(cache_dir, spec$effects), stringsAsFactors = FALSE, check.names = FALSE)
  overall <- summary_table[summary_table$group == "Overall", , drop = FALSE]
  if (nrow(overall) != 1L) stop("Canonical LOS summary must contain one Overall row: ", outcome_label)

  current <- kept[order(kept$row_order), c("study", "design", "yi", "vi_original"), drop = FALSE]
  canonical <- effects[order(effects$row_order), c("study", "design", "yi", "vi_original"), drop = FALSE]
  rownames(current) <- NULL
  rownames(canonical) <- NULL
  if (!identical(current$study, canonical$study) || !identical(current$design, canonical$design) ||
      !isTRUE(all.equal(current$yi, canonical$yi, tolerance = 1e-12)) ||
      !isTRUE(all.equal(current$vi_original, canonical$vi_original, tolerance = 1e-12))) {
    stop("LOS sensitivity main-analysis inputs do not match the canonical main analysis: ", outcome_label)
  }

  rct_weight <- sum(effects$weight_random[effects$design == "RCT"], na.rm = TRUE)
  data.frame(
    outcome = outcome_label,
    analysis = "Main analysis",
    scenario = "Main analysis",
    k = overall$k,
    k_rct = overall$k_rct,
    k_nrs = overall$k_nrs,
    effect = "MD days",
    estimate = overall$MD,
    lower_95 = overall$Lower_95,
    upper_95 = overall$Upper_95,
    effect_label = overall$MD_label,
    tau2 = overall$tau2,
    I2_percent = overall$I2_percent,
    p_heterogeneity = overall$p_heterogeneity,
    method = overall$method,
    note = "",
    rct_weight_percent = rct_weight,
    rct_weight_label = sprintf("%.1f%%", rct_weight),
    i2_label = paste0(sprintf("%.1f", overall$I2_percent), "%"),
    retained_studies = paste(effects$study[order(effects$row_order)], collapse = "; "),
    stringsAsFactors = FALSE
  )
}

fit_scenario <- function(dat, outcome_label, scenario) {
  kept <- scenario_data(dat, scenario)
  if (identical(scenario$filter, "main")) {
    return(canonical_los_main_row(kept, outcome_label))
  }
  fit <- fit_random(kept, outcome_label, scenario$scenario)
  estimate <- safe_model_value(fit$b)
  lower <- safe_model_value(fit$ci.lb)
  upper <- safe_model_value(fit$ci.ub)
  weights <- 1 / (kept$vi + safe_model_value(fit$tau2))
  rct_weight <- ifelse(sum(weights, na.rm = TRUE) > 0, 100 * sum(weights[kept$design == "RCT"], na.rm = TRUE) / sum(weights, na.rm = TRUE), NA_real_)
  data.frame(
    outcome = outcome_label,
    analysis = scenario$analysis,
    scenario = scenario$scenario,
    k = nrow(kept),
    k_rct = sum(kept$design == "RCT"),
    k_nrs = sum(kept$design == "NRS"),
    effect = "MD days",
    estimate = estimate,
    lower_95 = lower,
    upper_95 = upper,
    effect_label = fmt_effect(estimate, lower, upper, 2),
    tau2 = safe_model_value(fit$tau2),
    I2_percent = safe_model_value(fit$I2),
    p_heterogeneity = safe_model_value(fit$QEp),
    method = fit$method_used,
    note = scenario_note(dat, kept, scenario),
    rct_weight_percent = rct_weight,
    rct_weight_label = ifelse(is.na(rct_weight), "", sprintf("%.1f%%", rct_weight)),
    i2_label = ifelse(is.na(safe_model_value(fit$I2)), "", paste0(sprintf("%.1f", safe_model_value(fit$I2)), "%")),
    retained_studies = paste(kept$study, collapse = "; "),
    stringsAsFactors = FALSE
  )
}

empty_diagnostics <- function() {
  data.frame(
    outcome = character(),
    scenario = character(),
    step = integer(),
    study = character(),
    MD = numeric(),
    yi = numeric(),
    rstudent = numeric(),
    abs_rstudent = numeric(),
    cook_d = numeric(),
    dffits = numeric(),
    hat = numeric(),
    Potential_outlier = logical(),
    influential = logical(),
    bonferroni_outlier = logical(),
    bonferroni_cutoff = numeric(),
    eligible = logical(),
    classification = character(),
    stringsAsFactors = FALSE
  )
}

hospital <- read_los_file(file.path(script_dir, "source_csv", "hospital_los.csv"), "Hospital LOS")
icu <- read_los_file(file.path(script_dir, "source_csv", "ICU_los.csv"), "ICU LOS")
sa3_study_map <- list("Hospital LOS" = hospital$sa3_studies, "ICU LOS" = icu$sa3_studies)

all_los <- bind_rows(hospital$es, icu$es)
write.csv(
  all_los %>%
    transmute(
      outcome, study, design, mean_e, sd_e, n_e, mean_c, sd_c, n_c,
      yi, vi_original, se_original, MD, Lower_95, Upper_95,
      `special LOS construction` = special_los_construction,
      `SA3 potential outliers display` = sa3_potential_outliers_display,
      obvious_multicomponent_cointervention,
      `ICU denominator population only` = denominator_strictly_icu_population,
      `studies conducted in ICU only` = study_conducted_in_icu
    ),
  file.path(out_dir, "los_sensitivity_effect_sizes.csv"),
  row.names = FALSE
)
write.csv(
  bind_rows(hospital$excluded, icu$excluded),
  file.path(out_dir, "los_sensitivity_invalid_rows.csv"),
  row.names = FALSE
)
write.csv(
  all_los %>%
    transmute(
      outcome, study, design,
      `special LOS construction` = special_los_construction,
      `SA3 potential outliers display` = sa3_potential_outliers_display,
      obvious_multicomponent_cointervention,
      `ICU denominator population only` = denominator_strictly_icu_population,
      `studies conducted in ICU only` = study_conducted_in_icu,
      source_note, intervention_package_type, qualitative_setting
    ),
  file.path(out_dir, "los_sensitivity_exclusion_flags.csv"),
  row.names = FALSE
)

results <- list()
diagnostics_results <- list()
outcome_inputs <- list("Hospital LOS" = hospital$es, "ICU LOS" = icu$es)
for (outcome_name in names(outcome_inputs)) {
  dat <- outcome_inputs[[outcome_name]]
  specs <- scenario_specs
  if (outcome_name == "ICU LOS") specs <- c(specs, icu_scenario_specs)
  for (scenario in specs) {
    results[[length(results) + 1]] <- fit_scenario(dat, outcome_name, scenario)
    if (scenario$filter == "sa3") {
      sa3_dat <- scenario_data(dat, scenario)
      dx <- attr(sa3_dat, "diagnostics")
      if (!is.null(dx) && nrow(dx) > 0) {
        dx$outcome <- outcome_name
        dx$scenario <- scenario$scenario
        diagnostics_results[[length(diagnostics_results) + 1]] <- dx
      }
    }
  }
}

summary_table <- bind_rows(results)
write.csv(
  summary_table,
  file.path(out_dir, "los_sensitivity_final_results_summary.csv"),
  row.names = FALSE
)

diagnostics_table <- if (length(diagnostics_results) > 0) {
  bind_rows(diagnostics_results) %>% select(outcome, scenario, step, everything())
} else {
  empty_diagnostics()
}
write.csv(
  diagnostics_table,
  file.path(out_dir, "los_sensitivity_outlier_diagnostics.csv"),
  row.names = FALSE
)

message("Built length-of-stay sensitivity-analysis data in: ", out_dir)
