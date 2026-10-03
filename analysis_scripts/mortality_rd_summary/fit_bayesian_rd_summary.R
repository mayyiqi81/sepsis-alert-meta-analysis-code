# Bayesian mortality risk-difference summary analysis.

library(brms)
library(dplyr)
library(posterior)

script_arg <- c(grep("^--file=", commandArgs(FALSE), value = TRUE), grep("\\.R$", commandArgs(FALSE), value = TRUE, ignore.case = TRUE))
script_dir <- if (length(script_arg) > 0) dirname(normalizePath(sub("^--file=", "", script_arg[1]), mustWork = FALSE)) else getwd()
project_root <- normalizePath(file.path(script_dir, "..", ".."), winslash = "/", mustWork = TRUE)
source(file.path(project_root, "analysis_scripts", "common", "project_paths.R"))
options(sepsis_project_root = set_project_root_from_script(project_root))
script_dir <- project_root
workspace_dir <- script_dir
source(file.path(workspace_dir, "analysis_scripts", "common", "bayes_prior_scale_spec.R"))
source(file.path(workspace_dir, "analysis_scripts", "common", "bayes_mortality_functions.R"))
source(file.path(workspace_dir, "analysis_scripts", "common", "mortality_rd_functions.R"))

input_file <- Sys.getenv("MORTALITY_RD_INPUT_FILE", unset = file.path(workspace_dir, "source_csv", "mortality.csv"))
rd_inputs_file <- file.path(workspace_dir, "source_csv", "mortality_rd_inputs.csv")
out_dir <- Sys.getenv("BAYES_RD_SUMMARY_OUT_DIR", unset = "computed_data/mortality_rd_summary")
if (!grepl("^([A-Za-z]:[\\\\/]|/)", out_dir)) out_dir <- file.path(workspace_dir, out_dir)
fits_dir <- file.path(out_dir, "fits")

rd_prior_spec <- get_bayes_metareg_prior_spec("rd_standalone_current", response_scale = "rd")
RD_PRIOR_INTERCEPT_MEAN <- rd_prior_spec$intercept_prior_mean
RD_PRIOR_INTERCEPT_SD <- rd_prior_spec$intercept_prior_sd
RD_PRIOR_TAU_DF <- rd_prior_spec$tau_prior_df
RD_PRIOR_TAU_LOCATION <- rd_prior_spec$tau_prior_location
RD_PRIOR_TAU_SCALE <- rd_prior_spec$tau_prior_scale

settings <- bayes_meta_settings_from_env()
rd_analysis_specs <- data.frame(
  group_id = c("rct_only", "nrs_only", "overall"),
  scenario = c("RCT only", "NRS only", "Overall"),
  design_filter = c("RCT", "NRS", NA_character_),
  stringsAsFactors = FALSE
)

if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
if (!dir.exists(fits_dir)) dir.create(fits_dir, recursive = TRUE)

format_rd_pp <- function(rd, lower, upper) {
  sprintf("%.1f [%.1f, %.1f]", 100 * rd, 100 * lower, 100 * upper)
}

stop_for_invalid_rd_variance <- function(dat) {
  invalid <- dat[!is.finite(dat$vi_rd) | dat$vi_rd <= 0, , drop = FALSE]
  if (nrow(invalid) == 0) return(invisible(TRUE))

  details <- apply(
    invalid[, c("study", "events_e", "n_e", "events_c", "n_c", "vi_rd"), drop = FALSE],
    1,
    function(row) {
      paste0(
        row[["study"]],
        " (events_e=", row[["events_e"]],
        ", n_e=", row[["n_e"]],
        ", events_c=", row[["events_c"]],
        ", n_c=", row[["n_c"]],
        ", vi_rd=", row[["vi_rd"]],
        ")"
      )
    }
  )
  stop(
    "RD fit stopped because at least one study has vi_rd <= 0 or non-finite vi_rd. ",
    "No continuity correction was applied. Studies: ",
    paste(details, collapse = "; "),
    call. = FALSE
  )
}

prepare_mortality_rd_effect_sizes <- function(input_file) {
  raw <- read_mortality_csv(input_file, lower_names = TRUE)
  names(raw) <- tolower(trimws(sub("^\ufeff", "", names(raw))))
  if (length(names(raw)) > 0 && grepl("study$", names(raw)[1], ignore.case = TRUE)) {
    names(raw)[1] <- "study"
  }
  raw <- raw[, names(raw) != "", drop = FALSE]

  required_cols <- c("study", "events_e", "n_e", "events_c", "n_c", "design")
  optional_cols <- intersect(c("arr", "arr_lower_95", "arr_upper_95"), names(raw))
  missing_cols <- setdiff(required_cols, names(raw))
  if (length(missing_cols) > 0) {
    stop("Missing required mortality columns for RD analysis: ", paste(missing_cols, collapse = ", "))
  }

  dat <- raw[, c(required_cols, optional_cols), drop = FALSE]
  if (!"arr" %in% names(dat)) dat$arr <- NA_real_
  if (!"arr_lower_95" %in% names(dat)) dat$arr_lower_95 <- NA_real_
  if (!"arr_upper_95" %in% names(dat)) dat$arr_upper_95 <- NA_real_
  dat$row_order <- seq_len(nrow(dat))
  dat$study <- clean_mortality_value(dat$study)
  dat$design <- toupper(clean_mortality_value(dat$design))
  dat$events_e <- safe_mortality_num(dat$events_e)
  dat$n_e <- safe_mortality_num(dat$n_e)
  dat$events_c <- safe_mortality_num(dat$events_c)
  dat$n_c <- safe_mortality_num(dat$n_c)
  dat$aRR <- safe_mortality_num(dat$arr)
  dat$aRR_lower_95 <- safe_mortality_num(dat$arr_lower_95)
  dat$aRR_upper_95 <- safe_mortality_num(dat$arr_upper_95)
  dat$arr <- NULL
  dat$arr_lower_95 <- NULL
  dat$arr_upper_95 <- NULL

  validate_mortality_2x2(dat, context = basename(input_file))
  validate_adjusted_rr_fields(
    dat,
    estimate_col = "aRR",
    lower_col = "aRR_lower_95",
    upper_col = "aRR_upper_95",
    context = basename(input_file)
  )

  dat <- dat[
    !is.na(dat$study) &
      !is.na(dat$events_e) & !is.na(dat$n_e) &
      !is.na(dat$events_c) & !is.na(dat$n_c) &
      dat$design %in% c("RCT", "NRS"),
    ,
    drop = FALSE
  ]

  invalid_counts <- dat[
    !is.finite(dat$n_e) | !is.finite(dat$n_c) |
      dat$n_e <= 0 | dat$n_c <= 0 |
      dat$events_e < 0 | dat$events_c < 0 |
      dat$events_e > dat$n_e | dat$events_c > dat$n_c,
    ,
    drop = FALSE
  ]
  if (nrow(invalid_counts) > 0) {
    details <- apply(
      invalid_counts[, c("study", "events_e", "n_e", "events_c", "n_c"), drop = FALSE],
      1,
      function(row) {
        paste0(
          row[["study"]],
          " (events_e=", row[["events_e"]],
          ", n_e=", row[["n_e"]],
          ", events_c=", row[["events_c"]],
          ", n_c=", row[["n_c"]],
          ")"
        )
      }
    )
    stop("RD analysis found invalid event counts or denominators: ", paste(details, collapse = "; "), call. = FALSE)
  }

  dat <- add_mortality_rd_effects(dat, rd_inputs_file)

  stop_for_invalid_rd_variance(dat)
  dat
}

fit_bayesian_rd_meta <- function(dat, settings) {
  model_data <- dat
  model_data$study <- factor(model_data$study, levels = unique(model_data$study))

  rd_priors <- c(
    set_prior(sprintf("normal(%s, %s)", RD_PRIOR_INTERCEPT_MEAN, RD_PRIOR_INTERCEPT_SD), class = "Intercept"),
    set_prior(sprintf("student_t(%s, %s, %s)", RD_PRIOR_TAU_DF, RD_PRIOR_TAU_LOCATION, RD_PRIOR_TAU_SCALE), class = "sd")
  )

  fit <- brm(
    formula = bf(yi | se(se_rd, sigma = FALSE) ~ 1 + (1 | study)),
    data = model_data,
    prior = rd_priors,
    family = gaussian(),
    iter = settings$iter,
    warmup = settings$warmup,
    chains = settings$chains,
    cores = settings$cores,
    seed = settings$seed,
    refresh = settings$refresh,
    control = list(
      adapt_delta = settings$adapt_delta,
      max_treedepth = settings$max_treedepth
    )
  )
  fit$data <- model_data
  model_metadata <- make_bayes_model_cache_metadata(
    data = model_data,
    vi_col = "vi_rd",
    response_scale = "rd",
    prior_spec_id = "rd_standalone_current",
    uncertainty_spec_id = "supplied_sampling_variance_only_rd",
    settings = settings,
    study_col = "study",
    yi_col = "yi",
    formula_text = "yi | se(se_rd, sigma = FALSE) ~ 1 + (1 | study)"
  )
  attach_bayes_model_metadata(fit, model_metadata)
}

summarise_bayesian_rd_fit <- function(fit, dat, group_id, scenario, prob = 0.95) {
  alpha <- (1 - prob) / 2
  draws <- as_draws_df(fit)
  rd_draws <- as.numeric(draws$b_Intercept)
  tau_draws <- get_tau_draws(fit)
  tau <- median(tau_draws, na.rm = TRUE)
  tau2 <- median(tau_draws^2, na.rm = TRUE)
  diagnostics <- extract_bayes_diagnostics(fit)

  random_weights <- 1 / (dat$vi_rd + tau2)
  rct_weight_percent <- ifelse(
    sum(random_weights, na.rm = TRUE) > 0,
    100 * sum(random_weights[dat$design == "RCT"], na.rm = TRUE) / sum(random_weights, na.rm = TRUE),
    NA_real_
  )

  rd <- median(rd_draws, na.rm = TRUE)
  lower <- as.numeric(stats::quantile(rd_draws, alpha, na.rm = TRUE))
  upper <- as.numeric(stats::quantile(rd_draws, 1 - alpha, na.rm = TRUE))

  data.frame(
    group_id = group_id,
    analysis = "4. Absolute-effect sensitivity analysis",
    scenario = scenario,
    k = length(unique(dat$study)),
    k_rct = sum(dat$design == "RCT"),
    k_nrs = sum(dat$design == "NRS"),
    RD = rd,
    Lower_95_CrI = lower,
    Upper_95_CrI = upper,
    RD_label_pp = format_rd_pp(rd, lower, upper),
    tau = tau,
    tau2 = tau2,
    tau_pp = 100 * tau,
    tau_pp_label = sprintf("%.1f", 100 * tau),
    rhat_max = diagnostics$rhat_max,
    ess_bulk_min = diagnostics$ess_bulk_min,
    ess_tail_min = diagnostics$ess_tail_min,
    divergent_transitions = diagnostics$divergent_transitions,
    max_treedepth = diagnostics$max_treedepth,
    max_treedepth_hits = diagnostics$max_treedepth_hits,
    rct_weight_percent = rct_weight_percent,
    rct_weight_label = ifelse(is.na(rct_weight_percent), "", sprintf("%.1f%%", rct_weight_percent)),
    method = "Bayesian normal-normal random-effects model for risk difference",
    note = "RD used event-count estimates except for Adams 2022 (reported adjusted RD) and Arabi 2025 (adjusted RD derived from aRR and control risk, with delta-method variance incorporating adjusted-log-RR and control-risk sampling uncertainty under an independence assumption); no additional NRS uncertainty was incorporated.",
    stringsAsFactors = FALSE
  )
}

fit_rd_analysis <- function(spec, rd_effect_sizes, settings, fits_dir) {
  dat <- if (is.na(spec$design_filter)) {
    rd_effect_sizes
  } else {
    rd_effect_sizes[rd_effect_sizes$design == spec$design_filter, , drop = FALSE]
  }
  if (nrow(dat) == 0) {
    stop("No mortality RD rows available for analysis: ", spec$scenario, call. = FALSE)
  }

  fit_file <- file.path(fits_dir, paste0("mortality_bayesian_rd_", spec$group_id, ".brmsfit.rds"))
  model_data <- dat
  model_data$study <- factor(model_data$study, levels = unique(model_data$study))
  expected_formula <- "yi | se(se_rd, sigma = FALSE) ~ 1 + (1 | study)"
  expected_metadata <- make_bayes_model_cache_metadata(
    data = model_data,
    vi_col = "vi_rd",
    response_scale = "rd",
    prior_spec_id = "rd_standalone_current",
    uncertainty_spec_id = "supplied_sampling_variance_only_rd",
    settings = settings,
    study_col = "study",
    yi_col = "yi",
    formula_text = expected_formula
  )
  force_refit <- identical(Sys.getenv("BAYES_FORCE_REFIT", unset = ""), "1")
  fit <- NULL
  if (file.exists(fit_file) && !force_refit) {
    cached_fit <- readRDS(fit_file)
    if (isTRUE(assert_cached_bayes_model_matches(cached_fit, expected_metadata))) {
      fit <- cached_fit
    } else {
      fit <- restore_cached_bayes_model_data(
        cached_fit,
        expected_data = model_data,
        expected_metadata = expected_metadata,
        cached_se_col = "se_rd",
        expected_formula = expected_formula
      )
      if (!is.null(fit)) saveRDS(fit, fit_file)
    }
  }
  if (is.null(fit)) {
    fit <- fit_bayesian_rd_meta_with_quality_gate(dat, settings)
    saveRDS(fit, fit_file)
  }

  summary <- summarise_bayesian_rd_fit(
    fit = fit,
    dat = dat,
    group_id = spec$group_id,
    scenario = spec$scenario
  )
  diagnostics <- summary[, c(
    "group_id",
    "analysis",
    "scenario",
    "rhat_max",
    "ess_bulk_min",
    "ess_tail_min",
    "divergent_transitions",
    "max_treedepth",
    "max_treedepth_hits"
  ), drop = FALSE]

  list(summary = summary, diagnostics = diagnostics)
}

write_rd_metadata <- function(path, settings, input_file, analysis_specs) {
  metadata <- data.frame(
    key = c(
      "created_at",
      "input_file",
      "analysis_count",
      "analysis_order",
      "adjusted_rd_overrides",
      "model_formula",
      "pooled_rd_prior",
      "tau_prior",
      "iter",
      "warmup",
      "chains",
      "seed",
      "adapt_delta",
      "max_treedepth",
      "additional_nrs_uncertainty"
    ),
    value = c(
      format(Sys.time(), "%Y-%m-%d %H:%M:%S %z"),
      normalizePath(input_file, mustWork = FALSE),
      as.character(nrow(analysis_specs)),
      paste(analysis_specs$scenario, collapse = "; "),
      "Adjusted RD inputs are read from source_csv/mortality_rd_inputs.csv; Arabi 2025 delta-method variance incorporates adjusted-log-RR and control-risk sampling uncertainty under an independence assumption",
      "bf(yi | se(se_rd, sigma = FALSE) ~ 1 + (1 | study)); family = gaussian()",
      sprintf("normal(%s, %s) on RD proportion scale", RD_PRIOR_INTERCEPT_MEAN, RD_PRIOR_INTERCEPT_SD),
      sprintf("student_t(%s, %s, %s) on RD proportion scale", RD_PRIOR_TAU_DF, RD_PRIOR_TAU_LOCATION, RD_PRIOR_TAU_SCALE),
      as.character(settings$iter),
      as.character(settings$warmup),
      as.character(settings$chains),
      as.character(settings$seed),
      as.character(settings$adapt_delta),
      as.character(settings$max_treedepth),
      "none; no additional NRS uncertainty variance was incorporated"
    ),
    stringsAsFactors = FALSE
  )
  utils::write.csv(metadata, path, row.names = FALSE)
  metadata
}

rd_effect_sizes <- prepare_mortality_rd_effect_sizes(input_file)
utils::write.csv(
  rd_effect_sizes,
  file.path(out_dir, "mortality_bayesian_rd_effect_sizes.csv"),
  row.names = FALSE
)

legacy_fit_file <- file.path(fits_dir, "mortality_bayesian_rd.brmsfit.rds")
if (file.exists(legacy_fit_file)) {
  unlink(legacy_fit_file)
}

rd_fit_results <- lapply(
  seq_len(nrow(rd_analysis_specs)),
  function(i) fit_rd_analysis(rd_analysis_specs[i, , drop = FALSE], rd_effect_sizes, settings, fits_dir)
)

rd_summary <- dplyr::bind_rows(lapply(rd_fit_results, `[[`, "summary"))
utils::write.csv(
  rd_summary,
  file.path(out_dir, "mortality_bayesian_rd_summary.csv"),
  row.names = FALSE
)

rd_diagnostics <- dplyr::bind_rows(lapply(rd_fit_results, `[[`, "diagnostics"))
utils::write.csv(
  rd_diagnostics,
  file.path(out_dir, "mortality_bayesian_rd_diagnostics.csv"),
  row.names = FALSE
)

invisible(write_rd_metadata(
  file.path(out_dir, "mortality_bayesian_rd_metadata.csv"),
  settings,
  input_file,
  rd_analysis_specs
))

message("Saved Bayesian mortality RD summary analysis to: ", out_dir)
