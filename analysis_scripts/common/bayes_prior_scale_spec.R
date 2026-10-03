# Single source of truth for Bayesian meta-regression units and prior metadata.
#
# This file deliberately contains no plotting code.  It is the single source
# of truth for model response scales, canonical moderator units, priors, and
# uncertainty metadata used by the Bayesian analysis layer.

BAYES_PRIOR_SCALE_SPEC_VERSION <- "bayes_prior_scale_canonical_20260904"
# Intercept-only models retain their existing model specification when the
# continuous-moderator parameterisation is refactored. This separate version
# prevents an unrelated meta-regression change from invalidating unchanged
# intercept-only fits.
BAYES_INTERCEPT_MODEL_SPEC_VERSION <- "bayes_intercept_model_v1"

bayes_metareg_moderator_specs <- data.frame(
  moderator_id = c(
    "process_completion",
    "antibiotic_time",
    "baseline_mortality"
  ),
  raw_unit = c(
    "percentage_point",
    "hour",
    "proportion"
  ),
  raw_unit_label = c(
    "percentage points",
    "hours",
    "0-1 proportion"
  ),
  centre_raw = c(0, 0, 0.10),
  reported_increment_raw = c(10, 1, 0.10),
  # Compatibility alias for existing display metadata consumers.  Model
  # transforms use reported_increment_raw explicitly.
  reported_increment = c(10, 1, 0.10),
  reported_increment_label = c(
    "per 10-percentage-point increase",
    "per 1-h increase in intervention-minus-control time",
    "per 10-percentage-point increase"
  ),
  model_unit_description = c(
    "one model unit = 10 percentage points of process completion",
    "one model unit = 1 hour of intervention-minus-control antibiotic time",
    "one model unit = 0.10 baseline-mortality proportion; intercept is at 10% baseline mortality"
  ),
  canonical_model_transform = c(
    "(raw - 0) / 10",
    "(raw - 0) / 1",
    "(raw - 0.10) / 0.10"
  ),
  display_multiplier = c(1, 1, 100),
  stringsAsFactors = FALSE
)

bayes_response_scale_specs <- data.frame(
  response_scale = c("log_rr", "rd", "log_or"),
  response_label = c("log risk ratio", "risk difference", "log odds ratio"),
  response_unit = c("log-RR", "proportion", "log-OR"),
  response_coding = c(
    "yi = log(RR); known sampling SE supplied through se(..., sigma = FALSE)",
    "yi = risk difference on the 0-1 proportion scale; known sampling SE supplied through se(..., sigma = FALSE)",
    "yi = log(OR); known sampling SE supplied through se(..., sigma = FALSE)"
  ),
  display_unit = c("mortality RR", "mortality RD in percentage points", "odds ratio"),
  stringsAsFactors = FALSE
)

bayes_model_parameterisation_specs <- data.frame(
  parameterisation_id = c(
    "log_rr_random_intercept_meta",
    "rd_random_intercept_meta",
    "log_or_random_intercept_meta",
    "rd_random_intercept_standalone"
  ),
  response_scale = c("log_rr", "rd", "log_or", "rd"),
  formula = c(
    "yi | se(.bayes_se, sigma = FALSE) ~ 1 + (1 | .bayes_study)",
    "yi | se(.bayes_se, sigma = FALSE) ~ 1 + (1 | .bayes_study)",
    "yi | se(.bayes_se, sigma = FALSE) ~ 1 + (1 | .bayes_study)",
    "yi | se(se_rd, sigma = FALSE) ~ 1 + (1 | study)"
  ),
  random_effect = c(
    "study-level random intercept",
    "study-level random intercept",
    "study-level random intercept",
    "study-level random intercept"
  ),
  residual_sigma = c(
    "fixed to supplied SE; sigma = FALSE",
    "fixed to supplied SE; sigma = FALSE",
    "fixed to supplied SE; sigma = FALSE",
    "fixed to supplied SE; sigma = FALSE"
  ),
  categorical_parameterisation = c(
    "Generic random-intercept normal-normal model; population term is model-specific (1 for intercept-only, bayes_mod for continuous meta-regression)",
    "Generic random-intercept normal-normal model; population term is model-specific (1 for intercept-only, bayes_mod for continuous meta-regression)",
    "Generic random-intercept normal-normal model; population term is model-specific (1 for intercept-only, bayes_mod for continuous meta-regression)",
    "Generic random-intercept normal-normal model; population term is model-specific (1 for intercept-only, bayes_mod for continuous meta-regression)"
  ),
  model_role_description = c(
    "Intercept-only: population term = 1; continuous meta-regression: population term = bayes_mod",
    "Intercept-only: population term = 1; continuous meta-regression: population term = bayes_mod",
    "Intercept-only: population term = 1; continuous meta-regression: population term = bayes_mod",
    "Intercept-only: population term = 1; continuous meta-regression: population term = bayes_mod"
  ),
  stringsAsFactors = FALSE
)

bayes_uncertainty_specs <- data.frame(
  uncertainty_spec_id = c(
    "primary_nrs_log_rr_bias_sd_0.20",
    "bias_unadjusted_no_nrs_additional_uncertainty",
    "mortality_rd_no_nrs_additional_uncertainty",
    "strict_same_cohort_primary_nrs_bias_sd_0.20",
    "supplied_sampling_variance_only_log_rr",
    "supplied_sampling_variance_only_rd",
    "supplied_sampling_variance_only_log_or",
    "sensitivity_defined_variance_log_rr"
  ),
  response_scale = c("log_rr", "log_rr", "rd", "log_rr", "log_rr", "rd", "log_or", "log_rr"),
  description = c(
    "NRS vi = vi_original + 0.20^2; RCT vi unchanged",
    "No additional NRS variance; vi = vi_original",
    "No additional NRS uncertainty on the RD scale",
    "Same additional uncertainty specification as the primary mortality analysis",
    "Supplied log-RR sampling variance only; no additional uncertainty",
    "Supplied RD sampling variance only; no additional uncertainty",
    "Supplied log-OR sampling variance only; no additional uncertainty",
    "Variance supplied by the selected sensitivity scenario"
  ),
  stringsAsFactors = FALSE
)

# Canonical and preserved existing prior specifications. Intercept-only
# *_current IDs retain their existing numerical meaning. Legacy metareg IDs
# document the pre-canonical raw-unit fits and are not permitted for active
# fitting. Continuous moderator fits use the explicit *_metareg_canonical IDs.
bayes_metareg_prior_specs <- data.frame(
  prior_spec_id = c(
    "log_rr_shared_current",
    "log_rr_metareg_legacy_current",
    "log_rr_metareg_canonical",
    "rd_metareg_shared_current",
    "rd_metareg_legacy_current",
    "rd_metareg_canonical",
    "rd_standalone_current",
    "log_or_shared_current",
    "log_or_metareg_legacy_current",
    "log_or_metareg_canonical"
  ),
  response_scale = c("log_rr", "log_rr", "log_rr", "rd", "rd", "rd", "rd", "log_or", "log_or", "log_or"),
  intercept_prior = c(
    "normal(0, 1)", "normal(0, 1)", "normal(0, 1)",
    "normal(0, 1)", "normal(0, 1)", "normal(0, 0.10)",
    "normal(0, 0.10)", "normal(0, 1)", "normal(0, 1)",
    "normal(0, 1)"
  ),
  slope_prior = c(
    NA_character_, "normal(0, 1)", "normal(0, 0.5)",
    "normal(0, 1)", "normal(0, 1)", "normal(0, 0.10)",
    NA_character_, NA_character_, "normal(0, 1)", "normal(0, 0.5)"
  ),
  tau_prior = c(
    "student_t(3, 0, 0.5)", "student_t(3, 0, 0.5)", "student_t(3, 0, 0.5)",
    "student_t(3, 0, 0.5)", "student_t(3, 0, 0.5)", "student_t(3, 0, 0.05)",
    "student_t(3, 0, 0.05)", "student_t(3, 0, 0.5)", "student_t(3, 0, 0.5)",
    "student_t(3, 0, 0.5)"
  ),
  intercept_prior_mean = rep(0, 10),
  slope_prior_mean = c(NA_real_, 0, 0, 0, 0, 0, NA_real_, NA_real_, 0, 0),
  tau_prior_df = rep(3, 10),
  tau_prior_location = rep(0, 10),
  intercept_prior_sd = c(1, 1, 1, 1, 1, 0.10, 0.10, 1, 1, 1),
  slope_prior_sd = c(NA_real_, 1, 0.5, 1, 1, 0.10, NA_real_, NA_real_, 1, 0.5),
  tau_prior_scale = c(0.5, 0.5, 0.5, 0.5, 0.5, 0.05, 0.05, 0.5, 0.5, 0.5),
  prior_family = c(
    "relative_effect_intercept_only",
    "relative_effect_metareg_legacy",
    "relative_effect_metareg",
    "mortality_rd_metareg_legacy",
    "mortality_rd_metareg_legacy",
    "mortality_rd_metareg",
    "mortality_rd_standalone",
    "relative_effect_intercept_only",
    "relative_effect_metareg_legacy",
    "relative_effect_metareg"
  ),
  canonical = c(FALSE, FALSE, TRUE, FALSE, FALSE, TRUE, TRUE, FALSE, FALSE, TRUE),
  stringsAsFactors = FALSE
)

get_bayes_metareg_moderator_spec <- function(moderator_id) {
  if (length(moderator_id) != 1L || is.na(moderator_id) || !nzchar(moderator_id)) {
    stop("moderator_id must be one non-empty scalar string.")
  }
  rows <- bayes_metareg_moderator_specs[
    bayes_metareg_moderator_specs$moderator_id == moderator_id,
    ,
    drop = FALSE
  ]
  if (nrow(rows) != 1L) stop("Unknown Bayesian meta-regression moderator_id: ", moderator_id)
  as.list(rows[1, , drop = FALSE])
}

get_bayes_response_scale_spec <- function(response_scale) {
  if (length(response_scale) != 1L || is.na(response_scale) || !nzchar(response_scale)) {
    stop("response_scale must be one non-empty scalar string.")
  }
  rows <- bayes_response_scale_specs[
    bayes_response_scale_specs$response_scale == response_scale,
    ,
    drop = FALSE
  ]
  if (nrow(rows) != 1L) stop("Unknown Bayesian response_scale: ", response_scale)
  as.list(rows[1, , drop = FALSE])
}

get_bayes_parameterisation_spec <- function(parameterisation_id) {
  rows <- bayes_model_parameterisation_specs[
    bayes_model_parameterisation_specs$parameterisation_id == parameterisation_id,
    ,
    drop = FALSE
  ]
  if (nrow(rows) != 1L) stop("Unknown Bayesian parameterisation_id: ", parameterisation_id)
  as.list(rows[1, , drop = FALSE])
}

get_bayes_uncertainty_spec <- function(uncertainty_spec_id) {
  rows <- bayes_uncertainty_specs[
    bayes_uncertainty_specs$uncertainty_spec_id == uncertainty_spec_id,
    ,
    drop = FALSE
  ]
  if (nrow(rows) != 1L) stop("Unknown Bayesian uncertainty_spec_id: ", uncertainty_spec_id)
  as.list(rows[1, , drop = FALSE])
}

get_bayes_metareg_prior_spec <- function(prior_spec_id, response_scale = NULL) {
  if (length(prior_spec_id) != 1L || is.na(prior_spec_id) || !nzchar(prior_spec_id)) {
    stop("prior_spec_id must be one non-empty scalar string.")
  }
  rows <- bayes_metareg_prior_specs[
    bayes_metareg_prior_specs$prior_spec_id == prior_spec_id,
    ,
    drop = FALSE
  ]
  if (nrow(rows) != 1L) stop("Unknown Bayesian meta-regression prior_spec_id: ", prior_spec_id)
  if (!is.null(response_scale) && !identical(rows$response_scale[[1]], response_scale)) {
    stop("Prior specification ", prior_spec_id, " is not valid for response_scale ", response_scale)
  }
  as.list(rows[1, , drop = FALSE])
}

get_bayes_metareg_prior_spec_variant <- function(prior_spec_id,
                                                 response_scale,
                                                 prior_variant = "canonical") {
  if (length(prior_variant) != 1L || is.na(prior_variant) || !nzchar(prior_variant)) {
    stop("prior_variant must be one non-empty scalar string.")
  }
  base <- get_bayes_metareg_prior_spec(prior_spec_id, response_scale = response_scale)
  if (identical(prior_variant, "canonical")) return(base)
  if (!identical(prior_variant, "wide")) {
    stop("Unknown Bayesian prior variant: ", prior_variant)
  }
  if (!isTRUE(base$canonical) || is.na(base$slope_prior_sd)) {
    stop("Wide-prior sensitivity is defined only for canonical meta-regression specifications.")
  }
  if (response_scale %in% c("log_rr", "log_or")) {
    base$slope_prior_sd <- 1
    base$slope_prior <- "normal(0, 1)"
  } else if (identical(response_scale, "rd")) {
    base$intercept_prior_sd <- 0.20
    base$intercept_prior <- "normal(0, 0.20)"
    base$slope_prior_sd <- 0.20
    base$slope_prior <- "normal(0, 0.20)"
    base$tau_prior_scale <- 0.10
    base$tau_prior <- "student_t(3, 0, 0.10)"
  }
  base$prior_variant <- prior_variant
  base
}

scale_bayes_moderator <- function(x, moderator_id) {
  transform_bayes_moderator(x, moderator_id)
}

transform_bayes_moderator <- function(x, moderator_id) {
  spec <- get_bayes_metareg_moderator_spec(moderator_id)
  raw <- as.numeric(x)
  (raw - spec$centre_raw) / spec$reported_increment_raw
}

inverse_transform_bayes_moderator <- function(x, moderator_id) {
  spec <- get_bayes_metareg_moderator_spec(moderator_id)
  model <- as.numeric(x)
  spec$centre_raw + model * spec$reported_increment_raw
}

get_bayes_reported_increment <- function(moderator_id) {
  get_bayes_metareg_moderator_spec(moderator_id)$reported_increment_raw
}

get_bayes_reported_increment_label <- function(moderator_id,
                                               capitalise = FALSE) {
  label <- get_bayes_metareg_moderator_spec(moderator_id)$reported_increment_label
  if (isTRUE(capitalise)) {
    return(paste0(toupper(substr(label, 1, 1)), substr(label, 2, nchar(label))))
  }
  label
}

display_bayes_moderator <- function(x, moderator_id) {
  spec <- get_bayes_metareg_moderator_spec(moderator_id)
  as.numeric(x) * spec$display_multiplier
}

display_bayes_proportion_as_percent <- function(x) {
  # Display-only conversion.  This helper never changes a model response,
  # moderator, fit, cache, or statistical estimate.
  as.numeric(x) * 100
}

reported_increment_in_current_model_units <- function(moderator_id,
                                                       model_uses_canonical_units = TRUE) {
  spec <- get_bayes_metareg_moderator_spec(moderator_id)
  if (isTRUE(model_uses_canonical_units)) 1 else spec$reported_increment_raw
}

effective_slope_prior_sd <- function(prior_sd,
                                     moderator_id,
                                     model_uses_canonical_units = FALSE) {
  if (length(prior_sd) != 1L || !is.finite(prior_sd) || prior_sd < 0) {
    stop("prior_sd must be one finite non-negative number.")
  }
  prior_sd * reported_increment_in_current_model_units(
    moderator_id,
    model_uses_canonical_units = model_uses_canonical_units
  )
}

bayes_hash_object <- function(x) {
  if (requireNamespace("digest", quietly = TRUE)) {
    return(digest::digest(x, algo = "xxhash64"))
  }
  raw <- serialize(x, NULL, version = 2)
  if (length(raw) == 0L) return("empty")
  paste0("fallback-", length(raw), "-", sum(as.numeric(raw) * seq_along(raw)) %% 1000000007)
}

bayes_sampler_metadata <- function(settings) {
  required <- c("iter", "warmup", "chains", "seed", "adapt_delta", "max_treedepth")
  if (!is.list(settings) || !all(required %in% names(settings))) {
    stop("Bayesian sampler settings lack: ", paste(setdiff(required, names(settings)), collapse = ", "))
  }
  list(
    sampler_iter = as.integer(settings$iter),
    sampler_warmup = as.integer(settings$warmup),
    sampler_chains = as.integer(settings$chains),
    sampler_seed = as.integer(settings$seed),
    sampler_adapt_delta = as.numeric(settings$adapt_delta),
    sampler_max_treedepth = as.integer(settings$max_treedepth),
    sampler_policy_id = "initial_then_approved_retry_v1"
  )
}

bayes_sampler_metadata_names <- function() {
  c(
    "sampler_iter", "sampler_warmup", "sampler_chains", "sampler_seed",
    "sampler_adapt_delta", "sampler_max_treedepth", "sampler_policy_id"
  )
}

bayes_approved_retry_sampler_metadata <- function(expected) {
  retry <- expected
  retry$sampler_iter <- 8000L
  retry$sampler_warmup <- 4000L
  retry$sampler_adapt_delta <- 0.999
  retry$sampler_max_treedepth <- 18L
  retry
}

bayes_validated_legacy_high_control_sampler_metadata <- function(expected) {
  legacy <- expected
  legacy$sampler_adapt_delta <- 0.999
  legacy$sampler_max_treedepth <- 18L
  legacy
}

bayes_sampler_metadata_matches_policy <- function(actual, expected) {
  fields <- bayes_sampler_metadata_names()
  if (!all(fields %in% names(actual)) || !all(fields %in% names(expected))) return(FALSE)
  actual_values <- actual[fields]
  expected_values <- expected[fields]
  retry_values <- bayes_approved_retry_sampler_metadata(expected_values)
  legacy_values <- bayes_validated_legacy_high_control_sampler_metadata(expected_values)
  exact <- all(vapply(fields, function(name) identical(actual_values[[name]], expected_values[[name]]), logical(1)))
  retry <- all(vapply(fields, function(name) identical(actual_values[[name]], retry_values[[name]]), logical(1)))
  legacy <- all(vapply(fields, function(name) identical(actual_values[[name]], legacy_values[[name]]), logical(1)))
  exact || retry || legacy
}

cached_bayes_quality_gate_is_valid <- function(fit, actual_metadata) {
  gate <- attr(fit, "bayes_quality_gate")
  if (!is.list(gate) || !isTRUE(gate$passed) || !is.list(gate$settings) || !is.data.frame(gate$diagnostics)) {
    return(FALSE)
  }
  gate_sampler <- tryCatch(bayes_sampler_metadata(gate$settings), error = function(e) NULL)
  if (is.null(gate_sampler)) return(FALSE)
  fields <- bayes_sampler_metadata_names()
  if (!all(vapply(fields, function(name) identical(gate_sampler[[name]], actual_metadata[[name]]), logical(1)))) {
    return(FALSE)
  }
  required_diagnostics <- c(
    "rhat_max", "ess_bulk_min", "ess_tail_min", "divergent_transitions",
    "max_treedepth", "max_treedepth_hits"
  )
  if (nrow(gate$diagnostics) != 1L || !all(required_diagnostics %in% names(gate$diagnostics))) return(FALSE)
  diagnostics <- gate$diagnostics[1, required_diagnostics, drop = FALSE]
  isTRUE(
    is.finite(diagnostics$rhat_max) && diagnostics$rhat_max <= 1.01 &&
      is.finite(diagnostics$ess_bulk_min) && diagnostics$ess_bulk_min >= 400 &&
      is.finite(diagnostics$ess_tail_min) && diagnostics$ess_tail_min >= 400 &&
      is.finite(diagnostics$divergent_transitions) && diagnostics$divergent_transitions == 0 &&
      is.finite(diagnostics$max_treedepth) &&
      diagnostics$max_treedepth == actual_metadata$sampler_max_treedepth &&
      is.finite(diagnostics$max_treedepth_hits) && diagnostics$max_treedepth_hits == 0
  )
}

make_bayes_metareg_cache_metadata <- function(data,
                                              moderator_col,
                                              vi_col,
                                              moderator_id,
                                              response_scale,
                                              prior_spec_id,
                                              uncertainty_spec_id,
                                              settings,
                                              prior_variant = "canonical",
                                              study_col = "study") {
  moderator_spec <- get_bayes_metareg_moderator_spec(moderator_id)
  response_spec <- get_bayes_response_scale_spec(response_scale)
  prior_spec <- get_bayes_metareg_prior_spec_variant(
    prior_spec_id,
    response_scale = response_scale,
    prior_variant = prior_variant
  )
  uncertainty_spec <- get_bayes_uncertainty_spec(uncertainty_spec_id)
  if (!identical(uncertainty_spec$response_scale[[1]], response_scale)) {
    stop("Uncertainty specification ", uncertainty_spec_id, " is not valid for response_scale ", response_scale)
  }
  parameterisation_id <- if (identical(response_scale, "rd")) {
    "rd_random_intercept_meta"
  } else if (identical(response_scale, "log_rr")) {
    "log_rr_random_intercept_meta"
  } else {
    "log_or_random_intercept_meta"
  }
  required <- c(study_col, "yi", vi_col, moderator_col)
  missing_cols <- setdiff(required, names(data))
  if (length(missing_cols) > 0L) {
    stop("Cannot create Bayesian meta-regression cache metadata; missing: ", paste(missing_cols, collapse = ", "))
  }

  raw_moderator <- as.numeric(data[[moderator_col]])
  model_moderator <- scale_bayes_moderator(raw_moderator, moderator_id)
  if (any(!is.finite(raw_moderator)) || any(!is.finite(model_moderator))) {
    stop("Bayesian meta-regression moderator must be finite: ", moderator_col)
  }

  c(list(
    model_spec_version = BAYES_PRIOR_SCALE_SPEC_VERSION,
    response_scale = response_scale,
    response_unit = response_spec$response_unit,
    response_coding = response_spec$response_coding,
    parameterisation_id = parameterisation_id,
    moderator_id = moderator_id,
    antibiotic_time_direction = if (identical(moderator_id, "antibiotic_time")) "intervention_minus_control" else NA_character_,
    moderator_raw_unit = moderator_spec$raw_unit,
    moderator_centre_raw = moderator_spec$centre_raw,
    moderator_increment_raw = moderator_spec$reported_increment_raw,
    moderator_reported_increment = moderator_spec$reported_increment_raw,
    moderator_model_unit_description = moderator_spec$model_unit_description,
    moderator_transform = moderator_spec$canonical_model_transform,
    prior_spec_id = prior_spec_id,
    prior_variant = prior_variant,
    uncertainty_spec_id = uncertainty_spec_id,
    uncertainty_specification = uncertainty_spec$description,
    study_col = study_col,
    vi_col = vi_col,
    slope_prior = prior_spec$slope_prior,
    intercept_prior = prior_spec$intercept_prior,
    tau_prior = prior_spec$tau_prior,
    study_hash = bayes_hash_object(as.character(data[[study_col]])),
    yi_hash = bayes_hash_object(as.numeric(data$yi)),
    vi_hash = bayes_hash_object(as.numeric(data[[vi_col]])),
    moderator_raw_hash = bayes_hash_object(raw_moderator),
    moderator_model_hash = bayes_hash_object(model_moderator),
    # Compatibility aliases retained for older cache readers.
    moderator_hash = bayes_hash_object(raw_moderator),
    bayes_mod_hash = bayes_hash_object(model_moderator)
  ), bayes_sampler_metadata(settings))
}

make_bayes_model_cache_metadata <- function(data,
                                            vi_col,
                                            response_scale,
                                            prior_spec_id,
                                             uncertainty_spec_id,
                                             settings,
                                             study_col = "study",
                                            yi_col = "yi",
                                            formula_text = NA_character_) {
  response_spec <- get_bayes_response_scale_spec(response_scale)
  prior_spec <- get_bayes_metareg_prior_spec(prior_spec_id, response_scale = response_scale)
  uncertainty_spec <- get_bayes_uncertainty_spec(uncertainty_spec_id)
  if (!identical(uncertainty_spec$response_scale[[1]], response_scale)) {
    stop("Uncertainty specification ", uncertainty_spec_id, " is not valid for response_scale ", response_scale)
  }
  parameterisation_id <- if (identical(prior_spec_id, "rd_standalone_current")) {
    "rd_random_intercept_standalone"
  } else if (identical(response_scale, "rd")) {
    "rd_random_intercept_meta"
  } else if (identical(response_scale, "log_rr")) {
    "log_rr_random_intercept_meta"
  } else {
    "log_or_random_intercept_meta"
  }
  required <- c(study_col, yi_col, vi_col)
  missing_cols <- setdiff(required, names(data))
  if (length(missing_cols) > 0L) {
    stop("Cannot create Bayesian model cache metadata; missing: ", paste(missing_cols, collapse = ", "))
  }
  c(list(
    model_spec_version = BAYES_INTERCEPT_MODEL_SPEC_VERSION,
    response_scale = response_scale,
    response_unit = response_spec$response_unit,
    response_coding = response_spec$response_coding,
    parameterisation_id = parameterisation_id,
    prior_spec_id = prior_spec_id,
    intercept_prior = prior_spec$intercept_prior,
    slope_or_contrast_prior = prior_spec$slope_prior,
    tau_prior = prior_spec$tau_prior,
    residual_sigma = "fixed to supplied SE; sigma = FALSE",
    uncertainty_spec_id = uncertainty_spec_id,
    uncertainty_specification = uncertainty_spec$description,
    study_col = study_col,
    yi_col = yi_col,
    vi_col = vi_col,
    formula_text = formula_text,
    study_hash = bayes_hash_object(as.character(data[[study_col]])),
    yi_hash = bayes_hash_object(as.numeric(data[[yi_col]])),
    vi_hash = bayes_hash_object(as.numeric(data[[vi_col]]))
  ), bayes_sampler_metadata(settings))
}

attach_bayes_model_metadata <- function(fit, metadata) {
  attr(fit, "bayes_model_metadata") <- metadata
  fit
}

assert_cached_bayes_model_matches <- function(fit, expected_metadata) {
  actual <- attr(fit, "bayes_model_metadata")
  if (!is.list(actual)) return(FALSE)
  required <- c(
    "model_spec_version", "response_scale", "response_coding", "parameterisation_id",
    "prior_spec_id", "intercept_prior", "slope_or_contrast_prior", "tau_prior",
    "residual_sigma", "uncertainty_spec_id", "uncertainty_specification",
    "study_col", "yi_col", "vi_col", "study_hash", "yi_hash", "vi_hash"
  )
  if (!all(required %in% names(actual))) return(FALSE)
  if (!all(vapply(required, function(name) identical(actual[[name]], expected_metadata[[name]]), logical(1)))) {
    return(FALSE)
  }
  if (!bayes_sampler_metadata_matches_policy(actual, expected_metadata)) return(FALSE)
  if (!cached_bayes_quality_gate_is_valid(fit, actual)) return(FALSE)
  cached_data <- tryCatch(fit[["data"]], error = function(e) NULL)
  if (!is.data.frame(cached_data)) return(FALSE)
  required_data <- c(expected_metadata$study_col, expected_metadata$yi_col, expected_metadata$vi_col)
  if (!all(required_data %in% names(cached_data))) return(FALSE)
  if (!identical(
    bayes_hash_object(as.character(cached_data[[expected_metadata$study_col]])),
    expected_metadata$study_hash
  )) return(FALSE)
  if (!identical(
    bayes_hash_object(as.numeric(cached_data[[expected_metadata$yi_col]])),
    expected_metadata$yi_hash
  )) return(FALSE)
  if (!identical(
    bayes_hash_object(as.numeric(cached_data[[expected_metadata$vi_col]])),
    expected_metadata$vi_hash
  )) return(FALSE)
  TRUE
}

attach_bayes_metareg_metadata <- function(fit, metadata) {
  attr(fit, "bayes_metareg_metadata") <- metadata
  fit
}

assert_cached_metareg_fit_matches <- function(fit, expected_metadata) {
  actual <- attr(fit, "bayes_metareg_metadata")
  if (!is.list(actual)) return(FALSE)
  required <- c(
    "model_spec_version", "response_scale", "response_coding", "parameterisation_id",
    "moderator_id",
    "moderator_raw_unit", "moderator_centre_raw", "moderator_increment_raw",
    "moderator_reported_increment", "moderator_model_unit_description",
    "moderator_transform", "prior_spec_id", "prior_variant", "uncertainty_spec_id",
    "uncertainty_specification", "slope_prior",
    "intercept_prior", "tau_prior", "study_hash", "yi_hash",
    "vi_hash", "moderator_raw_hash", "moderator_model_hash",
    "moderator_hash", "bayes_mod_hash", "study_col", "vi_col"
  )
  if (!all(required %in% names(actual))) return(FALSE)
  if (!all(vapply(required, function(name) identical(actual[[name]], expected_metadata[[name]]), logical(1)))) {
    return(FALSE)
  }
  if (!bayes_sampler_metadata_matches_policy(actual, expected_metadata)) return(FALSE)
  if (!cached_bayes_quality_gate_is_valid(fit, actual)) return(FALSE)

  cached_data <- tryCatch(fit[["data"]], error = function(e) NULL)
  if (!is.data.frame(cached_data)) return(FALSE)
  required_data <- c(
    "bayes_mod_raw", "bayes_mod", "moderator_id", "response_scale",
    "prior_spec_id", "uncertainty_spec_id", expected_metadata$study_col,
    "yi", expected_metadata$vi_col
  )
  if (!all(required_data %in% names(cached_data))) return(FALSE)
  if (!identical(as.character(cached_data$moderator_id), rep(expected_metadata$moderator_id, nrow(cached_data)))) return(FALSE)
  if (!identical(as.character(cached_data$response_scale), rep(expected_metadata$response_scale, nrow(cached_data)))) return(FALSE)
  if (!identical(as.character(cached_data$prior_spec_id), rep(expected_metadata$prior_spec_id, nrow(cached_data)))) return(FALSE)
  if (!identical(as.character(cached_data$uncertainty_spec_id), rep(expected_metadata$uncertainty_spec_id, nrow(cached_data)))) return(FALSE)
  if (!identical(
    bayes_hash_object(as.character(cached_data[[expected_metadata$study_col]])),
    expected_metadata$study_hash
  )) return(FALSE)
  if (!identical(bayes_hash_object(as.numeric(cached_data$yi)), expected_metadata$yi_hash)) return(FALSE)
  if (!identical(
    bayes_hash_object(as.numeric(cached_data[[expected_metadata$vi_col]])),
    expected_metadata$vi_hash
  )) return(FALSE)
  if (!identical(bayes_hash_object(as.numeric(cached_data$bayes_mod_raw)), expected_metadata$moderator_raw_hash)) return(FALSE)
  if (!identical(bayes_hash_object(as.numeric(cached_data$bayes_mod)), expected_metadata$moderator_model_hash)) return(FALSE)
  TRUE
}

validate_bayes_metareg_spec <- function(moderator_id,
                                        response_scale,
                                        prior_spec_id,
                                        uncertainty_spec_id) {
  get_bayes_metareg_moderator_spec(moderator_id)
  get_bayes_response_scale_spec(response_scale)
  prior_spec <- get_bayes_metareg_prior_spec(prior_spec_id, response_scale = response_scale)
  if (!isTRUE(prior_spec$canonical) || is.na(prior_spec$slope_prior_sd)) {
    stop(
      "Active Bayesian meta-regression requires a canonical slope prior; received prior_spec_id=",
      prior_spec_id
    )
  }
  uncertainty_spec <- get_bayes_uncertainty_spec(uncertainty_spec_id)
  if (!identical(uncertainty_spec$response_scale[[1]], response_scale)) {
    stop("Uncertainty specification ", uncertainty_spec_id, " is not valid for response_scale ", response_scale)
  }
  invisible(TRUE)
}
