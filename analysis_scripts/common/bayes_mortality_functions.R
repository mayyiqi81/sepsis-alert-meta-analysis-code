# Shared Bayesian mortality meta-analysis helpers.

library(brms)
library(metafor)
library(posterior)

if (!exists("get_bayes_metareg_moderator_spec", mode = "function")) {
  option_base <- getOption("sepsis_code_dir", NULL)
  workspace_base <- if (exists("workspace_dir", inherits = TRUE)) {
    get("workspace_dir", inherits = TRUE)
  } else {
    NULL
  }
  spec_candidates <- unique(c(
    if (!is.null(option_base)) file.path(option_base, "analysis_scripts", "common", "bayes_prior_scale_spec.R") else character(),
    if (!is.null(workspace_base)) file.path(workspace_base, "analysis_scripts", "common", "bayes_prior_scale_spec.R") else character(),
    file.path(getwd(), "analysis_scripts", "common", "bayes_prior_scale_spec.R"),
    file.path(dirname(getwd()), "analysis_scripts", "common", "bayes_prior_scale_spec.R"),
    file.path(dirname(dirname(getwd())), "analysis_scripts", "common", "bayes_prior_scale_spec.R"),
    if (!is.null(option_base)) file.path(option_base, "bayes_prior_scale_spec.R") else character(),
    if (!is.null(workspace_base)) file.path(workspace_base, "bayes_prior_scale_spec.R") else character()
  ))
  spec_path <- spec_candidates[file.exists(spec_candidates)][1]
  if (is.na(spec_path)) {
    stop("Cannot locate bayes_prior_scale_spec.R before loading mortality helpers.")
  }
  source(spec_path, local = .GlobalEnv)
}

clean_mortality_value <- function(x) {
  x <- trimws(as.character(x))
  x[x %in% c("", "/", "NA", "Na", "na", "N/A")] <- NA_character_
  x
}

safe_mortality_num <- function(x) {
  suppressWarnings(as.numeric(gsub("%", "", trimws(as.character(x)), fixed = TRUE)))
}

validate_mortality_2x2 <- function(data,
                                   events_e_col = "events_e",
                                   n_e_col = "n_e",
                                   events_c_col = "events_c",
                                   n_c_col = "n_c",
                                   context = "mortality input") {
  study_col <- c("study", "Study")[c("study", "Study") %in% names(data)][1]
  design_col <- c("design", "Design")[c("design", "Design") %in% names(data)][1]
  if (is.na(study_col) || is.na(design_col)) {
    stop("Mortality validation requires study and design columns in ", context)
  }
  study <- trimws(as.character(data[[study_col]]))
  design <- toupper(trimws(as.character(data[[design_col]])))
  invalid_identity <- is.na(study) | !nzchar(study) | is.na(design) | !design %in% c("RCT", "NRS")
  if (any(invalid_identity)) {
    stop("Invalid study/design identity in ", context, " at rows: ", paste(which(invalid_identity), collapse = ", "))
  }
  if (anyDuplicated(study) > 0L) {
    stop("Mortality input must contain one row per study in ", context, "; duplicate study IDs found.")
  }
  required <- c(events_e_col, n_e_col, events_c_col, n_c_col)
  missing <- setdiff(required, names(data))
  if (length(missing) > 0L) {
    stop("Missing 2x2 count columns for ", context, ": ", paste(missing, collapse = ", "))
  }
  values <- lapply(data[required], function(x) suppressWarnings(as.numeric(x)))
  names(values) <- required
  finite_ok <- Reduce(`&`, lapply(values, is.finite))
  bounds_ok <- values[[events_e_col]] >= 0 & values[[events_e_col]] <= values[[n_e_col]] &
    values[[events_c_col]] >= 0 & values[[events_c_col]] <= values[[n_c_col]] &
    values[[n_e_col]] > 0 & values[[n_c_col]] > 0
  bad <- which(!(finite_ok & bounds_ok))
  if (length(bad) > 0L) {
    study_label <- if ("study" %in% names(data)) as.character(data$study[bad]) else as.character(bad)
    stop(
      "Invalid mortality 2x2 counts in ", context, " at rows: ",
      paste(paste0(bad, "[", study_label, "]"), collapse = ", ")
    )
  }
  invisible(TRUE)
}

validate_adjusted_rr_fields <- function(data,
                                        estimate_col = "aRR",
                                        lower_col = "aRR_lower_95",
                                        upper_col = "aRR_upper_95",
                                        context = "mortality input") {
  fields <- c(estimate_col, lower_col, upper_col)
  present <- fields %in% names(data)
  if (!any(present)) return(invisible(TRUE))
  if (!all(present)) {
    stop("Adjusted RR fields must be complete or absent in ", context, ": ", paste(fields, collapse = ", "))
  }
  values <- lapply(data[fields], function(x) suppressWarnings(as.numeric(x)))
  names(values) <- fields
  all_missing <- is.na(values[[estimate_col]]) & is.na(values[[lower_col]]) & is.na(values[[upper_col]])
  all_present <- !is.na(values[[estimate_col]]) & !is.na(values[[lower_col]]) & !is.na(values[[upper_col]])
  valid_present <- all_present &
    is.finite(values[[estimate_col]]) & is.finite(values[[lower_col]]) & is.finite(values[[upper_col]]) &
    values[[lower_col]] > 0 & values[[estimate_col]] > 0 & values[[upper_col]] > 0 &
    values[[lower_col]] <= values[[estimate_col]] & values[[estimate_col]] <= values[[upper_col]]
  bad <- which(!(all_missing | valid_present))
  if (length(bad) > 0L) {
    study_label <- if ("study" %in% names(data)) as.character(data$study[bad]) else as.character(bad)
    stop(
      "Adjusted RR triplet is partial, non-finite, or violates lower <= estimate <= upper in ",
      context, " at rows: ", paste(paste0(bad, "[", study_label, "]"), collapse = ", ")
    )
  }
  invisible(TRUE)
}

read_mortality_csv <- function(input_file, lower_names = FALSE) {
  raw <- utils::read.csv(
    input_file,
    stringsAsFactors = FALSE,
    check.names = FALSE,
    encoding = "UTF-8"
  )
  names(raw) <- trimws(sub("^\ufeff", "", names(raw)))
  raw <- raw[, names(raw) != "", drop = FALSE]
  if (lower_names) {
    names(raw) <- tolower(names(raw))
    if (length(names(raw)) > 0 && grepl("study$", names(raw)[1], ignore.case = TRUE)) {
      names(raw)[1] <- "study"
    }
  } else if (length(names(raw)) > 0 && grepl("study$", names(raw)[1], ignore.case = TRUE)) {
    names(raw)[1] <- "Study"
  }
  raw
}

prepare_mortality_effect_sizes <- function(input_file,
                                           nrs_bias_sd = 0.20,
                                           use_adjusted_rr = TRUE) {
  raw <- if (is.data.frame(input_file)) {
    input_file
  } else {
    read_mortality_csv(input_file, lower_names = TRUE)
  }
  names(raw) <- tolower(trimws(sub("^\ufeff", "", names(raw))))
  if (length(names(raw)) > 0 && grepl("study$", names(raw)[1], ignore.case = TRUE)) {
    names(raw)[1] <- "study"
  }
  raw <- raw[, names(raw) != "", drop = FALSE]

  required_cols <- c("study", "events_e", "n_e", "events_c", "n_c", "design")
  missing_cols <- setdiff(required_cols, names(raw))
  if (length(missing_cols) > 0) {
    stop("Missing required mortality columns: ", paste(missing_cols, collapse = ", "))
  }

  optional_cols <- intersect(c("arr", "arr_lower_95", "arr_upper_95", "arr_note"), names(raw))
  df <- raw[, c(required_cols, optional_cols), drop = FALSE]
  if (!"arr" %in% names(df)) df$arr <- NA_real_
  if (!"arr_lower_95" %in% names(df)) df$arr_lower_95 <- NA_real_
  if (!"arr_upper_95" %in% names(df)) df$arr_upper_95 <- NA_real_

  df$row_order <- seq_len(nrow(df))
  df$study <- clean_mortality_value(df$study)
  df$design_raw <- clean_mortality_value(df$design)
  df$design <- toupper(df$design_raw)
  df$events_e <- safe_mortality_num(df$events_e)
  df$n_e <- safe_mortality_num(df$n_e)
  df$events_c <- safe_mortality_num(df$events_c)
  df$n_c <- safe_mortality_num(df$n_c)
  df$aRR <- safe_mortality_num(df$arr)
  df$aRR_lower_95 <- safe_mortality_num(df$arr_lower_95)
  df$aRR_upper_95 <- safe_mortality_num(df$arr_upper_95)
  if ("arr_note" %in% names(df)) df$aRR_note <- clean_mortality_value(df$arr_note)
  df$arr <- NULL
  df$arr_lower_95 <- NULL
  df$arr_upper_95 <- NULL
  if ("arr_note" %in% names(df)) df$arr_note <- NULL

  validate_mortality_2x2(df, context = basename(if (is.character(input_file)) input_file else "mortality data"))
  validate_adjusted_rr_fields(df, context = basename(if (is.character(input_file)) input_file else "mortality data"))

  keep <- !is.na(df$study) &
    !is.na(df$events_e) & !is.na(df$n_e) &
    !is.na(df$events_c) & !is.na(df$n_c) &
    df$design %in% c("RCT", "NRS")
  df <- df[keep, , drop = FALSE]

  es <- metafor::escalc(
    measure = "RR",
    ai = events_e,
    bi = n_e - events_e,
    ci = events_c,
    di = n_c - events_c,
    data = df,
    add = 0.5,
    to = "all"
  )
  es$yi_crude <- es$yi
  es$vi_crude <- es$vi
  es$uses_adjusted_rr <- use_adjusted_rr &
    !is.na(es$aRR) & !is.na(es$aRR_lower_95) & !is.na(es$aRR_upper_95)
  adjusted_vi <- ((log(es$aRR_upper_95) - log(es$aRR_lower_95)) / (2 * 1.96))^2
  es$yi <- ifelse(es$uses_adjusted_rr, log(es$aRR), es$yi_crude)
  es$vi_original <- ifelse(es$uses_adjusted_rr, adjusted_vi, es$vi_crude)
  es$vi_bias_adjusted <- ifelse(es$design == "NRS", es$vi_original + nrs_bias_sd^2, es$vi_original)
  es$vi_adjusted <- es$vi_bias_adjusted
  es$se_original <- sqrt(es$vi_original)
  es$se_bias_adjusted <- sqrt(es$vi_bias_adjusted)
  es$se_adjusted <- es$se_bias_adjusted
  es$RR <- exp(es$yi)
  es$Lower_95 <- exp(es$yi - 1.96 * es$se_bias_adjusted)
  es$Upper_95 <- exp(es$yi + 1.96 * es$se_bias_adjusted)
  as.data.frame(es)
}

bayes_meta_settings <- function(iter = 4000,
                                warmup = floor(iter / 2),
                                chains = 4,
                                cores = min(chains, max(1, parallel::detectCores(logical = FALSE))),
                                seed = 20260722,
                                refresh = 0,
                                adapt_delta = 0.95,
                                max_treedepth = 12) {
  list(
    iter = iter,
    warmup = warmup,
    chains = chains,
    cores = cores,
    seed = seed,
    refresh = refresh,
    adapt_delta = adapt_delta,
    max_treedepth = max_treedepth
  )
}

bayes_meta_settings_from_env <- function(prefix = "BAYES_META") {
  get_int <- function(name, default) {
    value <- Sys.getenv(paste0(prefix, "_", name), unset = "")
    if (identical(value, "")) return(default)
    as.integer(value)
  }
  get_num <- function(name, default) {
    value <- Sys.getenv(paste0(prefix, "_", name), unset = "")
    if (identical(value, "")) return(default)
    as.numeric(value)
  }
  bayes_meta_settings(
    iter = get_int("ITER", 4000),
    warmup = get_int("WARMUP", 2000),
    chains = get_int("CHAINS", 4),
    cores = get_int("CORES", min(4, max(1, parallel::detectCores(logical = FALSE)))),
    seed = get_int("SEED", 20260722),
    refresh = get_int("REFRESH", 0),
    adapt_delta = get_num("ADAPT_DELTA", 0.95),
    max_treedepth = get_int("MAX_TREEDEPTH", 12)
  )
}

require_bayes_batch_enabled <- function(batch_name,
                                        env_var = "BAYES_ALLOW_NONCANONICAL",
                                        planned_script = "fit_canonical_main.R") {
  if (!identical(Sys.getenv(env_var, unset = ""), "1")) {
    stop(
      "This Bayesian batch is intentionally disabled for now: ", batch_name, ". ",
      "The current refactor only runs canonical main-analysis models via ", planned_script, ". ",
      "Refactor this batch into a separate fit script and plot script before long-chain use. ",
      "Set ", env_var, "=1 only for deliberate ad-hoc smoke testing."
    )
  }
  invisible(TRUE)
}

bayes_mortality_priors <- function(include_b = FALSE,
                                   include_intercept = TRUE,
                                   prior_spec_id = NULL,
                                   response_scale = NULL,
                                   prior_variant = "canonical") {
  if (is.null(prior_spec_id) || is.null(response_scale)) {
    stop("Bayesian model priors require explicit prior_spec_id and response_scale.")
  }
  prior_spec <- get_bayes_metareg_prior_spec_variant(
    prior_spec_id,
    response_scale = response_scale,
    prior_variant = prior_variant
  )
  if (include_b && is.na(prior_spec$slope_prior_sd)) {
    stop("The selected prior specification has no slope prior: ", prior_spec_id)
  }

  priors <- brms::prior_string(
    sprintf(
      "student_t(%s, %s, %s)",
      prior_spec$tau_prior_df,
      prior_spec$tau_prior_location,
      prior_spec$tau_prior_scale
    ),
    class = "sd"
  )
  if (include_intercept) {
    priors <- c(
      priors,
      brms::prior_string(
        sprintf("normal(%s, %s)", prior_spec$intercept_prior_mean, prior_spec$intercept_prior_sd),
        class = "Intercept"
      )
    )
  }
  if (include_b) {
    priors <- c(
      priors,
      brms::prior_string(
        sprintf("normal(%s, %s)", prior_spec$slope_prior_mean, prior_spec$slope_prior_sd),
        class = "b"
      )
    )
  }
  priors
}

prepare_bayes_mortality_data <- function(data, vi_col, study_col = "study") {
  required <- c("yi", vi_col, study_col)
  missing_cols <- setdiff(required, names(data))
  if (length(missing_cols) > 0) {
    stop("Missing required Bayesian mortality columns: ", paste(missing_cols, collapse = ", "))
  }
  out <- data
  out$.bayes_se <- sqrt(as.numeric(out[[vi_col]]))
  out$.bayes_study <- factor(out[[study_col]])
  if (any(!is.finite(out$yi)) || any(!is.finite(out$.bayes_se)) || any(out$.bayes_se <= 0)) {
    stop("Bayesian mortality model requires finite yi and positive finite standard errors.")
  }
  out
}

fit_bayes_mortality_meta <- function(data,
                                     vi_col,
                                     response_scale,
                                     prior_spec_id,
                                     uncertainty_spec_id,
                                     formula = NULL,
                                     settings = bayes_meta_settings(),
                                     study_col = "study",
                                     priors = NULL) {
  response_spec <- get_bayes_response_scale_spec(response_scale)
  prior_spec <- get_bayes_metareg_prior_spec(prior_spec_id, response_scale = response_scale)
  uncertainty_spec <- get_bayes_uncertainty_spec(uncertainty_spec_id)
  if (!identical(uncertainty_spec$response_scale[[1]], response_scale)) {
    stop("Uncertainty specification ", uncertainty_spec_id, " is not valid for response_scale ", response_scale)
  }
  model_data <- prepare_bayes_mortality_data(data, vi_col, study_col = study_col)
  model_formula <- if (is.null(formula)) {
    bf(yi | se(.bayes_se, sigma = FALSE) ~ 1 + (1 | .bayes_study))
  } else if (inherits(formula, "brmsformula")) {
    formula
  } else {
    bf(formula)
  }
  if (is.null(priors)) {
    formula_text <- paste(deparse(model_formula$formula), collapse = " ")
    has_zero_intercept <- grepl("~\\s*0\\s*\\+", formula_text) || grepl("~\\s*-1\\s*\\+", formula_text)
    has_population_terms <- !is.null(formula) && !identical(gsub("\\s+", "", formula_text), "yi|se(.bayes_se,sigma=FALSE)~1+(1|.bayes_study)")
    priors <- bayes_mortality_priors(
      include_b = has_population_terms,
      include_intercept = !has_zero_intercept,
      prior_spec_id = prior_spec_id,
      response_scale = response_scale
    )
  }

  fit <- brm(
    formula = model_formula,
    data = model_data,
    prior = priors,
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
    vi_col = vi_col,
    response_scale = response_scale,
    prior_spec_id = prior_spec_id,
    uncertainty_spec_id = uncertainty_spec_id,
    settings = settings,
    study_col = study_col,
    yi_col = "yi",
    formula_text = paste(deparse(model_formula$formula), collapse = " ")
  )
  attach_bayes_model_metadata(fit, model_metadata)
}

get_population_log_rr_draws <- function(fit, coefficient = "b_Intercept") {
  draws <- as_draws_df(fit)
  if (!coefficient %in% names(draws)) {
    stop("Coefficient not found in Bayesian fit draws: ", coefficient)
  }
  as.numeric(draws[[coefficient]])
}

get_tau_draws <- function(fit) {
  draws <- as_draws_df(fit)
  tau_cols <- grep("^sd_.*__Intercept$", names(draws), value = TRUE)
  if (length(tau_cols) == 0) return(rep(NA_real_, nrow(draws)))
  as.numeric(draws[[tau_cols[1]]])
}

summarise_bayes_mortality_fit <- function(fit,
                                          group,
                                          k_rct,
                                          k_nrs,
                                          coefficient = "b_Intercept",
                                          prob = 0.95,
                                          nrs_bias_sd = NA_real_) {
  alpha <- (1 - prob) / 2
  log_rr <- get_population_log_rr_draws(fit, coefficient = coefficient)
  rr <- exp(log_rr)
  tau <- get_tau_draws(fit)
  diagnostics <- extract_bayes_diagnostics(fit)

  data.frame(
    group = group,
    k = length(unique(fit$data$.bayes_study)),
    k_rct = k_rct,
    k_nrs = k_nrs,
    RR = median(rr, na.rm = TRUE),
    Lower_95_CrI = as.numeric(stats::quantile(rr, alpha, na.rm = TRUE)),
    Upper_95_CrI = as.numeric(stats::quantile(rr, 1 - alpha, na.rm = TRUE)),
    tau = median(tau, na.rm = TRUE),
    tau2 = median(tau^2, na.rm = TRUE),
    nrs_bias_sd = nrs_bias_sd,
    rhat_max = diagnostics$rhat_max,
    ess_bulk_min = diagnostics$ess_bulk_min,
    ess_tail_min = diagnostics$ess_tail_min,
    divergent_transitions = diagnostics$divergent_transitions,
    max_treedepth = diagnostics$max_treedepth,
    max_treedepth_hits = diagnostics$max_treedepth_hits,
    stringsAsFactors = FALSE
  )
}

infer_bayes_max_treedepth <- function(fit) {
  quality_gate <- attr(fit, "bayes_quality_gate")
  gate_depth <- if (is.list(quality_gate) && is.list(quality_gate$settings)) quality_gate$settings$max_treedepth else NULL
  if (length(gate_depth) == 1L && is.finite(gate_depth)) {
    return(as.integer(quality_gate$settings$max_treedepth))
  }
  for (attribute_name in c("bayes_metareg_metadata", "bayes_model_metadata")) {
    metadata <- attr(fit, attribute_name)
    metadata_depth <- if (is.list(metadata)) metadata$sampler_max_treedepth else NULL
    if (length(metadata_depth) == 1L && is.finite(metadata_depth)) {
      return(as.integer(metadata$sampler_max_treedepth))
    }
  }
  tryCatch({
    stan_fit <- fit[["fit"]]
    sampler_args <- attr(stan_fit@sim$samples[[1]], "args")
    as.integer(sampler_args$control$max_treedepth)
  }, error = function(e) NA_integer_)
}

extract_bayes_diagnostics <- function(fit, max_treedepth = NULL) {
  draws_summary <- summarise_draws(as_draws_df(fit))
  parameter_rows <- !grepl("^lp__$", draws_summary$variable)
  rhat_values <- draws_summary$rhat[parameter_rows]
  ess_bulk_values <- draws_summary$ess_bulk[parameter_rows]
  ess_tail_values <- draws_summary$ess_tail[parameter_rows]
  np <- nuts_params(fit)
  divergent <- if ("divergent__" %in% names(np)) {
    sum(np$divergent__, na.rm = TRUE)
  } else if (all(c("Parameter", "Value") %in% names(np)) && any(np$Parameter == "divergent__")) {
    sum(np$Value[np$Parameter == "divergent__"], na.rm = TRUE)
  } else {
    0
  }
  if (is.null(max_treedepth)) max_treedepth <- infer_bayes_max_treedepth(fit)
  max_treedepth <- as.integer(max_treedepth)[1]
  treedepth_values <- if ("treedepth__" %in% names(np)) {
    as.numeric(np$treedepth__)
  } else if (all(c("Parameter", "Value") %in% names(np)) && any(np$Parameter == "treedepth__")) {
    as.numeric(np$Value[np$Parameter == "treedepth__"])
  } else {
    numeric()
  }
  max_treedepth_hits <- if (is.finite(max_treedepth) && length(treedepth_values) > 0L) {
    sum(treedepth_values >= max_treedepth, na.rm = TRUE)
  } else {
    NA_real_
  }

  data.frame(
    rhat_max = suppressWarnings(max(rhat_values, na.rm = TRUE)),
    ess_bulk_min = suppressWarnings(min(ess_bulk_values, na.rm = TRUE)),
    ess_tail_min = suppressWarnings(min(ess_tail_values, na.rm = TRUE)),
    divergent_transitions = divergent,
    max_treedepth = max_treedepth,
    max_treedepth_hits = max_treedepth_hits,
    stringsAsFactors = FALSE
  )
}

write_bayes_diagnostics <- function(fit, path) {
  diagnostics <- extract_bayes_diagnostics(fit)
  utils::write.csv(diagnostics, path, row.names = FALSE)
  diagnostics
}

bayes_diagnostics_meet_gate <- function(diagnostics,
                                       rhat_max = 1.01,
                                       ess_bulk_min = 400,
                                       ess_tail_min = 400,
                                        divergent_transitions = 0,
                                        max_treedepth_hits = 0) {
  required <- c(
    "rhat_max", "ess_bulk_min", "ess_tail_min", "divergent_transitions",
    "max_treedepth", "max_treedepth_hits"
  )
  if (!is.data.frame(diagnostics) || nrow(diagnostics) != 1L ||
      !all(required %in% names(diagnostics))) {
    return(FALSE)
  }
  values <- diagnostics[1, required, drop = FALSE]
  isTRUE(
    is.finite(values$rhat_max) &&
      is.finite(values$ess_bulk_min) &&
      is.finite(values$ess_tail_min) &&
      is.finite(values$divergent_transitions) &&
      is.finite(values$max_treedepth) &&
      is.finite(values$max_treedepth_hits) &&
      values$rhat_max <= rhat_max &&
      values$ess_bulk_min >= ess_bulk_min &&
      values$ess_tail_min >= ess_tail_min &&
      values$divergent_transitions == divergent_transitions &&
      values$max_treedepth_hits == max_treedepth_hits
  )
}

upgrade_legacy_cached_bayes_fit <- function(fit,
                                            expected_metadata,
                                            metadata_attribute,
                                            matcher) {
  actual_metadata <- attr(fit, metadata_attribute)
  gate <- attr(fit, "bayes_quality_gate")
  if (!is.list(actual_metadata) || !is.list(gate) || !isTRUE(gate$passed) || !is.list(gate$settings)) {
    return(NULL)
  }

  actual_sampler <- tryCatch(bayes_sampler_metadata(gate$settings), error = function(e) NULL)
  if (is.null(actual_sampler) || !bayes_sampler_metadata_matches_policy(actual_sampler, expected_metadata)) {
    return(NULL)
  }
  for (name in names(actual_sampler)) actual_metadata[[name]] <- actual_sampler[[name]]
  diagnostics <- tryCatch(
    extract_bayes_diagnostics(fit, max_treedepth = gate$settings$max_treedepth),
    error = function(e) NULL
  )
  if (is.null(diagnostics) || !bayes_diagnostics_meet_gate(diagnostics)) return(NULL)

  gate$diagnostics <- diagnostics
  gate$passed <- TRUE
  attr(fit, metadata_attribute) <- actual_metadata
  attr(fit, "bayes_quality_gate") <- gate
  if (!isTRUE(matcher(fit, expected_metadata))) return(NULL)
  fit
}

extract_legacy_bayes_sampler_settings <- function(fit) {
  sampler_args <- tryCatch(
    lapply(fit[["fit"]]@sim$samples, function(chain) attr(chain, "args")),
    error = function(e) NULL
  )
  if (is.null(sampler_args) || length(sampler_args) == 0L ||
      any(!vapply(sampler_args, is.list, logical(1)))) {
    return(NULL)
  }

  extract_chain_value <- function(path) {
    vapply(sampler_args, function(args) {
      value <- args
      for (name in path) value <- value[[name]]
      as.numeric(value)[1]
    }, numeric(1))
  }
  iter <- extract_chain_value("iter")
  warmup <- extract_chain_value("warmup")
  seed <- extract_chain_value("random_seed")
  adapt_delta <- extract_chain_value(c("control", "adapt_delta"))
  max_treedepth <- extract_chain_value(c("control", "max_treedepth"))
  if (any(!is.finite(c(iter, warmup, seed, adapt_delta, max_treedepth))) ||
      length(unique(iter)) != 1L || length(unique(warmup)) != 1L ||
      length(unique(seed)) != 1L || length(unique(adapt_delta)) != 1L ||
      length(unique(max_treedepth)) != 1L) {
    return(NULL)
  }

  list(
    iter = as.integer(iter[1]),
    warmup = as.integer(warmup[1]),
    chains = as.integer(length(sampler_args)),
    cores = NA_integer_,
    seed = as.integer(seed[1]),
    refresh = NA_integer_,
    adapt_delta = as.numeric(adapt_delta[1]),
    max_treedepth = as.integer(max_treedepth[1])
  )
}

upgrade_unattributed_cached_intercept_fit <- function(fit,
                                                      expected_data,
                                                      expected_metadata) {
  if (!inherits(fit, "brmsfit") || !is.null(attr(fit, "bayes_model_metadata")) ||
      !is.null(attr(fit, "bayes_quality_gate"))) {
    return(NULL)
  }

  expected_formula <- "yi | se(.bayes_se, sigma = FALSE) ~ 1 + (1 | .bayes_study)"
  actual_formula <- tryCatch(paste(deparse(fit$formula$formula), collapse = " "), error = function(e) NA_character_)
  normalize_code <- function(x) gsub("\\s+", "", as.character(x))
  if (!identical(normalize_code(actual_formula), normalize_code(expected_formula)) ||
      !identical(normalize_code(expected_metadata$formula_text), normalize_code(expected_formula)) ||
      !identical(as.character(fit$family$family), "gaussian")) {
    return(NULL)
  }

  cached_data <- tryCatch(fit[["data"]], error = function(e) NULL)
  required_cached <- c("yi", ".bayes_se", ".bayes_study")
  required_expected <- c(expected_metadata$study_col, expected_metadata$yi_col, expected_metadata$vi_col)
  if (!is.data.frame(cached_data) || !is.data.frame(expected_data) ||
      !all(required_cached %in% names(cached_data)) ||
      !all(required_expected %in% names(expected_data)) ||
      nrow(cached_data) != nrow(expected_data)) {
    return(NULL)
  }
  same_numeric <- function(x, y) {
    isTRUE(all.equal(as.numeric(x), as.numeric(y), tolerance = 1e-12, check.attributes = FALSE))
  }
  if (!identical(as.character(cached_data$.bayes_study), as.character(expected_data[[expected_metadata$study_col]])) ||
      !same_numeric(cached_data$yi, expected_data[[expected_metadata$yi_col]]) ||
      !same_numeric(cached_data$.bayes_se^2, expected_data[[expected_metadata$vi_col]])) {
    return(NULL)
  }

  prior_table <- tryCatch(as.data.frame(fit$prior), error = function(e) NULL)
  if (!is.data.frame(prior_table) || !all(c("prior", "class", "source") %in% names(prior_table))) {
    return(NULL)
  }
  user_priors <- prior_table[prior_table$source == "user", c("prior", "class"), drop = FALSE]
  expected_priors <- data.frame(
    prior = c(expected_metadata$intercept_prior, expected_metadata$tau_prior),
    class = c("Intercept", "sd"),
    stringsAsFactors = FALSE
  )
  prior_key <- function(x) paste(normalize_code(x$prior), x$class, sep = "::")
  if (!setequal(prior_key(user_priors), prior_key(expected_priors))) return(NULL)

  sampler_settings <- extract_legacy_bayes_sampler_settings(fit)
  if (is.null(sampler_settings)) return(NULL)
  sampler_metadata <- bayes_sampler_metadata(sampler_settings)
  if (!bayes_sampler_metadata_matches_policy(sampler_metadata, expected_metadata)) return(NULL)

  diagnostics <- tryCatch(
    extract_bayes_diagnostics(fit, max_treedepth = sampler_settings$max_treedepth),
    error = function(e) NULL
  )
  if (is.null(diagnostics) || !bayes_diagnostics_meet_gate(diagnostics)) return(NULL)

  prepared_data <- prepare_bayes_mortality_data(
    expected_data,
    expected_metadata$vi_col,
    study_col = expected_metadata$study_col
  )
  for (name in names(sampler_metadata)) expected_metadata[[name]] <- sampler_metadata[[name]]
  fit$data <- prepared_data
  attr(fit, "bayes_model_metadata") <- expected_metadata
  sampler_label <- bayes_sampler_setting_label(sampler_settings)
  if (identical(sampler_label, "custom") &&
      identical(as.integer(sampler_settings$iter), 4000L) &&
      identical(as.integer(sampler_settings$warmup), 2000L) &&
      identical(as.integer(sampler_settings$chains), 4L) &&
      isTRUE(all.equal(as.numeric(sampler_settings$adapt_delta), 0.999)) &&
      identical(as.integer(sampler_settings$max_treedepth), 18L)) {
    sampler_label <- "validated_legacy_high_control"
  }
  attr(fit, "bayes_quality_gate") <- list(
    passed = TRUE,
    attempt = switch(
      sampler_label,
      enhanced = "retry",
      validated_legacy_high_control = "legacy_migration",
      "initial"
    ),
    settings = sampler_settings,
    sampler_setting_label = sampler_label,
    diagnostics = diagnostics,
    provenance = "strictly validated migration of an unattributed legacy fit"
  )
  if (!isTRUE(assert_cached_bayes_model_matches(fit, expected_metadata))) return(NULL)
  fit
}

restore_cached_bayes_model_data <- function(fit,
                                            expected_data,
                                            expected_metadata,
                                            cached_se_col,
                                            expected_formula) {
  actual_metadata <- attr(fit, "bayes_model_metadata")
  gate <- attr(fit, "bayes_quality_gate")
  if (!inherits(fit, "brmsfit") || !is.list(actual_metadata) || !is.list(gate)) return(NULL)

  sampler_names <- bayes_sampler_metadata_names()
  core_names <- setdiff(names(expected_metadata), sampler_names)
  if (!all(core_names %in% names(actual_metadata)) ||
      !all(vapply(core_names, function(name) identical(actual_metadata[[name]], expected_metadata[[name]]), logical(1))) ||
      !bayes_sampler_metadata_matches_policy(actual_metadata, expected_metadata) ||
      !cached_bayes_quality_gate_is_valid(fit, actual_metadata)) {
    return(NULL)
  }

  normalize_code <- function(x) gsub("\\s+", "", as.character(x))
  actual_formula <- tryCatch(paste(deparse(fit$formula$formula), collapse = " "), error = function(e) NA_character_)
  cached_data <- tryCatch(fit[["data"]], error = function(e) NULL)
  required_cached <- c(expected_metadata$yi_col, expected_metadata$study_col, cached_se_col)
  required_expected <- c(expected_metadata$yi_col, expected_metadata$study_col, expected_metadata$vi_col)
  if (!identical(normalize_code(actual_formula), normalize_code(expected_formula)) ||
      !is.data.frame(cached_data) || !is.data.frame(expected_data) ||
      !all(required_cached %in% names(cached_data)) || !all(required_expected %in% names(expected_data)) ||
      nrow(cached_data) != nrow(expected_data)) {
    return(NULL)
  }
  same_numeric <- function(x, y) {
    isTRUE(all.equal(as.numeric(x), as.numeric(y), tolerance = 1e-12, check.attributes = FALSE))
  }
  if (!identical(as.character(cached_data[[expected_metadata$study_col]]), as.character(expected_data[[expected_metadata$study_col]])) ||
      !same_numeric(cached_data[[expected_metadata$yi_col]], expected_data[[expected_metadata$yi_col]]) ||
      !same_numeric(cached_data[[cached_se_col]]^2, expected_data[[expected_metadata$vi_col]])) {
    return(NULL)
  }

  fit$data <- expected_data
  if (!isTRUE(assert_cached_bayes_model_matches(fit, expected_metadata))) return(NULL)
  fit
}

upgrade_legacy_cached_metareg_fit <- function(fit, expected_metadata) {
  upgrade_legacy_cached_bayes_fit(
    fit,
    expected_metadata,
    metadata_attribute = "bayes_metareg_metadata",
    matcher = assert_cached_metareg_fit_matches
  )
}

upgrade_legacy_cached_model_fit <- function(fit, expected_metadata) {
  upgrade_legacy_cached_bayes_fit(
    fit,
    expected_metadata,
    metadata_attribute = "bayes_model_metadata",
    matcher = assert_cached_bayes_model_matches
  )
}

bayes_retry_settings <- function(settings) {
  retry <- settings
  retry$iter <- 8000
  retry$warmup <- 4000
  retry$adapt_delta <- 0.999
  retry$max_treedepth <- 18
  retry
}

bayes_sampler_setting_label <- function(settings) {
  if (is.null(settings) || !is.list(settings)) return("unknown")
  if (identical(as.integer(settings$iter), 4000L) &&
      identical(as.integer(settings$warmup), 2000L) &&
      identical(as.integer(settings$chains), 4L) &&
      isTRUE(all.equal(as.numeric(settings$adapt_delta), 0.95)) &&
      identical(as.integer(settings$max_treedepth), 12L)) {
    return("standard")
  }
  if (identical(as.integer(settings$iter), 8000L) &&
      identical(as.integer(settings$warmup), 4000L) &&
      identical(as.integer(settings$chains), 4L) &&
      isTRUE(all.equal(as.numeric(settings$adapt_delta), 0.999)) &&
      identical(as.integer(settings$max_treedepth), 18L)) {
    return("enhanced")
  }
  "custom"
}

fit_reml_mortality_reference <- function(data, vi_col) {
  rma(
    yi = yi,
    vi = data[[vi_col]],
    data = data,
    method = "REML",
    control = list(maxiter = 10000, stepadj = 0.1)
  )
}

summarise_reml_mortality_reference <- function(fit, group, k_rct, k_nrs, nrs_bias_sd = NA_real_) {
  data.frame(
    group = group,
    k = fit$k,
    k_rct = k_rct,
    k_nrs = k_nrs,
    RR = exp(as.numeric(fit$b)[1]),
    Lower_95 = exp(as.numeric(fit$ci.lb)[1]),
    Upper_95 = exp(as.numeric(fit$ci.ub)[1]),
    tau2 = as.numeric(fit$tau2)[1],
    I2_percent = as.numeric(fit$I2)[1],
    p_heterogeneity = as.numeric(fit$QEp)[1],
    nrs_bias_sd = nrs_bias_sd,
    stringsAsFactors = FALSE
  )
}

format_bayes_rr <- function(rr, lower, upper) {
  sprintf("%.3f [%.3f, %.3f]", rr, lower, upper)
}

fit_bayes_mortality_metareg <- function(data,
                                        moderator_col,
                                        moderator_id,
                                        response_scale,
                                        prior_spec_id,
                                        uncertainty_spec_id,
                                        vi_col = "vi_adjusted",
                                        settings = bayes_meta_settings(),
                                        study_col = "study",
                                        priors = NULL,
                                        prior_variant = "canonical") {
  validate_bayes_metareg_spec(moderator_id, response_scale, prior_spec_id, uncertainty_spec_id)
  model_data <- prepare_bayes_mortality_data(data, vi_col, study_col = study_col)
  model_data$bayes_mod_raw <- as.numeric(model_data[[moderator_col]])
  model_data$bayes_mod <- scale_bayes_moderator(model_data$bayes_mod_raw, moderator_id)
  if (any(!is.finite(model_data$bayes_mod_raw)) || any(!is.finite(model_data$bayes_mod))) {
    stop("Bayesian mortality meta-regression requires a finite numeric moderator: ", moderator_col)
  }
  model_metadata <- make_bayes_metareg_cache_metadata(
    data = model_data,
    moderator_col = moderator_col,
    vi_col = vi_col,
    moderator_id = moderator_id,
    response_scale = response_scale,
    prior_spec_id = prior_spec_id,
    uncertainty_spec_id = uncertainty_spec_id,
    settings = settings,
    prior_variant = prior_variant,
    study_col = study_col
  )
  model_data$moderator_id <- rep(moderator_id, nrow(model_data))
  model_data$reported_increment <- rep(model_metadata$moderator_reported_increment, nrow(model_data))
  model_data$moderator_centre_raw <- rep(model_metadata$moderator_centre_raw, nrow(model_data))
  model_data$moderator_increment_raw <- rep(model_metadata$moderator_increment_raw, nrow(model_data))
  model_data$reported_increment_label <- rep(
    get_bayes_reported_increment_label(moderator_id),
    nrow(model_data)
  )
  model_data$model_unit_description <- rep(
    model_metadata$moderator_model_unit_description,
    nrow(model_data)
  )
  model_data$response_scale <- rep(response_scale, nrow(model_data))
  model_data$prior_spec_id <- rep(prior_spec_id, nrow(model_data))
  model_data$prior_variant <- rep(prior_variant, nrow(model_data))
  model_data$uncertainty_spec_id <- rep(uncertainty_spec_id, nrow(model_data))
  model_data$model_spec_version <- rep(model_metadata$model_spec_version, nrow(model_data))
  if (is.null(priors)) {
    priors <- bayes_mortality_priors(
      include_b = TRUE,
      include_intercept = TRUE,
      prior_spec_id = prior_spec_id,
      response_scale = response_scale,
      prior_variant = prior_variant
    )
  }
  fit <- brm(
    formula = bf(yi | se(.bayes_se, sigma = FALSE) ~ bayes_mod + (1 | .bayes_study)),
    data = model_data,
    prior = priors,
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
  # brms may retain only variables used by the formula in fit$data. Preserve
  # the complete prepared data so the canonical raw/model moderator, response,
  # uncertainty, and provenance columns remain auditable in the fit cache.
  fit$data <- model_data
  model_metadata$sampler_setting_label <- bayes_sampler_setting_label(settings)
  attach_bayes_metareg_metadata(fit, model_metadata)
}

fit_bayes_mortality_metareg_with_quality_gate <- function(
    data,
    moderator_col,
    moderator_id,
    response_scale,
    prior_spec_id,
    uncertainty_spec_id,
    vi_col = "vi_adjusted",
    settings = bayes_meta_settings(),
    study_col = "study",
    priors = NULL,
    prior_variant = "canonical") {
  fit <- fit_bayes_mortality_metareg(
    data = data,
    moderator_col = moderator_col,
    moderator_id = moderator_id,
    response_scale = response_scale,
    prior_spec_id = prior_spec_id,
    uncertainty_spec_id = uncertainty_spec_id,
    vi_col = vi_col,
    settings = settings,
    study_col = study_col,
    priors = priors,
    prior_variant = prior_variant
  )
  diagnostics <- extract_bayes_diagnostics(fit, max_treedepth = settings$max_treedepth)
  quality_gate_attempt <- "initial"
  if (!bayes_diagnostics_meet_gate(diagnostics)) {
    fit <- fit_bayes_mortality_metareg(
      data = data,
      moderator_col = moderator_col,
      moderator_id = moderator_id,
      response_scale = response_scale,
      prior_spec_id = prior_spec_id,
      uncertainty_spec_id = uncertainty_spec_id,
      vi_col = vi_col,
      settings = bayes_retry_settings(settings),
      study_col = study_col,
      priors = priors,
      prior_variant = prior_variant
    )
    retry_settings <- bayes_retry_settings(settings)
    diagnostics <- extract_bayes_diagnostics(fit, max_treedepth = retry_settings$max_treedepth)
    quality_gate_attempt <- "retry"
  }
  if (!bayes_diagnostics_meet_gate(diagnostics)) {
    stop(
      "Bayesian meta-regression failed the MCMC quality gate after the approved retry: ",
      "response_scale=", response_scale,
      ", moderator_id=", moderator_id,
      ", prior_spec_id=", prior_spec_id,
      ", rhat_max=", diagnostics$rhat_max,
      ", ess_bulk_min=", diagnostics$ess_bulk_min,
      ", ess_tail_min=", diagnostics$ess_tail_min,
      ", divergent_transitions=", diagnostics$divergent_transitions,
      ", max_treedepth_hits=", diagnostics$max_treedepth_hits
    )
  }
  attr(fit, "bayes_quality_gate") <- list(
    passed = TRUE,
    attempt = quality_gate_attempt,
    settings = if (identical(quality_gate_attempt, "retry")) {
      bayes_retry_settings(settings)
    } else {
      settings
    },
    sampler_setting_label = if (identical(quality_gate_attempt, "retry")) {
      bayes_sampler_setting_label(bayes_retry_settings(settings))
    } else {
      bayes_sampler_setting_label(settings)
    },
    diagnostics = diagnostics
  )
  fit
}

fit_bayes_mortality_meta_with_quality_gate <- function(
    data,
    vi_col,
    response_scale,
    prior_spec_id,
    uncertainty_spec_id,
    formula = NULL,
    settings = bayes_meta_settings(),
    study_col = "study",
    priors = NULL) {
  fit <- fit_bayes_mortality_meta(
    data = data,
    vi_col = vi_col,
    response_scale = response_scale,
    prior_spec_id = prior_spec_id,
    uncertainty_spec_id = uncertainty_spec_id,
    formula = formula,
    settings = settings,
    study_col = study_col,
    priors = priors
  )
  diagnostics <- extract_bayes_diagnostics(fit, max_treedepth = settings$max_treedepth)
  attempt <- "initial"
  if (!bayes_diagnostics_meet_gate(diagnostics)) {
    retry_settings <- bayes_retry_settings(settings)
    fit <- fit_bayes_mortality_meta(
      data = data,
      vi_col = vi_col,
      response_scale = response_scale,
      prior_spec_id = prior_spec_id,
      uncertainty_spec_id = uncertainty_spec_id,
      formula = formula,
      settings = retry_settings,
      study_col = study_col,
      priors = priors
    )
    diagnostics <- extract_bayes_diagnostics(fit, max_treedepth = retry_settings$max_treedepth)
    attempt <- "retry"
  } else {
    retry_settings <- settings
  }
  if (!bayes_diagnostics_meet_gate(diagnostics)) {
    stop(
      "Bayesian meta-analysis failed the MCMC quality gate after the approved retry: ",
      "response_scale=", response_scale,
      ", prior_spec_id=", prior_spec_id,
      ", rhat_max=", diagnostics$rhat_max,
      ", ess_bulk_min=", diagnostics$ess_bulk_min,
      ", ess_tail_min=", diagnostics$ess_tail_min,
      ", divergent_transitions=", diagnostics$divergent_transitions,
      ", max_treedepth_hits=", diagnostics$max_treedepth_hits
    )
  }
  attr(fit, "bayes_quality_gate") <- list(
    passed = TRUE,
    attempt = attempt,
    settings = if (identical(attempt, "retry")) retry_settings else settings,
    sampler_setting_label = if (identical(attempt, "retry")) {
      bayes_sampler_setting_label(retry_settings)
    } else {
      bayes_sampler_setting_label(settings)
    },
    diagnostics = diagnostics
  )
  fit
}

fit_bayesian_rd_meta_with_quality_gate <- function(
    data,
    settings = bayes_meta_settings(),
    fit_fun = NULL,
    ...) {
  if (is.null(fit_fun)) {
    if (!exists("fit_bayesian_rd_meta", mode = "function", inherits = TRUE)) {
      stop("fit_bayesian_rd_meta() must be available before using its quality gate.")
    }
    fit_fun <- get("fit_bayesian_rd_meta", mode = "function", inherits = TRUE)
  }
  fit <- fit_fun(data, settings, ...)
  diagnostics <- extract_bayes_diagnostics(fit, max_treedepth = settings$max_treedepth)
  attempt <- "initial"
  if (!bayes_diagnostics_meet_gate(diagnostics)) {
    retry_settings <- bayes_retry_settings(settings)
    fit <- fit_fun(data, retry_settings, ...)
    diagnostics <- extract_bayes_diagnostics(fit, max_treedepth = retry_settings$max_treedepth)
    attempt <- "retry"
  } else {
    retry_settings <- settings
  }
  if (!bayes_diagnostics_meet_gate(diagnostics)) {
    stop("Bayesian RD meta-analysis failed the approved MCMC quality gate after retry.")
  }
  attr(fit, "bayes_quality_gate") <- list(
    passed = TRUE,
    attempt = attempt,
    settings = if (identical(attempt, "retry")) retry_settings else settings,
    sampler_setting_label = if (identical(attempt, "retry")) {
      bayes_sampler_setting_label(retry_settings)
    } else {
      bayes_sampler_setting_label(settings)
    },
    diagnostics = diagnostics
  )
  fit
}

summarise_bayes_mortality_metareg <- function(fit, covariate_name) {
  draws <- as_draws_df(fit)
  beta <- as.numeric(draws$b_bayes_mod)
  metadata <- attr(fit, "bayes_metareg_metadata")
  response_scale <- if (is.list(metadata)) metadata$response_scale else NA_character_
  canonical_fit <- is.list(metadata) && identical(
    metadata$model_spec_version,
    BAYES_PRIOR_SCALE_SPEC_VERSION
  )
  response_display_scale <- if (identical(response_scale, "rd")) 100 else 1
  tau <- get_tau_draws(fit)
  diagnostics <- extract_bayes_diagnostics(fit)
  prob_positive <- mean(beta > 0, na.rm = TRUE)
  prob_negative <- mean(beta < 0, na.rm = TRUE)
  two_tail_probability <- 2 * min(prob_positive, prob_negative)

  data.frame(
    covariate = covariate_name,
    beta = median(beta, na.rm = TRUE),
    beta_lower = as.numeric(stats::quantile(beta, 0.025, na.rm = TRUE)),
    beta_upper = as.numeric(stats::quantile(beta, 0.975, na.rm = TRUE)),
    beta_model = if (canonical_fit) median(beta, na.rm = TRUE) else NA_real_,
    beta_lower_model = if (canonical_fit) as.numeric(stats::quantile(beta, 0.025, na.rm = TRUE)) else NA_real_,
    beta_upper_model = if (canonical_fit) as.numeric(stats::quantile(beta, 0.975, na.rm = TRUE)) else NA_real_,
    beta_reported = if (canonical_fit) median(beta, na.rm = TRUE) else NA_real_,
    beta_lower_reported = if (canonical_fit) as.numeric(stats::quantile(beta, 0.025, na.rm = TRUE)) else NA_real_,
    beta_upper_reported = if (canonical_fit) as.numeric(stats::quantile(beta, 0.975, na.rm = TRUE)) else NA_real_,
    beta_reported_unit = if (canonical_fit) median(beta, na.rm = TRUE) * response_display_scale else NA_real_,
    beta_lower_reported_unit = if (canonical_fit) as.numeric(stats::quantile(beta, 0.025, na.rm = TRUE)) * response_display_scale else NA_real_,
    beta_upper_reported_unit = if (canonical_fit) as.numeric(stats::quantile(beta, 0.975, na.rm = TRUE)) * response_display_scale else NA_real_,
    reported_increment_label = if (canonical_fit) {
      get_bayes_metareg_moderator_spec(metadata$moderator_id)$reported_increment_label
    } else {
      NA_character_
    },
    moderator_id = if (canonical_fit) metadata$moderator_id else NA_character_,
    prior_variant = if (canonical_fit) metadata$prior_variant else NA_character_,
    moderator_centre_raw = if (canonical_fit) metadata$moderator_centre_raw else NA_real_,
    moderator_increment_raw = if (canonical_fit) metadata$moderator_increment_raw else NA_real_,
    response_display_multiplier = response_display_scale,
    model_spec_version = if (canonical_fit) metadata$model_spec_version else NA_character_,
    sampler_iter = if (canonical_fit) metadata$sampler_iter else NA_integer_,
    sampler_warmup = if (canonical_fit) metadata$sampler_warmup else NA_integer_,
    sampler_chains = if (canonical_fit) metadata$sampler_chains else NA_integer_,
    sampler_seed = if (canonical_fit) metadata$sampler_seed else NA_integer_,
    sampler_adapt_delta = if (canonical_fit) metadata$sampler_adapt_delta else NA_real_,
    sampler_max_treedepth = if (canonical_fit) metadata$sampler_max_treedepth else NA_integer_,
    sampler_setting_label = if (canonical_fit) metadata$sampler_setting_label else NA_character_,
    prob_direction = max(prob_positive, prob_negative),
    p_value = two_tail_probability,
    tau = median(tau, na.rm = TRUE),
    tau2 = median(tau^2, na.rm = TRUE),
    rhat_max = diagnostics$rhat_max,
    ess_bulk_min = diagnostics$ess_bulk_min,
    ess_tail_min = diagnostics$ess_tail_min,
    divergent_transitions = diagnostics$divergent_transitions,
    max_treedepth = diagnostics$max_treedepth,
    max_treedepth_hits = diagnostics$max_treedepth_hits,
    method = "Bayesian random-effects meta-regression",
    stringsAsFactors = FALSE
  )
}

predict_bayes_mortality_metareg <- function(fit, x_grid, prob = 0.95) {
  alpha <- (1 - prob) / 2
  newdata <- data.frame(
    bayes_mod = x_grid,
    .bayes_se = 1,
    .bayes_study = fit$data$.bayes_study[1]
  )
  linpred <- posterior_linpred(fit, newdata = newdata, re_formula = NA, transform = FALSE)
  pred <- apply(linpred, 2, median, na.rm = TRUE)
  lower <- apply(linpred, 2, stats::quantile, probs = alpha, na.rm = TRUE)
  upper <- apply(linpred, 2, stats::quantile, probs = 1 - alpha, na.rm = TRUE)
  data.frame(
    x = x_grid,
    pred = pred,
    ci_lb = lower,
    ci_ub = upper,
    RR = exp(pred),
    lower_RR = exp(lower),
    upper_RR = exp(upper),
    stringsAsFactors = FALSE
  )
}

predict_bayes_mortality_metareg_raw <- function(fit,
                                                x_raw,
                                                moderator_id,
                                                prob = 0.95) {
  raw <- as.numeric(x_raw)
  if (any(!is.finite(raw))) stop("Raw prediction moderator values must be finite.")
  metadata <- attr(fit, "bayes_metareg_metadata")
  if (!is.list(metadata) || !identical(metadata$moderator_id, moderator_id)) {
    stop("Raw prediction moderator_id does not match the fitted model metadata.")
  }
  model <- transform_bayes_moderator(raw, moderator_id)
  out <- predict_bayes_mortality_metareg(fit, model, prob = prob)
  out$x_raw <- raw
  out$x_model <- model
  out$x <- raw
  out[, c("x", "x_raw", "x_model", "pred", "ci_lb", "ci_ub", "RR", "lower_RR", "upper_RR")]
}

canonical_bayes_group_specs <- function(es, nrs_bias_sd = 0.20) {
  list(
    rct_only = list(
      group_id = "rct_only",
      group = "RCT",
      data = es[es$design == "RCT", , drop = FALSE],
      vi_col = "vi_bias_adjusted",
      k_rct = sum(es$design == "RCT"),
      k_nrs = 0,
      nrs_bias_sd = nrs_bias_sd
    ),
    nrs_only = list(
      group_id = "nrs_only",
      group = sprintf("NRS with bias SD %.2f", nrs_bias_sd),
      data = es[es$design == "NRS", , drop = FALSE],
      vi_col = "vi_bias_adjusted",
      k_rct = 0,
      k_nrs = sum(es$design == "NRS"),
      nrs_bias_sd = nrs_bias_sd
    ),
    overall = list(
      group_id = "overall",
      group = sprintf("Overall with NRS bias SD %.2f", nrs_bias_sd),
      data = es,
      vi_col = "vi_bias_adjusted",
      k_rct = sum(es$design == "RCT"),
      k_nrs = sum(es$design == "NRS"),
      nrs_bias_sd = nrs_bias_sd
    )
  )
}

canonical_bayes_paths <- function(out_dir) {
  list(
    root = out_dir,
    fits = file.path(out_dir, "fits"),
    posterior_summary = file.path(out_dir, "canonical_mortality_posterior_summary.csv"),
    diagnostics = file.path(out_dir, "canonical_mortality_bayes_diagnostics.csv"),
    reml_reference = file.path(out_dir, "canonical_mortality_reml_reference.csv"),
    effect_sizes = file.path(out_dir, "canonical_mortality_effect_sizes.csv"),
    metadata = file.path(out_dir, "canonical_mortality_metadata.csv")
  )
}

ensure_canonical_bayes_dirs <- function(out_dir) {
  paths <- canonical_bayes_paths(out_dir)
  if (!dir.exists(paths$root)) dir.create(paths$root, recursive = TRUE)
  if (!dir.exists(paths$fits)) dir.create(paths$fits, recursive = TRUE)
  paths
}

canonical_file_md5 <- function(path) {
  if (!file.exists(path)) stop("Cannot hash missing canonical Bayesian artifact: ", path)
  unname(tools::md5sum(path))
}

canonical_fit_paths <- function(paths) {
  setNames(
    file.path(paths$fits, paste0(c("rct_only", "nrs_only", "overall"), ".brmsfit.rds")),
    c("rct_only", "nrs_only", "overall")
  )
}

write_canonical_bayes_metadata <- function(out_dir, settings, nrs_bias_sd, input_file) {
  paths <- ensure_canonical_bayes_dirs(out_dir)
  fit_paths <- canonical_fit_paths(paths)
  input_file <- normalizePath(input_file, winslash = "/", mustWork = TRUE)
  entries <- c(
    cache_schema_version = "canonical_bayes_manifest_v2",
    created_at = format(Sys.time(), "%Y-%m-%d %H:%M:%S %z"),
    input_file = input_file,
    input_md5 = canonical_file_md5(input_file),
    nrs_bias_sd = as.character(nrs_bias_sd),
    iter = as.character(settings$iter),
    warmup = as.character(settings$warmup),
    chains = as.character(settings$chains),
    seed = as.character(settings$seed),
    adapt_delta = as.character(settings$adapt_delta),
    max_treedepth = as.character(settings$max_treedepth),
    effect_sizes_md5 = canonical_file_md5(paths$effect_sizes),
    posterior_summary_md5 = canonical_file_md5(paths$posterior_summary),
    diagnostics_md5 = canonical_file_md5(paths$diagnostics),
    reml_reference_md5 = canonical_file_md5(paths$reml_reference),
    fit_rct_only_md5 = canonical_file_md5(fit_paths[["rct_only"]]),
    fit_nrs_only_md5 = canonical_file_md5(fit_paths[["nrs_only"]]),
    fit_overall_md5 = canonical_file_md5(fit_paths[["overall"]])
  )
  metadata <- data.frame(
    key = names(entries),
    value = unname(entries),
    stringsAsFactors = FALSE
  )
  utils::write.csv(metadata, paths$metadata, row.names = FALSE)
  metadata
}

fit_and_save_canonical_bayes_mortality <- function(es,
                                                   out_dir,
                                                   settings = bayes_meta_settings(),
                                                   nrs_bias_sd = 0.20,
                                                   input_file = "mortality.csv") {
  paths <- ensure_canonical_bayes_dirs(out_dir)
  utils::write.csv(es, paths$effect_sizes, row.names = FALSE)

  specs <- canonical_bayes_group_specs(es, nrs_bias_sd = nrs_bias_sd)
  posterior_rows <- list()
  diagnostic_rows <- list()
  reml_rows <- list()

  for (group_id in names(specs)) {
    spec <- specs[[group_id]]
    if (nrow(spec$data) < 2) {
      stop("Canonical Bayesian model requires at least two studies for ", group_id)
    }

    fit <- fit_bayes_mortality_meta_with_quality_gate(
      spec$data,
      spec$vi_col,
      response_scale = "log_rr",
      prior_spec_id = "log_rr_shared_current",
      uncertainty_spec_id = "primary_nrs_log_rr_bias_sd_0.20",
      settings = settings
    )
    saveRDS(fit, file.path(paths$fits, paste0(group_id, ".brmsfit.rds")))

    posterior_rows[[group_id]] <- summarise_bayes_mortality_fit(
      fit,
      group = spec$group,
      k_rct = spec$k_rct,
      k_nrs = spec$k_nrs,
      nrs_bias_sd = spec$nrs_bias_sd
    ) |>
      dplyr::mutate(group_id = group_id, .before = 1)

    diagnostic_rows[[group_id]] <- extract_bayes_diagnostics(fit) |>
      dplyr::mutate(group_id = group_id, group = spec$group, .before = 1)

    fit_reml <- fit_reml_mortality_reference(spec$data, spec$vi_col)
    reml_rows[[group_id]] <- summarise_reml_mortality_reference(
      fit_reml,
      spec$group,
      spec$k_rct,
      spec$k_nrs,
      spec$nrs_bias_sd
    ) |>
      dplyr::mutate(group_id = group_id, .before = 1)
  }

  posterior_summary <- dplyr::bind_rows(posterior_rows)
  diagnostics <- dplyr::bind_rows(diagnostic_rows)
  reml_reference <- dplyr::bind_rows(reml_rows)

  utils::write.csv(posterior_summary, paths$posterior_summary, row.names = FALSE)
  utils::write.csv(diagnostics, paths$diagnostics, row.names = FALSE)
  utils::write.csv(reml_reference, paths$reml_reference, row.names = FALSE)
  write_canonical_bayes_metadata(out_dir, settings, nrs_bias_sd, input_file)

  invisible(list(
    posterior_summary = posterior_summary,
    diagnostics = diagnostics,
    reml_reference = reml_reference,
    paths = paths
  ))
}

load_canonical_bayes_mortality <- function(out_dir,
                                           required_group_ids = c("rct_only", "nrs_only", "overall"),
                                           require_fits = TRUE,
                                           input_file = NULL,
                                           expected_settings = NULL,
                                           expected_nrs_bias_sd = NULL) {
  paths <- canonical_bayes_paths(out_dir)
  required_files <- c(paths$posterior_summary, paths$diagnostics, paths$reml_reference, paths$effect_sizes, paths$metadata)
  missing_files <- required_files[!file.exists(required_files)]
  if (length(missing_files) > 0) {
    stop(
      "Canonical Bayesian mortality results are missing. Run fit_canonical_main.R first. Missing: ",
      paste(missing_files, collapse = "; ")
    )
  }

  fit_files <- file.path(paths$fits, paste0(required_group_ids, ".brmsfit.rds"))
  if (require_fits) {
    missing_fit_files <- fit_files[!file.exists(fit_files)]
    if (length(missing_fit_files) > 0) {
      stop(
        "Canonical Bayesian brmsfit files are missing. Run fit_canonical_main.R first. Missing: ",
        paste(missing_fit_files, collapse = "; ")
      )
    }
  }

  metadata <- utils::read.csv(paths$metadata, stringsAsFactors = FALSE, check.names = FALSE)
  if (!all(c("key", "value") %in% names(metadata)) || anyDuplicated(metadata$key)) {
    stop("Canonical Bayesian metadata is malformed.")
  }
  manifest <- stats::setNames(as.character(metadata$value), metadata$key)
  required_manifest <- c(
    "cache_schema_version", "input_md5", "nrs_bias_sd", "iter", "warmup", "chains",
    "seed", "adapt_delta", "max_treedepth", "effect_sizes_md5", "posterior_summary_md5",
    "diagnostics_md5", "reml_reference_md5", "fit_rct_only_md5", "fit_nrs_only_md5",
    "fit_overall_md5"
  )
  missing_manifest <- setdiff(required_manifest, names(manifest))
  if (length(missing_manifest) > 0L || !identical(manifest[["cache_schema_version"]], "canonical_bayes_manifest_v2")) {
    stop("Canonical Bayesian metadata is stale or incomplete. Re-run fit_canonical_main.R.")
  }

  artifact_manifest <- c(
    effect_sizes_md5 = paths$effect_sizes,
    posterior_summary_md5 = paths$posterior_summary,
    diagnostics_md5 = paths$diagnostics,
    reml_reference_md5 = paths$reml_reference
  )
  for (key in names(artifact_manifest)) {
    if (!identical(canonical_file_md5(artifact_manifest[[key]]), manifest[[key]])) {
      stop("Canonical Bayesian artifact failed its integrity check: ", basename(artifact_manifest[[key]]))
    }
  }
  if (require_fits) {
    for (group_id in required_group_ids) {
      key <- paste0("fit_", group_id, "_md5")
      if (!identical(canonical_file_md5(file.path(paths$fits, paste0(group_id, ".brmsfit.rds"))), manifest[[key]])) {
        stop("Canonical Bayesian fit failed its integrity check: ", group_id)
      }
    }
  }

  if (!is.null(input_file)) {
    input_file <- normalizePath(input_file, winslash = "/", mustWork = TRUE)
    if (!identical(canonical_file_md5(input_file), manifest[["input_md5"]])) {
      stop("Canonical Bayesian mortality results are stale relative to the current mortality input.")
    }
  }
  if (!is.null(expected_nrs_bias_sd) &&
      !isTRUE(all.equal(as.numeric(manifest[["nrs_bias_sd"]]), as.numeric(expected_nrs_bias_sd), tolerance = 1e-12))) {
    stop("Canonical Bayesian NRS bias-SD setting does not match the requested analysis.")
  }

  manifest_settings <- bayes_meta_settings(
    iter = as.integer(manifest[["iter"]]),
    warmup = as.integer(manifest[["warmup"]]),
    chains = as.integer(manifest[["chains"]]),
    seed = as.integer(manifest[["seed"]]),
    adapt_delta = as.numeric(manifest[["adapt_delta"]]),
    max_treedepth = as.integer(manifest[["max_treedepth"]])
  )
  if (!is.null(expected_settings)) {
    expected_sampler <- bayes_sampler_metadata(expected_settings)
    manifest_sampler <- bayes_sampler_metadata(manifest_settings)
    if (!all(vapply(bayes_sampler_metadata_names(), function(name) {
      identical(expected_sampler[[name]], manifest_sampler[[name]])
    }, logical(1)))) {
      stop("Canonical Bayesian sampler settings do not match the requested analysis.")
    }
  }

  posterior_summary <- utils::read.csv(paths$posterior_summary, stringsAsFactors = FALSE, check.names = FALSE)
  diagnostics <- utils::read.csv(paths$diagnostics, stringsAsFactors = FALSE, check.names = FALSE)
  reml_reference <- utils::read.csv(paths$reml_reference, stringsAsFactors = FALSE, check.names = FALSE)
  effect_sizes <- utils::read.csv(paths$effect_sizes, stringsAsFactors = FALSE, check.names = FALSE)
  fit_effect_sizes <- effect_sizes
  if (!is.null(input_file)) {
    fit_effect_sizes <- prepare_mortality_effect_sizes(
      input_file,
      nrs_bias_sd = as.numeric(manifest[["nrs_bias_sd"]])
    )
    assert_canonical_effect_sizes_match(fit_effect_sizes, effect_sizes)
  }

  missing_groups <- setdiff(required_group_ids, posterior_summary$group_id)
  if (length(missing_groups) > 0) {
    stop("Canonical Bayesian posterior summary lacks groups: ", paste(missing_groups, collapse = ", "))
  }
  if (!all(vapply(seq_len(nrow(diagnostics)), function(i) {
    bayes_diagnostics_meet_gate(diagnostics[i, , drop = FALSE])
  }, logical(1)))) {
    stop("Canonical Bayesian diagnostics do not satisfy the production quality gate.")
  }

  if (require_fits) {
    if (is.null(input_file)) {
      stop("Validating canonical Bayesian fit metadata requires the current mortality input file.")
    }
    specs <- canonical_bayes_group_specs(fit_effect_sizes, nrs_bias_sd = as.numeric(manifest[["nrs_bias_sd"]]))
    for (group_id in required_group_ids) {
      fit <- readRDS(file.path(paths$fits, paste0(group_id, ".brmsfit.rds")))
      spec <- specs[[group_id]]
      expected_metadata <- make_bayes_model_cache_metadata(
        data = spec$data,
        vi_col = spec$vi_col,
        response_scale = "log_rr",
        prior_spec_id = "log_rr_shared_current",
        uncertainty_spec_id = "primary_nrs_log_rr_bias_sd_0.20",
        settings = manifest_settings,
        study_col = "study",
        yi_col = "yi",
        formula_text = "yi | se(.bayes_se, sigma = FALSE) ~ 1 + (1 | .bayes_study)"
      )
      if (!isTRUE(assert_cached_bayes_model_matches(fit, expected_metadata))) {
        stop("Canonical Bayesian fit metadata or quality gate is incompatible: ", group_id)
      }
    }
  }

  list(
    posterior_summary = posterior_summary,
    diagnostics = diagnostics,
    reml_reference = reml_reference,
    effect_sizes = effect_sizes,
    metadata = metadata,
    paths = paths
  )
}

canonical_summary_row <- function(canonical, group_id) {
  row <- canonical$posterior_summary[canonical$posterior_summary$group_id == group_id, , drop = FALSE]
  if (nrow(row) != 1) {
    stop("Expected exactly one canonical Bayesian summary row for group_id: ", group_id)
  }
  row
}

assert_canonical_effect_sizes_match <- function(current_es,
                                                canonical_es,
                                                tolerance = 1e-10,
                                                cols = c("study", "design", "yi", "vi_original", "vi_bias_adjusted")) {
  missing_current <- setdiff(cols, names(current_es))
  missing_canonical <- setdiff(cols, names(canonical_es))
  if (length(missing_current) > 0 || length(missing_canonical) > 0) {
    stop(
      "Cannot validate canonical Bayesian effect sizes. Missing current columns: ",
      paste(missing_current, collapse = ", "),
      "; missing canonical columns: ",
      paste(missing_canonical, collapse = ", ")
    )
  }

  current <- current_es[, cols, drop = FALSE]
  canonical <- canonical_es[, cols, drop = FALSE]
  if (nrow(current) != nrow(canonical)) {
    stop("Canonical Bayesian effect sizes are stale: row count differs from current mortality.csv.")
  }

  character_cols <- cols[vapply(current, function(x) is.character(x) || is.factor(x), logical(1))]
  for (col in character_cols) {
    if (!identical(as.character(current[[col]]), as.character(canonical[[col]]))) {
      stop("Canonical Bayesian effect sizes are stale: column differs: ", col)
    }
  }

  numeric_cols <- setdiff(cols, character_cols)
  for (col in numeric_cols) {
    diff_value <- max(abs(as.numeric(current[[col]]) - as.numeric(canonical[[col]])), na.rm = TRUE)
    if (!is.finite(diff_value) || diff_value > tolerance) {
      stop("Canonical Bayesian effect sizes are stale: numeric column differs: ", col)
    }
  }

  invisible(TRUE)
}
