# Fit cached Bayesian cumulative mortality models for publication overview.

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(metafor)
})

script_arg <- c(grep("^--file=", commandArgs(FALSE), value = TRUE), grep("\\.R$", commandArgs(FALSE), value = TRUE, ignore.case = TRUE))
script_dir <- if (length(script_arg) > 0) dirname(normalizePath(sub("^--file=", "", script_arg[1]), mustWork = FALSE)) else getwd()
project_root <- normalizePath(file.path(script_dir, "..", ".."), winslash = "/", mustWork = TRUE)
source(file.path(project_root, "analysis_scripts", "common", "project_paths.R"))
options(sepsis_project_root = set_project_root_from_script(project_root))
script_dir <- project_root
workspace_dir <- script_dir
source(file.path(workspace_dir, "analysis_scripts", "common", "bayes_prior_scale_spec.R"))
source(file.path(workspace_dir, "analysis_scripts", "common", "bayes_mortality_functions.R"))
source(file.path(workspace_dir, "analysis_scripts", "mortality_publication_cumulative", "cumulative_endpoint_reuse.R"))

MORTALITY_CSV <- file.path(workspace_dir, "source_csv", "mortality.csv")
CACHE_DIR <- Sys.getenv("BAYES_PUBLICATION_CUMULATIVE_DATA_DIR", unset = "computed_data/mortality_publication_cumulative")
if (!grepl("^([A-Za-z]:[\\\\/]|/)", CACHE_DIR)) CACHE_DIR <- file.path(workspace_dir, CACHE_DIR)
FITS_DIR <- file.path(CACHE_DIR, "fits")
dir.create(FITS_DIR, recursive = TRUE, showWarnings = FALSE)

BIAS_SD <- 0.20
CONTINUITY_CORRECTION <- 0.5
settings <- bayes_meta_settings_from_env()

read_mortality_data <- function(path) {
  if (!file.exists(path)) stop("Missing input: ", path)
  dat <- read_mortality_csv(path, lower_names = FALSE)
  required <- c("Study", "Year", "events_c", "n_c", "events_e", "n_e", "Design")
  missing <- setdiff(required, names(dat))
  if (length(missing) > 0) stop("mortality.csv is missing: ", paste(missing, collapse = ", "))
  optional <- intersect(c("aRR", "aRR_lower_95", "aRR_upper_95"), names(dat))
  out <- dat[, c(required, optional), drop = FALSE]
  out$Year <- as.integer(out$Year)
  out$events_c <- safe_mortality_num(out$events_c)
  out$n_c <- safe_mortality_num(out$n_c)
  out$events_e <- safe_mortality_num(out$events_e)
  out$n_e <- safe_mortality_num(out$n_e)
  out$Design <- toupper(trimws(out$Design))

  validate_mortality_2x2(out, context = basename(path))
  validate_adjusted_rr_fields(out, context = basename(path))
  if (any(!is.finite(out$Year))) {
    stop("Invalid or missing Year values in ", basename(path), "; no rows were silently dropped.")
  }

  out$Input_order <- seq_len(nrow(out))
  out$Publication_order <- out$Input_order
  out$Publication_order[out$Study == "Leisman 2024"] <- 202401
  out$Publication_order[out$Study == "Kijpaisalratana 2024"] <- 202402
  out$Publication_order[out$Study == "Persson 2024"] <- 202403
  out <- out[out$Design %in% c("RCT", "NRS"), , drop = FALSE]
  out[order(out$Year, out$Design != "RCT", out$Publication_order, out$Input_order), , drop = FALSE]
}

compute_rr_effects <- function(df, cc = 0.5) {
  effects <- prepare_mortality_effect_sizes(df, nrs_bias_sd = BIAS_SD) |>
    select(study, yi, vi_original, vi_bias_adjusted, uses_adjusted_rr)
  df |>
    left_join(effects, by = c("Study" = "study")) |>
    mutate(
      study = Study,
      Study_N = as.numeric(n_e) + as.numeric(n_c),
      vi = vi_original
    )
}

single_row_summary <- function(data, variance_col) {
  yi <- data$yi[1]
  se <- sqrt(data[[variance_col]][1])
  c(
    RR = exp(yi),
    Lower95 = exp(yi - qnorm(0.975) * se),
    Upper95 = exp(yi + qnorm(0.975) * se),
    tau2 = NA_real_,
    I2 = NA_real_
  )
}

fit_cached_cumulative <- function(data, model_id, variance_col, group_label, endpoint = NULL) {
  # Full-evidence fits remain owned by upstream even when BAYES_FORCE_REFIT=1.
  if (!is.null(endpoint)) {
    return(cumulative_endpoint_result(data, variance_col, endpoint, model_id, group_label))
  }
  if (nrow(data) == 1) {
    summary <- single_row_summary(data, variance_col)
    diagnostics <- data.frame(
      model_id = model_id,
      group = group_label,
      k = nrow(data),
      rhat_max = NA_real_,
      ess_bulk_min = NA_real_,
      ess_tail_min = NA_real_,
      divergent_transitions = NA_real_,
      fit_path = NA_character_,
      stringsAsFactors = FALSE
    )
    return(list(summary = summary, diagnostics = diagnostics))
  }

  model_data <- data
  model_data$study <- if ("study" %in% names(model_data)) as.character(model_data$study) else as.character(model_data$Study)
  fit_path <- file.path(FITS_DIR, paste0(model_id, ".brmsfit.rds"))
  uncertainty_spec_id <- if (identical(variance_col, "vi_bias_adjusted")) {
    "primary_nrs_log_rr_bias_sd_0.20"
  } else {
    "supplied_sampling_variance_only_log_rr"
  }
  expected_metadata <- make_bayes_model_cache_metadata(
    data = model_data,
    vi_col = variance_col,
    response_scale = "log_rr",
    prior_spec_id = "log_rr_shared_current",
    uncertainty_spec_id = uncertainty_spec_id,
    settings = settings,
    study_col = "study",
    yi_col = "yi",
    formula_text = "yi | se(.bayes_se, sigma = FALSE) ~ 1 + (1 | .bayes_study)"
  )
  force_refit <- identical(Sys.getenv("BAYES_FORCE_REFIT", unset = ""), "1")
  fit <- NULL
  if (file.exists(fit_path) && !force_refit) {
    cached_fit <- readRDS(fit_path)
    if (isTRUE(assert_cached_bayes_model_matches(cached_fit, expected_metadata))) {
      fit <- cached_fit
    } else {
      fit <- upgrade_unattributed_cached_intercept_fit(
        cached_fit,
        expected_data = model_data,
        expected_metadata = expected_metadata
      )
      if (!is.null(fit)) saveRDS(fit, fit_path)
    }
  }
  if (is.null(fit)) {
    fit <- fit_bayes_mortality_meta_with_quality_gate(
      model_data,
      variance_col,
      response_scale = "log_rr",
      prior_spec_id = "log_rr_shared_current",
      uncertainty_spec_id = uncertainty_spec_id,
      settings = settings
    )
    saveRDS(fit, fit_path)
  }
  posterior_summary <- summarise_bayes_mortality_fit(
    fit,
    group = group_label,
    k_rct = sum(model_data$Design == "RCT", na.rm = TRUE),
    k_nrs = sum(model_data$Design == "NRS", na.rm = TRUE),
    nrs_bias_sd = BIAS_SD
  )
  diagnostics <- extract_bayes_diagnostics(fit) |>
    mutate(
      model_id = model_id,
      group = group_label,
      k = nrow(model_data),
      fit_path = fit_path,
      .before = 1
    )
  list(
    summary = c(
      RR = posterior_summary$RR,
      Lower95 = posterior_summary$Lower_95_CrI,
      Upper95 = posterior_summary$Upper_95_CrI,
      tau2 = posterior_summary$tau2,
      I2 = NA_real_
    ),
    diagnostics = diagnostics
  )
}

study <- compute_rr_effects(read_mortality_data(MORTALITY_CSV), CONTINUITY_CORRECTION)
endpoints <- load_publication_endpoints(workspace_dir, MORTALITY_CSV, settings, BIAS_SD)
endpoints_only <- identical(Sys.getenv("BAYES_PUBLICATION_CUMULATIVE_ENDPOINTS_ONLY", unset = ""), "1")
if (endpoints_only) {
  read_cache <- function(name) read.csv(file.path(CACHE_DIR, name), stringsAsFactors = FALSE, check.names = FALSE)
  refreshed <- refresh_publication_endpoints(
    read_cache("row_2_overall_cumulative_effect_data.csv"),
    read_cache("row_2_design_specific_effect_data.csv"),
    read_cache("publication_cumulative_bayes_diagnostics.csv"), study, endpoints)
  overall_results <- refreshed$overall
  design_results <- refreshed$design
  diagnostics <- refreshed$diagnostics
} else {
diagnostic_rows <- list()
overall_rows <- list()
design_rows <- list()

for (year in sort(unique(study$Year))) {
  cumulative <- study[study$Year <= year, , drop = FALSE]
  raw <- fit_cached_cumulative(cumulative, paste0("overall_raw_", year), "vi_original", paste0("Overall raw <= ", year),
    endpoint = if (year == max(study$Year)) endpoints$overall_raw else NULL)
  adjusted <- fit_cached_cumulative(cumulative, paste0("overall_bias_adjusted_", year), "vi_bias_adjusted", paste0("Overall bias-adjusted <= ", year),
    endpoint = if (year == max(study$Year)) endpoints$overall_adjusted else NULL)
  diagnostic_rows[[length(diagnostic_rows) + 1]] <- raw$diagnostics
  diagnostic_rows[[length(diagnostic_rows) + 1]] <- adjusted$diagnostics
  overall_rows[[length(overall_rows) + 1]] <- data.frame(
    Year = year,
    Cumulative_studies = nrow(cumulative),
    Raw_RR = raw$summary[["RR"]],
    Raw_Lower95 = raw$summary[["Lower95"]],
    Raw_Upper95 = raw$summary[["Upper95"]],
    BiasAdjusted_RR = adjusted$summary[["RR"]],
    BiasAdjusted_Lower95 = adjusted$summary[["Lower95"]],
    BiasAdjusted_Upper95 = adjusted$summary[["Upper95"]]
  )
}

for (design in c("RCT", "NRS")) {
  design_data <- study[study$Design == design, , drop = FALSE]
  for (year in sort(unique(design_data$Year))) {
    cumulative <- design_data[design_data$Year <= year, , drop = FALSE]
    fit <- fit_cached_cumulative(cumulative, paste0("design_", tolower(design), "_", year), "vi", paste0(design, " <= ", year),
      endpoint = if (year == max(design_data$Year)) endpoints[[if (design == "RCT") "rct" else "nrs_raw"]] else NULL)
    diagnostic_rows[[length(diagnostic_rows) + 1]] <- fit$diagnostics
    design_rows[[length(design_rows) + 1]] <- data.frame(
      Design = design,
      Year = year,
      Cumulative_studies = nrow(cumulative),
      RR = fit$summary[["RR"]],
      Lower95 = fit$summary[["Lower95"]],
      Upper95 = fit$summary[["Upper95"]]
    )
  }
}

overall_results <- bind_rows(overall_rows)
design_results <- bind_rows(design_rows)
diagnostics <- bind_rows(diagnostic_rows)
}

write.csv(overall_results, file.path(CACHE_DIR, "row_2_overall_cumulative_effect_data.csv"), row.names = FALSE)
write.csv(design_results, file.path(CACHE_DIR, "row_2_design_specific_effect_data.csv"), row.names = FALSE)
write.csv(diagnostics, file.path(CACHE_DIR, "publication_cumulative_bayes_diagnostics.csv"), row.names = FALSE)
if (!endpoints_only) write.csv(
  data.frame(
    key = c("created_at", "iter", "warmup", "chains", "seed", "adapt_delta", "max_treedepth"),
    value = c(
      format(Sys.time(), "%Y-%m-%d %H:%M:%S %z"),
      settings$iter,
      settings$warmup,
      settings$chains,
      settings$seed,
      settings$adapt_delta,
      settings$max_treedepth
    )
  ),
  file.path(CACHE_DIR, "publication_cumulative_cache_metadata.csv"),
  row.names = FALSE
)

message("Saved publication cumulative Bayesian cache to: ", CACHE_DIR)
