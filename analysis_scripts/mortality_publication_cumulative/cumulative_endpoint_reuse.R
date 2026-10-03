# Full-evidence endpoints belong to the main/sensitivity analyses, not new fits.
load_publication_endpoints <- function(workspace_dir, input_file, settings, bias_sd = 0.20) {
  canonical_dir <- Sys.getenv("BAYES_CANONICAL_OUT_DIR", unset = "computed_data/mortality_main/canonical")
  if (!grepl("^([A-Za-z]:[\\\\/]|/)", canonical_dir)) canonical_dir <- file.path(workspace_dir, canonical_dir)
  canonical <- load_canonical_bayes_mortality(canonical_dir, require_fits = FALSE,
    input_file = input_file, expected_settings = settings, expected_nrs_bias_sd = bias_sd)
  es <- prepare_mortality_effect_sizes(input_file, nrs_bias_sd = bias_sd)
  assert_canonical_effect_sizes_match(es, canonical$effect_sizes)
  sensitivity_dir <- Sys.getenv("MORTALITY_SENSITIVITY_OUT_DIR", unset = "computed_data/mortality_sensitivity")
  if (!grepl("^([A-Za-z]:[\\\\/]|/)", sensitivity_dir)) sensitivity_dir <- file.path(workspace_dir, sensitivity_dir)
  summary_file <- file.path(sensitivity_dir, "mortality_sensitivity_final_results_summary.csv")
  if (!file.exists(summary_file)) stop("Missing mortality sensitivity summary; build upstream results first.")
  sensitivity <- read.csv(summary_file, stringsAsFactors = FALSE, check.names = FALSE)
  endpoints <- list()
  for (key in c("overall_adjusted", "overall_raw", "rct", "nrs_raw")) {
    primary <- key %in% c("overall_adjusted", "rct")
    data <- if (key == "rct") es[es$design == "RCT", , drop = FALSE] else
      if (key == "nrs_raw") es[es$design == "NRS", , drop = FALSE] else es
    data$endpoint_variance <- if (key == "overall_adjusted") data$vi_bias_adjusted else data$vi_original
    if (primary) {
      group_id <- if (key == "rct") "rct_only" else "overall"
      row <- canonical_summary_row(canonical, group_id)
      source_file <- canonical$paths$posterior_summary
      source_row <- paste0("group_id=", group_id)
      fit_file <- file.path(canonical_dir, "fits", paste0(group_id, ".brmsfit.rds"))
    } else {
      scenario <- if (key == "nrs_raw") "NRS only" else "Total"
      row <- sensitivity[sensitivity$analysis == "No down-weighting reference" &
        sensitivity$scenario == scenario, , drop = FALSE]
      if (nrow(row) != 1L) stop("Expected one no-down-weighting source row for ", key)
      source_file <- summary_file
      source_row <- paste0("analysis=No down-weighting reference; scenario=", scenario)
      fit_file <- file.path(sensitivity_dir, "fits", paste0("no_down_weighting_reference__",
        if (key == "nrs_raw") "nrs_only" else "total", ".brmsfit.rds"))
      if (!file.exists(fit_file)) stop("Missing upstream sensitivity fit: ", fit_file)
      # Check saved studies, variance, prior, sampler and quality metadata.
      # A read-only posterior-summary check detects a stale CSV; no MCMC is run.
      data$vi_adjusted_current <- data$vi_original
      metadata <- make_bayes_model_cache_metadata(data, "vi_adjusted_current", response_scale = "log_rr",
        prior_spec_id = "log_rr_shared_current", uncertainty_spec_id = "sensitivity_defined_variance_log_rr",
        settings = settings, study_col = "study", yi_col = "yi",
        formula_text = "yi | se(.bayes_se, sigma = FALSE) ~ 1 + (1 | .bayes_study)")
      fit <- readRDS(fit_file)
      if (!isTRUE(assert_cached_bayes_model_matches(fit, metadata))) {
        stop("Upstream sensitivity fit is stale or incompatible: ", key, "; rebuild upstream first.")
      }
      check <- summarise_bayes_mortality_fit(fit, group = key,
        k_rct = sum(data$design == "RCT"), k_nrs = sum(data$design == "NRS"), nrs_bias_sd = bias_sd)
      for (field in c("RR", "Lower_95_CrI", "Upper_95_CrI", "tau", "tau2")) {
        if (!isTRUE(all.equal(row[[field]], check[[field]], tolerance = 1e-12))) {
          stop("Upstream sensitivity summary does not match its saved fit: ", key, " / ", field)
        }
      }
    }
    if (row$k != nrow(data) || row$k_rct != sum(data$design == "RCT") ||
        row$k_nrs != sum(data$design == "NRS") || !bayes_diagnostics_meet_gate(row)) {
      stop("Invalid upstream endpoint study counts or diagnostics: ", key)
    }
    endpoints[[key]] <- list(data = data, row = row, fit_file = fit_file,
      source_file = source_file, source_row = source_row, source_md5 = canonical_file_md5(source_file),
      fit_source = if (primary) "canonical_main" else "mortality_sensitivity_cache")
  }
  endpoints
}

cumulative_endpoint_result <- function(data, variance_col, endpoint, model_id, group_label) {
  if (anyDuplicated(data$study) || nrow(data) != nrow(endpoint$data) ||
      !setequal(data$study, endpoint$data$study)) stop("Endpoint reuse requires the same complete study set: ", model_id)
  actual <- data[order(data$study), , drop = FALSE]
  expected <- endpoint$data[order(endpoint$data$study), , drop = FALSE]
  if (!identical(as.character(actual$Design), as.character(expected$design)) ||
      !isTRUE(all.equal(as.numeric(actual$yi), as.numeric(expected$yi), tolerance = 1e-12)) ||
      !isTRUE(all.equal(as.numeric(actual[[variance_col]]), as.numeric(expected$endpoint_variance), tolerance = 1e-12))) {
    stop("Endpoint reuse requires identical design, effects and variances: ", model_id)
  }
  row <- endpoint$row
  fields <- c("rhat_max", "ess_bulk_min", "ess_tail_min", "divergent_transitions", "max_treedepth", "max_treedepth_hits")
  diagnostics <- data.frame(model_id = model_id, group = group_label, k = row$k,
    row[, fields, drop = FALSE], fit_path = endpoint$fit_file, fit_source = endpoint$fit_source,
    summary_file = endpoint$source_file, summary_row = endpoint$source_row, summary_md5 = endpoint$source_md5,
    stringsAsFactors = FALSE, check.names = FALSE)
  list(summary = c(RR = row$RR, Lower95 = row$Lower_95_CrI, Upper95 = row$Upper_95_CrI,
    tau2 = row$tau2, I2 = NA_real_), diagnostics = diagnostics)
}

refresh_publication_endpoints <- function(overall, design, diagnostics, study, endpoints) {
  # Update existing final rows only; missing/stale historical grids fail closed.
  if (anyDuplicated(overall$Year) || !setequal(overall$Year, unique(study$Year))) stop("Cumulative overall year grid is stale.")
  expected_design <- unique(study[, c("Design", "Year")])
  if (anyDuplicated(paste(design$Design, design$Year)) ||
      !setequal(paste(design$Design, design$Year), paste(expected_design$Design, expected_design$Year))) {
    stop("Cumulative design year grid is stale.")
  }
  if (anyDuplicated(diagnostics$model_id)) stop("Duplicate cumulative diagnostic model IDs.")
  for (i in seq_len(nrow(overall))) {
    if (overall$Cumulative_studies[i] != sum(study$Year <= overall$Year[i])) stop("Cumulative overall counts are stale.")
  }
  for (i in seq_len(nrow(design))) {
    if (design$Cumulative_studies[i] != sum(study$Design == design$Design[i] & study$Year <= design$Year[i])) {
      stop("Cumulative design counts are stale.")
    }
  }
  for (key in names(endpoints)) {
    is_overall <- grepl("^overall", key)
    subset <- if (is_overall) study else study[study$Design == if (key == "rct") "RCT" else "NRS", , drop = FALSE]
    year <- max(subset$Year)
    model_id <- paste0(switch(key, overall_adjusted = "overall_bias_adjusted_", overall_raw = "overall_raw_",
      rct = "design_rct_", nrs_raw = "design_nrs_"), year)
    ix <- which(diagnostics$model_id == model_id)
    if (length(ix) != 1L) stop("Missing cumulative endpoint diagnostic: ", model_id)
    result <- cumulative_endpoint_result(subset,
      if (key == "overall_adjusted") "vi_bias_adjusted" else "vi_original",
      endpoints[[key]], model_id, diagnostics$group[ix])
    if (is_overall) {
      row <- which.max(overall$Year)
      prefix <- if (key == "overall_adjusted") "BiasAdjusted" else "Raw"
      overall[row, paste0(prefix, c("_RR", "_Lower95", "_Upper95"))] <- as.list(result$summary[c("RR", "Lower95", "Upper95")])
    } else {
      row <- which(design$Design == subset$Design[1] & design$Year == year)
      design[row, c("RR", "Lower95", "Upper95")] <- as.list(result$summary[c("RR", "Lower95", "Upper95")])
    }
    for (field in setdiff(names(result$diagnostics), names(diagnostics))) diagnostics[[field]] <- NA_character_
    diagnostics[ix, names(result$diagnostics)] <- result$diagnostics
  }
  list(overall = overall, design = design, diagnostics = diagnostics)
}
