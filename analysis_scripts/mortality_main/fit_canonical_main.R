# Fit canonical Bayesian mortality main-analysis models.
# This script runs MCMC and saves reusable brmsfit/posterior artifacts.

library(metafor)
library(dplyr)

script_arg <- c(grep("^--file=", commandArgs(FALSE), value = TRUE), grep("\\.R$", commandArgs(FALSE), value = TRUE, ignore.case = TRUE))
script_dir <- if (length(script_arg) > 0) dirname(normalizePath(sub("^--file=", "", script_arg[1]), mustWork = FALSE)) else getwd()
project_root <- normalizePath(file.path(script_dir, "..", ".."), winslash = "/", mustWork = TRUE)
source(file.path(project_root, "analysis_scripts", "common", "project_paths.R"))
options(sepsis_project_root = set_project_root_from_script(project_root))
script_dir <- project_root
workspace_dir <- script_dir
source(file.path(workspace_dir, "analysis_scripts", "common", "bayes_prior_scale_spec.R"))
source(file.path(workspace_dir, "analysis_scripts", "common", "bayes_mortality_functions.R"))

input_file <- file.path(workspace_dir, "source_csv", "mortality.csv")
out_dir <- Sys.getenv(
  "BAYES_CANONICAL_OUT_DIR",
  unset = "computed_data/mortality_main/canonical"
)
if (!grepl("^([A-Za-z]:[\\\\/]|/)", out_dir)) out_dir <- file.path(workspace_dir, out_dir)

nrs_bias_sd <- 0.20
bayes_settings <- bayes_meta_settings_from_env()

clean_value <- function(x) {
  x <- trimws(as.character(x))
  x[x %in% c("", "/", "NA", "Na", "na", "N/A")] <- NA_character_
  x
}

safe_num <- function(x) {
  suppressWarnings(as.numeric(gsub("%", "", trimws(as.character(x)), fixed = TRUE)))
}

es <- prepare_mortality_effect_sizes(input_file, nrs_bias_sd = nrs_bias_sd)

message("Preparing canonical Bayesian mortality models:")
message("  RCT only: ", sum(es$design == "RCT"), " studies")
message("  NRS only: ", sum(es$design == "NRS"), " studies")
message("  Overall: ", nrow(es), " studies")
message("  Output: ", out_dir)

force_refit <- identical(Sys.getenv("BAYES_FORCE_REFIT", unset = ""), "1")
canonical_cache <- NULL
if (!force_refit && dir.exists(out_dir)) {
  canonical_cache <- tryCatch(
    load_canonical_bayes_mortality(
      out_dir,
      require_fits = TRUE,
      input_file = input_file,
      expected_settings = bayes_settings,
      expected_nrs_bias_sd = nrs_bias_sd
    ),
    error = function(e) NULL
  )
  if (!is.null(canonical_cache)) {
    cache_matches <- tryCatch({
      assert_canonical_effect_sizes_match(es, canonical_cache$effect_sizes)
      TRUE
    }, error = function(e) FALSE)
    if (!cache_matches) canonical_cache <- NULL
  }
}

if (is.null(canonical_cache)) {
  message("Canonical cache is absent, stale, or invalid; fitting three models.")
  fit_and_save_canonical_bayes_mortality(
    es = es,
    out_dir = out_dir,
    settings = bayes_settings,
    nrs_bias_sd = nrs_bias_sd,
    input_file = input_file
  )
} else {
  message("Reused validated canonical Bayesian mortality cache.")
}

message("Saved canonical Bayesian mortality artifacts to: ", out_dir)
