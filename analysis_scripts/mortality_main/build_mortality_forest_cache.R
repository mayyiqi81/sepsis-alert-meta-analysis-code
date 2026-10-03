# Mortality forest-plot analysis cache builder.

library(metafor)
library(dplyr)

script_arg <- c(grep("^--file=", commandArgs(FALSE), value = TRUE), grep("\\.R$", commandArgs(FALSE), value = TRUE, ignore.case = TRUE))
script_dir <- if (length(script_arg) > 0) dirname(normalizePath(sub("^--file=", "", script_arg[1]), mustWork = FALSE)) else getwd()
project_root <- normalizePath(file.path(script_dir, "..", ".."), winslash = "/", mustWork = TRUE)
source(file.path(project_root, "analysis_scripts", "common", "bayes_prior_scale_spec.R"))
source(file.path(project_root, "analysis_scripts", "common", "bayes_mortality_functions.R"))
input_file <- file.path(project_root, "source_csv", "mortality.csv")
data_out_dir <- file.path(project_root, "computed_data", "mortality_main")
canonical_out_dir <- file.path(data_out_dir, "canonical")
out_dir <- data_out_dir
nrs_bias_sd <- 0.20
nrs_bias_sd_label <- sprintf("%.2f", nrs_bias_sd)
if (!dir.exists(data_out_dir)) dir.create(data_out_dir, recursive = TRUE)

clean_value <- function(x) {
  x <- trimws(as.character(x))
  x[x %in% c("", "/", "NA", "Na", "na", "N/A")] <- NA_character_
  x
}

fmt_rr <- function(rr, lower, upper) {
  sprintf("%.3f [%.3f, %.3f]", rr, lower, upper)
}

fit_rma <- function(data, vi_col) {
  rma(
    yi = yi,
    vi = data[[vi_col]],
    data = data,
    method = "REML",
    control = list(maxiter = 10000, stepadj = 0.1)
  )
}

format_output_3dp <- function(dat) {
  count_cols <- c("k", "k_rct", "k_nrs", "Q_df")
  out <- dat
  for (nm in names(out)) {
    if (is.numeric(out[[nm]])) {
      if (tolower(nm) %in% tolower(count_cols)) {
        out[[nm]] <- ifelse(is.na(out[[nm]]), "", sprintf("%.0f", out[[nm]]))
      } else {
        out[[nm]] <- ifelse(is.na(out[[nm]]), "", sprintf("%.3f", out[[nm]]))
      }
    }
  }
  out
}

output_path <- function(filename) {
  file.path(out_dir, filename)
}

write_output_csv <- function(dat, filename) {
  write.csv(dat, file.path(data_out_dir, filename), row.names = FALSE)
}

df_raw <- read.csv(input_file, stringsAsFactors = FALSE, check.names = FALSE)
names(df_raw) <- tolower(trimws(names(df_raw)))
if (length(names(df_raw)) > 0 && grepl("study$", names(df_raw)[1], ignore.case = TRUE)) names(df_raw)[1] <- "study"
df_raw <- df_raw[, names(df_raw) != "", drop = FALSE]

required_cols <- c("study", "events_e", "n_e", "events_c", "n_c", "design")
missing_cols <- setdiff(required_cols, names(df_raw))
if (length(missing_cols) > 0) {
  stop("Missing required columns: ", paste(missing_cols, collapse = ", "))
}

df_all <- df_raw[, required_cols] %>%
  mutate(
    row_order = row_number(),
    study = clean_value(study),
    design_raw = clean_value(design),
    design = toupper(design_raw),
    events_e = as.numeric(events_e),
    n_e = as.numeric(n_e),
    events_c = as.numeric(events_c),
    n_c = as.numeric(n_c)
  )

excluded <- df_all %>%
  filter(
    is.na(study) | is.na(events_e) | is.na(n_e) | is.na(events_c) | is.na(n_c) |
      !design %in% c("RCT", "NRS")
  )

df <- df_all %>%
  filter(
    !is.na(study),
    !is.na(events_e), !is.na(n_e),
    !is.na(events_c), !is.na(n_c),
    design %in% c("RCT", "NRS")
  )

es <- prepare_mortality_effect_sizes(input_file, nrs_bias_sd = nrs_bias_sd) %>%
  mutate(rr_label = fmt_rr(RR, Lower_95, Upper_95))

canonical <- load_canonical_bayes_mortality(
  canonical_out_dir,
  require_fits = FALSE,
  input_file = input_file,
  expected_nrs_bias_sd = nrs_bias_sd
)
assert_canonical_effect_sizes_match(es, canonical$effect_sizes)

fit_all_het <- fit_rma(es, "vi_original")
fit_rct_het <- fit_rma(es %>% filter(design == "RCT"), "vi_original")
fit_nrs_het <- fit_rma(es %>% filter(design == "NRS"), "vi_original")

fit_all_het_bias <- fit_rma(es, "vi_bias_adjusted")
fit_nrs_het_bias <- fit_rma(es %>% filter(design == "NRS"), "vi_bias_adjusted")

summary_rct <- canonical_summary_row(canonical, "rct_only")
summary_nrs <- canonical_summary_row(canonical, "nrs_only")
summary_all <- canonical_summary_row(canonical, "overall")

es <- es %>%
  mutate(
    weight_random = 100 * (1 / (vi_bias_adjusted + summary_all$tau2)) /
      sum(1 / (vi_bias_adjusted + summary_all$tau2)),
    weight_label = sprintf("%.1f", weight_random)
  )

summary_table <- bind_rows(
  summary_rct,
  summary_nrs,
  summary_all
)

heterogeneity_table <- bind_rows(
  summary_rct %>% mutate(variance_basis = "Bayesian posterior"),
  summary_nrs %>% mutate(variance_basis = "Bayesian posterior"),
  summary_all %>% mutate(variance_basis = "Bayesian posterior")
) %>%
  select(group, variance_basis, k, k_rct, k_nrs, tau, tau2, rhat_max, ess_bulk_min, ess_tail_min, divergent_transitions, nrs_bias_sd)

reml_reference_table <- bind_rows(
  summarise_reml_mortality_reference(fit_rct_het, "RCT original vi", sum(es$design == "RCT"), 0, nrs_bias_sd),
  summarise_reml_mortality_reference(fit_nrs_het, paste0("NRS original vi with bias SD ", nrs_bias_sd_label), 0, sum(es$design == "NRS"), nrs_bias_sd),
  summarise_reml_mortality_reference(fit_nrs_het_bias, paste0("NRS vi + bias variance ", nrs_bias_sd_label), 0, sum(es$design == "NRS"), nrs_bias_sd),
  summarise_reml_mortality_reference(fit_all_het, paste0("Overall original vi with bias SD ", nrs_bias_sd_label), sum(es$design == "RCT"), sum(es$design == "NRS"), nrs_bias_sd),
  summarise_reml_mortality_reference(fit_all_het_bias, paste0("Overall vi + bias variance ", nrs_bias_sd_label), sum(es$design == "RCT"), sum(es$design == "NRS"), nrs_bias_sd)
)

write_output_csv(
  es %>%
    select(
      study, design, events_e, n_e, events_c, n_c,
      yi, vi_original, vi_bias_adjusted, se_original, se_bias_adjusted,
      RR, Lower_95, Upper_95
    ),
  paste0("mortality_meta_analysis_nrs_bias_sd_", nrs_bias_sd_label, "_effect_sizes.csv")
)
write_output_csv(
  format_output_3dp(summary_table),
  paste0("mortality_meta_analysis_nrs_bias_sd_", nrs_bias_sd_label, "_summary.csv")
)
write_output_csv(
  format_output_3dp(heterogeneity_table),
  paste0("mortality_meta_analysis_nrs_bias_sd_", nrs_bias_sd_label, "_heterogeneity.csv")
)
write_output_csv(
  format_output_3dp(reml_reference_table),
  paste0("mortality_meta_analysis_nrs_bias_sd_", nrs_bias_sd_label, "_reml_reference.csv")
)
write_output_csv(
  format_output_3dp(canonical$diagnostics),
  paste0("mortality_meta_analysis_nrs_bias_sd_", nrs_bias_sd_label, "_bayes_diagnostics.csv")
)
write_output_csv(
  excluded,
  paste0("mortality_meta_analysis_nrs_bias_sd_", nrs_bias_sd_label, "_excluded.csv")
)

write_output_csv(
  es |> arrange(row_order),
  paste0("mortality_meta_analysis_nrs_bias_sd_", nrs_bias_sd_label, "_plot_data.csv")
)
message("Wrote mortality forest-plot analysis cache to: ", data_out_dir)
