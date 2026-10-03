# Build the non-Bayesian plotting caches for Figure 4: publication trends and TSA.

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
})

script_arg <- c(
  grep("^--file=", commandArgs(FALSE), value = TRUE),
  grep("\\.R$", commandArgs(FALSE), value = TRUE, ignore.case = TRUE)
)
script_dir <- if (length(script_arg) > 0) {
  dirname(normalizePath(sub("^--file=", "", script_arg[1]), mustWork = FALSE))
} else {
  getwd()
}
project_root <- normalizePath(file.path(script_dir, "..", ".."), winslash = "/", mustWork = TRUE)
source(file.path(project_root, "analysis_scripts", "common", "bayes_prior_scale_spec.R"))
source(file.path(project_root, "analysis_scripts", "common", "bayes_mortality_functions.R"))
source(file.path(project_root, "analysis_scripts", "common", "tsa_functions.R"))

publication_csv <- file.path(project_root, "source_csv", "narrative_analysis_combined_conclusions.csv")
mortality_csv <- file.path(project_root, "source_csv", "mortality.csv")
cache_dir <- Sys.getenv(
  "BAYES_PUBLICATION_CUMULATIVE_DATA_DIR",
  unset = file.path(project_root, "computed_data", "mortality_publication_cumulative")
)
if (!grepl("^([A-Za-z]:[\\\\/]|/)", cache_dir)) cache_dir <- file.path(project_root, cache_dir)
dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)

alpha <- 0.05
power <- 0.80
rrr <- 0.20
bias_sd <- 0.20

read_publication_data <- function(path) {
  if (!file.exists(path)) stop("Missing input: ", path)
  dat <- read.csv(path, stringsAsFactors = FALSE, check.names = FALSE, encoding = "UTF-8")
  names(dat) <- trimws(sub("^\ufeff", "", names(dat)))
  framework_col <- "Sepsis framework"
  algorithm_col <- "Algorithm class"
  required <- c("Study", "Year", "Design", framework_col, algorithm_col)
  missing <- setdiff(required, names(dat))
  if (length(missing) > 0) {
    stop("narrative_analysis_combined_conclusions.csv is missing: ", paste(missing, collapse = ", "))
  }
  dat$Year <- as.integer(dat$Year)
  dat <- dat[!is.na(dat$Year), , drop = FALSE]
  dat$Design <- toupper(trimws(dat$Design))
  dat$Framework <- trimws(as.character(dat[[framework_col]]))
  dat$Algorithm <- trimws(as.character(dat[[algorithm_col]]))
  dat$Algorithm <- ifelse(tolower(dat$Algorithm) == "rule", "rule", dat$Algorithm)
  dat
}

build_cumulative_data <- function(data, years, group_col, levels, panel, total_levels = NULL) {
  group_sym <- rlang::sym(group_col)
  grouped <- data |>
    mutate(.group = trimws(as.character(!!group_sym))) |>
    filter(.group %in% levels) |>
    count(Year, .group, name = "Annual_count") |>
    complete(Year = years, .group = levels, fill = list(Annual_count = 0)) |>
    arrange(.group, Year) |>
    group_by(.group) |>
    mutate(Cumulative_count = cumsum(Annual_count)) |>
    ungroup() |>
    transmute(Panel = panel, Year, Series = .group, Annual_count, Cumulative_count)

  total_source <- data |>
    mutate(.group = trimws(as.character(!!group_sym)))
  if (!is.null(total_levels)) total_source <- total_source |> filter(.group %in% total_levels)

  total <- total_source |>
    count(Year, name = "Annual_count") |>
    complete(Year = years, fill = list(Annual_count = 0)) |>
    arrange(Year) |>
    mutate(Panel = panel, Series = "Total", Cumulative_count = cumsum(Annual_count)) |>
    select(Panel, Year, Series, Annual_count, Cumulative_count)

  bind_rows(grouped, total)
}

publication <- read_publication_data(publication_csv)
publication_years <- seq(min(publication$Year), max(publication$Year), by = 1)
publication_trends <- bind_rows(
  build_cumulative_data(publication, publication_years, "Design", c("RCT", "NRS"), "Study design"),
  build_cumulative_data(
    publication, publication_years, "Framework", c("2", "3"),
    "Sepsis framework", total_levels = c("2", "3")
  ),
  build_cumulative_data(publication, publication_years, "Algorithm", c("rule", "ML"), "Algorithm class")
)

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
  if (any(!is.finite(out$Year))) stop("Invalid or missing Year values in ", basename(path))

  out$Input_order <- seq_len(nrow(out))
  out$Publication_order <- out$Input_order
  out$Publication_order[out$Study == "Leisman 2024"] <- 202401
  # Explicit within-2024 TSA order requested for Figure 4.
  out$Publication_order[out$Study == "Persson 2024"] <- 202402
  out$Publication_order[out$Study == "Kijpaisalratana 2024"] <- 202403
  out <- out[out$Design %in% c("RCT", "NRS"), , drop = FALSE]
  out[order(out$Year, out$Design != "RCT", out$Publication_order, out$Input_order), , drop = FALSE]
}

compute_rr_effects <- function(df) {
  effects <- prepare_mortality_effect_sizes(df, nrs_bias_sd = bias_sd) |>
    select(study, yi, vi_original, uses_adjusted_rr)
  out <- df |>
    left_join(effects, by = c("Study" = "study"))
  out$logRR <- out$yi
  out$Variance <- out$vi_original
  out$SE <- sqrt(out$Variance)
  out$RR <- exp(out$logRR)
  out$Lower95 <- exp(out$logRR - 1.96 * out$SE)
  out$Upper95 <- exp(out$logRR + 1.96 * out$SE)
  out$Study_N <- as.numeric(out$n_e) + as.numeric(out$n_c)
  out
}

study <- compute_rr_effects(read_mortality_data(mortality_csv))
rct <- study[study$Design == "RCT", , drop = FALSE]
if (nrow(rct) == 0) stop("No RCT rows found in mortality.csv")
meta <- fit_tsa_reml(rct, yi_col = "logRR", vi_col = "Variance", context = "Final mortality TSA")
tsa_cumulative <- cumulative_reml_z(rct, yi_col = "logRR", vi_col = "Variance", n_col = "Study_N")
cer <- sum(as.numeric(rct$events_c)) / sum(as.numeric(rct$n_c))
d_squared <- tsa_diversity_d2(rct$Variance, as.numeric(meta$tau2)[1])
ris <- calculate_tsa_information_size(cer, rrr, alpha, power, d_squared)
tsa_cumulative$Information_fraction <- pmin(tsa_cumulative$Cumulative_N / ris$DARIS[[1]], 1)
tsa_boundary <- calibrate_two_sided_obrien_fleming(tsa_cumulative$Information_fraction, alpha)
tsa_boundary_constant <- attr(tsa_boundary, "boundary_constant")
tsa_cumulative$TSA_upper_boundary <- as.numeric(tsa_boundary)
tsa_cumulative$TSA_lower_boundary <- -tsa_cumulative$TSA_upper_boundary

tsa_summary <- data.frame(
  Metric = c(
    "RCT studies", "RCT total sample size", "DARIS", "Control event rate",
    "Assumed RRR", "D-squared", "REML tau-squared", "O'Brien-Fleming boundary constant"
  ),
  Value = c(
    nrow(rct), sum(rct$Study_N), ris$DARIS[[1]], cer, rrr, d_squared,
    as.numeric(meta$tau2)[1], tsa_boundary_constant
  ),
  stringsAsFactors = FALSE
)

write.csv(publication_trends, file.path(cache_dir, "row_1_publication_trends_data.csv"), row.names = FALSE)
write.csv(tsa_cumulative, file.path(cache_dir, "row_3_tsa_cumulative_data.csv"), row.names = FALSE)
write.csv(tsa_summary, file.path(cache_dir, "row_3_tsa_summary.csv"), row.names = FALSE)

message("Saved Figure 4 publication-trend and TSA caches to: ", cache_dir)
