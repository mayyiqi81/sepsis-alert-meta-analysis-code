# Process-outcome subgroup analyses (calculation only).

library(metafor)
library(dplyr)

script_arg <- c(grep("^--file=", commandArgs(FALSE), value = TRUE), grep("\\.R$", commandArgs(FALSE), value = TRUE, ignore.case = TRUE))
script_dir <- if (length(script_arg) > 0) dirname(normalizePath(sub("^--file=", "", script_arg[1]), mustWork = FALSE)) else getwd()
project_root <- normalizePath(file.path(script_dir, "..", ".."), winslash = "/", mustWork = TRUE)
script_dir <- project_root
out_dir <- file.path(project_root, "computed_data", "process_subgroup")
fit_control <- list(maxiter = 10000, stepadj = 0.1)
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

clean_value <- function(x) {
  x <- trimws(as.character(x))
  x[x %in% c("", "/", "NA", "Na", "na", "N/A")] <- NA_character_
  x
}
safe_num <- function(x) suppressWarnings(as.numeric(trimws(as.character(x))))

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

safe_model_value <- function(x) ifelse(is.null(x) || length(x) == 0, NA_real_, as.numeric(x)[1])
fmt_effect <- function(x, lower, upper, digits = 3) sprintf(paste0("%.", digits, "f [%.", digits, "f, %.", digits, "f]"), x, lower, upper)
fmt_p <- function(x) {
  if (is.na(x)) return("")
  if (x < 0.001) return("<0.001")
  sprintf("%.3f", x)
}
fmt_p_subgroup_label <- function(p_subgroup) {
  p_subgroup_text <- fmt_p(p_subgroup)
  if (p_subgroup_text == "") return("")
  paste0("P for subgroup difference = ", p_subgroup_text)
}
fmt_heterogeneity <- function(i2, p) {
  if (is.na(i2) && is.na(p)) return("")
  paste0("I2=", sprintf("%.0f", i2), "%; p=", fmt_p(p))
}
fmt_i2 <- function(i2) {
  if (is.na(i2)) return("")
  paste0(sprintf("%.1f", i2), "%")
}

fit_random <- function(dat, context_label = "unknown subgroup") {
  if (nrow(dat) < 1) stop("No studies retained for subgroup: ", context_label)
  if (nrow(dat) == 1) {
    yi <- dat$yi[1]; se <- sqrt(dat$vi[1])
    return(list(k = 1, b = yi, ci.lb = yi - 1.96 * se, ci.ub = yi + 1.96 * se, tau2 = NA_real_, I2 = NA_real_, QE = NA_real_, QEp = NA_real_, p = 1, method_used = "single study"))
  }
  fit <- try(rma(yi = yi, vi = vi, data = dat, method = "REML", control = fit_control), silent = TRUE)
  if (inherits(fit, "try-error")) {
    stop(
      "Process subgroup REML model failed to converge; ",
      "no alternative heterogeneity estimator was substituted. ",
      "context = ", context_label, "; k = ", nrow(dat),
      "; retained studies = ", paste(dat$study, collapse = "; ")
    )
  }
  fit$method_used <- "REML"
  fit
}

subgroup_p <- function(dat, group_col) {
  groups <- dat[[group_col]]
  if (length(unique(groups)) < 2 || nrow(dat) <= length(unique(groups))) return(NA_real_)
  fit <- try(rma(yi = yi, vi = vi, mods = ~ factor(groups), data = dat, method = "REML", control = fit_control), silent = TRUE)
  if (inherits(fit, "try-error")) return(NA_real_)
  safe_model_value(fit$QMp)
}

collapse_setting <- function(x) {
  dplyr::case_when(
    x %in% c("Wards", "ICU") ~ "Wards/ICU",
    x == "ED" ~ "ED",
    TRUE ~ NA_character_
  )
}

read_binary <- function(file, outcome) {
  raw <- read.csv(file, stringsAsFactors = FALSE, check.names = FALSE, encoding = "UTF-8")
  names(raw) <- trimws(sub("^\ufeff", "", names(raw)))
  dat <- raw %>%
    transmute(
      outcome = outcome,
      study = clean_value(Study),
      design = toupper(clean_value(Design)),
      events_e = safe_num(`Intervention events`),
      n_e = safe_num(`Intervention total`),
      events_c = safe_num(`Control events`),
      n_c = safe_num(`Control total`),
      baseline_subgroup = clean_value(`Control baseline subgroup`),
      setting = collapse_setting(clean_value(`Deployment setting`)),
      adjusted_rd = safe_num(`Adjusted RD, percentage points`),
      adjusted_lower_95 = safe_num(`Adjusted lower 95, percentage points`),
      adjusted_upper_95 = safe_num(`Adjusted upper 95, percentage points`),
      effect_type = "RD percentage points"
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
      yi = selected_rd / 100,
      vi = (selected_se / 100)^2
    ) %>%
    filter(!is.na(yi), !is.na(vi), vi > 0) %>%
    as.data.frame()
}

read_time <- function(file, outcome) {
  raw <- read.csv(file, stringsAsFactors = FALSE, check.names = FALSE, encoding = "UTF-8")
  names(raw) <- trimws(sub("^\ufeff", "", names(raw)))
  raw %>%
    transmute(
      outcome = outcome,
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
      baseline_subgroup = clean_value(`Control baseline subgroup`),
      setting = collapse_setting(clean_value(`Deployment setting`)),
      time_origin = clean_value(`Antibiotic time origin`),
      effect_type = "MD hours"
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

specs <- list(
  list(file = "process_lactate_measurement.csv", outcome = "Lactate measurement", type = "rd"),
  list(file = "process_blood_culture.csv", outcome = "Blood culture", type = "rd"),
  list(file = "process_antibiotic_proportion.csv", outcome = "Antibiotic proportion", type = "rd"),
  list(file = "process_antibiotic_time.csv", outcome = "Antibiotic time", type = "md"),
  list(file = "process_fluid_resuscitation.csv", outcome = "Fluid resuscitation", type = "rd"),
  list(file = "process_sepsis_bundle.csv", outcome = "Sepsis bundle completion", type = "rd")
)
specs <- lapply(specs, function(spec) {
  spec$file <- file.path(script_dir, "source_csv", spec$file)
  spec
})

display_outcome_order <- c(
  "Lactate measurement",
  "Blood culture",
  "Antibiotic proportion",
  "Fluid resuscitation",
  "Sepsis bundle completion",
  "Antibiotic time"
)

display_outcome_labels <- c(
  "Lactate measurement" = "Lactate measurement",
  "Blood culture" = "Blood culture",
  "Antibiotic proportion" = "Timely antibiotics",
  "Antibiotic time" = "Time to antibiotics",
  "Fluid resuscitation" = "Fluid resuscitation",
  "Sepsis bundle completion" = "Sepsis bundle"
)

display_outcome_label <- function(outcome) {
  if (outcome %in% names(display_outcome_labels)) {
    return(unname(display_outcome_labels[[outcome]]))
  }
  outcome
}

all_dat <- list()
for (spec in specs) {
  dat <- if (spec$type == "rd") read_binary(spec$file, spec$outcome) else read_time(spec$file, spec$outcome)
  all_dat[[length(all_dat) + 1]] <- dat
}
all_dat <- bind_rows(all_dat)

format_group_order <- function(group_type, effect_type) {
  if (group_type == "Control baseline subgroup" && effect_type == "RD percentage points") {
    c("Low baseline (<33.3%)", "Moderate baseline (33.3% to <66.7%)", "High baseline (>=66.7%)")
  } else if (group_type == "Control baseline subgroup") {
    c("Short baseline time (<2 h)", "Moderate baseline time (2 to <6 h)", "Long baseline time (>=6 h)", "NR")
  } else if (group_type == "Time-origin subgroup") {
    c("Arrival, triage, or presentation", "Clinical recognition or sepsis time origin", "Alert or screen trigger", "Other or unclear")
  } else {
    c("ED", "Wards/ICU")
  }
}

summarise_group <- function(dat, outcome, effect_type, subgroup_type, group_name, group_col, subgroup_p_value) {
  fit <- fit_random(dat, paste(outcome, subgroup_type, group_name, sep = " / "))
  tau2 <- safe_model_value(fit$tau2)
  weights <- 1 / (dat$vi + tau2)
  subgroup_weight <- sum(weights, na.rm = TRUE)
  if (effect_type == "RD percentage points") {
    estimate <- 100 * safe_model_value(fit$b); lower <- 100 * safe_model_value(fit$ci.lb); upper <- 100 * safe_model_value(fit$ci.ub)
    label <- fmt_effect(estimate, lower, upper, 2)
  } else {
    estimate <- safe_model_value(fit$b); lower <- safe_model_value(fit$ci.lb); upper <- safe_model_value(fit$ci.ub)
    label <- fmt_effect(estimate, lower, upper, 2)
  }
  data.frame(
    outcome = outcome,
    subgroup_type = subgroup_type,
    subgroup = group_name,
    k = nrow(dat),
    k_rct = sum(dat$design == "RCT"),
    k_nrs = sum(dat$design == "NRS"),
    effect = effect_type,
    estimate = estimate,
    lower_95 = lower,
    upper_95 = upper,
    effect_label = label,
    tau2 = safe_model_value(fit$tau2),
    I2_percent = safe_model_value(fit$I2),
    p_heterogeneity = safe_model_value(fit$QEp),
    heterogeneity_label = fmt_heterogeneity(safe_model_value(fit$I2), safe_model_value(fit$QEp)),
    i2_label = fmt_i2(safe_model_value(fit$I2)),
    p_subgroup = subgroup_p_value,
    method = fit$method_used,
    subgroup_weight = subgroup_weight,
    studies = paste(dat$study, collapse = "; "),
    stringsAsFactors = FALSE
  )
}

results <- list()
for (outcome_name in unique(all_dat$outcome)) {
  odat <- all_dat %>% filter(outcome == outcome_name)
  effect_type <- unique(odat$effect_type)[1]
  subgroup_types <- c("Control baseline subgroup", "Deployment setting")
  if (outcome_name == "Antibiotic time") subgroup_types <- c(subgroup_types, "Time-origin subgroup")
  for (subgroup_type in subgroup_types) {
    group_col <- dplyr::case_when(
      subgroup_type == "Deployment setting" ~ "setting",
      subgroup_type == "Time-origin subgroup" ~ "time_origin",
      TRUE ~ "baseline_subgroup"
    )
    dat <- odat %>% filter(!is.na(.data[[group_col]]), .data[[group_col]] != "NR")
    p_int <- subgroup_p(dat, group_col)
    orders <- format_group_order(subgroup_type, effect_type)
    group_names <- intersect(orders, unique(dat[[group_col]]))
    group_names <- c(group_names, setdiff(unique(dat[[group_col]]), group_names))
    for (group_name in group_names) {
      gdat <- dat %>% filter(.data[[group_col]] == group_name)
      results[[length(results) + 1]] <- summarise_group(gdat, outcome_name, effect_type, subgroup_type, group_name, group_col, p_int)
    }
  }
}
summary_table <- bind_rows(results) %>%
  arrange(match(outcome, display_outcome_order))
write.csv(summary_table, file.path(out_dir, "process_outcomes_subgroup_summary.csv"), row.names = FALSE)

pvals <- summary_table %>%
  group_by(outcome, subgroup_type) %>%
  summarise(
    p_subgroup = first(p_subgroup),
    p_subgroup_formatted = fmt_p(first(p_subgroup)),
    p_subgroup_label = fmt_p_subgroup_label(first(p_subgroup)),
    groups = paste(subgroup, collapse = "; "),
    .groups = "drop"
  )
pvals <- pvals %>% arrange(match(outcome, display_outcome_order))
write.csv(pvals, file.path(out_dir, "process_outcomes_subgroup_pvalues.csv"), row.names = FALSE)


