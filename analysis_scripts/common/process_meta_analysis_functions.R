# Shared calculations for the six primary process-outcome meta-analyses.

library(metafor)
library(dplyr)

process_titles <- c(
  "Lactate measurement completion proportion",
  "Blood culture collection completion proportion",
  "Timely antibiotic administration proportion",
  "Fluid resuscitation completion proportion",
  "Overall sepsis bundle completion proportion",
  "Time to antibiotic administration"
)

kim_2021_fluid_exclusion_sentence <- paste(
  "Kim 2021 was not meta-analysed because both groups had 100% event rates,",
  "resulting in a non-informative effect estimate."
)

clean_value <- function(x) {
  x <- trimws(as.character(x))
  x[x %in% c("", "/", "NA", "Na", "na", "N/A")] <- NA_character_
  x
}

safe_num <- function(x) {
  suppressWarnings(as.numeric(trimws(as.character(x))))
}

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
      selected_lower_95 = ifelse(use_adjusted_rd, adjusted_lower_95, selected_rd - 1.96 * selected_se),
      selected_upper_95 = ifelse(use_adjusted_rd, adjusted_upper_95, selected_rd + 1.96 * selected_se),
      effect_basis = ifelse(use_adjusted_rd, "Adjusted RD from 4.3 double-dagger-marked row", "Crude RD from 2x2 counts")
    )
}

safe_model_value <- function(x) {
  ifelse(is.null(x) || length(x) == 0, NA_real_, as.numeric(x)[1])
}

fmt_effect_2dp <- function(x, lower, upper) {
  sprintf("%.2f [%.2f, %.2f]", x, lower, upper)
}

fmt_p <- function(p) {
  ifelse(is.na(p), "", ifelse(p < 0.001, "<0.001", sprintf("%.3f", p)))
}

format_output_2dp <- function(dat) {
  count_cols <- c("k", "n_e", "n_c", "events_e", "events_c", "q_df")
  out <- dat
  for (nm in names(out)) {
    if (is.numeric(out[[nm]])) {
      if (tolower(nm) %in% count_cols) {
        out[[nm]] <- ifelse(is.na(out[[nm]]), "", sprintf("%.0f", out[[nm]]))
      } else {
        out[[nm]] <- ifelse(is.na(out[[nm]]), "", sprintf("%.2f", out[[nm]]))
      }
    }
  }
  out
}

parse_iqr_sd <- function(x) {
  vals <- regmatches(as.character(x), gregexpr("[0-9]+\\.?[0-9]*", as.character(x)))[[1]]
  vals <- as.numeric(vals)
  if (length(vals) < 2) return(NA_real_)
  (vals[2] - vals[1]) / 1.35
}

disp_to_sd <- function(value, type, n) {
  type <- toupper(trimws(as.character(type)))
  if (type == "SD") return(safe_num(value))
  if (type == "SE") return(safe_num(value) * sqrt(safe_num(n)))
  if (type == "IQR") return(parse_iqr_sd(value))
  safe_num(value)
}

fit_random <- function(es) {
  if (nrow(es) == 1) {
    yi <- es$yi[1]
    se <- sqrt(es$vi[1])
    return(list(
      k = 1,
      b = yi,
      ci.lb = yi - 1.96 * se,
      ci.ub = yi + 1.96 * se,
      tau2 = NA_real_,
      I2 = NA_real_,
      QE = NA_real_,
      QEp = NA_real_,
      H2 = NA_real_,
      p = 1,
      method_used = "single study"
    ))
  }

  fit <- try(
    rma(yi = yi, vi = vi, data = es, method = "REML", control = list(maxiter = 1000)),
    silent = TRUE
  )
  if (inherits(fit, "try-error")) {
    fit <- try(
      rma(yi = yi, vi = vi, data = es, method = "REML", control = list(maxiter = 10000)),
      silent = TRUE
    )
    if (inherits(fit, "try-error")) {
      stop(
        "Process outcome REML model failed to converge; ",
        "no alternative heterogeneity estimator was substituted. ",
        "outcome = ", paste(unique(es$outcome), collapse = "; "),
        "; k = ", nrow(es), "; retained studies = ",
        paste(es$study, collapse = "; ")
      )
    }
  }
  fit$method_used <- "REML"
  fit
}

fit_random_if_possible <- function(es) {
  if (nrow(es) < 2) return(NULL)
  fit_random(es)
}

fmt_footer_p <- function(p) {
  if (is.na(p)) {
    "p = NA"
  } else if (p < 0.001) {
    "p < 0.001"
  } else {
    paste0("p = ", sprintf("%.3f", p))
  }
}

fmt_footer_heterogeneity <- function(fit) {
  if (is.null(fit)) {
    "  Heterogeneity: not estimable (k < 2)"
  } else {
    paste0(
      "  \u03C4\u00B2 = ", sprintf("%.3f", safe_model_value(fit$tau2)),
      "; \u03C7\u00B2 = ", sprintf("%.2f", safe_model_value(fit$QE)),
      ", df = ", sprintf("%.0f", safe_model_value(fit$k - fit$p)),
      " (", fmt_footer_p(safe_model_value(fit$QEp)), ")",
      "; I\u00B2 = ", sprintf("%.1f", safe_model_value(fit$I2)), "%"
    )
  }
}

blank_cols <- function(n_cols) {
  cols <- as.list(rep("", n_cols))
  names(cols) <- paste0("col", seq_len(n_cols))
  cols
}

header_row <- function(label, n_cols) {
  data.frame(
    row_type = "header",
    label = label,
    blank_cols(n_cols),
    mean = NA_real_,
    lower = NA_real_,
    upper = NA_real_,
    effect_label = "",
    stringsAsFactors = FALSE
  )
}

heterogeneity_row <- function(fit, n_cols) {
  data.frame(
    row_type = "heterogeneity",
    label = fmt_footer_heterogeneity(fit),
    blank_cols(n_cols),
    mean = NA_real_,
    lower = NA_real_,
    upper = NA_real_,
    effect_label = "",
    stringsAsFactors = FALSE
  )
}

note_row <- function(label, n_cols) {
  data.frame(
    row_type = "note",
    label = label,
    blank_cols(n_cols),
    mean = NA_real_,
    lower = NA_real_,
    upper = NA_real_,
    effect_label = "",
    stringsAsFactors = FALSE
  )
}

summary_from_fit <- function(fit, analysis, effect_name, transformed = FALSE) {
  est <- safe_model_value(fit$b)
  lower <- safe_model_value(fit$ci.lb)
  upper <- safe_model_value(fit$ci.ub)
  if (transformed) {
    est_out <- exp(est)
    lower_out <- exp(lower)
    upper_out <- exp(upper)
  } else {
    est_out <- est
    lower_out <- lower
    upper_out <- upper
  }
  data.frame(
    analysis = analysis,
    k = safe_model_value(fit$k),
    effect = effect_name,
    estimate = est_out,
    lower_95 = lower_out,
    upper_95 = upper_out,
    tau2 = safe_model_value(fit$tau2),
    I2_percent = safe_model_value(fit$I2),
    Q = safe_model_value(fit$QE),
    Q_df = safe_model_value(fit$k - fit$p),
    Q_p_value = safe_model_value(fit$QEp),
    H2 = safe_model_value(fit$H2),
    method = fit$method_used,
    stringsAsFactors = FALSE
  )
}

write_common_outputs <- function(out_dir, prefix, es, fit, summary_table) {
  write.csv(es, file.path(out_dir, paste0(prefix, "_effect_sizes.csv")), row.names = FALSE)
  write.csv(summary_table, file.path(out_dir, paste0(prefix, "_summary.csv")), row.names = FALSE, na = "")
}

process_scope_summary <- function(dat, scope, effect_name, fit = NULL) {
  if (nrow(dat) == 0) return(NULL)
  if (is.null(fit) && nrow(dat) >= 2) fit <- fit_random(dat)
  if (is.null(fit)) {
    return(data.frame(
      analysis_scope = scope,
      k = nrow(dat),
      effect = effect_name,
      estimate = dat$yi[[1]],
      lower_95 = dat$yi[[1]] - 1.96 * sqrt(dat$vi[[1]]),
      upper_95 = dat$yi[[1]] + 1.96 * sqrt(dat$vi[[1]]),
      tau2 = NA_real_,
      I2_percent = NA_real_,
      Q = NA_real_,
      Q_df = NA_real_,
      Q_p_value = NA_real_,
      H2 = NA_real_,
      method = "single study",
      stringsAsFactors = FALSE
    ))
  }
  out <- summary_from_fit(fit, scope, effect_name, transformed = FALSE)
  names(out)[names(out) == "analysis"] <- "analysis_scope"
  out
}

write_process_plot_summary <- function(out_dir, prefix, es, overall_fit, effect_name) {
  plot_summary <- bind_rows(
    process_scope_summary(es, "Overall", effect_name, fit = overall_fit),
    process_scope_summary(filter(es, design == "RCT"), "RCT", effect_name),
    process_scope_summary(filter(es, design == "NRS"), "NRS", effect_name)
  )
  write.csv(plot_summary, file.path(out_dir, paste0(prefix, "_plot_summary.csv")), row.names = FALSE, na = "")
  invisible(plot_summary)
}

build_binary_process_main <- function(
  index,
  csv_file,
  folder_slug,
  title,
  out_dir = NULL
) {
  base_dir <- getOption("sepsis_code_dir", getwd())
  if (is.null(out_dir)) {
    out_dir <- file.path(base_dir, "computed_data", "process_main")
  }
  if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

  raw <- read.csv(csv_file, check.names = FALSE, encoding = "UTF-8")
  names(raw) <- trimws(sub("^\ufeff", "", names(raw)))
  all_rows <- raw %>%
    transmute(
      study = clean_value(Study),
      design = toupper(clean_value(Design)),
      events_e = safe_num(`Intervention events`),
      n_e = safe_num(`Intervention total`),
      events_c = safe_num(`Control events`),
      n_c = safe_num(`Control total`),
      adjusted_rd = safe_num(`Adjusted RD, percentage points`),
      adjusted_lower_95 = safe_num(`Adjusted lower 95, percentage points`),
      adjusted_upper_95 = safe_num(`Adjusted upper 95, percentage points`),
      bi = n_e - events_e,
      di = n_c - events_c
    ) %>%
    add_binary_process_rd_effect()

  exclusions <- all_rows %>%
    mutate(
      reason = case_when(
        study == "Kim 2021" &
          is.finite(events_e) & is.finite(n_e) & is.finite(events_c) & is.finite(n_c) &
          events_e == n_e & events_c == n_c ~ kim_2021_fluid_exclusion_sentence,
        is.na(study) ~ "Not meta-analysed because the study identifier was missing.",
        !design %in% c("RCT", "NRS") ~ "Not meta-analysed because study design was not classified as RCT or NRS.",
        is.na(events_e) | is.na(n_e) | is.na(events_c) | is.na(n_c) ~
          "Not meta-analysed because event counts or denominators were missing.",
        n_e <= 0 | n_c <= 0 | events_e < 0 | events_e > n_e | events_c < 0 | events_c > n_c ~
          "Not meta-analysed because event counts or denominators were invalid.",
        !is.finite(selected_rd) | !is.finite(selected_se) | selected_se <= 0 ~
          "Not meta-analysed because the effect estimate was non-informative.",
        TRUE ~ NA_character_
      )
    ) %>%
    filter(!is.na(reason)) %>%
    transmute(study, design, events_e, n_e, events_c, n_c, reason)

  es <- all_rows %>%
    filter(
      !is.na(study),
      design %in% c("RCT", "NRS"),
      !is.na(events_e), !is.na(n_e), !is.na(events_c), !is.na(n_c),
      n_e > 0, n_c > 0,
      events_e >= 0, events_e <= n_e,
      events_c >= 0, events_c <= n_c
    ) %>%
    filter(
      !is.na(selected_rd),
      !is.na(selected_se),
      selected_se > 0
    ) %>%
    transmute(
      study = study,
      design = design,
      events_e = events_e,
      n_e = n_e,
      events_c = events_c,
      n_c = n_c,
      adjusted_rd = adjusted_rd,
      adjusted_lower_95 = adjusted_lower_95,
      adjusted_upper_95 = adjusted_upper_95,
      crude_rd = crude_rd,
      crude_se = crude_se,
      adjusted_se = adjusted_se,
      use_adjusted_rd = use_adjusted_rd,
      selected_rd = selected_rd,
      selected_se = selected_se,
      selected_lower_95 = selected_lower_95,
      selected_upper_95 = selected_upper_95,
      effect_basis = effect_basis,
      yi = selected_rd,
      vi = selected_se^2,
      se = selected_se,
      input_type = ifelse(use_adjusted_rd, "adjusted RD+SE", "crude count-based RD")
    )
  validate_process_2x2 <- function(x) {
    if (exists("validate_mortality_2x2", mode = "function", inherits = TRUE)) {
      validate_mortality_2x2(x, context = basename(csv_file))
    }
    invisible(TRUE)
  }
  validate_process_2x2(es)
  fit <- fit_random(es)
  summary_table <- summary_from_fit(fit, title, "RD", transformed = FALSE)

  prefix <- paste0("3.", index, "_", folder_slug)
  write_common_outputs(out_dir, prefix, es, fit, summary_table)
  write.csv(
    exclusions,
    file.path(out_dir, paste0(prefix, "_exclusions.csv")),
    row.names = FALSE,
    na = ""
  )
  plot_summary <- write_process_plot_summary(out_dir, prefix, es, fit, "RD")
  invisible(list(summary = summary_table, plot_summary = plot_summary, effect_sizes = es))
}

build_antibiotic_time_main <- function(out_dir = NULL) {
  base_dir <- getOption("sepsis_code_dir", getwd())
  if (is.null(out_dir)) {
    out_dir <- file.path(base_dir, "computed_data", "process_main")
  }
  if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

  raw <- read.csv(file.path(base_dir, "source_csv", "process_antibiotic_time.csv"), check.names = FALSE, encoding = "UTF-8")
  names(raw) <- trimws(sub("^\ufeff", "", names(raw)))
  dat <- raw %>%
    transmute(
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
      true_antibiotic_administration = clean_value(`True antibiotic administration`),
      antibiotic_time_origin = clean_value(`Antibiotic time origin`)
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
      yi_direct = ifelse(use_adjusted_md, adjusted_md, crude_md),
      se_direct = ifelse(use_adjusted_md, adjusted_se, crude_se)
    )

  es <- dat %>%
    filter(
      !is.na(study),
      design %in% c("RCT", "NRS"),
      !is.na(yi_direct), !is.na(se_direct), se_direct > 0
    ) %>%
    transmute(
      study = study,
      design = design,
      mean_e = mean_e,
      sd_e = sd_e,
      n_e = n_e,
      mean_c = mean_c,
      sd_c = sd_c,
      n_c = n_c,
      adjusted_md = adjusted_md,
      adjusted_lower_95 = adjusted_lower_95,
      adjusted_upper_95 = adjusted_upper_95,
      crude_md = crude_md,
      crude_se = crude_se,
      adjusted_se = adjusted_se,
      use_adjusted_md = use_adjusted_md,
      selected_md = yi_direct,
      selected_se = se_direct,
      true_antibiotic_administration = true_antibiotic_administration,
      antibiotic_time_origin = antibiotic_time_origin,
      yi = yi_direct,
      vi = se_direct^2,
      input_type = "selected direct MD+SE"
    )

  if (exists("validate_mortality_2x2", mode = "function", inherits = TRUE)) {
    validate_mortality_2x2(es, context = "process_antibiotic_time.csv")
  }

  fit <- fit_random(es)
  summary_table <- summary_from_fit(fit, process_titles[6], "MD", transformed = FALSE)
  prefix <- "3.6_time_to_antibiotic_administration_meta_analysis"
  write_common_outputs(out_dir, prefix, es, fit, summary_table)
  plot_summary <- write_process_plot_summary(out_dir, prefix, es, fit, "MD")
  invisible(list(summary = summary_table, plot_summary = plot_summary, effect_sizes = es))
}
