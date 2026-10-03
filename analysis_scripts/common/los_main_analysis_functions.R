suppressPackageStartupMessages({
  library(dplyr)
  library(metafor)
})

los_clean_value <- function(x) {
  x <- trimws(as.character(x))
  x[x %in% c("", "/", "NA", "Na", "na", "N/A")] <- NA_character_
  x
}

los_safe_num <- function(x) {
  suppressWarnings(as.numeric(trimws(as.character(x))))
}

los_model_value <- function(x) {
  ifelse(is.null(x) || length(x) == 0, NA_real_, as.numeric(x)[1])
}

fit_los_reml <- function(dat, outcome_label) {
  fit <- try(
    metafor::rma(
      yi = yi,
      vi = vi_original,
      data = dat,
      method = "REML",
      control = list(maxiter = 1000)
    ),
    silent = TRUE
  )
  if (inherits(fit, "try-error")) {
    fit <- try(
      metafor::rma(
        yi = yi,
        vi = vi_original,
        data = dat,
        method = "REML",
        control = list(maxiter = 10000)
      ),
      silent = TRUE
    )
  }
  if (inherits(fit, "try-error")) {
    stop(
      outcome_label,
      " REML model failed to converge; no alternative heterogeneity estimator was substituted. ",
      "k = ", nrow(dat), "; retained studies = ", paste(dat$study, collapse = "; ")
    )
  }
  fit$method_used <- "REML"
  fit
}

summarise_los_reml <- function(fit, group, dat) {
  estimate <- los_model_value(fit$b)
  lower <- los_model_value(fit$ci.lb)
  upper <- los_model_value(fit$ci.ub)
  data.frame(
    group = group,
    k = fit$k,
    k_rct = sum(dat$design == "RCT"),
    k_nrs = sum(dat$design == "NRS"),
    MD = estimate,
    Lower_95 = lower,
    Upper_95 = upper,
    MD_label = sprintf("%.2f [%.2f, %.2f]", estimate, lower, upper),
    tau2 = los_model_value(fit$tau2),
    I2_percent = los_model_value(fit$I2),
    Q = los_model_value(fit$QE),
    Q_df = los_model_value(fit$k - fit$p),
    p_heterogeneity = los_model_value(fit$QEp),
    method = fit$method_used,
    stringsAsFactors = FALSE
  )
}

format_los_summary_3dp <- function(dat) {
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

build_los_main_analysis <- function(input_file,
                                    output_dir,
                                    outcome_label,
                                    effect_file,
                                    summary_file,
                                    excluded_file) {
  if (!file.exists(input_file)) stop("Missing LOS source CSV: ", input_file)
  if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

  raw <- read.csv(
    input_file,
    stringsAsFactors = FALSE,
    check.names = FALSE,
    encoding = "UTF-8"
  )
  names(raw) <- trimws(sub("^\ufeff", "", names(raw)))
  raw <- raw[, names(raw) != "", drop = FALSE]

  required_cols <- c(
    "Study", "Control mean", "Control SD", "Control n",
    "Intervention mean", "Intervention SD", "Intervention n", "Design",
    "Adjusted effect estimate, days", "Adjusted lower 95, days", "Adjusted upper 95, days"
  )
  missing_cols <- setdiff(required_cols, names(raw))
  if (length(missing_cols) > 0) {
    stop("Missing required columns: ", paste(missing_cols, collapse = ", "))
  }

  all_rows <- raw %>%
    transmute(
      row_order = row_number(),
      study = los_clean_value(Study),
      design_raw = los_clean_value(Design),
      design = toupper(design_raw),
      mean_c = los_safe_num(.data[["Control mean"]]),
      sd_c = los_safe_num(.data[["Control SD"]]),
      n_c = los_safe_num(.data[["Control n"]]),
      mean_e = los_safe_num(.data[["Intervention mean"]]),
      sd_e = los_safe_num(.data[["Intervention SD"]]),
      n_e = los_safe_num(.data[["Intervention n"]]),
      adjusted_yi = los_safe_num(.data[["Adjusted effect estimate, days"]]),
      adjusted_lower = los_safe_num(.data[["Adjusted lower 95, days"]]),
      adjusted_upper = los_safe_num(.data[["Adjusted upper 95, days"]]),
      adjusted_se = (adjusted_upper - adjusted_lower) / (2 * 1.96),
      eligible_adjusted = !is.na(adjusted_yi) & !is.na(adjusted_lower) &
        !is.na(adjusted_upper) & !is.na(adjusted_se) & adjusted_se > 0 &
        adjusted_upper > adjusted_lower,
      eligible_group = !is.na(mean_e) & !is.na(sd_e) & !is.na(n_e) &
        !is.na(mean_c) & !is.na(sd_c) & !is.na(n_c) &
        sd_e > 0 & sd_c > 0 & n_e > 1 & n_c > 1
    )

  excluded <- all_rows %>%
    filter(is.na(study) | !design %in% c("RCT", "NRS") | !(eligible_adjusted | eligible_group))

  group_rows <- all_rows %>%
    filter(
      !is.na(study),
      design %in% c("RCT", "NRS"),
      eligible_group,
      !eligible_adjusted
    )

  group_effects <- metafor::escalc(
    measure = "MD",
    m1i = mean_e,
    sd1i = sd_e,
    n1i = n_e,
    m2i = mean_c,
    sd2i = sd_c,
    n2i = n_c,
    data = group_rows
  ) %>%
    mutate(
      vi_original = vi,
      se_original = sqrt(vi_original),
      MD = yi,
      Lower_95 = yi - 1.96 * se_original,
      Upper_95 = yi + 1.96 * se_original,
      md_label = sprintf("%.2f [%.2f, %.2f]", MD, Lower_95, Upper_95),
      input_type = "group mean+SD"
    )

  adjusted_effects <- all_rows %>%
    filter(!is.na(study), design %in% c("RCT", "NRS"), eligible_adjusted) %>%
    transmute(
      row_order,
      study,
      design_raw,
      design,
      mean_c,
      sd_c,
      n_c,
      mean_e,
      sd_e,
      n_e,
      yi = adjusted_yi,
      vi = adjusted_se^2,
      vi_original = vi,
      se_original = adjusted_se,
      MD = adjusted_yi,
      Lower_95 = adjusted_lower,
      Upper_95 = adjusted_upper,
      md_label = sprintf("%.2f [%.2f, %.2f]", MD, Lower_95, Upper_95),
      input_type = "adjusted direct effect"
    )

  effects <- bind_rows(as.data.frame(group_effects), adjusted_effects)
  fit_all <- fit_los_reml(effects, outcome_label)
  fit_rct <- fit_los_reml(effects %>% filter(design == "RCT"), outcome_label)
  fit_nrs <- fit_los_reml(effects %>% filter(design == "NRS"), outcome_label)

  tau2_all <- los_model_value(fit_all$tau2)
  effects <- effects %>%
    mutate(
      weight_random = 100 * (1 / (vi_original + tau2_all)) /
        sum(1 / (vi_original + tau2_all)),
      weight_label = sprintf("%.1f", weight_random)
    )

  summary_table <- bind_rows(
    summarise_los_reml(fit_rct, "RCT", effects %>% filter(design == "RCT")),
    summarise_los_reml(fit_nrs, "NRS", effects %>% filter(design == "NRS")),
    summarise_los_reml(fit_all, "Overall", effects)
  )

  effect_output <- effects %>%
    select(
      study, design, mean_e, sd_e, n_e, mean_c, sd_c, n_c,
      yi, vi_original, se_original, MD, Lower_95, Upper_95,
      row_order, input_type, weight_random, weight_label, md_label
    )

  write.csv(effect_output, file.path(output_dir, effect_file), row.names = FALSE)
  write.csv(summary_table, file.path(output_dir, summary_file), row.names = FALSE, na = "")
  write.csv(excluded, file.path(output_dir, excluded_file), row.names = FALSE)

  invisible(list(effects = effect_output, summary = summary_table, excluded = excluded))
}
