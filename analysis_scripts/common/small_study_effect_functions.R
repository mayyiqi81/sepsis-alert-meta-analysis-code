# Shared small-study-effect diagnostics used by the evidence-map analysis.

small_study_model_value <- function(x) {
  if (is.null(x) || length(x) == 0L) return(NA_real_)
  as.numeric(x)[1]
}

publication_bias_tests <- function(es, vi_col = "vi") {
  required <- c("yi", vi_col)
  missing <- setdiff(required, names(es))
  if (length(missing) > 0L) {
    stop("Small-study-effect data lack: ", paste(missing, collapse = ", "))
  }

  test_dat <- data.frame(
    yi = as.numeric(es$yi),
    vi = as.numeric(es[[vi_col]])
  )
  test_dat <- test_dat[
    is.finite(test_dat$yi) & is.finite(test_dat$vi) & test_dat$vi > 0,
    ,
    drop = FALSE
  ]
  if (nrow(test_dat) < 3L) {
    return(data.frame(
      egger_intercept = NA_real_,
      egger_p = NA_real_,
      begg_tau = NA_real_,
      begg_p = NA_real_
    ))
  }

  test_dat$se <- sqrt(test_dat$vi)
  test_dat$z_value <- test_dat$yi / test_dat$se
  test_dat$precision <- 1 / test_dat$se

  egger_fit <- try(stats::lm(z_value ~ precision, data = test_dat), silent = TRUE)
  if (inherits(egger_fit, "try-error")) {
    egger_intercept <- NA_real_
    egger_p <- NA_real_
  } else {
    egger_summary <- summary(egger_fit)
    egger_intercept <- small_study_model_value(stats::coef(egger_fit)[1])
    egger_p <- suppressWarnings(as.numeric(stats::coef(egger_summary)[1, "Pr(>|t|)"]))
  }

  # Begg-Mazumdar uses the rank correlation implementation and variance
  # correction provided by metafor.
  begg_fit <- try(metafor::ranktest(test_dat$yi, vi = test_dat$vi), silent = TRUE)
  if (inherits(begg_fit, "try-error")) {
    begg_tau <- NA_real_
    begg_p <- NA_real_
  } else {
    begg_tau <- small_study_model_value(begg_fit$tau)
    begg_p <- small_study_model_value(begg_fit$pval)
  }

  data.frame(
    egger_intercept = egger_intercept,
    egger_p = egger_p,
    begg_tau = begg_tau,
    begg_p = begg_p
  )
}
