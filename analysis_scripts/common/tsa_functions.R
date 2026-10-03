# Shared calculations for the mortality trial sequential analysis.

fit_tsa_reml <- function(data, yi_col = "yi", vi_col = "vi", context = "TSA") {
  yi <- as.numeric(data[[yi_col]])
  vi <- as.numeric(data[[vi_col]])
  if (length(yi) < 2L) {
    return(list(
      b = yi[[1]],
      se = sqrt(vi[[1]]),
      tau2 = NA_real_,
      method_used = "single study"
    ))
  }

  fit <- try(
    metafor::rma(yi = yi, vi = vi, method = "REML", control = list(maxiter = 1000)),
    silent = TRUE
  )
  if (inherits(fit, "try-error")) {
    fit <- try(
      metafor::rma(yi = yi, vi = vi, method = "REML", control = list(maxiter = 10000)),
      silent = TRUE
    )
  }
  if (inherits(fit, "try-error")) {
    stop(context, " REML model failed to converge; no alternative heterogeneity estimator was substituted.")
  }
  fit$method_used <- "REML"
  fit
}

cumulative_reml_z <- function(data, yi_col = "yi", vi_col = "vi", n_col = "n") {
  required <- c(yi_col, vi_col, n_col)
  missing <- setdiff(required, names(data))
  if (length(missing) > 0L) stop("Cumulative TSA data lack: ", paste(missing, collapse = ", "))

  out <- data
  out$Cumulative_N <- cumsum(as.numeric(out[[n_col]]))
  out$Cumulative_logRR_REML <- NA_real_
  out$Cumulative_SE_REML <- NA_real_
  out$Cumulative_Z <- NA_real_
  out$Cumulative_tau2_REML <- NA_real_
  out$Cumulative_method <- NA_character_

  for (i in seq_len(nrow(out))) {
    fit <- fit_tsa_reml(
      out[seq_len(i), , drop = FALSE],
      yi_col = yi_col,
      vi_col = vi_col,
      context = paste0("Cumulative TSA at study ", i)
    )
    estimate <- as.numeric(fit$b)[1]
    se <- as.numeric(fit$se)[1]
    out$Cumulative_logRR_REML[[i]] <- estimate
    out$Cumulative_SE_REML[[i]] <- se
    out$Cumulative_Z[[i]] <- estimate / se
    out$Cumulative_tau2_REML[[i]] <- as.numeric(fit$tau2)[1]
    out$Cumulative_method[[i]] <- fit$method_used
  }
  out
}

tsa_diversity_d2 <- function(vi, tau2) {
  vi <- as.numeric(vi)
  if (any(!is.finite(vi)) || any(vi <= 0)) stop("TSA sampling variances must be positive and finite.")
  if (!is.finite(tau2) || tau2 < 0) stop("TSA REML tau-squared must be finite and non-negative.")
  variance_fixed <- 1 / sum(1 / vi)
  variance_random <- 1 / sum(1 / (vi + tau2))
  pmin(pmax(1 - variance_fixed / variance_random, 0), 1 - 1e-12)
}

calculate_tsa_information_size <- function(cer, rrr, alpha, power, d_squared) {
  if (!is.finite(cer) || cer <= 0 || cer >= 1) stop("TSA control event rate must be between zero and one.")
  if (!is.finite(rrr) || rrr <= 0 || rrr >= 1) stop("TSA relative risk reduction must be between zero and one.")
  if (!is.finite(d_squared) || d_squared < 0 || d_squared >= 1) stop("TSA D-squared must be in [0, 1).")

  p_control <- cer
  p_intervention <- cer * (1 - rrr)
  p_bar <- (p_control + p_intervention) / 2
  z_alpha <- stats::qnorm(1 - alpha / 2)
  z_beta <- stats::qnorm(power)
  numerator <- (
    z_alpha * sqrt(2 * p_bar * (1 - p_bar)) +
      z_beta * sqrt(p_control * (1 - p_control) + p_intervention * (1 - p_intervention))
  )^2
  risk_difference <- p_control - p_intervention
  ris <- 2 * numerator / risk_difference^2
  data.frame(
    CER = cer,
    RRR = rrr,
    RIS = ris,
    D_squared = d_squared,
    DARIS = ris / (1 - d_squared),
    stringsAsFactors = FALSE
  )
}

tsa_canonical_correlation <- function(information_fraction) {
  timing <- sort(unique(as.numeric(information_fraction)))
  if (length(timing) == 0L || any(!is.finite(timing)) || any(timing <= 0) || any(timing > 1)) {
    stop("TSA information fractions must be finite values in (0, 1].")
  }
  outer(timing, timing, function(x, y) sqrt(pmin(x, y) / pmax(x, y)))
}

calibrate_two_sided_obrien_fleming <- function(information_fraction, alpha = 0.05) {
  timing_input <- pmin(pmax(as.numeric(information_fraction), 1e-12), 1)
  timing <- sort(unique(timing_input))
  corr <- tsa_canonical_correlation(timing)
  coverage <- function(constant) {
    boundary <- constant / sqrt(timing)
    algorithm <- if (length(timing) <= 5L) {
      mvtnorm::Miwa()
    } else {
      mvtnorm::GenzBretz(maxpts = 1000000L, abseps = 1e-7, releps = 0)
    }
    had_seed <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
    if (had_seed) old_seed <- get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
    on.exit({
      if (had_seed) {
        assign(".Random.seed", old_seed, envir = .GlobalEnv)
      } else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
        rm(".Random.seed", envir = .GlobalEnv)
      }
    }, add = TRUE)
    set.seed(20260905)
    as.numeric(mvtnorm::pmvnorm(
      lower = -boundary,
      upper = boundary,
      corr = corr,
      algorithm = algorithm
    ))
  }
  objective <- function(constant) coverage(constant) - (1 - alpha)
  lower <- 0.5
  upper <- 8
  while (objective(lower) > 0 && lower > 1e-6) lower <- lower / 2
  while (objective(upper) < 0 && upper < 64) upper <- upper * 2
  if (objective(lower) > 0 || objective(upper) < 0) {
    stop("Could not bracket the calibrated O'Brien-Fleming boundary constant.")
  }
  constant <- stats::uniroot(
    objective,
    interval = c(lower, upper),
    tol = if (length(timing) <= 5L) 1e-10 else 1e-7
  )$root
  structure(
    constant / sqrt(timing_input),
    boundary_constant = constant,
    alpha = alpha,
    information_fraction = timing_input
  )
}
