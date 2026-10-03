# Shared mortality risk-difference construction used by all RD analyses.

read_mortality_rd_inputs <- function(path) {
  if (!file.exists(path)) stop("Missing mortality RD input file: ", path)
  inputs <- utils::read.csv(path, stringsAsFactors = FALSE, check.names = FALSE, encoding = "UTF-8")
  names(inputs) <- trimws(sub("^\ufeff", "", names(inputs)))
  required <- c(
    "Study", "rd_method", "rd_percentage_points",
    "lower_95_percentage_points", "upper_95_percentage_points", "source_note"
  )
  missing <- setdiff(required, names(inputs))
  if (length(missing) > 0L) stop("Mortality RD input file lacks: ", paste(missing, collapse = ", "))
  if (anyDuplicated(inputs$Study)) stop("Mortality RD input file contains duplicate study rows.")
  inputs
}

add_mortality_rd_effects <- function(data, rd_inputs_file) {
  dat <- as.data.frame(data)
  if (!"study" %in% names(dat) && "Study" %in% names(dat)) dat$study <- dat$Study
  if (!"design" %in% names(dat) && "Design" %in% names(dat)) dat$design <- dat$Design
  if (!"aRR" %in% names(dat) && "arr" %in% names(dat)) dat$aRR <- dat$arr
  if (!"aRR_lower_95" %in% names(dat) && "arr_lower_95" %in% names(dat)) dat$aRR_lower_95 <- dat$arr_lower_95
  if (!"aRR_upper_95" %in% names(dat) && "arr_upper_95" %in% names(dat)) dat$aRR_upper_95 <- dat$arr_upper_95
  for (alias in c("Study", "Design", "arr", "arr_lower_95", "arr_upper_95")) {
    canonical <- switch(
      alias,
      Study = "study",
      Design = "design",
      arr = "aRR",
      arr_lower_95 = "aRR_lower_95",
      arr_upper_95 = "aRR_upper_95"
    )
    if (alias %in% names(dat) && canonical %in% names(dat) && !identical(alias, canonical)) dat[[alias]] <- NULL
  }

  required <- c("study", "design", "events_e", "n_e", "events_c", "n_c")
  missing <- setdiff(required, names(dat))
  if (length(missing) > 0L) stop("Mortality RD data lack: ", paste(missing, collapse = ", "))
  for (name in c("events_e", "n_e", "events_c", "n_c", "aRR", "aRR_lower_95", "aRR_upper_95")) {
    if (!name %in% names(dat)) dat[[name]] <- NA_real_
    dat[[name]] <- suppressWarnings(as.numeric(trimws(as.character(dat[[name]]))))
  }
  dat$study <- trimws(as.character(dat$study))
  dat$design <- toupper(trimws(as.character(dat$design)))

  invalid <- !is.finite(dat$n_e) | !is.finite(dat$n_c) | dat$n_e <= 0 | dat$n_c <= 0 |
    !is.finite(dat$events_e) | !is.finite(dat$events_c) |
    dat$events_e < 0 | dat$events_e > dat$n_e | dat$events_c < 0 | dat$events_c > dat$n_c
  if (any(invalid)) {
    stop("Mortality RD construction found invalid event counts or denominators: ", paste(dat$study[invalid], collapse = "; "))
  }

  event_rd <- metafor::escalc(
    measure = "RD",
    ai = dat$events_e,
    bi = dat$n_e - dat$events_e,
    ci = dat$events_c,
    di = dat$n_c - dat$events_c,
    add = 0,
    to = "none"
  )
  dat$p_e <- dat$events_e / dat$n_e
  dat$p_c <- dat$events_c / dat$n_c
  dat$yi_rd_event_counts <- as.numeric(event_rd$yi)
  dat$vi_rd_event_counts <- as.numeric(event_rd$vi)
  dat$yi_rd <- dat$yi_rd_event_counts
  dat$vi_rd <- dat$vi_rd_event_counts
  dat$rd_source <- "event_counts"
  dat$rd_note <- "Risk difference from study-level event counts."

  inputs <- read_mortality_rd_inputs(rd_inputs_file)
  unknown <- setdiff(inputs$Study, dat$study)
  if (length(unknown) > 0L) stop("Mortality RD inputs refer to unknown studies: ", paste(unknown, collapse = "; "))

  for (i in seq_len(nrow(inputs))) {
    row <- inputs[i, , drop = FALSE]
    index <- which(dat$study == row$Study)
    if (length(index) != 1L) stop("Expected one mortality row for RD input: ", row$Study)

    if (identical(row$rd_method, "reported_adjusted_rd")) {
      values <- suppressWarnings(as.numeric(unlist(row[1, c(
        "rd_percentage_points", "lower_95_percentage_points", "upper_95_percentage_points"
      ), drop = FALSE], use.names = FALSE))) / 100
      if (any(!is.finite(values)) || values[[3]] <= values[[2]]) {
        stop("Invalid reported adjusted RD input for: ", row$Study)
      }
      dat$yi_rd[index] <- values[[1]]
      dat$vi_rd[index] <- ((values[[3]] - values[[2]]) / (2 * 1.96))^2
      dat$rd_source[index] <- "adjusted_rd_reported"
      dat$rd_note[index] <- row$source_note
    } else if (identical(row$rd_method, "derive_from_adjusted_rr_and_control_risk")) {
      rr_values <- unlist(dat[index, c("aRR", "aRR_lower_95", "aRR_upper_95"), drop = FALSE], use.names = FALSE)
      if (any(!is.finite(rr_values)) || any(rr_values <= 0) || rr_values[[3]] <= rr_values[[2]]) {
        stop("Cannot derive adjusted mortality RD because adjusted RR fields are invalid for: ", row$Study)
      }
      log_rr_variance <- ((log(rr_values[[3]]) - log(rr_values[[2]])) / (2 * 1.96))^2
      control_risk_variance <- dat$p_c[index] * (1 - dat$p_c[index]) / dat$n_c[index]
      dat$yi_rd[index] <- dat$p_c[index] * (rr_values[[1]] - 1)
      dat$vi_rd[index] <-
        (dat$p_c[index] * rr_values[[1]])^2 * log_rr_variance +
        (rr_values[[1]] - 1)^2 * control_risk_variance
      dat$rd_source[index] <- "adjusted_rd_from_arr"
      dat$rd_note[index] <- row$source_note
    } else {
      stop("Unknown mortality RD construction method for ", row$Study, ": ", row$rd_method)
    }
  }

  dat$se_rd <- sqrt(dat$vi_rd)
  dat$rd_lower_95 <- dat$yi_rd - 1.96 * dat$se_rd
  dat$rd_upper_95 <- dat$yi_rd + 1.96 * dat$se_rd
  dat$RD <- dat$yi_rd
  dat$yi <- dat$yi_rd

  invalid_variance <- !is.finite(dat$vi_rd) | dat$vi_rd <= 0
  if (any(invalid_variance)) {
    stop(
      "Mortality RD construction produced non-positive or non-finite variance for: ",
      paste(dat$study[invalid_variance], collapse = "; "),
      ". No continuity correction was applied."
    )
  }
  dat
}
