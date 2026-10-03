# Canonical convention for all antibiotic-time differences: intervention minus
# control. Faster treatment is negative; this does not recode event-time HRs.
assert_antibiotic_time_source_direction <- function(data) {
  col <- 'Antibiotic time MD, hours'
  stopifnot(all(c('Study',col) %in% names(data)))
  # Independent anchors spanning shortened, prolonged and source-reported
  # reduction contrasts. These catch accidental restoration of the old axis.
  expected <- c('Hayden 2016'=-50.6/60,'Ostermayer 2026'=0.8,
                'Kurczewski 2015'=-1.3,'Adams 2022'=-1.85)
  idx <- match(names(expected),data$Study)
  if(anyNA(idx) || any(abs(as.numeric(data[[col]][idx])-expected)>1e-7))
    stop('Antibiotic-time source direction mismatch: require intervention minus control, negative = faster treatment. Audit source_csv/mortality.csv before fitting.')
  invisible(TRUE)
}

assert_antibiotic_time_plot_direction <- function(data, project_root) {
  src <- read.csv(file.path(project_root,'source_csv','mortality.csv'),check.names=FALSE)
  assert_antibiotic_time_source_direction(src)
  rows <- data[data$covariate=='antibiotic_time',,drop=FALSE]
  idx <- match(rows$study,src$Study)
  expected <- as.numeric(src[['Antibiotic time MD, hours']][idx])
  if(!nrow(rows) || anyNA(idx) || any(!is.finite(expected)) ||
     any(abs(rows$covariate_value-expected)>1e-7))
    stop('Stale antibiotic-time plotting cache: moderator values must match the intervention-minus-control source.')
  invisible(TRUE)
}
