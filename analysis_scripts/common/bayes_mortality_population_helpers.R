# Shared display-layer population definitions for mortality figures.
#
# These helpers mirror the population taxonomy used by the mortality analyses.
# They do not read data, fit models, or write outputs. Both the current
# *-targeted encoding and the equivalent parenthetical legacy encoding are
# accepted so callers cannot silently invent a competing regex.

population_scope_text <- function(x) {
  tolower(trimws(ifelse(is.na(x), "", as.character(x))))
}

is_all_screened <- function(x) {
  grepl("(^|[;|, ])all screened($|[;|, ])", population_scope_text(x))
}

is_infection_targeted <- function(x) {
  scope <- population_scope_text(x)
  grepl("infection-targeted|\\(infection\\)", scope)
}

is_severe_sepsis_targeted <- function(x) {
  scope <- population_scope_text(x)
  grepl("\\(severe sepsis\\)|severe sepsis-targeted", scope)
}

is_sepsis3_targeted <- function(x) {
  scope <- population_scope_text(x)
  grepl("\\(sepsis-3\\)|sepsis-3-targeted", scope)
}

is_organ_dysfunction_targeted <- function(x) {
  is_severe_sepsis_targeted(x) | is_sepsis3_targeted(x)
}
