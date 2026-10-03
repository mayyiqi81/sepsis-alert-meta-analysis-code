# Run registered computed-data producers in dependency order.
#
# Usage:
#   Rscript analysis_scripts/run_all_analysis.R
#   Rscript analysis_scripts/run_all_analysis.R --skip-heavy
#   Rscript analysis_scripts/run_all_analysis.R --domain=process_main,evidence_map
#   Rscript analysis_scripts/run_all_analysis.R --dry-run

args <- commandArgs(trailingOnly = TRUE)
script_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
script_path <- if (length(script_arg) > 0) sub("^--file=", "", script_arg[1]) else "analysis_scripts/run_all_analysis.R"
project_root <- normalizePath(file.path(dirname(script_path), ".."), winslash = "/", mustWork = TRUE)

source(file.path(project_root, "analysis_scripts", "pipeline_registry.R"))
registry <- analysis_pipeline_registry(project_root)

domain_arg <- grep("^--domain=", args, value = TRUE)
skip_heavy <- "--skip-heavy" %in% args
dry_run <- "--dry-run" %in% args
list_only <- "--list" %in% args

if (length(domain_arg) > 1) stop("Use at most one --domain= argument.")
if (length(domain_arg) == 1) {
  requested <- trimws(strsplit(sub("^--domain=", "", domain_arg), ",", fixed = TRUE)[[1]])
  unknown <- setdiff(requested, registry$domain)
  if (length(unknown) > 0) stop("Unknown analysis domain(s): ", paste(unknown, collapse = ", "))
  registry <- registry[registry$domain %in% requested, , drop = FALSE]
}
if (skip_heavy) registry <- registry[!registry$heavy, , drop = FALSE]
registry <- registry[order(registry$stage, registry$domain), , drop = FALSE]

if (list_only || dry_run) {
  print(registry, row.names = FALSE)
  quit(save = "no", status = 0)
}

rscript <- file.path(R.home("bin"), if (.Platform$OS.type == "windows") "Rscript.exe" else "Rscript")
for (i in seq_len(nrow(registry))) {
  domain <- registry$domain[i]
  script <- file.path(project_root, registry$script[i])
  message("[", i, "/", nrow(registry), "] Running ", domain, ": ", registry$script[i])
  status <- system2(rscript, shQuote(script), stdout = "", stderr = "")
  if (!identical(as.integer(status), 0L)) {
    stop("Analysis domain failed with status ", status, ": ", domain)
  }
}

message("Completed ", nrow(registry), " registered analysis domain(s).")
