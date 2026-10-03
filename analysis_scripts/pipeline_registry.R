# Canonical producer registry for every computed_data domain.

analysis_pipeline_registry <- function(project_root = getwd()) {
  data.frame(
    domain = c(
      "los_hospital_main",
      "los_icu_main",
      "process_main",
      "mortality_main",
      "los_sensitivity",
      "process_sensitivity",
      "process_subgroup",
      "mortality_rd_summary",
      "mortality_sensitivity",
      "mortality_subgroup",
      "mortality_baseline_distribution",
      "mortality_baseline_process",
      "mortality_baseline_rd",
      "mortality_process_metareg",
      "mortality_process_rd_metareg",
      "mortality_publication_cumulative",
      "evidence_map",
      "narrative_analysis"
    ),
    script = c(
      "analysis_scripts/los_hospital_main/build_hospital_los_main.R",
      "analysis_scripts/los_icu_main/build_icu_los_main.R",
      "analysis_scripts/process_main/build_process_main.R",
      "analysis_scripts/mortality_main/build_mortality_main.R",
      "analysis_scripts/los_sensitivity/build_los_sensitivity.R",
      "analysis_scripts/process_sensitivity/build_process_sensitivity.R",
      "analysis_scripts/process_subgroup/build_process_subgroup.R",
      "analysis_scripts/mortality_rd_summary/fit_bayesian_rd_summary.R",
      "analysis_scripts/mortality_sensitivity/fit_sensitivity_models.R",
      "analysis_scripts/mortality_subgroup/fit_subgroup_analysis.R",
      "analysis_scripts/mortality_baseline_distribution/build_baseline_mortality_distribution.R",
      "analysis_scripts/mortality_baseline_process/build_baseline_process_analysis.R",
      "analysis_scripts/mortality_baseline_rd/fit_baseline_mortality_rd_metareg.R",
      "analysis_scripts/mortality_process_metareg/fit_process_metareg.R",
      "analysis_scripts/mortality_process_rd_metareg/fit_process_rd_metareg.R",
      "analysis_scripts/mortality_publication_cumulative/build_mortality_publication_cumulative.R",
      "analysis_scripts/evidence_map/build_evidence_map.R",
      "analysis_scripts/narrative_analysis/build_narrative_analysis.R"
    ),
    stage = c(
      10L, 10L, 10L, 10L,
      20L, 20L, 20L,
      30L, 30L, 30L, 30L, 30L, 30L, 30L, 30L, 30L,
      40L, 40L
    ),
    heavy = c(
      FALSE, FALSE, FALSE, TRUE,
      FALSE, FALSE, FALSE,
      TRUE, TRUE, TRUE, FALSE, FALSE, TRUE, TRUE, TRUE, TRUE,
      FALSE, FALSE
    ),
    stringsAsFactors = FALSE
  )
}
