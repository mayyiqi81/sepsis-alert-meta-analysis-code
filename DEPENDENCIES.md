# Dependencies

## Software dependencies

Manuscript-specified environment: R 4.4.2; brms 2.23.0; metafor 4.8.0; mvtnorm 1.3.3.

Additional non-base packages referenced in the source: `cowplot`, `digest`, `dplyr`, `forestplot`, `ggplot2`, `ggtext`, `gridExtra`, `posterior`, `rlang`, `tidyr`. No additional package versions are asserted. R-supplied packages referenced include `grid`, `parallel`, `stats`, `tools`, `utils`. A working Stan backend and compilation environment appropriate to brms are needed for Bayesian fitting.

## Internal code dependencies

All executable code dependencies reside within `analysis_scripts/`. The project-path helper is located at `analysis_scripts/common/project_paths.R`; it recognises this package by its analysis tree and registry. The tables below were checked against parsed R calls.

| Module | Direct shared/helper or module code dependencies | Input/output and upstream relationship |
| --- | --- | --- |
| `common/` | `common/bayes_prior_scale_spec.R` | Shared functions take module-supplied inputs, configuration and fitted objects. |
| `mortality_main/` | `common/bayes_mortality_functions.R`; `common/bayes_prior_scale_spec.R`; `common/project_paths.R`; `mortality_main/build_mortality_forest_cache.R`; `mortality_main/fit_canonical_main.R` | The wrapper fits/validates canonical mortality models before constructing forest summaries. |
| `mortality_sensitivity/` | `common/bayes_mortality_functions.R`; `common/bayes_mortality_population_helpers.R`; `common/bayes_prior_scale_spec.R`; `common/project_paths.R` | Canonical mortality summaries and fits from mortality_main. |
| `mortality_subgroup/` | `common/bayes_mortality_functions.R`; `common/bayes_mortality_population_helpers.R`; `common/bayes_prior_scale_spec.R`; `common/project_paths.R` | Fits subgroup models from mortality inputs and reuses compatible local subgroup caches. |
| `mortality_publication_cumulative/` | `common/bayes_mortality_functions.R`; `common/bayes_prior_scale_spec.R`; `common/project_paths.R`; `common/tsa_functions.R`; `mortality_publication_cumulative/build_publication_trends_and_tsa.R`; `mortality_publication_cumulative/cumulative_endpoint_reuse.R`; `mortality_publication_cumulative/fit_publication_cumulative_cache.R` | Canonical mortality summaries plus mortality_sensitivity summaries and saved fits for full-evidence endpoints. |
| `mortality_rd_summary/` | `common/bayes_mortality_functions.R`; `common/bayes_prior_scale_spec.R`; `common/mortality_rd_functions.R`; `common/project_paths.R` | Mortality and supporting adjusted RD inputs; compatible local RD fit caches. |
| `los_hospital_main/` | `common/los_main_analysis_functions.R` | Hospital LOS study-level input. |
| `los_icu_main/` | `common/los_main_analysis_functions.R` | ICU LOS study-level input. |
| `los_sensitivity/` | Self-contained R code; package functions | Hospital/ICU LOS inputs and their main-analysis summaries/effect sizes. |
| `process_main/` | `common/process_meta_analysis_functions.R` | Five proportional-process inputs and antibiotic-administration-time input. |
| `process_sensitivity/` | Self-contained R code; package functions | Process inputs and process_main output summaries/effect sizes. |
| `process_subgroup/` | Self-contained R code; package functions | Process inputs and subgroup fields. |
| `mortality_process_rd_metareg/` | `common/antibiotic_time_direction_contract.R`; `common/bayes_mortality_functions.R`; `common/bayes_mortality_population_helpers.R`; `common/bayes_prior_scale_spec.R`; `common/project_paths.R` | Canonical mortality summaries/fits; mortality records containing process covariates and baseline mortality. |
| `mortality_process_metareg/` | `common/antibiotic_time_direction_contract.R`; `common/bayes_mortality_functions.R`; `common/bayes_mortality_population_helpers.R`; `common/bayes_prior_scale_spec.R`; `common/mortality_rd_functions.R`; `common/project_paths.R` | Canonical mortality summaries/fits; mortality and supporting adjusted RD inputs. |
| `mortality_baseline_process/` | `common/bayes_prior_scale_spec.R`; `common/project_paths.R` | Process inputs plus process_main antibiotic-time effect-size output. |
| `mortality_baseline_rd/` | `common/bayes_mortality_functions.R`; `common/bayes_mortality_population_helpers.R`; `common/bayes_prior_scale_spec.R`; `common/mortality_rd_functions.R`; `common/project_paths.R` | Mortality and supporting adjusted RD inputs; compatible local meta-regression caches. |
| `mortality_baseline_distribution/` | `common/bayes_mortality_population_helpers.R` | Mortality input with population definitions. |
| `evidence_map/` | `common/project_paths.R`; `common/small_study_effect_functions.R`; `evidence_map/build_evidence_map_overview.R`; `evidence_map/build_small_study_effect_diagnostics.R` | Mortality, hospital LOS, ICU LOS and process main outputs, plus GRADE certainty input; diagnostics precede overview assembly. |
| `narrative_analysis/` | `common/project_paths.R`; `narrative_analysis/data.R`; `narrative_analysis/definitions.R` | Narrative classifications and antibiotic-time baseline validation input; definitions.R precedes data.R. |

The root runner sources `pipeline_registry.R` and invokes its returned entry points. An initial run should produce primary results, then mortality sensitivity results, before cumulative mortality results. The explicit sequence in the README supplies this ordering while preserving the original registry stages and analysis specifications.

## Analysis input relationship

Study-level analysis inputs are documented in the Supplementary Appendix of the associated article and are referenced by the analysis scripts through their corresponding input structures.

| Input structure under `source_csv/` | Use |
| --- | --- |
| `mortality.csv` | Mortality effects, study design, subgroup fields, baseline mortality and process covariates |
| `mortality_rd_inputs.csv` | Supporting adjusted mortality risk-difference inputs |
| `hospital_los.csv`; `ICU_los.csv` | Hospital and ICU length-of-stay synthesis |
| `process_lactate_measurement.csv`; `process_blood_culture.csv`; `process_antibiotic_proportion.csv`; `process_fluid_resuscitation.csv`; `process_sepsis_bundle.csv` | Proportional process outcomes |
| `process_antibiotic_time.csv` | Time to antibiotic administration |
| `narrative_analysis_combined_conclusions.csv` | Study-level narrative outcomes, baseline classes and characteristics; publication accumulation |
| `narrative_analysis_antibiotic_time_audit.csv` | Analysis input used to validate study identities, antibiotic-time outcome status and baseline classes against the narrative input |
| `grade_certainty.csv` | Outcome-level certainty used in evidence-map summaries |

Input column names, selection rules, effect definitions and validation assertions are retained in the scripts. The antibiotic-time validation input above is part of the narrative analysis input structure. All paths are resolved within the release root or through the existing environment-variable overrides.

## Per-script path and cache inventory

The following flags refer to direct calls or path strings; dependencies inherited through shared helpers are covered by the module table above. `computed_data/` contains outputs and, for downstream modules, upstream analysis inputs. RDS objects are written/read by the fitting and caching code.

| Script | Input path structure | Derived-output path structure | RDS use | Environment variables |
| --- | --- | --- | --- | --- |
| `analysis_scripts/common/antibiotic_time_direction_contract.R` | source_csv/ | Module-supplied paths or no direct output | No direct use | None directly |
| `analysis_scripts/common/bayes_mortality_functions.R` | Module-supplied or upstream objects | Module-supplied paths or no direct output | Yes | `BAYES_META_*` (dynamic sampler controls); helper accepts a batch-control variable |
| `analysis_scripts/common/bayes_mortality_population_helpers.R` | Module-supplied or upstream objects | Module-supplied paths or no direct output | No direct use | None directly |
| `analysis_scripts/common/bayes_prior_scale_spec.R` | Module-supplied or upstream objects | Module-supplied paths or no direct output | No direct use | None directly |
| `analysis_scripts/common/los_main_analysis_functions.R` | Module-supplied or upstream objects | Module-supplied paths or no direct output | No direct use | None directly |
| `analysis_scripts/common/mortality_rd_functions.R` | Module-supplied or upstream objects | Module-supplied paths or no direct output | No direct use | None directly |
| `analysis_scripts/common/process_meta_analysis_functions.R` | source_csv/ | computed_data/ | No direct use | None directly |
| `analysis_scripts/common/project_paths.R` | source_csv/ | computed_data/ | No direct use | None directly |
| `analysis_scripts/common/small_study_effect_functions.R` | Module-supplied or upstream objects | Module-supplied paths or no direct output | No direct use | None directly |
| `analysis_scripts/common/tsa_functions.R` | Module-supplied or upstream objects | Module-supplied paths or no direct output | No direct use | None directly |
| `analysis_scripts/evidence_map/build_evidence_map.R` | Module-supplied or upstream objects | Module-supplied paths or no direct output | No direct use | None directly |
| `analysis_scripts/evidence_map/build_evidence_map_overview.R` | source_csv/ | computed_data/ | No direct use | `BAYES_EVIDENCE_MAP_DATA_DIR` |
| `analysis_scripts/evidence_map/build_small_study_effect_diagnostics.R` | Module-supplied or upstream objects | computed_data/ | No direct use | `BAYES_EVIDENCE_MAP_DATA_DIR` |
| `analysis_scripts/los_hospital_main/build_hospital_los_main.R` | source_csv/ | computed_data/ | No direct use | None directly |
| `analysis_scripts/los_icu_main/build_icu_los_main.R` | source_csv/ | computed_data/ | No direct use | None directly |
| `analysis_scripts/los_sensitivity/build_los_sensitivity.R` | source_csv/ | computed_data/ | No direct use | None directly |
| `analysis_scripts/mortality_baseline_distribution/build_baseline_mortality_distribution.R` | source_csv/ | computed_data/ | No direct use | None directly |
| `analysis_scripts/mortality_baseline_process/build_baseline_process_analysis.R` | source_csv/ | computed_data/ | Yes | None directly |
| `analysis_scripts/mortality_baseline_rd/fit_baseline_mortality_rd_metareg.R` | source_csv/ | computed_data/ | Yes | `BAYES_FORCE_REFIT`; `BAYES_PRIOR_VARIANT`; `S6E_DATA_DIR` |
| `analysis_scripts/mortality_main/build_mortality_forest_cache.R` | source_csv/ | computed_data/ | No direct use | None directly |
| `analysis_scripts/mortality_main/build_mortality_main.R` | Module-supplied or upstream objects | Module-supplied paths or no direct output | No direct use | None directly |
| `analysis_scripts/mortality_main/fit_canonical_main.R` | source_csv/ | computed_data/ | No direct use | `BAYES_CANONICAL_OUT_DIR`; `BAYES_FORCE_REFIT` |
| `analysis_scripts/mortality_process_metareg/fit_process_metareg.R` | source_csv/ | computed_data/ | Yes | `BAYES_CANONICAL_OUT_DIR`; `BAYES_FORCE_REFIT`; `BAYES_PRIOR_VARIANT`; `MORTALITY_PROCESS_METAREG_DATA_DIR` |
| `analysis_scripts/mortality_process_rd_metareg/fit_process_rd_metareg.R` | source_csv/ | computed_data/ | Yes | `BAYES_CANONICAL_OUT_DIR`; `BAYES_FORCE_REFIT`; `BAYES_PRIOR_VARIANT`; `MORTALITY_PROCESS_RD_METAREG_DATA_DIR` |
| `analysis_scripts/mortality_publication_cumulative/build_mortality_publication_cumulative.R` | Module-supplied or upstream objects | Module-supplied paths or no direct output | No direct use | None directly |
| `analysis_scripts/mortality_publication_cumulative/build_publication_trends_and_tsa.R` | source_csv/ | computed_data/ | No direct use | `BAYES_PUBLICATION_CUMULATIVE_DATA_DIR` |
| `analysis_scripts/mortality_publication_cumulative/cumulative_endpoint_reuse.R` | Module-supplied or upstream objects | computed_data/ | Yes | `BAYES_CANONICAL_OUT_DIR`; `MORTALITY_SENSITIVITY_OUT_DIR` |
| `analysis_scripts/mortality_publication_cumulative/fit_publication_cumulative_cache.R` | source_csv/ | computed_data/ | Yes | `BAYES_FORCE_REFIT`; `BAYES_PUBLICATION_CUMULATIVE_DATA_DIR`; `BAYES_PUBLICATION_CUMULATIVE_ENDPOINTS_ONLY` |
| `analysis_scripts/mortality_rd_summary/fit_bayesian_rd_summary.R` | source_csv/ | computed_data/ | Yes | `BAYES_FORCE_REFIT`; `BAYES_RD_SUMMARY_OUT_DIR`; `MORTALITY_RD_INPUT_FILE` |
| `analysis_scripts/mortality_sensitivity/fit_sensitivity_models.R` | source_csv/ | computed_data/ | Yes | `BAYES_CANONICAL_OUT_DIR`; `BAYES_FORCE_REFIT`; `MORTALITY_SENSITIVITY_OUT_DIR` |
| `analysis_scripts/mortality_subgroup/fit_subgroup_analysis.R` | source_csv/ | computed_data/ | Yes | `BAYES_FORCE_REFIT`; `MORTALITY_SUBGROUP_DATA_DIR` |
| `analysis_scripts/narrative_analysis/build_narrative_analysis.R` | source_csv/ | computed_data/ | Yes | None directly |
| `analysis_scripts/narrative_analysis/data.R` | source_csv/ | Module-supplied paths or no direct output | No direct use | None directly |
| `analysis_scripts/narrative_analysis/definitions.R` | Module-supplied or upstream objects | Module-supplied paths or no direct output | No direct use | None directly |
| `analysis_scripts/pipeline_registry.R` | Module-supplied or upstream objects | computed_data/ | No direct use | None directly |
| `analysis_scripts/process_main/build_process_main.R` | source_csv/ | computed_data/ | No direct use | None directly |
| `analysis_scripts/process_sensitivity/build_process_sensitivity.R` | source_csv/ | computed_data/ | No direct use | None directly |
| `analysis_scripts/process_subgroup/build_process_subgroup.R` | source_csv/ | computed_data/ | No direct use | None directly |
| `analysis_scripts/run_all_analysis.R` | Module-supplied or upstream objects | Module-supplied paths or no direct output | No direct use | None directly |

## Environment-variable controls

- Path/output selectors: `BAYES_CANONICAL_OUT_DIR`, `BAYES_EVIDENCE_MAP_DATA_DIR`, `BAYES_PUBLICATION_CUMULATIVE_DATA_DIR`, `BAYES_RD_SUMMARY_OUT_DIR`, `MORTALITY_PROCESS_METAREG_DATA_DIR`, `MORTALITY_PROCESS_RD_METAREG_DATA_DIR`, `MORTALITY_SENSITIVITY_OUT_DIR`, `MORTALITY_SUBGROUP_DATA_DIR`, `S6E_DATA_DIR`, and `MORTALITY_RD_INPUT_FILE`. Consult each script for its default path and resolution logic.
- Cache controls: `BAYES_FORCE_REFIT` requests refitting where supported; `BAYES_PUBLICATION_CUMULATIVE_ENDPOINTS_ONLY` selects cumulative endpoint refresh.
- Prior selection: `BAYES_PRIOR_VARIANT` selects existing prespecified variants; its canonical default is retained.
- Shared sampler controls: `BAYES_META_ITER`, `BAYES_META_WARMUP`, `BAYES_META_CHAINS`, `BAYES_META_CORES`, `BAYES_META_SEED`, `BAYES_META_REFRESH`, `BAYES_META_ADAPT_DELTA`, and `BAYES_META_MAX_TREEDEPTH`. Script defaults are preserved.
- The shared helper defines a `BAYES_ALLOW_NONCANONICAL` batch guard; no retained entry point calls that guard.

## Computational considerations

Bayesian models may require substantial computation time. Some scripts support cached fits and validate input, model, prior and sampling metadata before reuse. Primary fitting scripts fit models when compatible caches are absent; downstream consumers require their upstream producers to have completed. Cumulative analyses reuse full-evidence endpoints from primary and sensitivity analyses. TSA boundary calibration and exploratory meta-regression can require additional runtime. The packaging checks parse scripts and inspect dependencies; they do not fit statistical models.
