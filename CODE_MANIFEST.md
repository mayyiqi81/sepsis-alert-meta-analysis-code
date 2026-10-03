# Code manifest

| Analysis area | Script/module | Purpose |
| --- | --- | --- |
| Shared functions | `analysis_scripts/common/` | Bayesian models and priors, population and effect definitions, LOS/process synthesis, small-study diagnostics, TSA, and project paths |
| Mortality primary | `analysis_scripts/mortality_main/` | Canonical Bayesian random-effects mortality synthesis and derived forest summaries |
| Mortality sensitivity | `analysis_scripts/mortality_sensitivity/` | Prespecified additional-uncertainty, power-prior-type, exclusion, adjusted-estimate and population analyses |
| Mortality subgroup | `analysis_scripts/mortality_subgroup/` | Mortality subgroup estimates and between-group comparisons |
| Cumulative mortality and TSA | `analysis_scripts/mortality_publication_cumulative/` | Publication accumulation, Bayesian cumulative mortality synthesis, and RCT trial sequential analysis |
| Mortality RD | `analysis_scripts/mortality_rd_summary/` | Supporting Bayesian mortality risk-difference synthesis |
| Hospital LOS | `analysis_scripts/los_hospital_main/` | Random-effects mean-difference synthesis for hospital length of stay |
| ICU LOS | `analysis_scripts/los_icu_main/` | Random-effects mean-difference synthesis for ICU length of stay |
| LOS sensitivity | `analysis_scripts/los_sensitivity/` | Hospital and ICU length-of-stay sensitivity analyses |
| Process main analyses | `analysis_scripts/process_main/` | Five proportional care-process outcomes and time to antibiotic administration |
| Process sensitivity | `analysis_scripts/process_sensitivity/` | Prespecified sensitivity analyses for proportional and time process outcomes |
| Process subgroup | `analysis_scripts/process_subgroup/` | Process-outcome subgroup synthesis and between-group comparisons |
| Process–mortality meta-regression | `analysis_scripts/mortality_process_rd_metareg/` | Primary univariable mortality meta-regression using process absolute changes/time differences and baseline mortality |
| Meta-regression sensitivity | `analysis_scripts/mortality_process_metareg/` | Supporting and sensitivity analyses of study-level process–mortality associations |
| Baseline process analyses | `analysis_scripts/mortality_baseline_process/` | Association between baseline process levels and process improvement |
| Baseline mortality analyses | `analysis_scripts/mortality_baseline_rd/` | Bayesian association between baseline mortality and mortality risk difference |
| Baseline mortality distribution | `analysis_scripts/mortality_baseline_distribution/` | Descriptive baseline mortality distributions by analytic population |
| Evidence map | `analysis_scripts/evidence_map/` | Combine upstream outcome summaries, GRADE certainty and small-study-effect diagnostics |
| Narrative analysis | `analysis_scripts/narrative_analysis/` | Study-level narrative classifications and process, population and alert-characteristic summaries |
| Orchestration | `analysis_scripts/pipeline_registry.R`; `analysis_scripts/run_all_analysis.R` | Register and invoke the 18 analysis domains |

Within `mortality_publication_cumulative/`, `build_publication_trends_and_tsa.R` implements publication accumulation and RCT TSA; `fit_publication_cumulative_cache.R` fits cumulative Bayesian models, and `cumulative_endpoint_reuse.R` validates and reuses upstream endpoints.
