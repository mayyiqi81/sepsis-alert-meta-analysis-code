# Effectiveness of clinically implemented automated sepsis alerting systems: a systematic review and meta-analysis

## Scope

This repository contains the analysis code used for the systematic review and meta-analysis *Effectiveness of clinically implemented automated sepsis alerting systems: a systematic review and meta-analysis*. The extracted study-level data underlying the review and evidence syntheses are provided in the Supplementary Appendix of the associated article. Analysis code is publicly available in this versioned repository.

## Repository structure

`analysis_scripts/` contains the analysis modules and their shared code. `pipeline_registry.R` records the module entry points; `run_all_analysis.R` selects and invokes them. See [CODE_MANIFEST.md](CODE_MANIFEST.md) for module purposes and [DEPENDENCIES.md](DEPENDENCIES.md) for code, input and upstream-output relationships.

| Module | Analysis |
| --- | --- |
| `analysis_scripts/common/` | Bayesian models and priors, population and effect definitions, LOS/process synthesis, small-study diagnostics, TSA, and project paths |
| `analysis_scripts/mortality_main/` | Canonical Bayesian random-effects mortality synthesis and derived forest summaries |
| `analysis_scripts/mortality_sensitivity/` | Prespecified additional-uncertainty, power-prior-type, exclusion, adjusted-estimate and population analyses |
| `analysis_scripts/mortality_subgroup/` | Mortality subgroup estimates and between-group comparisons |
| `analysis_scripts/mortality_publication_cumulative/` | Publication accumulation, Bayesian cumulative mortality synthesis, and RCT trial sequential analysis |
| `analysis_scripts/mortality_rd_summary/` | Supporting Bayesian mortality risk-difference synthesis |
| `analysis_scripts/los_hospital_main/` | Random-effects mean-difference synthesis for hospital length of stay |
| `analysis_scripts/los_icu_main/` | Random-effects mean-difference synthesis for ICU length of stay |
| `analysis_scripts/los_sensitivity/` | Hospital and ICU length-of-stay sensitivity analyses |
| `analysis_scripts/process_main/` | Five proportional care-process outcomes and time to antibiotic administration |
| `analysis_scripts/process_sensitivity/` | Prespecified sensitivity analyses for proportional and time process outcomes |
| `analysis_scripts/process_subgroup/` | Process-outcome subgroup synthesis and between-group comparisons |
| `analysis_scripts/mortality_process_rd_metareg/` | Primary univariable mortality meta-regression using process absolute changes/time differences and baseline mortality |
| `analysis_scripts/mortality_process_metareg/` | Supporting and sensitivity analyses of study-level process–mortality associations |
| `analysis_scripts/mortality_baseline_process/` | Association between baseline process levels and process improvement |
| `analysis_scripts/mortality_baseline_rd/` | Bayesian association between baseline mortality and mortality risk difference |
| `analysis_scripts/mortality_baseline_distribution/` | Descriptive baseline mortality distributions by analytic population |
| `analysis_scripts/evidence_map/` | Combine upstream outcome summaries, GRADE certainty and small-study-effect diagnostics |
| `analysis_scripts/narrative_analysis/` | Study-level narrative classifications and process, population and alert-characteristic summaries |

## Analysis framework

- **Mortality:** Bayesian random-effects meta-analysis on the log risk-ratio scale, with separate RCT and NRS syntheses and an integrated synthesis. The integrated primary analysis incorporates prespecified additional uncertainty for NRS. Alternative additional-uncertainty settings and power-prior-type down-weighting are separate sensitivity strategies. Supporting analyses include subgroups, cumulative mortality analyses and RCT trial sequential analysis.
- **Length of stay:** random-effects meta-analysis of mean differences for hospital and ICU length of stay, with sensitivity analyses.
- **Care processes:** risk differences for proportional outcomes and mean differences for time to antibiotic administration, with sensitivity and subgroup analyses.
- **Exploratory analyses:** study-level process–mortality meta-regression, baseline mortality analyses, baseline process analyses and related descriptive summaries. These analyses assess associations.

## Software

The manuscript reports R 4.4.2, brms 2.23.0, metafor 4.8.0 and mvtnorm 1.3.3. Additional R packages are loaded within individual analysis scripts; their names are listed in [DEPENDENCIES.md](DEPENDENCIES.md). Bayesian fitting uses the Stan environment configured for brms.

## Analysis inputs

The study-level inputs corresponding to these analyses are documented in the Supplementary Appendix, including mortality, length-of-stay, process-outcome, adjusted-effect, subgroup, and explanatory-analysis data.

The Supplementary Appendix contains study characteristics (Table S1A), comparator and alert characteristics (Tables S2A–S2D), extraction and synthesis rules (Supplementary outcome synthesis 1A–1C and 2A–2C), mortality evidence (Tables S3A–S3C), length-of-stay evidence (Tables S4A–S4D), process-outcome evidence (Tables S5A–S5I), supporting analyses (Tables S6A–S6B), risk-of-bias assessments (S7) and GRADE assessments (S8).

## Running the analyses

Paths are relative to this repository root. The scripts reference study-level inputs through the `source_csv/` input structure and write derived summaries and model objects under `computed_data/`. The input names and upstream dependencies are documented in [DEPENDENCIES.md](DEPENDENCIES.md).

```text
Study-level analysis inputs
    -> analysis_scripts/
    -> model fitting / sensitivity / subgroup analyses
    -> derived analysis summaries
```

To inspect the registered entry points without fitting models, run from the repository root:

```sh
Rscript analysis_scripts/run_all_analysis.R --list
Rscript analysis_scripts/run_all_analysis.R --dry-run
```

For an initial analysis run with inputs arranged in the corresponding input structure, invoke modules in dependency order. In particular, cumulative analyses reuse both primary and sensitivity mortality results. The following sequence resolves these dependencies explicitly:

```sh
Rscript analysis_scripts/run_all_analysis.R --domain=mortality_main,los_hospital_main,los_icu_main,process_main
Rscript analysis_scripts/run_all_analysis.R --domain=mortality_sensitivity
Rscript analysis_scripts/run_all_analysis.R --domain=los_sensitivity,process_sensitivity,process_subgroup,mortality_subgroup,mortality_rd_summary,mortality_baseline_distribution,mortality_baseline_process,mortality_baseline_rd,mortality_process_metareg,mortality_process_rd_metareg
Rscript analysis_scripts/run_all_analysis.R --domain=mortality_publication_cumulative
Rscript analysis_scripts/run_all_analysis.R --domain=evidence_map,narrative_analysis
```

Some Bayesian analyses support cached model objects to reduce repeated computation. When validated cached model objects are absent, the corresponding fitting scripts fit models from the analysis inputs. Modules that consume upstream fits should be run after their producers. Bayesian fitting, meta-regression and TSA can require substantial computation time. Environment-variable controls are listed in [DEPENDENCIES.md](DEPENDENCIES.md); model specifications remain as defined in the scripts.

## Versioning

Versioned releases correspond to defined manuscript stages. The v1.0.0-submission release represents the analysis code accompanying the initial journal submission. Subsequent peer-review revisions, if needed, will be preserved as separate versioned releases.

## Citation

Please cite the associated systematic review and meta-analysis when using this code. Full bibliographic details will be added after publication.
