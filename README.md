# ICU Extubation Outcomes: Caregiver Volume, Failed Extubation Rate, and Patient Survival

A clinical data science project investigating whether caregiver-level failed
extubation rates predict patient outcomes in the ICU, built on MIMIC-IV v3.1.

---

## Research Question

Extubation — removing a patient from mechanical ventilation — is a high-stakes
clinical decision. Failed extubation (reintubation within 72 hours) is associated
with significantly higher mortality, longer ICU stays, and increased risk of
ventilator-associated pneumonia. The decision of *when* and *how* to extubate
involves substantial provider judgment.

This project asks: **does a caregiver's historical failed extubation rate predict
their patients' 12-month survival?** We also ask whether inter-provider variation
in FE rates can be explained by patient case mix, and whether any apparent
survival effect survives propensity score adjustment for measured confounders.

---

## Key Findings

**The main finding is a null result, and understanding why is the contribution.**

Initial analyses using the full dataset — including extubation events inferred
from ventilation records when explicit documentation was absent — suggested a
significant association between caregiver FE rate and patient survival (PSM OR =
0.733, p < 0.001). Subsequent sensitivity analyses revealed this was driven by
a data quality artifact rather than a true caregiver performance effect.

**Data quality analysis:**
- 87.8% of patients in the high-FE-rate caregiver group (top quartile, FE rate
  ≥ 32%) had extubation events classified as *inferred* from ventilation records,
  versus 39.4% in the low-FE-rate group.
- Inferred events are reconstructed algorithmically from ventilation segment
  boundaries and are not explicitly documented procedure events. They generate
  spurious failed extubation classifications that inflate caregiver FE rates.
- Caregivers with predominantly inferred events also have patients with worse
  baseline outcomes — creating an apparent but non-causal association.

**Sensitivity analysis restricted to explicitly documented events:**
- FE rate distribution becomes clinically plausible: most caregivers 10–19%,
  tail reaching ~45% rather than 100%.
- PSM OR attenuates from 0.733 to 0.810 (p = 0.021).
- Splitting by caregiver volume: survival difference is −8.9pp in low-volume
  caregivers and −1.3pp in high-volume caregivers — concentrated exactly where
  small-sample rate instability would predict it.
- Volume-weighted correlation between caregiver FE rate and patient survival
  at the caregiver level: r = −0.134, **p = 0.124** (not significant).

**Mixed effects models:**
- Mixed effects logistic regression (glmmTMB, frequentist; rstanarm, Bayesian)
  with random intercepts for caregiver (520 groups) and ICU unit (13 groups)
  confirms the null result with proper uncertainty quantification.
- ICU unit random effect variance: 0.264 (glmmTMB); posterior mean 0.4,
  95% CI (0.1, 1.2) (rstanarm) — unit-level factors explain meaningful
  survival variation independently of patient severity.
- Caregiver random effect variance: ~0 boundary estimate (glmmTMB); posterior
  mean 0.0, 97.5% = 0.1 (rstanarm) — no residual caregiver-level variation
  once unit and patient covariates are controlled.
- `caregiver_fe_rate` fixed effect: coef −0.083, p = 0.832 (glmmTMB);
  posterior mean 0.0, 95% CI (−0.1, 0.1) (rstanarm). Consistent null across
  both frameworks.
- Bayesian shrinkage: the rstanarm posterior provides honest uncertainty
  quantification where glmmTMB hits the boundary of the parameter space,
  confirming the null is not a numerical artifact.

**Conclusion:** Caregiver FE rate as derivable from MIMIC-IV procedure events
is not a reliable predictor of patient 12-month survival once data quality is
adequately controlled. The apparent effect in the full dataset reflects inferred
event misclassification and small-sample caregiver rate instability rather than
true provider performance differences.

**Case mix analysis:**
- LASSO feature selection across 12 ICD chapters retained minimal signal.
- LDA topic modeling on 388 CCSR diagnosis codes produced 8 clinically coherent
  condition profiles, but cluster separation was weak (maximum variance
  explained ~15%), reflecting the multimorbid, overlapping nature of ICU
  presentations.
- Case mix does not explain inter-provider FE rate variation.

---

## Limitations and Open Questions

**Inferred event imputation:** The caregiver ID fallback cascade (recorded
procedureevents caregiver → nearest chartevents caregiver within ±15 minutes
→ most frequent caregiver for that ICU stay) assigns caregivers to events where
the true responsible provider is unknown. This introduces systematic noise into
caregiver-level metrics. Improving this imputation — or restricting entirely to
explicitly documented events — is the most important methodological next step.

**Explicit-only sample size:** Restricting to explicit events substantially
reduces the analyzable cohort, limiting statistical power for caregiver-level
analyses. Only 133 caregivers with more than 5 explicit extubations are
available for volume-weighted analyses.

**Outcome attribution:** `survival_12mo` is derived from date of death in
MIMIC-IV regardless of cause. Deaths unrelated to the index extubation are
counted equally, attenuating any caregiver effect on extubation-attributable
mortality.

**Unmeasured confounding:** PSM controls only for observed covariates. Shift
patterns, time of day of extubation, and attending physician supervision are
not captured in MIMIC-IV. ICU unit is now included as a random effect in the
mixed effects models (Section 10), partially addressing unit-level confounding.

**External validity:** MIMIC-IV is a single-center dataset (Beth Israel
Deaconess Medical Center). Provider-level metrics may not generalize.

---

## Project Structure

```
.
├── sql/
│   └── ICU_Last_Extubations.sql            # BigQuery pipeline: MIMIC-IV → analysis dataset
├── r/
│   ├── MIMIC_analysis_clean.R              # Primary analysis: cohort construction,
│   │                                       #   logistic regression, caregiver analysis,
│   │                                       #   LASSO, random forest, 5-fold CV
│   ├── psm_analysis.R                      # Propensity score analysis: CCSR clustering,
│   │                                       #   PSM, IPW, data quality sensitivity analyses
│   └── generate_synthetic_extubation_data.R  # Synthetic cohort generator (see below)
├── data/
│   ├── synthetic_last_extubations.csv      # 800-patient synthetic cohort (Tableau input)
│   └── synthetic_caregiver_rate.csv        # 80-caregiver summary (Tableau input)
├── figures/
│   └── psm_love_plot.png                   # Covariate balance before and after PSM
├── tableau/
│   └── tableau_dashboard_spec.md           # Full build specification for Tableau dashboards
└── README.md
```

**Note:** The `data/` directory contains only synthetic data (see below).
No MIMIC-IV patient records are included in this repository.

---

## Data

### Source: MIMIC-IV v3.1

The analysis was conducted on MIMIC-IV, a large de-identified EHR dataset from
the Beth Israel Deaconess Medical Center ICU, accessed via PhysioNet and Google
BigQuery. MIMIC-IV contains data for over 65,000 ICU patients admitted between
2008 and 2022.

Access requires credentialing through PhysioNet and agreement to the MIMIC-IV
Data Use Agreement. The DUA prohibits redistribution of the data or derivatives
that could be used to re-identify patients. **No MIMIC-IV data appears in this
repository.**

### Synthetic cohort

The CSV files in `data/` are fully synthetic. They were generated by
`generate_synthetic_extubation_data.R` using distributional parameters drawn
exclusively from published literature — no MIMIC records were used as input.

Key parameter sources:

| Variable | Source | Value used |
|---|---|---|
| Age distribution | Johnson et al., *Scientific Data* 2016 (MIMIC-III) | Median 65.8, IQR 52.8–77.8 |
| Failed extubation rate | Fernandez et al., *Respiratory Care* 2024 | 15.4% (72h window) |
| Hospital mortality | Johnson et al. 2016 | ~11.5% |
| Charlson comorbidity index | MIMIC-Sepsis, arXiv 2025 | Median 5 |
| SOFA score | Published MV cohorts | Mean ~8 |
| Ventilation hours | Published MV cohorts | Median ~96h |
| Vasopressor use | Published MV cohorts | ~30% of patients |

The synthetic data is appropriate for portfolio demonstration and public
visualization. It is not appropriate for clinical inference.

### SQL pipeline

`ICU_Last_Extubations.sql` is a single BigQuery CTE chain that:

1. Combines explicit intubation/extubation procedure events (itemids 224385,
   227194, 225468, 225477) with inferred events from ventilation segments
   (≥4h InvasiveVent), deduplicating inferred events within 30 minutes of
   explicit events
2. Classifies each intubation as `surgical`, `medical-respiratory`, or
   `medical-non-respiratory` based on proximate surgical procedure codes and
   respiratory ICD chapter flags
3. Selects the *last* extubation per patient (one row per patient)
4. Joins patient demographics, Charlson comorbidity index, norepinephrine
   equivalent dose, first-day SOFA score, hospital mortality, primary ICD
   chapter, and an array of CCSR categories (up to seq_num 10)
5. Assigns caregiver ID via a three-level fallback: recorded procedureevents
   caregiver → nearest chartevents caregiver within ±15 minutes → most frequent
   caregiver for that ICU stay
6. Computes `caregiver_fe_rate` and `caregiver_n` at query time
7. Applies final filters: vent_hours ∈ (1, 1000), norepinephrine < 1
   mcg/kg/min, failed_extubations < 10

---

## Analysis

### Primary analysis (`MIMIC_analysis_clean.R`)

| Section | Content |
|---|---|
| 0 | Libraries and setup |
| 1 | Data ingestion (SQL output CSVs → R) |
| 2 | Cohort construction: filtering, centering, `cleandata` / `supercleandata` |
| 3 | Patient-level logistic regression for 12-month survival: three GLM variants (logreg1–3) plus mixed effects models — frequentist (glmmTMB) and Bayesian (rstanarm::stan_glmer) — with random intercepts for caregiver and ICU unit |
| 4 | Caregiver-level analysis: volume vs. FE rate, volume vs. survival, volume vs. hospital mortality |
| 5 | Diagnosis case mix: LASSO on wide ICD chapter × caregiver matrix; ICD indicator construction |
| 6 | Full interaction model (`logreg_icd`) with 5-fold cross-validation; classification metrics |
| 7 | Random forest on caregiver case mix; variable importance (%IncMSE, IncNodePurity) |
| 8 | RSBI-based logistic regression (separate cohort with pre-extubation vitals) |
| 9 | Appendix: patient timeline visualization (data quality checking) |

### Propensity score and sensitivity analysis (`psm_analysis.R`)

| Section | Content |
|---|---|
| 0 | Libraries and setup |
| 1 | Data ingestion: `last_extubations_clean.csv`, `patient_ccsr_long.csv` |
| 2 | CCSR document-term matrix construction |
| 3 | LDA topic modeling (K = 5, 8, 10, 12); perplexity comparison; topic interpretation |
| 4a | Topic proportion extraction (soft cluster membership) |
| 4b | K-means clustering on binary CCSR matrix; variance explained comparison |
| 4c | Revised covariate strategy: high-prevalence CCSR binary flags (≥10% prevalence) |
| 5 | PSM dataset construction; treatment variable definition (top vs. bottom quartile) |
| 6 | Propensity score matching: nearest neighbor 1:1, caliper = 0.2 SD |
| 7 | Treatment effect estimation: unadjusted and doubly-adjusted outcome models |
| 8 | IPW sensitivity analysis; covariate balance love plot; effect modification test |
| 9 | Data quality sensitivity: explicit events only; volume-stratified analysis; weighted correlation |

### Results summary

| Analysis | OR | 95% CI | p | Note |
|---|---|---|---|---|
| PSM — full dataset | 0.733 | 0.643–0.836 | <0.001 | Driven by inferred event artifact |
| IPW — full dataset | 0.789 | 0.735–0.846 | <0.001 | Driven by inferred event artifact |
| PSM — explicit events only | 0.810 | 0.678–0.968 | 0.021 | Concentrated in low-volume caregivers |
| Weighted correlation (explicit, caregiver-level) | r = −0.134 | — | 0.124 | **Not significant** |
| glmmTMB — caregiver FE rate fixed effect | coef −0.083 | — | 0.832 | RE variance ~0 (boundary) |
| rstanarm — caregiver FE rate posterior | mean 0.0 | −0.1 to 0.1 | — | **Not significant; Rhat = 1.0** |

---

## Tableau Dashboard

The interactive dashboard is published at:
**[ICU Extubation Outcomes Dashboard](https://public.tableau.com/views/ExtubationsSimulated/Dashboard1-CohortOverview)**

Three dashboards built on a synthetic cohort:

**1. Patient Cohort Overview** — age distribution, intubation type breakdown,
ICD chapter volume and FE rates, acuity scatter (SOFA vs. Charlson).

**2. Provider Volume & Outcomes** — caregiver-level scatter of volume vs. FE
rate and volume vs. 12-month survival rate, with trend lines and volume tier
box plots.

**3. Risk Score Simulator** — parameter sliders drive a logistic regression
equation implemented as Tableau calculated fields, showing estimated 12-month
survival probability and factor contributions.

> The simulator uses stylized logistic coefficients consistent with published
> effect directions. It does not reproduce fitted coefficients from the MIMIC
> analysis and must not be used for clinical decision-making.

Full build instructions are in `tableau/tableau_dashboard_spec.md`.

---

## Reproducing the Synthetic Data

```r
install.packages(c("MASS", "dplyr", "tidyr", "purrr", "readr",
                   "tibble", "stringr", "forcats"))
source("r/generate_synthetic_extubation_data.R")
# Outputs: synthetic_last_extubations.csv, synthetic_caregiver_rate.csv
```

The generator is fully deterministic given `set.seed(237)`.

---

## Reproducing the Full Analysis

Full reproduction requires MIMIC-IV access:

1. Apply for access at [physionet.org](https://physionet.org/content/mimiciv/3.1/)
2. Run `sql/ICU_Last_Extubations.sql` on BigQuery
3. Export result to `data/last_extubations_clean.csv`
4. Run the CCSR long query (see `r/psm_analysis.R` header) and export to
   `data/patient_ccsr_long.csv`
5. Run `r/MIMIC_analysis_clean.R`
6. Run `r/psm_analysis.R`

The SQL query references a derived ventilation table and CCSR mapping table
stored in a private BigQuery project. The ventilation table follows the
[MIMIC-IV-derived](https://github.com/MIT-LCP/mimic-iv) pipeline. The CCSR
mapping table uses the publicly available
[AHRQ CCSR](https://hcup-us.ahrq.gov/toolssoftware/ccsr/ccs_refined.jsp)
crosswalk for ICD-10-CM codes.

---

## Dependencies

**R:** `tidyverse`, `lubridate`, `splines`, `glmnet`, `randomForest`, `caret`,
`MASS`, `MatchIt`, `cobalt`, `WeightIt`, `topicmodels`, `tidytext`, `Matrix`,
`slam`, `weights`, `glmmTMB`, `rstanarm`

**SQL:** Google BigQuery with MIMIC-IV v3.1 access

**Visualization:** Tableau Desktop 2020.1+

---

## Acknowledgments

MIMIC-IV was developed by the MIT Laboratory for Computational Physiology and
is made available through PhysioNet. This project was developed in collaboration
with Andrew Hersh, MD, who contributed clinical domain expertise, motivated the
caregiver-level research questions, and proposed the propensity score matching
design.

Johnson, A., Bulgarelli, L., Shen, L. et al. MIMIC-IV, a freely accessible
electronic health record dataset. *Sci Data* 10, 1 (2023).
https://doi.org/10.1038/s41597-022-01899-x

---

## License

Code in this repository (SQL, R) is released under the MIT License.
Synthetic data files are released under CC0 (public domain).
No MIMIC-IV data is included; use of MIMIC-IV is governed by the
[PhysioNet Credentialed Health Data License](https://physionet.org/content/mimiciv/view-license/3.1/).
