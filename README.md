# ICU Extubation Outcomes: Caregiver Failed Extubation Rate and Patient Survival

A clinical data science project investigating whether caregiver-level failed
extubation rates predict patient 12-month survival in the ICU, using MIMIC-IV
v3.1 and a multi-method analytical approach including propensity score matching,
mixed effects logistic regression, and caregiver-level partial correlation.

---

## Research Question

Extubation — removing a patient from mechanical ventilation — is a high-stakes
clinical decision. Failed extubation (reintubation within 72 hours) is associated
with significantly higher mortality, longer ICU stays, and increased risk of
ventilator-associated pneumonia. This project asks: **does a caregiver's
historical failed extubation rate predict their patients' 12-month survival,
after accounting for patient case mix, illness severity, and ICU unit?**

A secondary question motivated by clinical intuition: is there an optimal
FE rate range, where rates that are too low reflect overly conservative practice
and rates that are too high reflect excessive aggressiveness?

---

## Key Findings

**Cohort:** 12,662 patients with a documented last ICU extubation event in
MIMIC-IV v3.1 — one of the largest published ICU extubation cohorts derived
from this dataset. Of these, 4,387 had explicitly documented procedure events
with a genuine caregiver assignment and form the analytical cohort (271
caregivers, 12 ICU units).

**Propensity score matching** comparing patients of high-FE-rate caregivers
(top quartile, ≥5.3%) to patients of low-FE-rate caregivers (bottom quartile,
≤2.9%), after matching on illness severity (SOFA, Charlson comorbidity index,
vasopressor dose, ventilation hours), intubation type, and 21 CCSR diagnosis
flags:

- 12-month survival: **63.1%** (high FE caregivers) vs **78.0%** (low FE
  caregivers) — a 14.9 percentage point absolute difference
- PSM odds ratio (doubly adjusted): **OR = 0.571 (95% CI 0.465–0.700)**
- IPW sensitivity analysis: **OR = 0.699 (95% CI 0.616–0.792)**
- No effect modification by intubation type (p = 0.241)

**Mixed effects logistic regression** with random intercepts for caregiver
(271 groups) and ICU unit (12 groups):

- ICU unit random effect variance: 0.265 (glmmTMB), posterior mean 0.41
  (rstanarm) — ICU unit explains meaningful survival variation independently
  of patient severity
- Caregiver random effect variance: ~0 — no residual caregiver-level variation
  once unit and patient covariates are controlled
- caregiver_fe_rate fixed effect: coef −2.80, p = 0.097 (glmmTMB); posterior
  mean −0.067, 95% CI (−0.148, 0.017), P(effect < 0) = 94.4% (rstanarm)

**Caregiver-level partial correlation** (134 caregivers with ≥5 sample
patients):

- FE rate independently predicts survival controlling for volume:
  **r = −0.363, p < 0.001**
- Volume does not predict survival controlling for FE rate:
  r = −0.121, p = 0.167
- FE rate is the signal; caregiver volume is not an independent predictor

**Functional form** of the FE rate — survival relationship (caregivers with
FE rate > 0, N = 111):

- Asymptotic exponential model: P(survival) = b + (1−b) × exp(−k × FE rate)
- **b = 0.590 (95% CI 0.521–0.633)** — survival floor. Approximately 59% of
  patients survive 12 months regardless of caregiver FE rate, representing
  irreducible mortality driven by underlying illness severity.
- **k = 47.2 (95% CI 31.2–71.0)** — decay rate. Most of the survival
  difference occurs in the 0–8% FE rate range.
- The remaining ~41% of patients represent outcomes potentially sensitive to
  caregiver extubation practice — the upper bound on what quality improvement
  could achieve in this population.
- **There is no evidence of a sweet spot.** Survival decreases monotonically
  with FE rate; the relationship is best described by diminishing marginal harm
  rather than a U-shaped optimum.

**Physiological trajectory analysis** using all explicit extubation events
(6,063 events, 5,505 patients) with per-event reintubation within 72h as
outcome:

- RSBI and P/F ratio in the 72h before extubation do not predict per-event
  extubation failure beyond the level at the moment of extubation
- RSBI level at extubation: significant (p = 0.032), mean 52.3 (failed) vs
  48.7 (succeeded)
- RSBI slope 72h before extubation: not significant (p = 0.634)
- P/F ratio level and slope: not significant (p = 0.505, p = 0.726)
- The two groups are physiologically indistinguishable throughout the 72h
  pre-extubation window — consistent with clinicians already conditioning on
  a plateau criterion before extubating

---

## Data Quality and Cohort Restriction

MIMIC-IV extubation events are either explicitly documented procedure events
or algorithmically inferred from ventilation segment boundaries. Caregiver
identity is often missing or unreliable for inferred events, and imputation
approaches tested so far have not produced reliable caregiver assignments for
this subset.

As a result, the project operates on two cohorts: a broad descriptive cohort
of 12,662 patients that characterizes the full population of ICU extubations,
and an analytical cohort of 4,387 patients restricted to explicitly documented
events with a genuine caregiver assignment. The analytical cohort is smaller
but still large enough to support the mixed effects and PSM analyses, and
the restriction produces a clinically plausible caregiver FE rate distribution
(median 4.3%, IQR 2.9–5.3%). Improving caregiver assignment for inferred
events remains an open methodological problem.

---

## Project Structure

```
.
├── sql/
│   ├── ICU_Last_Extubations.sql        # BigQuery pipeline: MIMIC-IV → last extubation per patient
│   └── all_extubations.sql             # BigQuery pipeline: all explicit extubation events
│                                       #   with per-event failed_extubation_flag outcome
│                                       #   and leave-one-out caregiver FE rate
├── r/
│   ├── 00_run_all.R                    # Master script: sources 01, 02, 03 in order
│   ├── 01_cohort_and_descriptive.R     # Cohort construction, Table 1, descriptive plots
│   ├── 02_logistic_and_mixed_effects.R # glmmTMB and rstanarm mixed effects models
│   ├── 03_psm_and_sensitivity.R        # PSM, IPW, partial correlation, functional form,
│   │                                   #   asymptotic exponential model (nls_fit)
│   ├── 04_trajectory_analysis.R        # Per-event RSBI and P/F trajectory analysis
│   │                                   #   requires bigrquery and BigQuery access
│   └── generate_synthetic_extubation_data.R  # Synthetic cohort generator (Tableau)
├── data/
│   ├── synthetic_last_extubations.csv  # 800-patient synthetic cohort (Tableau input)
│   └── synthetic_caregiver_rate.csv    # 80-caregiver summary (Tableau input)
├── figures/
│   ├── psm_love_plot.png               # Covariate balance before and after PSM
│   ├── fe_rate_survival_exponential.png # Asymptotic exponential fit with survival floor
│   └── preextubation_trajectories_by_event_outcome.png  # RSBI and P/F ratio in 72h
│                                                        #   before extubation by outcome
├── tableau/
│   └── tableau_dashboard_spec.md       # Full build specification for Tableau dashboards
└── README.md
```

**Note:** The `data/` directory contains only synthetic data. No MIMIC-IV
patient records are included in this repository.

---

## Analytical Pipeline

### Script 01: Cohort Construction and Descriptive Analysis

Constructs two cohorts from `last_extubations.csv`:

**`patient_cohort`** (N = 12,662) — all patients with a documented last
extubation and complete covariates. Includes both explicit and inferred events.
Used for population description and Table 1.

**`explicit_extubations`** (N = 4,387) — restricted to explicitly documented
procedure events with a genuine caregiver assignment (caregiver_fe_rate not
NULL) and caregiver_n > 10. Used for all caregiver-level analyses.

Produces Table 1 stratified by tube event source and by FE rate quartile,
caregiver-level summary statistics, and population visualizations. Saves
`cohorts.RData` for downstream scripts.

### Script 02: Mixed Effects Models

Fits two models on `explicit_extubations`:

**glmmTMB** (frequentist): random intercepts for caregiver_id and
first_careunit; centered continuous predictors.

**rstanarm** (Bayesian): same structure; scaled predictors (unit variance)
to resolve Stan initialization; weakly informative priors normal(0, 2.5) on
fixed effects, decov(scale = 1) on random effect covariance. Runtime ~30–40
minutes. Saves `mixed_effects_models.RData`.

### Script 03: PSM and Sensitivity Analysis

**CCSR binary flags:** 21 CCSR codes with ≥10% prevalence in the explicit
cohort (RSP012 excluded — near-universal in ventilated patients) serve as
PSM covariates. LDA topic modeling and k-means clustering were explored but
showed weak patient differentiation (short-document problem with median ~4
CCSR codes per patient); direct binary flags are used instead.

**PSM:** method = "quick", ratio 1:1, caliper 0.2 SD on logit propensity
score. Treatment defined as top vs bottom FE rate quartile (≥5.3% vs ≤2.9%).

**IPW:** WeightIt package, ATE estimand, as sensitivity analysis.

**Partial correlation:** ppcor package, caregiver-level (N = 134 caregivers
with ≥5 sample patients), controlling for caregiver volume.

**Functional form:** weighted asymptotic exponential model P(survival) = b +
(1−b)×exp(−k×FE rate) at the caregiver level (nls_fit); bootstrap CI on
parameters (2,000 resamples). Quadratic model retained for comparison.

Saves `psm_results.RData` and two figures.

### SQL: all_extubations.sql

Produces one row per explicit extubation event (not just the last per patient)
for trajectory analysis. Key differences from `ICU_Last_Extubations.sql`:

- All explicit extubation events per patient, not just the last
- `failed_extubation_flag`: 1 if this specific extubation was followed by
  reintubation within 72h (both events explicit)
- `caregiver_fe_rate_loo`: leave-one-out caregiver FE rate — computed from
  all other events by this caregiver, avoiding circularity
- No diagnosis joins — lean output for trajectory analysis

Output: 6,063 events, 5,505 patients, 110 caregivers, 13 ICU units,
overall FE rate 8.3%.

### Script 04: Trajectory Analysis

Pulls RSBI and P/F ratio measurements from BigQuery via `bigrquery` for all
explicit extubation events, fits per-event slopes, and tests whether trajectory
predicts failure beyond level. Requires active Google Cloud authentication.

**RSBI:** 43,510 measurements across 5,046 events (median 6 per event).
Slope p = 0.634; level p = 0.032.

**P/F ratio:** 28,101 measurements across 5,327 events (median 4 per event).
Slope p = 0.726; level p = 0.505.

Produces `figures/preextubation_trajectories_by_event_outcome.png`.

---

## Data

### Source: MIMIC-IV v3.1

MIMIC-IV is a large de-identified EHR dataset from Beth Israel Deaconess
Medical Center, accessed via PhysioNet and Google BigQuery. Access requires
credentialing through PhysioNet and agreement to the MIMIC-IV Data Use
Agreement. **No MIMIC-IV data appears in this repository.**

### SQL Pipeline

`ICU_Last_Extubations.sql` is a BigQuery CTE chain producing one row per
patient (last extubation event). Key steps:

1. Combines explicit procedure events (itemids 224385, 227194, 225468, 225477)
   with inferred events from ventilation segments (≥4h InvasiveVent);
   deduplicates inferred events within 30 min of explicit events
2. Caregiver assignment cascade (explicit events only):
   recorded → other procedureevents ±30 min → chartevents ±15 min →
   inputevents ±5 min → mode caregiver ±2h window → NULL for inferred events
3. `caregiver_fe_rate` computed from explicit events only; NULL if no
   caregiver assignment
4. Joins: Charlson comorbidity index, norepinephrine equivalent dose, first-day
   SOFA, hospital mortality, primary ICD chapter, CCSR code arrays (seq_num
   ≤10), first_careunit from icustays
5. `caregiver_imputation_source` column tracks which tier assigned each
   caregiver_id

**Note on `failed_extubations`:** The `failed_extubations` column in the
`last_extubations` output is a patient lifetime count of reintubation pairs
across all extubation events during the hospital stay — not a flag for whether
the last extubation specifically failed. It is a patient complexity marker.
The per-event `failed_extubation_flag` (available in `all_extubations.sql`)
is the correct variable for extubation outcome analyses.

`all_extubations.sql` produces one row per explicit extubation event with
`failed_extubation_flag` and a leave-one-out caregiver FE rate. See the SQL
file header for details.

### Synthetic Cohort

The CSV files in `data/` are fully synthetic, generated by
`generate_synthetic_extubation_data.R` from published literature parameters.
Used only for the Tableau Public dashboard.

| Variable | Source | Value used |
|---|---|---|
| Age distribution | Johnson et al., *Scientific Data* 2016 | Median 65.8, IQR 52.8–77.8 |
| Failed extubation rate | Fernandez et al., *Respiratory Care* 2024 | 15.4% |
| Charlson comorbidity index | MIMIC-Sepsis, arXiv 2025 | Median 5 |
| SOFA score | Published MV cohorts | Mean ~8 |
| Ventilation hours | Published MV cohorts | Median ~96h |

---

## Reproducing the Analysis

Full reproduction requires MIMIC-IV access:

1. Apply at [physionet.org](https://physionet.org/content/mimiciv/3.1/)
2. Run `sql/ICU_Last_Extubations.sql` on BigQuery; save as
   `data/last_extubations.csv`
3. Run the CCSR long query (see `r/03_psm_and_sensitivity.R` header); save
   as `data/patient_ccsr_long.csv`
4. From the `r/` directory: `source("00_run_all.R")`

The SQL references a derived ventilation table and CCSR mapping table in a
private BigQuery project. The ventilation table follows the
[MIMIC-IV-derived](https://github.com/MIT-LCP/mimic-iv) pipeline. The CCSR
mapping uses the [AHRQ CCSR](https://hcup-us.ahrq.gov/toolssoftware/ccsr/ccs_refined.jsp)
crosswalk for ICD-10-CM.

### Reproducing the Synthetic Data

```r
install.packages(c("MASS", "dplyr", "tidyr", "purrr", "readr",
                   "tibble", "stringr", "forcats"))
source("r/generate_synthetic_extubation_data.R")
```

Fully deterministic given `set.seed(237)`.

---

## Dependencies

**R:** `dplyr`, `tidyr`, `readr`, `tibble`, `ggplot2`, `lubridate`, `splines`,
`tableone`, `glmmTMB`, `rstanarm`, `MatchIt`, `cobalt`, `WeightIt`, `ppcor`,
`MASS`, `bigrquery`, `patchwork`

**SQL:** Google BigQuery with MIMIC-IV v3.1 access

**Visualization:** Tableau Desktop 2020.1+

---

## Tableau Dashboard

Published at: **[ICU Extubation Outcomes Dashboard](https://public.tableau.com/app/profile/zachary.treisman/viz/ExtubationsSimulated/Dashboard1-CohortOverview)**

Three interactive dashboards built on the synthetic cohort:

**1. Patient Cohort Overview** — age distribution, intubation type breakdown,
ICD chapter volume, acuity scatter (SOFA vs. Charlson).

**2. Provider Volume & Outcomes** — caregiver-level scatter of volume vs. FE
rate and vs. 12-month survival, with trend lines and volume tier box plots.

**3. Risk Score Simulator** — parameter sliders drive a logistic regression
equation implemented as Tableau calculated fields; factor contribution chart
shows deviation from population-average patient.

> The simulator uses stylized coefficients consistent with published effect
> directions. It must not be used for clinical decision-making.

Full build instructions in `tableau/tableau_dashboard_spec.md`.

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

Code (SQL, R) is released under the MIT License.
Synthetic data files are released under CC0 (public domain).
No MIMIC-IV data is included; use of MIMIC-IV is governed by the
[PhysioNet Credentialed Health Data License](https://physionet.org/content/mimiciv/view-license/3.1/).
