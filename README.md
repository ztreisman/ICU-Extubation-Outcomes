# ICU Extubation Outcomes: Caregiver Failed Extubation Rate and Patient Mortality

A clinical data science project investigating whether caregiver-level failed
extubation rates predict patient 12-month mortality in the ICU, using MIMIC-IV
v3.1 and a multi-method analytical approach including marginal structural
modeling, mixed effects logistic regression, caregiver-level partial correlation,
and double/debiased machine learning.

---

## Research Question

Extubation — removing a patient from mechanical ventilation — is a high-stakes
clinical decision. Failed extubation (reintubation within 72 hours) is associated
with significantly higher mortality, longer ICU stays, and increased risk of
ventilator-associated pneumonia. This project asks: **does a caregiver's
historical failed extubation rate predict their patients' 12-month mortality,
after accounting for patient case mix, illness severity, and ICU unit?**

A secondary question motivated by clinical intuition: is there an optimal
FE rate range, where rates that are too low reflect overly conservative practice
and rates that are too high reflect excessive risk-taking?

---

## Key Findings

**Cohort:** 12,662 patients with a documented last ICU extubation event in
MIMIC-IV v3.1. Of these, 3,426 form the analytical cohort: explicitly documented
procedure events with a caregiver performing ≥20 unit-specific extubations
(`caregiver_unit_n >= 20`, 134 caregivers). `caregiver_fe_rate` is computed
within each (caregiver_id, first_careunit) pair (median 3.0%, IQR 0–6.1%),
benchmarking caregivers against peers extubating similar patient populations.
12-month mortality in the analytical cohort is ~30%. Likely comfort extubations
(in-hospital death within 3 days of extubation) account for approximately
4–5% of patients.

**Marginal structural model (MSM)** stratified by ICU unit group — CVICU
(post-cardiac surgery), Medical (MICU + MICU/SICU), Surgical/Trauma (SICU +
TSICU) — using inverse probability of treatment weighting with continuous
treatment and stabilized weights. OR per 5 percentage-point increase in FE rate:

| Group | N | Caregivers | OR per 5 pp | 95% CI | p |
|---|---|---|---|---|---|
| CVICU | 1,435 | 102 | **1.432** | 0.931–2.200 | 0.102 |
| Medical MICU+MICU/SICU | 1,001 | 52 | 0.898 | 0.750–1.074 | 0.239 |
| Surgical/Trauma SICU+TSICU | 809 | 53 | 0.893 | 0.627–1.270 | 0.529 |

The CVICU point estimate (OR 1.43) is directionally consistent with a harmful
effect of higher FE rate on mortality but does not reach conventional significance
(p = 0.102). Medical and Surgical/Trauma groups are null. Denominator R² (OLS)
is 0.057 / 0.042 / 0.069, confirming near-random assignment of patients to
caregivers within units after case-mix adjustment.

**CVICU comfort sensitivity:** Excluding comfort extubations, MSM OR shifts
from 1.432 to 1.398 (95% CI 0.914–2.138, p = 0.122). The finding is robust.

**Mixed effects logistic regression** with random intercepts for caregiver
and ICU unit:

- ICU unit random effect variance: 0.271 (glmmTMB), posterior mean 0.43
  (rstanarm) — ICU unit explains meaningful variation independently of patient
  severity
- Caregiver random effect variance: ~0. No residual caregiver-level variation
  once unit and patient covariates are controlled
- caregiver_fe_rate fixed effect: coef −0.554, p = 0.539 (glmmTMB); posterior
  mean −0.022, 95% CI (−0.102, 0.059) (rstanarm)

**Caregiver-level partial correlation** (134 caregivers, controlling for
caregiver volume):

- Global (N=134): **r = +0.216, p = 0.012** — higher FE rate associated with
  higher mortality
- **CVICU (N=79): r = +0.261, p = 0.021** — signal holds within the largest
  unit alone
- Volume ~ mortality controlling for FE rate: **r = +0.282, p = 0.001** —
  both FE rate and volume independently carry information
- Sensitivity (n_patients ≥ 10, N=70): r = +0.432, p < 0.001
- Sensitivity (n_patients ≥ 20, N=42): r = +0.485, p = 0.001

The signal strengthens substantially with more established caregivers,
arguing against a small-sample artifact.

**Functional form** of the FE rate — mortality relationship:

- Asymptotic exponential model: P(mortality) = m_max × (1 − exp(−k × FE rate))
  fit at the caregiver level (N = 90 caregivers with FE rate > 0)
- **m_max = 0.357 (95% CI 0.321–0.397)** — mortality ceiling. Approximately
  36% of patients die within 12 months at high FE rates; ~64% survive regardless.
- **k = 76.4 (95% CI 50.8–124.0)** — decay rate. The mortality penalty
  concentrates sharply in the 0–5% FE rate range.
- **CVICU unit-specific fit** (N=44 caregivers): m_max = 0.318 (0.263–0.398)
- **Sweet spot hypothesis unsupported.** Mortality increases monotonically with
  FE rate quartile within CVICU; there is no detectable floor on the benefit of
  lower FE rates.

**Deep learning pipeline** (Python, `python/`) trained on the analytical cohort
(N ≈ 3,437) with the same features and random seed as the R pipeline:

| Model | AUC-ROC | AUC-PR | Brier |
|---|---|---|---|
| Logistic regression | 0.829 | 0.663 | 0.149 |
| Gradient boosting | 0.852 | 0.674 | 0.141 |
| MLP with entity embeddings (PyTorch) | 0.841 | 0.670 | 0.166 |

Top SHAP features (mean |SHAP|): charlson, sofa, first_careunit, vent_hours,
primary_diagnosis, anchor_age, norepinephrine, caregiver_fe_rate (8th),
caregiver_n (9th). caregiver_fe_rate direction: positive (higher FE rate →
higher predicted mortality), consistent across R and Python.

**Double/Debiased ML** (Chernozhukov et al. 2018) using 5-fold cross-fitting
with gradient boosting nuisance models:

*Global* (all units pooled, N ≈ 3,437):
- θ̂ = −0.333, 95% CI [−1.010, +0.344], p = 0.335; IQR effect −1.4 pp

*Unit-stratified*:

| Group | N | θ̂ | 95% CI | p | IQR effect |
|---|---|---|---|---|---|
| CVICU | 1,435 | **+1.062** | −0.550 to +2.673 | 0.197 | +2.3 pp |
| Medical | 1,001 | −0.542 | −1.980 to +0.896 | 0.460 | −2.2 pp |
| Surgical/Trauma | 809 | −0.696 | −2.163 to +0.771 | 0.352 | −2.2 pp |

CVICU θ̂ direction (+1.06, higher FE rate → higher mortality) is consistent
with the MSM and partial correlation results.

**Physiological trajectory analysis** using all explicit extubation events
(6,063 events, 5,505 patients) with per-event reintubation within 72h as
outcome:

- RSBI and P/F ratio in the 72h before extubation do not predict per-event
  extubation failure beyond the level at the moment of extubation
- RSBI level at extubation: significant (p = 0.032), mean 52.3 (failed) vs
  48.7 (succeeded); slope not significant (p = 0.634)
- P/F ratio level and slope: not significant (p = 0.505, p = 0.726)
- Consistent with clinicians already conditioning on a plateau criterion before
  extubating

---

## Cohort Characteristics

### Cohort Flow

| Stage | N patients | N caregivers |
|---|---|---|
| patient_cohort (complete covariates) | 12,662 | 570 |
| explicit_extubations (caregiver_n ≥ 20) | 4,026 | 166 |
| analysis_base (caregiver_unit_n ≥ 20) | 3,426 | 134 |
| — CVICU | 1,435 | 102 |
| — Medical (MICU + MICU/SICU) | 1,001 | 52 |
| — Surgical/Trauma (SICU + TSICU) | 809 | 53 |
| — Other (descriptive only) | 181 | — |

The volume threshold of ≥20 unit-specific extubations is justified by the
false-zero problem: a caregiver with a true 5% FE rate has a 77% probability
of appearing 0% at n = 5; this drops to 36% at n = 20.

### Table 1: Patient Characteristics by Extubation Event Source

|  | Overall (N=12,662) | Explicit (N=4,836) | Inferred (N=7,826) | SMD |
|---|---|---|---|---|
| Age, mean (SD) | 62.7 (16.2) | 63.1 (15.9) | 62.5 (16.4) | 0.034 |
| Female, % | 40.7 | 38.5 | 42.1 | 0.073 |
| Charlson index, mean (SD) | 5.57 (3.25) | 5.45 (3.25) | 5.64 (3.26) | 0.059 |
| SOFA score, mean (SD) | 6.91 (3.86) | 6.39 (3.54) | 7.22 (4.02) | 0.220 |
| Norepinephrine equivalent, mean (SD) | 0.17 (0.23) | 0.14 (0.19) | 0.19 (0.24) | 0.242 |
| Ventilation hours, mean (SD) | 94.3 (134.1) | 64.6 (94.7) | 112.6 (150.6) | 0.382 |
| Intubation type (SMD 0.406) | | | | |
| — Surgical | 1.4% | 0.9% | 1.6% | |
| — Medical-respiratory | 74.7% | 64.4% | 81.1% | |
| — Medical-non-respiratory | 23.9% | 34.7% | 17.3% | |
| Primary diagnosis (SMD 0.379) | | | | |
| — Circulatory | 32.9% | 40.2% | 28.4% | |
| — Injury/Poisoning | 16.0% | 14.4% | 16.9% | |
| — Infectious | 14.8% | 10.0% | 17.7% | |
| — Digestive | 9.9% | 11.6% | 8.9% | |
| — Respiratory | 9.2% | 7.6% | 10.3% | |
| — Neoplasms | 5.8% | 5.9% | 5.8% | |
| — Other | 11.4% | 10.3% | 12.3% | |
| 12-month mortality, % | 45.2 | 32.3 | 53.2 | 0.432 |
| Hospital mortality, % | 27.3 | 12.4 | 36.5 | 0.585 |

Explicit events are enriched for circulatory diagnoses (40.2% vs 28.4%) and
have markedly lower 12-month mortality (32.3% vs 53.2%), reflecting the large
CVICU contribution to the explicit cohort. SOFA, vent hours, and intubation
type differ substantially (SMD > 0.2), justifying the analytical cohort
restriction.

### Unit-Level Summary

| ICU Unit | N (analysis_base) | N caregivers | Median FE rate | 12-mo mortality |
|---|---|---|---|---|
| Cardiac Vascular ICU (CVICU) | 1,435 | 102 | ~1.3% | ~14% |
| Medical ICU (MICU) | ~600 | ~40 | ~5–6% | ~50% |
| Medical/Surgical ICU (MICU/SICU) | ~400 | ~15 | ~5–7% | ~50% |
| Surgical ICU (SICU) | ~400 | ~30 | ~4–5% | ~30% |
| Trauma SICU (TSICU) | ~400 | ~25 | ~5–6% | ~30% |

Unit FE rates range from ~1% (CVICU) to ~7% (MICU), primarily reflecting
patient population rather than caregiver skill. `caregiver_fe_rate` is computed
within each unit, so CVICU caregivers are benchmarked only against other CVICU
caregivers.

---

## Data Quality and Cohort Restriction

MIMIC-IV extubation events are either explicitly documented procedure events
or algorithmically inferred from ventilation segment boundaries. Caregiver
identity is often missing or unreliable for inferred events. The project operates
on a broad descriptive `patient_cohort` (12,662 patients) and an `analysis_base`
(3,426 patients) restricted to explicitly documented events with a caregiver
performing ≥20 unit-specific extubations.

**Key data notes:**

- `caregiver_id` in MIMIC-IV identifies the documenting user, not necessarily
  the performing provider
- `dod` is only recorded for deaths within 1 year of hospital stay;
  `!is.na(dod)` correctly identifies 12-month mortality in both R and Python
- 8,068 patients are lost at "no computable vent duration" — a structural
  mismatch between procedureevents and the derived ventilation table, not a
  data quality error

---

## Project Structure

```
.
├── sql/
│   ├── ICU_Last_Extubations.sql        # BigQuery pipeline: MIMIC-IV → last extubation per patient
│   └── all_extubations.sql             # All explicit extubation events with per-event
│                                       #   failed_extubation_flag and leave-one-out FE rate
├── r/
│   ├── 00_run_all.R                    # Orchestration: sources 01–08 in order
│   ├── 01_cohort.R                     # Cohort construction, caregiver_summary, analysis_base;
│   │                                   #   saves cohorts.RData
│   ├── 02_characterization.R           # Patient/caregiver descriptives, CVICU strata
│   ├── 03_nls_partial_cor.R            # NLS mortality ceiling model, caregiver partial
│   │                                   #   correlations; saves nls_results.RData
│   ├── 04_mixed_effects.R              # glmmTMB and rstanarm mixed effects models;
│   │                                   #   saves mixed_effects_models.RData
│   ├── 05_msm.R                        # Continuous treatment MSM (primary causal analysis);
│   │                                   #   saves msm_results.RData
│   ├── 06_trajectory.R                 # Per-event RSBI and P/F trajectory analysis
│   │                                   #   (requires bigrquery and BigQuery access)
│   ├── 07_data_quality.R               # Data quality investigations
│   ├── 08_build_json.R                 # Aggregates RData files → analysis_outputs.json
│   ├── 09_consort.R                    # CONSORT diagram figure
│   └── generate_synthetic_extubation_data.R  # Synthetic cohort generator (Tableau)
├── python/
│   ├── 00_run_all.py                   # Master script: runs 02 → 03 → 04 → 05 → 06
│   ├── 02_baselines.py                 # Logistic regression and gradient boosting
│   ├── 03_train_mlp.py                 # MLP with entity embeddings (PyTorch)
│   ├── 04_shap_interpret.py            # SHAP feature importance (KernelExplainer)
│   ├── 05_dml_caregiver_fe.py          # Global double/debiased ML causal estimate
│   ├── 06_dml_cvicu.py                 # Unit-stratified DML (CVICU / Medical / Surg-Trauma)
│   ├── config.py                       # Shared hyperparameters and paths
│   ├── data_utils.py                   # Data loading, CCSR flags, feature encoding
│   ├── model.py                        # MLPWithEmbeddings architecture
│   ├── figures/                        # Output figures (gitignored, reproducible)
│   └── models/                         # Saved model checkpoints (gitignored)
├── data/
│   ├── analysis_outputs.json           # Summary statistics for document verification
│   ├── synthetic_last_extubations.csv  # 800-patient synthetic cohort (Tableau input)
│   └── synthetic_caregiver_rate.csv    # 80-caregiver summary (Tableau input)
├── figures/                            # All manuscript figures
│   ├── consort_diagram.png
│   ├── caregiver_vol_fe.png
│   ├── unit_sofa_mortality.png
│   ├── cvicu_strata_covariates.png
│   ├── msm_balance_cvicu.png
│   ├── msm_balance_medical__micu___micu_sicu_.png
│   ├── msm_balance_surgical_trauma__sicu___tsicu_.png
│   ├── msm_forest.png
│   ├── msm_secondary_forest.png
│   ├── fe_rate_mortality_exponential.png
│   ├── nls_stratified.png
│   └── mortality_by_quartile.png
├── doc/
│   ├── discussion_document.tex         # Working note (LaTeX source)
│   └── discussion_document.pdf         # Compiled discussion document
├── tableau/
│   └── tableau_dashboard_spec.md       # Full build specification for Tableau dashboards
└── README.md
```

**Note:** The `data/` directory contains only synthetic data and the
analysis_outputs.json summary. No MIMIC-IV patient records are included.

---

## Analytical Pipeline

### Script 01: Cohort Construction (`01_cohort.R`)

Constructs cohorts from `last_extubations.csv`:

**`patient_cohort`** (N = 12,662) — all patients with a documented last
extubation and complete covariates. Includes explicit and inferred events.

**`explicit_extubations`** — restricted to explicit procedure events with
`caregiver_n >= 20` (total volume across all units).

**`analysis_base`** — further restricted to `caregiver_unit_n >= 20`
(unit-specific volume). Primary analytical cohort: N = 3,426, 134 caregivers.

`caregiver_summary` is built from `analysis_base`. All downstream scripts
depend on `cohorts.RData`, which also includes `psm_ccsr_codes` and
`unit_groups`.

### Script 03: NLS and Partial Correlation (`03_nls_partial_cor.R`)

**Asymptotic exponential model:** P(mortality) = m_max × (1 − exp(−k × FE rate))
fit at the caregiver level (N = 90 caregivers with FE rate > 0). Bootstrap CIs
(2,000 resamples). Unit-stratified fits for CVICU (N=44) and TSICU (N=12).

**Partial correlations:** ppcor package, 134 caregivers, controlling for
caregiver volume. Global r = +0.216 (p = 0.012); CVICU r = +0.261 (p = 0.021).

Saves `nls_results.RData`.

### Script 04: Mixed Effects Models (`04_mixed_effects.R`)

**glmmTMB** (frequentist): random intercepts for caregiver_id and
first_careunit; centered continuous predictors. Mortality outcome.

**rstanarm** (Bayesian): same structure; scaled predictors, weakly informative
priors normal(0, 2.5) on fixed effects. Runtime ~10–30 minutes.

Saves `mixed_effects_models.RData`.

### Script 05: Marginal Structural Model (`05_msm.R`)

Continuous treatment MSM using stabilized inverse probability weights.
Denominator model: OLS of `caregiver_fe_rate` on case-mix covariates within
each unit group. Numerator: marginal distribution (caregiver-level mean and SD).
Weights truncated at the 1st and 99th percentiles. Outcome: weighted logistic
regression of 12-month mortality on `caregiver_fe_rate`.

Covariate balance assessed via weighted SMDs; balance plots saved to `figures/`.

Saves `msm_results.RData`.

### SQL: all_extubations.sql

One row per explicit extubation event for trajectory analysis:

- `failed_extubation_flag`: 1 if this extubation was followed by reintubation
  within 72h (both events explicit)
- `caregiver_fe_rate_loo`: leave-one-out caregiver FE rate — avoids circularity
- Output: 6,063 events, 5,505 patients, 110 caregivers, 13 ICU units,
  overall FE rate 8.3%

### Script 06: Trajectory Analysis (`06_trajectory.R`)

RSBI and P/F ratio measurements from BigQuery for all explicit extubation events.
Requires active Google Cloud authentication.

**RSBI:** 43,510 measurements across 5,046 events (median 6 per event).
Level p = 0.032; slope p = 0.634.

**P/F ratio:** 28,101 measurements across 5,327 events. Level and slope not
significant.

---

## Python Deep Learning Pipeline

Runs on the same analytical cohort (N ≈ 3,437) with the same random seed
(237), CCSR threshold (≥10%), and cohort filters as the R pipeline.
Target: `mortality_12mo` (`dod.notna()`).

### Scripts 02–03: Baseline and MLP Models

**Script 02** fits logistic regression and gradient boosting using 8 continuous
features, 4 categorical features, and 21 CCSR flags. Train/val/test split
70/15/15% stratified on mortality.

**Script 03** trains an MLP with entity embeddings (PyTorch): categorical
features passed through learned embedding layers. Architecture: three hidden
layers (256 → 128 → 64), dropout 0.30, Adam with ReduceLROnPlateau and early
stopping (patience = 20).

### Script 04: SHAP Interpretability

KernelExplainer on 200 held-out test patients with a 100-sample k-means
background set. Top feature: charlson; first_careunit and primary_diagnosis
together capture the unit and case-mix effects. caregiver_fe_rate ranks 8th,
quantitatively consistent with the glmmTMB finding that patient severity
dominates.

### Scripts 05–06: Double/Debiased ML

5-fold cross-fit DML estimating the average partial effect of `caregiver_fe_rate`
on P(12-month mortality). Nuisance models: GradientBoostingClassifier for
E[Y|W] and GradientBoostingRegressor for E[D|W].

**Script 05** (global, all units pooled): θ̂ = −0.333, p = 0.335.

**Script 06** (unit-stratified): CVICU θ̂ = +1.062 (p = 0.197), consistent
with the MSM direction. Medical and Surgical/Trauma null.

Run the full pipeline from the `python/` directory:

```bash
python 00_run_all.py                        # all steps
python 00_run_all.py --skip-shap            # skip SHAP (~5 min)
python 00_run_all.py --skip-shap --skip-dml # baselines + MLP only
```

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
3. `caregiver_fe_rate` computed within each `(caregiver_id, first_careunit)`
   pair from explicit events only. `caregiver_unit_n` is the unit-specific
   extubation count. Patients without a unit-specific rate are excluded from
   the analytical cohort via `!is.na(caregiver_unit_n)`.
4. Joins: Charlson comorbidity index, norepinephrine equivalent dose, first-day
   SOFA, hospital mortality, primary ICD chapter, CCSR code arrays (seq_num
   ≤10), first_careunit from icustays
5. `caregiver_imputation_source` column tracks which tier assigned each
   caregiver_id

**Note on `failed_extubations`:** The `failed_extubations` column in the
`last_extubations` output is a patient lifetime count of reintubation pairs
across all extubation events during the hospital stay — not a flag for whether
the last extubation specifically failed. It is a patient complexity marker.
The per-event `failed_extubation_flag` (in `all_extubations.sql`) is the
correct variable for extubation outcome analyses.

### Synthetic Cohort

The CSV files in `data/` are fully synthetic, generated by
`generate_synthetic_extubation_data.R` from published literature parameters.
Used only for the Tableau Public dashboard. Fully deterministic given
`set.seed(237)`.

---

## Reproducing the Analysis

Full reproduction requires MIMIC-IV access:

1. Apply at [physionet.org](https://physionet.org/content/mimiciv/3.1/)
2. Run `sql/ICU_Last_Extubations.sql` on BigQuery; save as
   `data/last_extubations.csv`
3. Run the CCSR long query (see `r/01_cohort.R` header); save as
   `data/patient_ccsr_long.csv`
4. From the `r/` directory: `source("00_run_all.R")`

The SQL references a derived ventilation table and CCSR mapping table in a
private BigQuery project. The ventilation table follows the
[MIMIC-IV-derived](https://github.com/MIT-LCP/mimic-iv) pipeline. The CCSR
mapping uses the [AHRQ CCSR](https://hcup-us.ahrq.gov/toolssoftware/ccsr/ccs_refined.jsp)
crosswalk for ICD-10-CM.

---

## Dependencies

**R:** `dplyr`, `tidyr`, `readr`, `tibble`, `ggplot2`, `lubridate`, `splines`,
`tableone`, `glmmTMB`, `rstanarm`, `WeightIt`, `cobalt`, `ppcor`, `MASS`,
`bigrquery`, `patchwork`, `DiagrammeR`

**Python:** `torch`, `numpy>=1.26,<2`, `pandas>=2.1`, `scikit-learn`, `scipy`,
`shap`, `matplotlib`

**SQL:** Google BigQuery with MIMIC-IV v3.1 access

**Visualization:** Tableau Desktop 2020.1+

---

## Tableau Dashboard

Published at: **[ICU Extubation Outcomes Dashboard](https://public.tableau.com/app/profile/zachary.treisman/viz/ExtubationsSimulated/Dashboard1-CohortOverview)**

Three interactive dashboards built on the synthetic cohort:

**1. Patient Cohort Overview** — age distribution, intubation type breakdown,
ICD chapter volume, acuity scatter (SOFA vs. Charlson).

**2. Provider Volume & Outcomes** — caregiver-level scatter of volume vs. FE
rate and vs. 12-month mortality, with trend lines and volume tier box plots.

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
with Andrew Hersh, MD, who contributed clinical domain expertise and motivated
the caregiver-level research questions.

Johnson, A., Bulgarelli, L., Shen, L. et al. MIMIC-IV, a freely accessible
electronic health record dataset. *Sci Data* 10, 1 (2023).
https://doi.org/10.1038/s41597-022-01899-x

---

## License

Code (SQL, R, Python) is released under the MIT License.
Synthetic data files are released under CC0 (public domain).
No MIMIC-IV data is included; use of MIMIC-IV is governed by the
[PhysioNet Credentialed Health Data License](https://physionet.org/content/mimiciv/view-license/3.1/).
