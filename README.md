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
and rates that are too high reflect excessive desire to extubate?

---

## Key Findings

**Cohort:** 12,662 patients with a documented last ICU extubation event in
MIMIC-IV v3.1 — one of the largest published ICU extubation cohorts derived
from this dataset. Of these, 4,387 had explicitly documented procedure events
with a genuine caregiver assignment and form the analytical cohort (271
caregivers, 12 ICU units). `caregiver_fe_rate` is computed within each ICU
unit (median 3.0%, IQR 0–6.1%), so caregivers are benchmarked against peers
extubating similar patient populations rather than across-unit comparisons.
30-day mortality in the analytical cohort is 11.7%. Likely comfort extubations
(in-hospital death within 3 days of extubation) account for 198 patients (4.5%).

**Propensity score matching** stratified by ICU unit group — CVICU
(post-cardiac surgery), Medical (MICU + MICU/SICU), Surgical/Trauma
(SICU + TSICU) — comparing top vs bottom quartile FE rate caregivers within
each group (caregivers with ≥5 unit-specific extubations), matched on illness
severity, intubation type, and 21 CCSR diagnosis flags. **The FE rate →
survival association is specific to CVICU and absent elsewhere:**

Primary outcome (12-month survival):

| Group | Survival high vs low FE | PSM OR (adj.) | IPW OR |
|---|---|---|---|
| CVICU (N=797 matched) | 84.4% vs 90.1% | **0.581 (0.356–0.937)**, p=0.027 | **0.656 (0.487–0.879)**, p=0.005 |
| Medical MICU+MICU/SICU (N=607) | 49.4% vs 49.2% | 1.091 (0.754–1.580), p=0.645 | 1.157 (0.930–1.439), p=0.190 |
| Surgical/Trauma SICU+TSICU (N=524) | 72.9% vs 68.2% | 1.202 (0.766–1.891), p=0.423 | 1.171 (0.907–1.511), p=0.226 |

Secondary outcomes (PSM-matched, unadjusted OR for hospital mortality and 30-day mortality):

| Group | Hospital mortality OR | 30-day mortality OR |
|---|---|---|
| CVICU | 1.89 (0.886–4.28), p=0.11 | 1.69 (0.672–4.58), p=0.28 |
| Medical | 1.09 (0.734–1.61), p=0.68 | 1.13 (0.764–1.67), p=0.54 |
| Surgical/Trauma | **0.557 (0.321–0.949)**, p=0.034 | **0.546 (0.304–0.958)**, p=0.038 |

The earlier cross-unit PSM result (OR = 0.646) was driven by CVICU patient
composition. In Medical ICUs, all outcomes are null. In Surgical/Trauma ICUs,
12-month survival is null but short-term mortality endpoints are reversed —
high-FE-rate caregivers are associated with lower hospital and 30-day mortality
(OR ≈ 0.55). One interpretation: in surgical populations, patients who fail
extubation are reintubated and bridge to hospital discharge (lower acute
mortality), while long-term prognosis is governed by the underlying condition
rather than extubation management. In CVICU, the short-term mortality
endpoints are null (hospital OR=1.89, p=0.11; 30-day OR=1.69, p=0.28) despite
the significant 12-month survival signal — consistent with post-cardiac-surgery
failed extubations causing prolonged ICU stays and delayed downstream mortality
rather than early in-hospital death.

**Comfort extubation sensitivity (CVICU):** Excluding the 11 CVICU patients
likely extubated as part of comfort/withdrawal care (in-hospital death ≤3 days,
0.7% of CVICU cases), PSM OR shifts from 0.581 to 0.619 (p=0.055) and IPW OR
from 0.656 to 0.673 (p=0.010). The finding is robust to comfort care exclusion.

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

- Global (134 caregivers): r = −0.288, p < 0.001 — but this reflects CVICU's
  dominance in the caregiver sample, not a universal effect
- **CVICU (73 caregivers): r = −0.342, p = 0.003** — the primary signal;
  holds within the single largest unit, controlling for volume
- MICU (30 caregivers): r = −0.298, p = 0.116 — directionally consistent,
  underpowered
- SICU (8 caregivers): r = +0.646, p = 0.117 — reversed direction; N too
  small to interpret, consistent with the Surgical/Trauma PSM null
- Caregiver volume is not an independent predictor (r = −0.109, p = 0.211)
- **Volume threshold sensitivity:** restricting to caregivers with ≥10 sample
  patients (N=88): r = −0.255, p = 0.017; ≥20 patients (N=52): r = −0.354,
  p = 0.011. The signal strengthens with more established caregivers, arguing
  against a sample-size artifact.

**Functional form** of the FE rate — survival relationship (caregivers with
FE rate > 0, N = 111):

- Asymptotic exponential model: P(survival) = b + (1−b) × exp(−k × FE rate)
  (caregivers with FE rate > 0, N = 78)
- **b = 0.624 (95% CI 0.586–0.658)** — survival floor. Approximately 62% of
  patients survive 12 months regardless of caregiver FE rate, representing
  irreducible mortality driven by underlying illness severity.
- **k = 77.5 (95% CI 52.7–122.4)** — decay rate. The survival penalty
  concentrates sharply in the 0–5% FE rate range.
- The remaining ~38% of patients represent outcomes potentially sensitive to
  caregiver extubation practice — the upper bound on what quality improvement
  could achieve in this population.
- **Sweet spot hypothesis unsupported.** Within-group quartile analysis shows
  Q1 (lowest FE rate) has the best outcomes in CVICU (90.9% vs 84.4% in Q4),
  declining monotonically through Q3. Within CVICU, lower FE rate is better with 
  no detectable floor.

**Deep learning pipeline** (Python, `python/`) trained on the same 4,387-patient
explicit cohort with the same features and random seed as the R pipeline:

- Logistic regression: AUC-ROC 0.806, AUC-PR 0.935, Brier 0.141
- Gradient boosting: AUC-ROC 0.828, AUC-PR 0.943, Brier 0.132
- MLP with entity embeddings (PyTorch): AUC-ROC 0.841, AUC-PR 0.948, Brier 0.171

Top SHAP features (mean |SHAP|): charlson (0.119), first_careunit (0.069),
primary_diagnosis (0.042), sofa (0.038), anchor_age (0.035).
caregiver_fe_rate ranks 10th — large predictive improvement from patient severity
features, small independent contribution from caregiver, consistent with the R results.

**Double/Debiased ML** (Chernozhukov et al. 2018) for the causal effect of
caregiver_fe_rate on survival, using 5-fold cross-fitting with gradient boosting
nuisance models:

- θ̂ = −0.172 (change in P(survival) per unit of FE rate), SE = 0.311, p = 0.581
- IQR survival difference (FE rate 2.9% → 5.3%): −0.4 pp
- Approximate log-odds equivalent: −0.977 (vs. glmmTMB fixed effect −2.80, p = 0.097)
- Treatment R² ~0.40–0.46: confounders explain ~40% of variation in FE rate; DML
  uses the remaining unexplained variation as the identifying signal

DML is more conservative than glmmTMB, likely because the flexible GBM nuisance
model absorbs nonlinear confounding and because DML does not incorporate caregiver
random effects or unit stratification. The group-stratified PSM subsequently
showed the signal is confined to CVICU; DML's cross-unit null is consistent with
that finding — the CVICU signal is diluted by the null Medical and Surgical/Trauma
populations in a pooled causal estimate.

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

## Cohort Characteristics

### Cohort Flow

| Step | N patients |
|---|---|
| Raw SQL output (MIMIC-IV v3.1) | 12,668 |
| patient_cohort (complete covariates) | 12,662 |
| — of which explicit procedure events | 4,836 |
| — of which algorithmically inferred events | 7,826 |
| explicit_extubations (analytical cohort) | 4,387 |
| — Unique caregivers | 271 |
| — Unique ICU units | 12 |

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
| Failed extubations, mean (SD) | 0.10 (0.34) | 0.09 (0.32) | 0.10 (0.34) | 0.026 |
| 12-month survival, % | 54.8 | 67.7 | 46.8 | 0.432 |
| Hospital mortality, % | 27.3 | 12.4 | 36.5 | 0.585 |

Explicit events are enriched for circulatory diagnoses (40.2% vs 28.4%) and
have markedly better survival (67.7% vs 46.8%), reflecting the large CVICU
contribution to the explicit cohort (see unit breakdown below). SOFA, vent
hours, and intubation type differ substantially (SMD > 0.2), justifying the
analytical cohort restriction.

### Table 1b: Analytical Cohort by Caregiver FE Rate Quartile

FE rate quartile thresholds (unit-specific rates): Q1 ≤ 0%, Q2 0–3.2%,
Q3 3.2–5.3%, Q4 ≥ 5.3%.

|  | Q1 lowest (N=1,103) | Q2 (N=1,143) | Q3 (N=1,046) | Q4 highest (N=1,095) | p |
|---|---|---|---|---|---|
| Age, mean (SD) | 64.1 (15.7) | 65.5 (13.6) | 62.0 (15.8) | 61.0 (17.3) | <0.001 |
| Female, % | 36.3 | 35.8 | 40.1 | 39.8 | 0.066 |
| Charlson, mean (SD) | 5.49 (3.23) | 5.31 (3.02) | 5.48 (3.32) | 5.49 (3.40) | 0.470 |
| SOFA, mean (SD) | 6.40 (3.66) | 6.07 (3.19) | 6.47 (3.59) | 6.75 (3.66) | <0.001 |
| Vent hours, mean (SD) | 56.2 (92.9) | 46.5 (81.1) | 77.1 (107.3) | 79.9 (94.5) | <0.001 |
| Primary diagnosis (SMD 0.705) | | | | | <0.001 |
| — Circulatory | 48.9% | 68.2% | 29.3% | 17.7% | |
| — Infectious | 8.5% | 5.2% | 9.8% | 15.2% | |
| — Injury/Poisoning | 13.4% | 8.3% | 18.3% | 18.5% | |
| — Respiratory | 6.3% | 3.6% | 8.4% | 11.1% | |
| — Digestive | 8.7% | 4.1% | 15.1% | 17.2% | |
| 12-month survival, % | 67.9 | **78.2** | 64.5 | 60.3 | <0.001 |
| Hospital mortality, % | 13.5 | 6.3 | 14.0 | 15.3 | <0.001 |

The large diagnosis SMD (0.705) across quartiles reflects CVICU domination
of the low-FE-rate groups (circulatory: 48.9% in Q1, 68.2% in Q2) versus
more medical/trauma mix in Q3–Q4. Q2 shows higher survival than Q1 (78.2%
vs 67.9%), and the quadratic functional form model is now significant
(p = 0.037), though this pattern may partly reflect case-mix differences
between units.

### Unit-Level Summary (Analytical Cohort)

| ICU Unit | N patients | N caregivers | Mean FE rate | Median FE rate | 12-mo survival | Mean SOFA | Mean Charlson |
|---|---|---|---|---|---|---|---|
| Cardiac Vascular ICU (CVICU) | 1,543 | 148 | 1.7% | 1.5% | 86.5% | — | — |
| Coronary Care Unit (CCU) | 268 | 72 | 2.9% | 0.0% | 43.3% | — | — |
| Medical/Surgical ICU | 465 | 90 | 4.6% | 3.2% | 44.5% | — | — |
| Surgical ICU (SICU) | 573 | 116 | 5.0% | 4.3% | 66.7% | — | — |
| Neuro Surgical ICU | 53 | 36 | 5.2% | 0.0% | 66.0% | — | — |
| Medical ICU (MICU) | 850 | 129 | 5.5% | 6.2% | 52.0% | — | — |
| Trauma SICU (TSICU) | 554 | 118 | 5.6% | 5.3% | 72.9% | — | — |
| Neuro Intermediate | 51 | 31 | 8.6% | 0.0% | 68.6% | — | — |
| Neuro Stepdown | 25 | 20 | 8.7% | 0.0% | 72.0% | — | — |
| Other (≤3 patients each) | 5 | 5 | — | — | — | — | — |

Unit FE rates range from 1.7% (CVICU, post-cardiac-surgery patients) to 8.7%
(Neuro Stepdown), a five-fold difference driven primarily by patient population
rather than caregiver skill. `caregiver_fe_rate` is now computed within each
unit (see SQL pipeline), so CVICU caregivers are compared only to other CVICU
caregivers, not against MICU caregivers extubating fundamentally different patients.

Natural analytical groupings suggested by this breakdown:
- **Cardiac** (CVICU + CCU): post-cardiac-surgery and acute coronary patients; low-FE, high-acuity split
- **Medical/Mixed** (MICU + MICU/SICU): general medical ICU patients; highest FE rates
- **Surgical/Trauma** (SICU + TSICU): surgical and trauma patients; intermediate FE rates, better survival
- **Neuro** (Neuro Surgical ICU + Neuro Intermediate + Neuro Stepdown): small N, borderline for stratified analysis

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
the restriction produces a clinically plausible within-unit caregiver FE rate
distribution (median 3.0%, IQR 0–6.1%). Improving caregiver assignment for
inferred events remains an open methodological problem.

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
├── python/
│   ├── 00_run_all.py                   # Master script: runs 02 → 03 → 04 → 05
│   ├── 02_baselines.py                 # Logistic regression and gradient boosting
│   ├── 03_train_mlp.py                 # MLP with entity embeddings (PyTorch)
│   ├── 04_shap_interpret.py            # SHAP feature importance (KernelExplainer)
│   ├── 05_dml_caregiver_fe.py          # Double/Debiased ML causal estimate
│   ├── config.py                       # Shared hyperparameters and paths
│   ├── data_utils.py                   # Data loading, CCSR flags, feature encoding
│   ├── model.py                        # MLPWithEmbeddings architecture
│   ├── figures/                        # Output figures
│   └── models/                         # Saved model checkpoints
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

**PSM:** Group-stratified — CVICU, Medical (MICU + MICU/SICU), and
Surgical/Trauma (SICU + TSICU). Within each group: method = "quick", ratio
1:1, caliper 0.2 SD on logit propensity score; exact matching within subunit
for multi-unit groups. Treatment defined as top vs bottom quartile FE rate
within the group, restricted to caregivers with ≥5 unit-specific extubations
(`caregiver_unit_n >= 5`). CCSR codes with zero variance within a group are
dropped from the PS model. Outcome model on matched data without MatchIt
subclass weights (subclass weights re-introduce stratum-size imbalance and
invert the estimate in some settings).

**IPW:** WeightIt package, ATE estimand, run within each group alongside PSM.

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

## Python Deep Learning Pipeline

Runs on the same explicit cohort (N = 4,387) and uses the same random seed
(237), CCSR threshold (≥10%), and cohort filters as the R pipeline to keep
results directly comparable.

### Scripts 02–03: Baseline and MLP Models

**Script 02** fits logistic regression and gradient boosting classifiers using
the same 8 continuous features, 4 categorical features, and 21 CCSR flags.
Train/val/test split is 70/15/15% stratified on `survival_12mo`.

**Script 03** trains an MLP with entity embeddings (PyTorch): categorical
features (gender, intubation type, primary diagnosis, first careunit) are
passed through learned embedding layers rather than one-hot encoded, allowing
the model to capture similarity structure across categories. Architecture:
three hidden layers (256 → 128 → 64), dropout 0.30, Adam optimizer with
ReduceLROnPlateau scheduling and early stopping (patience = 20).

### Script 04: SHAP Interpretability

KernelExplainer on 200 held-out test patients with a 100-sample k-means
background set. Produces summary, bar, and dependence plots for
`caregiver_fe_rate`. Top feature by mean |SHAP|: charlson (0.119); first
careunit (0.069) and primary diagnosis (0.042) together show the ICU unit
and case-mix effects the R random effects model captures. caregiver_fe_rate
ranks 10th (0.007) — quantitatively consistent with the glmmTMB finding that
patient severity dominates.

### Script 05: Double/Debiased ML

5-fold cross-fit DML estimating the average partial effect of caregiver_fe_rate
on P(12-month survival). Nuisance models: GradientBoostingClassifier for
E[Y|W] and GradientBoostingRegressor for E[D|W], where W is all confounders
excluding FE rate.

Key diagnostics:
- Outcome nuisance AUC per fold: 0.806–0.828 (nuisance model adequate)
- Treatment nuisance R² per fold: 0.36–0.46 (substantial unexplained variation
  in FE rate remains as the identifying signal)

Result: θ̂ = −0.172, SE = 0.311, p = 0.581. IQR effect (FE rate 2.9% → 5.3%):
−0.4 pp. The flexible nonparametric confounder model does not find a significant
residual effect, consistent with the overall finding that the caregiver signal
is real but small and sensitive to modeling assumptions.

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
   pair from explicit events only; falls back to global caregiver rate, then
   overall mean. `caregiver_unit_n` is the unit-specific extubation count.
   Patients without a unit-specific rate (cross-unit imputation artifacts)
   are excluded from the analytical cohort via `!is.na(caregiver_unit_n)`.
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
