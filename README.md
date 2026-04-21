

# ICU Extubation Outcomes: Does Provider Volume Predict Patient Survival?

A clinical data science project investigating the relationship between
caregiver extubation volume and patient outcomes in the ICU, built on
MIMIC-IV and published as an interactive Tableau Public dashboard.

------------------------------------------------------------------------

## Research Question

Extubation — removing a patient from mechanical ventilation — is a
high-stakes clinical decision. Failed extubation (reintubation within 72
hours) is associated with significantly higher mortality, longer ICU
stays, and increased risk of ventilator-associated pneumonia. The
decision of *when* and *how* to extubate involves substantial provider
judgment.

This project asks: **does the volume of extubations a caregiver has
performed predict their patients' outcomes?** A volume-outcome
relationship is well documented in surgical procedures, but less studied
in ventilator weaning. We also ask whether inter-provider variation in
failed extubation rates can be explained by differences in patient case
mix — and whether the survival effect of caregiver FE rate persists
after adjusting for measured confounders.

------------------------------------------------------------------------

## Key Findings

**Descriptive:** - Low-volume caregivers show **higher and more
variable** failed extubation rates than high-volume caregivers,
consistent with a learning curve effect. - The volume-outcome signal
persists after adjusting for intubation type, illness severity (SOFA,
Charlson), and vasopressor use in logistic regression.

**Case mix analysis:** - We explored whether inter-provider variation in
FE rates could be explained by differences in patient diagnosis
profiles. LASSO feature selection across 12 ICD chapters retained
minimal signal. LDA topic modeling on 388 CCSR diagnosis codes produced
8 clinically coherent condition profiles, but cluster separation was
weak (maximum variance explained \~15%), reflecting the multimorbid,
overlapping nature of ICU presentations. **Case mix does not explain
inter-provider FE rate variation.**

**Propensity score analysis:** - To assess whether the caregiver FE rate
effect on 12-month survival survives adjustment for measured
confounders, we conducted propensity score matching (PSM) comparing
patients of high-FE-rate caregivers (top quartile, FE rate ≥ 32%) to
patients of low-FE-rate caregivers (bottom quartile, FE rate ≤ 12%). -
After matching on age, illness severity (SOFA, Charlson, norepinephrine,
ventilation hours), intubation type, and 21 high-prevalence CCSR
diagnosis flags, **2,353 matched pairs** achieved excellent covariate
balance (all standardized mean differences \< 0.10 post-match). -
**12-month survival: 50.8% (high-FE caregivers) vs. 56.8% (low-FE
caregivers)** — a 6 percentage point absolute difference. - PSM odds
ratio (doubly adjusted): **OR = 0.733 (95% CI: 0.643–0.836)** - IPW
sensitivity analysis: **OR = 0.789 (95% CI: 0.735–0.846)** - No effect
modification by intubation type (interaction p = 0.623) - The caregiver
FE rate effect on 12-month survival **is not explained by the patient
factors we can measure**.

------------------------------------------------------------------------

## Limitations

This is an observational study and causal claims are not warranted.

**Unmeasured confounding:** PSM controls only for observed covariates.
Shift patterns, time of day of extubation, ICU unit type, attending
physician supervision, and institutional protocols are not captured in
MIMIC-IV and may confound the caregiver-outcome relationship.

**Treatment definition:** `caregiver_fe_rate` is computed from the full
patient cohort and is itself an outcome-adjacent measure. Caregivers
with high FE rates may have high rates because their patients are
systematically harder to extubate in ways not fully captured by
available covariates.

**Outcome attribution:** `survival_12mo` is derived from date of death
in MIMIC-IV regardless of cause. Deaths unrelated to the index
extubation are counted equally, which attenuates the estimated caregiver
effect on extubation-attributable mortality.

**Survival rates:** 12-month survival of 50–57% in the matched cohort is
lower than general MIMIC population figures, reflecting the restriction
to ventilated patients with sufficient caregiver volume data
(`caregiver_n > 10`) and the top/bottom quartile selection.

**CCSR clustering:** LDA topic modeling on diagnosis codes was initially
explored as a semi-supervised approach to identify patient condition
profiles. The short-document problem (median \~4 CCSR codes per patient)
limited topic differentiation. We pivoted to using high-prevalence CCSR
binary flags directly as PSM covariates, which is more transparent and
interpretable.

------------------------------------------------------------------------

## Project Structure

```         
.
├── sql/
│   └── ICU_Last_Extubations.sql            # BigQuery pipeline: MIMIC-IV → analysis dataset
├── r/
│   ├── MIMIC_analysis_clean.R              # Primary analysis: cohort construction,
│   │                                       #   logistic regression, caregiver analysis,
│   │                                       #   LASSO, random forest, 5-fold CV
│   ├── psm_analysis.R                      # Propensity score analysis: LDA topic modeling,
│   │                                       #   CCSR binary matrix, PSM, IPW sensitivity,
│   │                                       #   covariate balance diagnostics
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

**Note:** The `data/` directory contains only synthetic data (see
below). No MIMIC-IV patient records are included in this repository.

------------------------------------------------------------------------

## Data

### Source: MIMIC-IV v3.1

The analysis was conducted on MIMIC-IV, a large de-identified EHR
dataset from the Beth Israel Deaconess Medical Center ICU, accessed via
PhysioNet and Google BigQuery. MIMIC-IV contains data for over 65,000
ICU patients admitted between 2008 and 2022.

Access requires credentialing through PhysioNet and agreement to the
MIMIC-IV Data Use Agreement. The DUA prohibits redistribution of the
data or derivatives that could be used to re-identify patients. **No
MIMIC-IV data appears in this repository.**

### Synthetic cohort

The CSV files in `data/` are fully synthetic. They were generated by
`generate_synthetic_extubation_data.R` using distributional parameters
drawn exclusively from published literature — no MIMIC records were used
as input.

The generator preserves: - Realistic marginal distributions for each
variable - Clinically plausible correlations among predictors - Outcome
rates consistent with published ICU extubation cohorts - The
caregiver-level structure (volume, FE rate) central to the analysis

Key parameter sources:

| Variable | Source | Value used |
|------------------------|------------------------|------------------------|
| Age distribution | Johnson et al., *Scientific Data* 2016 (MIMIC-III) | Median 65.8, IQR 52.8–77.8 |
| Failed extubation rate | Fernandez et al., *Respiratory Care* 2024 | 15.4% (72h window) |
| Hospital mortality | Johnson et al. 2016 | \~11.5% |
| Charlson comorbidity index | MIMIC-Sepsis, arXiv 2025 | Median 5 |
| SOFA score | Published MV cohorts | Mean \~8 |
| Ventilation hours | Published MV cohorts | Median \~96h |
| Vasopressor use | Published MV cohorts | \~30% of patients |

The synthetic data is appropriate for portfolio demonstration and public
visualization. It is not appropriate for clinical inference.

### SQL pipeline

`ICU_Last_Extubations.sql` is a single BigQuery CTE chain that:

1.  Combines explicit intubation/extubation procedure events (itemids
    224385, 227194, 225468, 225477) with inferred events from
    ventilation segments (≥4h InvasiveVent), deduplicating inferred
    events within 30 minutes of explicit events
2.  Classifies each intubation as `surgical`, `medical-respiratory`, or
    `medical-non-respiratory` based on proximate surgical procedure
    codes and respiratory ICD chapter flags
3.  Selects the *last* extubation per patient (one row per patient)
4.  Joins patient demographics, Charlson comorbidity index,
    norepinephrine equivalent dose, first-day SOFA score, hospital
    mortality, primary ICD chapter, and an array of CCSR categories (up
    to seq_num 10)
5.  Assigns caregiver ID via a three-level fallback: recorded
    procedureevents caregiver → nearest chartevents caregiver within ±15
    minutes → most frequent caregiver for that ICU stay
6.  Computes `caregiver_fe_rate` and `caregiver_n` at query time
7.  Applies final filters: vent_hours ∈ (1, 1000), norepinephrine \< 1
    mcg/kg/min, failed_extubations \< 10

------------------------------------------------------------------------

## Analysis

### Primary analysis (`MIMIC_analysis_clean.R`)

| Section | Content |
|------------------------------------|------------------------------------|
| 0 | Libraries and setup |
| 1 | Data ingestion (SQL output CSVs → R) |
| 2 | Cohort construction: filtering, centering, `cleandata` / `supercleandata` |
| 3 | Patient-level logistic regression for 12-month survival (three model variants) |
| 4 | Caregiver-level analysis: volume vs. FE rate, volume vs. survival, volume vs. hospital mortality |
| 5 | Diagnosis case mix: LASSO on wide ICD chapter × caregiver matrix; ICD indicator construction |
| 6 | Full interaction model (`logreg_icd`) with 5-fold cross-validation; classification metrics |
| 7 | Random forest on caregiver case mix; variable importance (%IncMSE, IncNodePurity) |
| 8 | RSBI-based logistic regression (separate cohort with pre-extubation vitals) |
| 9 | Appendix: patient timeline visualization (data quality checking) |

### Propensity score analysis (`psm_analysis.R`)

| Section | Content |
|------------------------------------|------------------------------------|
| 0 | Libraries and setup |
| 1 | Data ingestion: `last_extubations_clean.csv`, `patient_ccsr_long.csv` |
| 2 | CCSR document-term matrix construction |
| 3 | LDA topic modeling (K = 5, 8, 10, 12); perplexity comparison; topic interpretation |
| 4a | Topic proportion extraction (soft cluster membership) |
| 4b | K-means clustering on binary CCSR matrix; variance explained comparison |
| 4c | Revised covariate strategy: high-prevalence CCSR binary flags (≥10% prevalence) |
| 5 | PSM analysis dataset construction; treatment variable definition (top vs. bottom quartile) |
| 6 | Propensity score matching: nearest neighbor 1:1, caliper = 0.2 SD |
| 7 | Treatment effect estimation: unadjusted and doubly-adjusted outcome models |
| 8 | Sensitivity and validity checks: IPW, covariate balance love plot, effect modification test |
| 9 | Output: results summary, matched dataset, love plot |

### PSM results summary

| Analysis              | OR    | 95% CI      | N     |
|-----------------------|-------|-------------|-------|
| PSM — unadjusted      | 0.786 | 0.700–0.881 | 4,706 |
| PSM — doubly adjusted | 0.733 | 0.643–0.836 | 4,706 |
| IPW — full sample     | 0.789 | 0.735–0.846 | 6,191 |

All p-values \< 0.001. Covariate balance plot:
`figures/psm_love_plot.png`.

------------------------------------------------------------------------

## Tableau Dashboard

The interactive dashboard is published at: **[Tableau Public link — add
after publishing]**

Three dashboards:

**1. Patient Cohort Overview** Age distribution, intubation type
breakdown, ICD chapter volume and FE rates, and an acuity scatter (SOFA
vs. Charlson). Establishes the population before the modeling story
begins.

**2. Provider Volume & Outcomes** Caregiver-level scatter of volume vs.
failed extubation rate (polynomial trend, confidence bands) and volume
vs. 12-month survival rate. Box plots showing FE rate spread within
volume tiers surface the narrowing variance at high volume — the visual
center of the project.

**3. Risk Score Simulator** Parameter sliders for Charlson index, SOFA
score, vasopressor dose, ventilation hours, intubation type, and
caregiver FE rate drive a logistic regression equation implemented as
Tableau calculated fields. The predicted 12-month survival probability
updates in real time. A factor contribution chart shows how each input
deviates from the population-average patient.

> The simulator uses stylized logistic coefficients consistent with
> published effect directions. It does not reproduce fitted coefficients
> from the MIMIC analysis and must not be used for clinical
> decision-making.

Full build instructions, calculated field formulas, layout
specifications, and parameter configurations are in
`tableau/tableau_dashboard_spec.md`.

------------------------------------------------------------------------

## Reproducing the Synthetic Data

``` r
# Install dependencies
install.packages(c("MASS", "dplyr", "tidyr", "purrr", "readr",
                   "tibble", "stringr", "forcats"))

# Run generator (set working directory to r/)
source("generate_synthetic_extubation_data.R")

# Outputs written to working directory:
#   synthetic_last_extubations.csv   (800 rows)
#   synthetic_caregiver_rate.csv     (80 rows)
```

The generator is fully deterministic given `set.seed(237)`. Validation
output confirms simulated rates against published targets before writing
files.

------------------------------------------------------------------------

## Reproducing the Full Analysis

Full reproduction requires MIMIC-IV access:

1.  Apply for access at
    [physionet.org](https://physionet.org/content/mimiciv/3.1/)
2.  Run `sql/ICU_Last_Extubations.sql` on BigQuery
    (`physionet-data.mimiciv_3_1_*` datasets)
3.  Export query result to `r/../data/last_extubations_clean.csv`
4.  Run the CCSR long query (see `psm_analysis.R` header) and export to
    `r/../data/patient_ccsr_long.csv`
5.  Run `r/MIMIC_analysis_clean.R` from the `r/` directory
6.  Run `r/psm_analysis.R` from the `r/` directory

The SQL query also references a derived ventilation table
(`Failed_Extubation_Rate.ventilation`) and a CCSR mapping table
(`Failed_Extubation_Rate.ccsr_mappings`) stored in a private BigQuery
project. The ventilation table is derived from MIMIC-IV's
`mimiciv_derived` schema following the
[MIMIC-IV-derived](https://github.com/MIT-LCP/mimic-iv) pipeline. The
CCSR mapping table uses the publicly available [AHRQ
CCSR](https://hcup-us.ahrq.gov/toolssoftware/ccsr/ccs_refined.jsp)
crosswalk for ICD-10-CM codes.

------------------------------------------------------------------------

## Dependencies

**R packages:** `tidyverse`, `lubridate`, `splines`, `glmnet`,
`randomForest`, `caret`, `MASS`, `MatchIt`, `cobalt`, `WeightIt`,
`topicmodels`, `tidytext`, `Matrix`, `slam`

**SQL:** Google BigQuery with MIMIC-IV v3.1 access
(`physionet-data.mimiciv_3_1_hosp`, `physionet-data.mimiciv_3_1_icu`,
`physionet-data.mimiciv_3_1_derived`)

**Visualization:** Tableau Desktop 2020.1+ (publish to Tableau Public)

------------------------------------------------------------------------

## Acknowledgments

MIMIC-IV was developed by the MIT Laboratory for Computational
Physiology and is made available through PhysioNet. This project was
developed in collaboration with Andrew Hersh, MD, who contributed
clinical domain expertise, motivated the caregiver-level research
questions, and proposed the propensity score matching design to assess
the causal plausibility of the volume-outcome relationship.

Johnson, A., Bulgarelli, L., Shen, L. et al. MIMIC-IV, a freely
accessible electronic health record dataset. *Sci Data* 10, 1 (2023).
<https://doi.org/10.1038/s41597-022-01899-x>

------------------------------------------------------------------------

## License

Code in this repository (SQL, R) is released under the MIT License.
Synthetic data files are released under CC0 (public domain). No MIMIC-IV
data is included; use of MIMIC-IV is governed by the [PhysioNet
Credentialed Health Data
License](https://physionet.org/content/mimiciv/view-license/3.1/).
