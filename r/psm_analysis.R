################################################################################
##
##  ICU Last Extubations — Propensity Score and Data Quality Analysis
##
##  Companion script to MIMIC_analysis_clean.R
##
##  Research question:
##    Does caregiver failed extubation (FE) rate predict patient 12-month
##    survival after adjustment for measured confounders?
##
##  Summary of findings:
##    Initial PSM on the full dataset (including inferred tube events) found
##    OR = 0.733 (p < 0.001). Sensitivity analyses revealed this was driven
##    by a data quality artifact: caregivers in the high-FE-rate group had
##    87.8% of their events classified as inferred (algorithmically reconstructed
##    from ventilation records), vs. 39.4% in the low-FE-rate group. Inferred
##    events generate spurious failed extubation classifications.
##
##    Restricting to explicitly documented procedure events only, and weighting
##    by caregiver volume, the association is not statistically significant
##    (weighted r = -0.134, p = 0.124). The apparent effect in the full dataset
##    reflects inferred event misclassification and small-sample caregiver rate
##    instability rather than true provider performance differences.
##
##  Input files (require MIMIC-IV access — not included in repository):
##    last_extubations_clean.csv   — patient-level cohort from BigQuery
##    patient_ccsr_long.csv        — long-format CCSR codes per patient
##
##  Outputs:
##    psm_love_plot.png            — covariate balance before/after PSM
##    psm_results_summary.csv      — results table (see Section 9)
##
##  CCSR long query (run in BigQuery, export as patient_ccsr_long.csv):
##
##    SELECT
##      e.subject_id,
##      e.hadm_id,
##      e.caregiver_fe_rate,
##      e.caregiver_n,
##      e.failed_extubations,
##      code AS ccsr_code
##    FROM `[your-project].Failed_Extubation_Rate.last_extubations` e,
##    UNNEST(e.ccsr_codes) AS code
##    WHERE e.caregiver_n > 10
##      AND e.ccsr_codes IS NOT NULL
##    ORDER BY e.subject_id, code
##
################################################################################


## ── 0. LIBRARIES AND SETUP ───────────────────────────────────────────────────

library(MASS)         # mvrnorm — load before dplyr to avoid select() conflict
library(dplyr)
library(tidyr)
library(purrr)
library(readr)
library(tibble)
library(stringr)
library(forcats)
library(ggplot2)
library(MatchIt)      # propensity score matching
library(cobalt)       # balance diagnostics and love plots
library(WeightIt)     # IPW sensitivity analysis
library(topicmodels)  # LDA topic modeling
library(tidytext)     # cast_dtm for document-term matrix
library(Matrix)       # sparse matrix support
library(weights)      # wtd.cor for volume-weighted correlation

set.seed(237)

# setwd("~/Projects/Extubate/r")   # adjust as needed


## ── 1. DATA INGESTION ────────────────────────────────────────────────────────

ccsr_long <- read_csv(
  "../data/patient_ccsr_long.csv",
  col_types = cols(
    subject_id   = col_character(),
    hadm_id      = col_character(),
    ccsr_code    = col_character()
  )
)

patients <- read_csv(
  "../data/last_extubations_clean.csv",
  col_types = cols(
    subject_id      = col_character(),
    hadm_id         = col_character(),
    gender          = col_factor(),
    intubation_type = col_factor()
  )
) %>%
  mutate(survival_12mo = is.na(dod))

cat("Patients:", nrow(patients), "\n")
cat("CCSR rows:", nrow(ccsr_long), "\n")
cat("Unique CCSR codes:", n_distinct(ccsr_long$ccsr_code), "\n")


## ── 2. CCSR DOCUMENT-TERM MATRIX ─────────────────────────────────────────────
##
##  Each patient is a "document", each CCSR code is a "term".
##  Used as input to LDA topic modeling (Section 3).

ccsr_counts <- ccsr_long %>%
  count(subject_id, ccsr_code) %>%
  rename(count = n)

ccsr_dtm <- ccsr_counts %>%
  cast_dtm(document = subject_id,
           term     = ccsr_code,
           value    = count)

cat("DTM dimensions:", dim(ccsr_dtm), "\n")


## ── 3. LDA TOPIC MODELING ────────────────────────────────────────────────────
##
##  Fit LDA for K = 5, 8, 10, 12 and compare perplexity.
##  K = 8 selected based on elbow in perplexity and clinical interpretability.
##
##  Perplexity results:
##    K=5:  90.65
##    K=8:  86.80
##    K=10: 85.70
##    K=12: 84.50
##
##  Note: topic proportions showed weak differentiation (SD ~0.02-0.03 per
##  topic, max loading ~0.17). The short-document problem — median 4.2 CCSR
##  codes per patient — limits LDA performance. K-means on the binary CCSR
##  matrix also showed weak separation (max variance explained ~16% at K=10).
##  Both approaches were superseded by direct CCSR binary flags (Section 4c).

control_lda <- list(seed = 237, burnin = 100, iter = 500, thin = 10)

cat("Fitting LDA K=8 (primary model)...\n")
lda_8 <- LDA(ccsr_dtm, k = 8, method = "Gibbs", control = control_lda)

cat("Perplexity K=8:", perplexity(lda_8, newdata = ccsr_dtm), "\n")

# Top terms per topic — clinical interpretation
cat("\nTop 10 CCSR codes per topic:\n")
print(terms(lda_8, k = 10))

##  Topic clinical interpretation (K=8):
##    Topic 1: Cardiac — heart failure / vascular (CIR019, CIR008, CIR009)
##    Topic 2: Neurologic — CNS with cardiac (NVS020, CIR021, NVS008)
##    Topic 3: Coronary artery disease / metabolic (CIR011, CIR007, END010)
##    Topic 4: Trauma / injury with complications (BLD004, RSP016, INJ037)
##    Topic 5: Respiratory failure / metabolic (RSP012, END011, RSP010)
##    Topic 6: Sepsis / infectious / digestive (INF002, RSP012, END008)
##    Topic 7: Respiratory infection / sepsis (RSP012, INF002, GEN002)
##    Topic 8: Renal / metabolic / digestive (GEN002, BLD006, END011)


## ── 4a. TOPIC PROPORTION EXTRACTION ─────────────────────────────────────────
##
##  Extract patient-level topic proportions from LDA posterior.
##  Used for initial PSM covariate exploration; superseded by Section 4c.

topic_proportions <- posterior(lda_8)$topics %>%
  as.data.frame() %>%
  tibble::rownames_to_column("subject_id") %>%
  rename(
    topic_cardiac_hf       = `1`,
    topic_neurologic       = `2`,
    topic_cad_metabolic    = `3`,
    topic_trauma           = `4`,
    topic_respiratory      = `5`,
    topic_sepsis_abdominal = `6`,
    topic_respiratory_infx = `7`,
    topic_renal_metabolic  = `8`
  )

patients_topics <- patients %>%
  left_join(topic_proportions, by = "subject_id")


## ── 4b. K-MEANS CLUSTERING (EXPLORATORY) ─────────────────────────────────────
##
##  Binary patient × CCSR matrix, k-means for K = 6, 8, 10.
##  Variance explained: K=6: 12.6%, K=8: 14.6%, K=10: 16.1%.
##  Weak separation consistent with multimorbid ICU population.
##  Not used as primary PSM covariate — see Section 4c.

ccsr_binary <- ccsr_long %>%
  distinct(subject_id, ccsr_code) %>%
  mutate(present = 1) %>%
  pivot_wider(
    id_cols     = subject_id,
    names_from  = ccsr_code,
    values_from = present,
    values_fill = 0
  )

cat("Binary matrix dimensions:", dim(ccsr_binary), "\n")

set.seed(123)
km_8 <- kmeans(ccsr_binary[, -1], centers = 8, nstart = 25, iter.max = 100)
cat("K-means K=8 variance explained:",
    round(km_8$betweenss / km_8$totss, 3), "\n")


## ── 4c. CCSR BINARY FLAGS (PRIMARY COVARIATE STRATEGY) ───────────────────────
##
##  Given weak clustering, use high-prevalence CCSR codes directly as binary
##  flags. Codes present in >= 10% of patients (22 codes identified).
##  RSP012 (respiratory failure, 66.7% prevalence) excluded — near-universal
##  in ventilated patients, near-zero variance.
##
##  CCSR codes and prevalence:
##    RSP012: 66.7%  GEN002: 41.3%  END011: 38.1%  INF002: 32.0%
##    CIR019: 24.5%  BLD004: 24.2%  RSP002: 23.0%  NVS020: 21.2%
##    BLD006: 19.2%  CIR017: 17.0%  RSP010: 16.4%  INJ030: 16.0%
##    SYM003: 15.6%  CIR008: 14.9%  CIR009: 14.2%  INF003: 13.6%
##    END008: 13.6%  RSP011: 12.8%  FAC025: 12.0%  CIR011: 11.8%
##    CIR007: 11.5%  GEN004: 10.2%

ccsr_prevalence <- ccsr_long %>%
  distinct(subject_id, ccsr_code) %>%
  count(ccsr_code, name = "n_patients") %>%
  mutate(prevalence = n_patients / n_distinct(ccsr_long$subject_id)) %>%
  arrange(desc(prevalence))

psm_ccsr_codes <- ccsr_prevalence %>%
  filter(prevalence >= 0.10, ccsr_code != "RSP012") %>%
  pull(ccsr_code)

cat("CCSR binary flags used as PSM covariates:", length(psm_ccsr_codes), "\n")

ccsr_flags <- ccsr_long %>%
  distinct(subject_id, ccsr_code) %>%
  filter(ccsr_code %in% psm_ccsr_codes) %>%
  mutate(present = 1) %>%
  pivot_wider(
    id_cols     = subject_id,
    names_from  = ccsr_code,
    values_from = present,
    values_fill = 0
  )

analysis_data <- patients_topics %>%
  left_join(ccsr_flags, by = "subject_id") %>%
  mutate(across(all_of(psm_ccsr_codes), ~ replace_na(.x, 0))) %>%
  mutate(survival_12mo = is.na(dod))

cat("Analysis dataset N:", nrow(analysis_data), "\n")


## ── 5. PSM DATASET CONSTRUCTION ──────────────────────────────────────────────
##
##  Treatment variable: caregiver_fe_rate above 75th percentile = "high"
##                      caregiver_fe_rate below 25th percentile = "low"
##  Middle 50% excluded from PSM to sharpen the contrast.
##
##  Full dataset quartiles:
##    Q25 = 0.120 (12% FE rate)
##    Q75 = 0.319 (32% FE rate)
##
##  Note: the high-FE group (>32%) is dominated by caregivers with inferred
##  events (87.8% inferred vs 39.4% in low group) — see Section 9 for details.

analysis_data <- analysis_data %>%
  mutate(
    fe_q25 = quantile(caregiver_fe_rate, 0.25, na.rm = TRUE),
    fe_q75 = quantile(caregiver_fe_rate, 0.75, na.rm = TRUE),
    fe_group = case_when(
      caregiver_fe_rate >= fe_q75 ~ "high",
      caregiver_fe_rate <= fe_q25 ~ "low",
      TRUE                        ~ "middle"
    )
  )

cat("FE rate Q25:", round(analysis_data$fe_q25[1], 3), "\n")
cat("FE rate Q75:", round(analysis_data$fe_q75[1], 3), "\n")
cat("Group sizes:\n")
print(table(analysis_data$fe_group))

psm_data <- analysis_data %>%
  filter(fe_group != "middle") %>%
  mutate(
    treated        = as.integer(fe_group == "high"),
    gender_m       = as.integer(gender == "M"),
    intub_surgical = as.integer(intubation_type == "surgical"),
    intub_med_resp = as.integer(intubation_type == "medical-respiratory")
  ) %>%
  filter(
    !is.na(charlson),
    !is.na(sofa),
    !is.na(norepinephrine),
    !is.na(vent_hours)
  )

cat("PSM dataset N:", nrow(psm_data), "\n")

ccsr_formula_terms <- paste(psm_ccsr_codes, collapse = " + ")

psm_formula <- as.formula(paste(
  "treated ~ anchor_age + gender_m + charlson + sofa + norepinephrine +
   vent_hours + intub_surgical + intub_med_resp +",
  ccsr_formula_terms
))


## ── 6. PROPENSITY SCORE MATCHING ─────────────────────────────────────────────
##
##  Nearest neighbor 1:1 matching within caliper = 0.2 SD of logit PS.
##  Matched on: age, charlson, sofa, norepinephrine, vent_hours,
##              intubation type, 21 CCSR binary flags.
##
##  Results: 2,353 matched pairs; all SMDs < 0.10 post-match.

set.seed(237)
match_out <- matchit(
  psm_formula,
  data        = psm_data,
  method      = "nearest",
  ratio       = 1,
  caliper     = 0.2,
  std.caliper = TRUE,
  distance    = "logit"
)

summary(match_out)


## ── 7. TREATMENT EFFECT ESTIMATION ───────────────────────────────────────────
##
##  Primary outcome: survival_12mo (12-month survival, derived from dod)
##
##  Full dataset PSM results:
##    Survival — high FE caregivers: 50.8%
##    Survival — low FE caregivers:  56.8%
##    Unadjusted OR: 0.786 (0.700-0.881), p < 0.001
##    Doubly adjusted OR: 0.733 (0.643-0.836), p < 0.001
##
##  Note: these results are driven by inferred event artifact — see Section 9.

matched_data <- match.data(match_out)

cat("Matched N:", nrow(matched_data), "\n")
cat("Survival — high FE:",
    round(mean(matched_data$survival_12mo[matched_data$treated == 1],
               na.rm = TRUE), 3), "\n")
cat("Survival — low FE:",
    round(mean(matched_data$survival_12mo[matched_data$treated == 0],
               na.rm = TRUE), 3), "\n")

# Unadjusted outcome model
outcome_model <- glm(
  survival_12mo ~ treated,
  data   = matched_data,
  family = binomial,
  weights = weights
)
summary(outcome_model)

or <- exp(coef(outcome_model)["treated"])
ci <- exp(confint(outcome_model)["treated", ])
cat("\nOR (unadjusted):", round(or, 3), "\n")
cat("95% CI:", round(ci, 3), "\n")
cat("p-value:", round(summary(outcome_model)$coefficients["treated", "Pr(>|z|)"], 4), "\n")

# Doubly adjusted outcome model
outcome_model_adj <- glm(
  survival_12mo ~ treated + anchor_age + charlson + sofa +
    norepinephrine + vent_hours + intub_surgical + intub_med_resp,
  data   = matched_data,
  family = binomial,
  weights = weights
)

or_adj <- exp(coef(outcome_model_adj)["treated"])
ci_adj <- exp(confint(outcome_model_adj)["treated", ])
cat("\nOR (doubly adjusted):", round(or_adj, 3), "\n")
cat("95% CI:", round(ci_adj, 3), "\n")
cat("p-value:", round(summary(outcome_model_adj)$coefficients["treated", "Pr(>|z|)"], 4), "\n")


## ── 8. IPW SENSITIVITY AND BALANCE DIAGNOSTICS ───────────────────────────────
##
##  IPW full-dataset result: OR = 0.789 (0.735-0.846), p < 0.001
##  Consistent with PSM — both driven by inferred event artifact.
##
##  Effect modification by intubation type: p = 0.623 (no interaction)

# IPW
set.seed(237)
ipw_out <- weightit(
  psm_formula,
  data     = psm_data,
  method   = "ps",
  estimand = "ATE"
)
summary(ipw_out)

ipw_model <- glm(
  survival_12mo ~ treated,
  data    = psm_data,
  family  = binomial,
  weights = ipw_out$weights
)

or_ipw <- exp(coef(ipw_model)["treated"])
ci_ipw <- exp(confint(ipw_model)["treated", ])
cat("\nIPW OR:", round(or_ipw, 3), "\n")
cat("IPW 95% CI:", round(ci_ipw, 3), "\n")
cat("IPW p-value:", round(summary(ipw_model)$coefficients["treated", "Pr(>|z|)"], 4), "\n")

# Effect modification
outcome_interaction <- glm(
  survival_12mo ~ treated * intub_med_resp,
  data   = matched_data,
  family = binomial,
  weights = weights
)
cat("\nInteraction p-value (treated × medical-respiratory):",
    round(summary(outcome_interaction)$coefficients[
      "treated:intub_med_resp", "Pr(>|z|)"], 4), "\n")

# Love plot — covariate balance
p_love <- love.plot(
  match_out,
  threshold  = 0.1,
  title      = "Covariate Balance Before and After PSM",
  colors     = c("#D95F02", "#1D9E8E"),
  stars      = "std",
  var.order  = "unadjusted",
  abs        = TRUE
)
print(p_love)
ggsave("../figures/psm_love_plot.png", plot = p_love,
       width = 10, height = 12, dpi = 150)
cat("Love plot saved.\n")


## ── 9. DATA QUALITY SENSITIVITY ANALYSES ─────────────────────────────────────
##
##  Key finding: the PSM result in Sections 6-8 is driven by inferred events.
##  This section documents the diagnosis and the evidence.

## 9a. Inferred event prevalence by FE rate group
cat("\n── DATA QUALITY ANALYSIS ────────────────────────────────────────────────\n")

analysis_data %>%
  mutate(fe_tier = case_when(
    caregiver_fe_rate < 0.15  ~ "Low (<15%)",
    caregiver_fe_rate < 0.32  ~ "Middle (15-32%)",
    TRUE                      ~ "High (>32%)"
  )) %>%
  group_by(fe_tier) %>%
  summarise(
    n_patients          = n(),
    n_caregivers        = n_distinct(caregiver_id),
    median_caregiver_n  = median(caregiver_n),
    pct_inferred        = mean(tube_event_source == "inferred"),
    mean_sofa           = mean(sofa, na.rm = TRUE),
    survival_12mo       = mean(survival_12mo, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  print()

##  Result:
##    Low (<15%):    39.4% inferred, median caregiver_n = 979
##    Middle:        63.5% inferred, median caregiver_n = 159
##    High (>32%):   87.8% inferred, median caregiver_n = 78
##
##  The high-FE group is dominated by caregivers whose extubation events were
##  inferred from ventilation records, not explicitly documented. These caregivers
##  also have much lower extubation volume, suggesting their FE rates are driven
##  by both misclassification and small-sample instability.

## 9b. FE rate distribution — explicit events only
cat("\nFE rate distribution (explicit events only):\n")
analysis_data %>%
  filter(tube_event_source == "explicit") %>%
  pull(caregiver_fe_rate) %>%
  summary() %>%
  print()

## 9c. Survival by FE category — explicit events only
cat("\nSurvival by FE category (explicit events only):\n")
analysis_data %>%
  filter(tube_event_source == "explicit") %>%
  mutate(fe_cat = case_when(
    caregiver_fe_rate < 0.05  ~ "1. <5%",
    caregiver_fe_rate < 0.15  ~ "2. 5-15%",
    caregiver_fe_rate < 0.25  ~ "3. 15-25%",
    TRUE                      ~ "4. >25%"
  )) %>%
  group_by(fe_cat) %>%
  summarise(
    n_patients    = n(),
    n_caregivers  = n_distinct(caregiver_id),
    mean_fe_rate  = mean(caregiver_fe_rate),
    survival_12mo = mean(survival_12mo, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  print()

##  Result: above 5% FE rate, survival is essentially flat (67-69%).
##  No monotonic gradient in explicitly documented data.

## 9d. PSM restricted to explicit events only
cat("\nPSM — explicit events only:\n")

psm_data_explicit <- analysis_data %>%
  filter(tube_event_source == "explicit") %>%
  mutate(
    fe_q25 = quantile(caregiver_fe_rate, 0.25, na.rm = TRUE),
    fe_q75 = quantile(caregiver_fe_rate, 0.75, na.rm = TRUE),
    fe_group = case_when(
      caregiver_fe_rate >= fe_q75 ~ "high",
      caregiver_fe_rate <= fe_q25 ~ "low",
      TRUE                        ~ "middle"
    )
  ) %>%
  filter(fe_group != "middle") %>%
  mutate(
    treated        = as.integer(fe_group == "high"),
    gender_m       = as.integer(gender == "M"),
    intub_surgical = as.integer(intubation_type == "surgical"),
    intub_med_resp = as.integer(intubation_type == "medical-respiratory")
  ) %>%
  filter(!is.na(charlson), !is.na(sofa),
         !is.na(norepinephrine), !is.na(vent_hours))

cat("Explicit-only PSM N:", nrow(psm_data_explicit), "\n")
cat("Q25:", round(psm_data_explicit$fe_q25[1], 3), "\n")
cat("Q75:", round(psm_data_explicit$fe_q75[1], 3), "\n")

set.seed(237)
match_out_explicit <- matchit(
  psm_formula,
  data        = psm_data_explicit,
  method      = "nearest",
  ratio       = 1,
  caliper     = 0.2,
  std.caliper = TRUE,
  distance    = "logit"
)

matched_explicit <- match.data(match_out_explicit)

cat("\nSurvival — high FE (explicit):",
    round(mean(matched_explicit$survival_12mo[matched_explicit$treated == 1],
               na.rm = TRUE), 3), "\n")
cat("Survival — low FE (explicit):",
    round(mean(matched_explicit$survival_12mo[matched_explicit$treated == 0],
               na.rm = TRUE), 3), "\n")

outcome_explicit <- glm(
  survival_12mo ~ treated,
  data   = matched_explicit,
  family = binomial,
  weights = weights
)
or_exp <- exp(coef(outcome_explicit)["treated"])
ci_exp <- exp(confint(outcome_explicit)["treated", ])
cat("OR (explicit-only):", round(or_exp, 3), "\n")
cat("95% CI:", round(ci_exp, 3), "\n")
cat("p-value:", round(summary(outcome_explicit)$coefficients["treated", "Pr(>|z|)"], 4), "\n")

##  Result: OR = 0.810 (0.678-0.968), p = 0.021. Attenuated but still
##  significant. Volume stratification (Section 9e) explains the residual.

## 9e. Volume-stratified analysis (explicit matched cohort)
cat("\nSurvival difference by caregiver volume tier (explicit matched cohort):\n")

matched_explicit %>%
  mutate(volume_group = ifelse(caregiver_n < median(caregiver_n),
                               "Lower volume", "Higher volume")) %>%
  group_by(volume_group, treated) %>%
  summarise(
    survival = mean(survival_12mo, na.rm = TRUE),
    n        = n(),
    .groups  = "drop"
  ) %>%
  pivot_wider(names_from  = treated,
              values_from = c(survival, n),
              names_prefix = "treated_") %>%
  mutate(survival_diff = survival_treated_1 - survival_treated_0) %>%
  print()

##  Result:
##    Higher-volume caregivers: survival diff = -1.3pp (no meaningful effect)
##    Lower-volume caregivers:  survival diff = -8.9pp (small-sample noise)
##
##  The residual effect in explicit-only PSM is concentrated in low-volume
##  caregivers whose FE rates are most unstable.

## 9f. Volume-weighted caregiver-level correlation (primary conclusion)
cat("\nVolume-weighted caregiver-level correlation (explicit events, n > 5):\n")

caregivers_explicit <- analysis_data %>%
  filter(tube_event_source == "explicit") %>%
  group_by(caregiver_id) %>%
  summarise(
    survival         = mean(survival_12mo, na.rm = TRUE),
    caregiver_fe_rate = mean(caregiver_fe_rate, na.rm = TRUE),
    n                = n(),
    .groups          = "drop"
  ) %>%
  filter(n > 5)

cat("Caregivers included:", nrow(caregivers_explicit), "\n")

# Volume-weighted correlation
wtd_cor <- wtd.cor(caregivers_explicit$caregiver_fe_rate,
                   caregivers_explicit$survival,
                   weight = caregivers_explicit$n)
cat("Weighted correlation r:", round(wtd_cor[1, "correlation"], 3), "\n")
cat("p-value:", round(wtd_cor[1, "p.value"], 4), "\n")

# Volume-weighted linear model
lm_weighted <- lm(survival ~ caregiver_fe_rate,
                  data    = caregivers_explicit,
                  weights = n)
cat("\nWeighted regression:\n")
summary(lm_weighted)

##  Result: r = -0.134, p = 0.124. Not statistically significant.
##  This is the primary conclusion of the sensitivity analysis.

# Scatter plot: caregiver FE rate vs survival (explicit, volume-weighted)
ggplot(caregivers_explicit, aes(caregiver_fe_rate, survival, size = n)) +
  geom_point(alpha = 0.5) +
  geom_smooth(aes(weight = n), method = "lm",
              se = TRUE, color = "#2E75B6", linewidth = 0.8) +
  geom_smooth(method = "lm", se = FALSE,
              color = "#D95F02", linetype = "dashed", linewidth = 0.8) +
  annotate("text", x = 0.35, y = 0.90,
           label = "Unweighted trend", color = "#D95F02", size = 3.5) +
  annotate("text", x = 0.35, y = 0.83,
           label = "Volume-weighted trend", color = "#2E75B6", size = 3.5) +
  labs(
    title    = "12-Mo Survival by Caregiver FE Rate — Explicit Events Only",
    subtitle = "Bubble size = patients in sample. Trend lines nearly identical (r = -0.134, p = 0.124).",
    x        = "Caregiver FE Rate",
    y        = "12-Month Survival Rate"
  ) +
  theme_minimal()


## ── 10. RESULTS SUMMARY ───────────────────────────────────────────────────────

cat("\n── FINAL RESULTS SUMMARY ────────────────────────────────────────────────\n")
cat("Full dataset PSM (unadjusted):          OR = 0.786 (0.700-0.881), p < 0.001\n")
cat("Full dataset PSM (doubly adjusted):     OR = 0.733 (0.643-0.836), p < 0.001\n")
cat("Full dataset IPW:                       OR = 0.789 (0.735-0.846), p < 0.001\n")
cat("  → Above results driven by inferred event artifact\n\n")
cat("Explicit events PSM:                    OR = 0.810 (0.678-0.968), p = 0.021\n")
cat("  → Concentrated in low-volume caregivers (small-sample noise)\n\n")
cat("Volume-weighted correlation (explicit): r = -0.134, p = 0.124\n")
cat("  → NOT SIGNIFICANT — primary conclusion\n\n")
cat("Conclusion: caregiver FE rate is not a reliable predictor of patient\n")
cat("12-month survival in high-quality, volume-weighted explicit-event data.\n")
cat("The apparent effect reflects inferred event misclassification and\n")
cat("small-sample caregiver rate instability.\n")
cat("─────────────────────────────────────────────────────────────────────────\n")

results_summary <- tibble(
  Analysis = c(
    "PSM — full dataset, unadjusted",
    "PSM — full dataset, doubly adjusted",
    "IPW — full dataset",
    "PSM — explicit events only",
    "Volume-weighted correlation (explicit, caregiver-level)"
  ),
  Estimate = c("OR 0.786", "OR 0.733", "OR 0.789", "OR 0.810", "r = -0.134"),
  CI       = c("0.700-0.881", "0.643-0.836", "0.735-0.846", "0.678-0.968", "—"),
  p_value  = c("<0.001", "<0.001", "<0.001", "0.021", "0.124"),
  Note     = c(
    "Driven by inferred event artifact",
    "Driven by inferred event artifact",
    "Driven by inferred event artifact",
    "Concentrated in low-volume caregivers",
    "Primary conclusion — NOT significant"
  )
)

print(results_summary)
write_csv(results_summary, "../data/psm_results_summary.csv")
cat("Results summary written.\n")
