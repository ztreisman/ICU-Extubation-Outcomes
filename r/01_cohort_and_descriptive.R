################################################################################
##
##  ICU Extubation Outcomes — Script 01: Cohort Construction and Descriptive Analysis
##
##  Data source: MIMIC-IV v3.1 (BigQuery), extracted via ICU_Last_Extubations.sql
##
##  This script constructs two analytical cohorts and produces all descriptive
##  statistics, Table 1, and population-level visualizations.
##
##  Cohorts:
##    patient_cohort       — all patients with a documented last extubation event
##                           meeting basic quality filters. Used for population
##                           description. One of the largest ICU extubation
##                           cohorts derived from MIMIC-IV.
##
##    explicit_extubations — subset of patient_cohort restricted to explicitly
##                           documented procedure events with a genuine caregiver
##                           assignment (caregiver_fe_rate not missing). Used for
##                           all caregiver-level analyses in Scripts 02 and 03.
##
##  SQL pipeline summary (ICU_Last_Extubations.sql):
##    - Combines explicit procedure events (itemids 224385, 227194, 225468, 225477)
##      with inferred events from ventilation segments (>= 4h InvasiveVent)
##    - Deduplicates inferred events within 30 min of explicit events
##    - Classifies intubation type: surgical / medical-respiratory /
##      medical-non-respiratory
##    - Selects last extubation per patient (one row per patient)
##    - Caregiver assignment cascade (explicit events only):
##        1. Directly recorded caregiver_id
##        2. Other procedureevents within ±30 min
##        3. chartevents within ±15 min (nearest-neighbor, bug fixed from v1)
##        4. inputevents/outputevents within ±5 min
##        5. Mode caregiver within ±2h window
##        6. NULL for inferred events (no imputation)
##    - caregiver_fe_rate computed from explicit events only; NULL if no
##      caregiver assignment or caregiver_n <= 0
##    - Joins: Charlson comorbidity index, norepinephrine equivalent dose,
##             first-day SOFA, hospital mortality, ICD chapter, CCSR arrays,
##             first_careunit from icustays
##
##  Output objects (passed to Scripts 02 and 03 via save/load or run_all.R):
##    patient_cohort, explicit_extubations
##
################################################################################


## ── 0. LIBRARIES ─────────────────────────────────────────────────────────────

library(MASS)       # load before dplyr to avoid select() conflict
library(dplyr)
library(tidyr)
library(readr)
library(tibble)
library(stringr)
library(forcats)
library(ggplot2)
library(lubridate)
library(splines)
library(tableone)   # Table 1 construction

set.seed(237)

# setwd("~/Projects/Extubate/r")   # adjust as needed


## ── 1. DATA INGESTION ────────────────────────────────────────────────────────

last_extubations <- read_csv(
  "../data/last_extubations.csv",
  col_types = cols(
    subject_id                  = col_character(),
    hadm_id                     = col_character(),
    caregiver_id                = col_character(),
    gender                      = col_factor(),
    intubation_type             = col_factor(
      levels = c("surgical", "medical-respiratory", "medical-non-respiratory")
    ),
    tube_event_source           = col_factor(levels = c("explicit", "inferred")),
    caregiver_imputation_source = col_factor(),
    first_careunit              = col_factor(),
    event_time                  = col_datetime(),
    dod                         = col_date()
  )
) %>%
  mutate(survival_12mo = is.na(dod))

cat("Raw rows from SQL:", nrow(last_extubations), "\n")
cat("Unique patients:", n_distinct(last_extubations$subject_id), "\n")


## ── 2. COHORT CONSTRUCTION ───────────────────────────────────────────────────

## patient_cohort: broad cohort for population description
##   All patients with a valid ventilation duration and plausible
##   clinical values. Includes both explicit and inferred events.
##   No caregiver-level filters applied.

patient_cohort <- last_extubations %>%
  filter(
    !is.na(vent_hours),
    !is.na(sofa),
    !is.na(charlson),
    !is.na(intubation_type),
    !is.na(primary_diagnosis)
  )

cat("\npatient_cohort N:", nrow(patient_cohort), "\n")
cat("Explicit events:", sum(patient_cohort$tube_event_source == "explicit"), "\n")
cat("Inferred events:", sum(patient_cohort$tube_event_source == "inferred"), "\n")
cat("12-month survival rate:", round(mean(patient_cohort$survival_12mo, na.rm = TRUE), 3), "\n")
cat("Failed extubation rate:", round(mean(patient_cohort$failed_extubations > 0, na.rm = TRUE), 3), "\n")

## explicit_extubations: analytical cohort for caregiver-level analysis
##   - Explicit tube events only (directly documented procedure events)
##   - Genuine caregiver assignment (caregiver_fe_rate not missing)
##   - caregiver_n > 10: FE rate estimated from > 10 extubations
##   - Complete covariates for modeling

explicit_extubations <- patient_cohort %>%
  filter(
    tube_event_source == "explicit",
    !is.na(caregiver_fe_rate),
    caregiver_n > 10
  ) %>%
  mutate(
    # Centered predictors for mixed effects models
    charlson_centered       = charlson       - mean(charlson,       na.rm = TRUE),
    sofa_centered           = sofa           - mean(sofa,           na.rm = TRUE),
    norepinephrine_centered = norepinephrine - mean(norepinephrine, na.rm = TRUE),
    vent_hours_centered     = vent_hours     - mean(vent_hours,     na.rm = TRUE),
    # Scaled predictors for Bayesian models
    charlson_scaled         = as.numeric(scale(charlson)),
    sofa_scaled             = as.numeric(scale(sofa)),
    norepinephrine_scaled   = as.numeric(scale(norepinephrine)),
    vent_hours_scaled       = as.numeric(scale(vent_hours)),
    caregiver_fe_scaled     = as.numeric(scale(caregiver_fe_rate))
  )

cat("\nexplicit_extubations N:", nrow(explicit_extubations), "\n")
cat("Unique caregivers:", n_distinct(explicit_extubations$caregiver_id), "\n")
cat("Unique ICU units:", n_distinct(explicit_extubations$first_careunit), "\n")
cat("12-month survival rate:", round(mean(explicit_extubations$survival_12mo, na.rm = TRUE), 3), "\n")
cat("Failed extubation rate:", round(mean(explicit_extubations$failed_extubations > 0, na.rm = TRUE), 3), "\n")
cat("Caregiver FE rate — median:", round(median(explicit_extubations$caregiver_fe_rate), 3),
    "IQR:", round(quantile(explicit_extubations$caregiver_fe_rate, 0.25), 3),
    "—", round(quantile(explicit_extubations$caregiver_fe_rate, 0.75), 3), "\n")


## ── 3. TABLE 1: PATIENT CHARACTERISTICS ──────────────────────────────────────
##
##  Stratified by tube_event_source (explicit vs inferred) to characterize
##  the two populations and justify the analytical cohort restriction.

table1_vars <- c(
  "anchor_age", "gender", "charlson", "sofa", "norepinephrine",
  "vent_hours", "intubation_type", "primary_diagnosis",
  "failed_extubations", "survival_12mo", "hospital_expire_flag"
)

table1_cat <- c(
  "gender", "intubation_type", "primary_diagnosis",
  "survival_12mo", "hospital_expire_flag"
)

tab1 <- CreateTableOne(
  vars       = table1_vars,
  strata     = "tube_event_source",
  data       = patient_cohort,
  factorVars = table1_cat,
  addOverall = TRUE
)
print(tab1, showAllLevels = FALSE, smd = TRUE)

## Table 1b: explicit_extubations stratified by FE rate quartile
patient_cohort_explicit <- explicit_extubations %>%
  mutate(
    fe_quartile = case_when(
      caregiver_fe_rate <= quantile(caregiver_fe_rate, 0.25) ~ "Q1 (lowest)",
      caregiver_fe_rate <= quantile(caregiver_fe_rate, 0.50) ~ "Q2",
      caregiver_fe_rate <= quantile(caregiver_fe_rate, 0.75) ~ "Q3",
      TRUE                                                   ~ "Q4 (highest)"
    ) %>% factor(levels = c("Q1 (lowest)", "Q2", "Q3", "Q4 (highest)"))
  )

tab1b <- CreateTableOne(
  vars       = table1_vars,
  strata     = "fe_quartile",
  data       = patient_cohort_explicit,
  factorVars = table1_cat,
  addOverall = TRUE
)
print(tab1b, showAllLevels = FALSE, smd = TRUE)


## ── 4. CAREGIVER-LEVEL DESCRIPTIVE ───────────────────────────────────────────
##
##  One row per caregiver. Used for all caregiver-level visualizations
##  and the partial correlation analysis in Script 03.

caregiver_summary <- explicit_extubations %>%
  group_by(caregiver_id) %>%
  summarise(
    caregiver_fe_rate  = mean(caregiver_fe_rate),
    caregiver_n        = mean(caregiver_n),
    survival_rate      = mean(survival_12mo, na.rm = TRUE),
    hospital_expire    = mean(hospital_expire_flag, na.rm = TRUE),
    n_patients_sample  = n(),
    first_careunit     = first(first_careunit),
    .groups            = "drop"
  ) %>%
  filter(n_patients_sample > 5)

cat("\nCaregivers with > 5 sample patients:", nrow(caregiver_summary), "\n")
cat("Caregiver FE rate range:",
    round(range(caregiver_summary$caregiver_fe_rate), 3), "\n")
cat("Caregiver volume range:",
    round(range(caregiver_summary$caregiver_n)), "\n")


## ── 5. POPULATION VISUALIZATIONS ─────────────────────────────────────────────

## 5a. Caregiver FE rate distribution
ggplot(caregiver_summary, aes(caregiver_fe_rate)) +
  geom_histogram(bins = 30, fill = "#2E75B6", color = "white") +
  labs(
    title = "Distribution of Caregiver Failed Extubation Rates",
    subtitle = "Explicit events only; caregivers with > 10 total extubations",
    x = "Caregiver FE Rate",
    y = "Number of Caregivers"
  ) +
  theme_minimal()

## 5b. Caregiver volume vs FE rate
ggplot(caregiver_summary, aes(caregiver_n, caregiver_fe_rate)) +
  geom_point(aes(size = n_patients_sample), alpha = 0.5) +
  geom_smooth(method = "loess", se = TRUE, color = "#2E75B6") +
  labs(
    title   = "Caregiver Volume vs. Failed Extubation Rate",
    subtitle = "Bubble size = patients in analytical sample",
    x       = "Caregiver Total Extubation Volume",
    y       = "Caregiver FE Rate",
    size    = "N patients"
  ) +
  theme_minimal()

## 5c. FE rate vs 12-month survival (caregiver level, FE rate > 0)
caregivers_nonzero <- caregiver_summary %>%
  filter(caregiver_fe_rate > 0)

# Fit competing functional forms (used in Script 03 and paper)
lm_linear    <- lm(survival_rate ~ caregiver_fe_rate,
                   data = caregivers_nonzero, weights = n_patients_sample)
lm_quadratic <- lm(survival_rate ~ caregiver_fe_rate + I(caregiver_fe_rate^2),
                   data = caregivers_nonzero, weights = n_patients_sample)
lm_log       <- lm(log(survival_rate) ~ caregiver_fe_rate,
                   data = caregivers_nonzero, weights = n_patients_sample)

cat("\nModel comparison (AIC):\n")
print(AIC(lm_linear, lm_quadratic))

cat("\nQuadratic vs linear ANOVA:\n")
print(anova(lm_linear, lm_quadratic))

fe_seq  <- seq(0.001, max(caregivers_nonzero$caregiver_fe_rate), by = 0.001)
pred_df <- data.frame(caregiver_fe_rate = fe_seq) %>%
  mutate(
    linear      = predict(lm_linear,    .),
    quadratic   = predict(lm_quadratic, .),
    exponential = exp(predict(lm_log,   .))
  )

ggplot() +
  geom_point(
    data  = caregivers_nonzero,
    aes(caregiver_fe_rate, survival_rate, size = n_patients_sample),
    alpha = 0.5
  ) +
  geom_line(data = pred_df,
            aes(caregiver_fe_rate, linear),
            color = "#2E75B6", linewidth = 0.8) +
  geom_line(data = pred_df,
            aes(caregiver_fe_rate, quadratic),
            color = "#D95F02", linetype = "dashed", linewidth = 0.8) +
  geom_line(data = pred_df,
            aes(caregiver_fe_rate, exponential),
            color = "#1D9E8E", linetype = "dotted", linewidth = 0.8) +
  labs(
    title    = "Caregiver FE Rate vs. 12-Month Patient Survival",
    subtitle = "Caregivers with FE rate > 0. Linear (blue), quadratic (orange), exponential (teal).",
    x        = "Caregiver Failed Extubation Rate",
    y        = "12-Month Patient Survival Rate",
    size     = "N patients"
  ) +
  theme_minimal()

## 5d. Age distribution
ggplot(patient_cohort, aes(anchor_age, fill = tube_event_source)) +
  geom_histogram(bins = 30, position = "dodge", color = "white") +
  scale_fill_manual(values = c("#1D9E8E", "#D95F02"),
                    labels = c("Explicit", "Inferred")) +
  labs(
    title = "Age Distribution by Tube Event Source",
    x     = "Patient Age",
    y     = "Count",
    fill  = "Event Source"
  ) +
  theme_minimal()

## 5e. ICD chapter breakdown
patient_cohort %>%
  count(primary_diagnosis, sort = TRUE) %>%
  mutate(
    pct               = n / sum(n) * 100,
    primary_diagnosis = fct_reorder(primary_diagnosis, n)
  ) %>%
  slice_head(n = 12) %>%
  ggplot(aes(primary_diagnosis, pct)) +
  geom_col(fill = "#2E75B6") +
  coord_flip() +
  labs(
    title = "Primary ICD Chapter Distribution",
    x     = NULL,
    y     = "Percentage of Patients"
  ) +
  theme_minimal()


## ── 6. COHORT FLOW SUMMARY ───────────────────────────────────────────────────

cat("\n── COHORT FLOW ──────────────────────────────────────────────────────────\n")
cat("Raw SQL output:                        ", nrow(last_extubations), "patients\n")
cat("patient_cohort (complete covariates):  ", nrow(patient_cohort), "patients\n")
cat("  of which explicit events:            ",
    sum(patient_cohort$tube_event_source == "explicit"), "\n")
cat("  of which inferred events:            ",
    sum(patient_cohort$tube_event_source == "inferred"), "\n")
cat("explicit_extubations (analytical):     ", nrow(explicit_extubations), "patients\n")
cat("  Unique caregivers:                   ",
    n_distinct(explicit_extubations$caregiver_id), "\n")
cat("  Unique ICU units:                    ",
    n_distinct(explicit_extubations$first_careunit), "\n")
cat("─────────────────────────────────────────────────────────────────────────\n")

## Save cohorts for Scripts 02 and 03
save(patient_cohort, explicit_extubations, caregiver_summary,
     caregivers_nonzero, lm_linear, lm_quadratic, lm_log,
     file = "../data/cohorts.RData")
cat("Cohorts saved to ../data/cohorts.RData\n")
