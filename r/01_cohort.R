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
##  Output objects (saved to cohorts.RData for Scripts 02–05):
##    patient_cohort, explicit_extubations, analysis_base,
##    psm_ccsr_codes, unit_groups, caregiver_summary
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

#setwd("D:/Projects/Extubate/ICU-Extubation-Outcomes/r")   # adjust as needed
setwd("/shared/Projects/Extubate/ICU-Extubation-Outcomes/r")   # adjust as needed


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
    dod                         = col_date(),
    caregiver_unit_n            = col_double()
  )
) %>%
  mutate(
    survival_12mo  = is.na(dod),
    mortality_12mo = 1L - as.integer(is.na(dod)),
    days_to_death  = as.numeric(as.Date(dod) - as.Date(event_time)),
    died_30d       = as.integer(!is.na(dod) & days_to_death <= 30),
    likely_comfort = hospital_expire_flag == 1 & !is.na(dod) & days_to_death <= 3
  )

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
##   - caregiver_n >= 20: caregiver has >= 20 total extubations (quality threshold)
##   - !is.na(caregiver_unit_n): unit-specific FE rate is genuinely available
##     (not a fallback to global rate from a cross-unit imputed assignment)
##   - Complete covariates for modeling

explicit_extubations <- patient_cohort %>%
  filter(
    tube_event_source == "explicit",
    !is.na(caregiver_unit_n),
    caregiver_n >= 20
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

icu_units <- c(
  "Cardiac Vascular Intensive Care Unit (CVICU)",
  "Coronary Care Unit (CCU)",
  "Medical/Surgical Intensive Care Unit (MICU/SICU)",
  "Surgical Intensive Care Unit (SICU)",
  "Neuro Surgical Intensive Care Unit",
  "Medical Intensive Care Unit (MICU)",
  "Trauma SICU (TSICU)",
  "Neuro Intermediate",
  "Neuro Stepdown"
)

explicit_extubations <- explicit_extubations %>%
  filter(first_careunit %in% icu_units)

cat("\nexplicit_extubations N:", nrow(explicit_extubations), "\n")
cat("Unique caregivers:", n_distinct(explicit_extubations$caregiver_id), "\n")
cat("Unique ICU units:", n_distinct(explicit_extubations$first_careunit), "\n")
cat("12-month mortality rate:", round(mean(explicit_extubations$mortality_12mo, na.rm = TRUE), 3), "\n")
cat("Failed extubation rate:", round(mean(explicit_extubations$failed_extubations > 0, na.rm = TRUE), 3), "\n")
cat("Caregiver FE rate — median:", round(median(explicit_extubations$caregiver_fe_rate), 3),
    "IQR:", round(quantile(explicit_extubations$caregiver_fe_rate, 0.25), 3),
    "—", round(quantile(explicit_extubations$caregiver_fe_rate, 0.75), 3), "\n")
cat("30-day mortality rate:", round(mean(explicit_extubations$died_30d, na.rm = TRUE), 3), "\n")
cat("Likely comfort extubations (in-hosp death <= 3d from extubation):",
    sum(explicit_extubations$likely_comfort),
    sprintf("(%.1f%%)\n", mean(explicit_extubations$likely_comfort) * 100))


## ── 3. ANALYSIS BASE WITH CCSR FLAGS ─────────────────────────────────────────
##
##  Joins CCSR diagnostic flags onto explicit_extubations and applies the
##  caregiver_unit_n >= 20 volume threshold. Placed before the descriptive
##  sections so that caregiver_summary and cohort flow counts reference a
##  freshly constructed analysis_base rather than stale in-memory objects.
##  Saved into cohorts.RData so that scripts 02–05 all share a single,
##  consistently constructed analytic cohort without needing to reload
##  patient_ccsr_long.csv.
##
##  CCSR selection: >= 10% prevalence in the analytic cohort, excluding RSP012
##  (respiratory failure, ~47% prevalence — near-universal, near-zero variance).

ccsr_long <- read_csv(
  "../data/patient_ccsr_long.csv",
  col_types = cols(
    subject_id = col_character(),
    hadm_id    = col_character(),
    ccsr_code  = col_character()
  )
)

cat("\nCCSR rows:", nrow(ccsr_long), "\n")
cat("Unique CCSR codes:", n_distinct(ccsr_long$ccsr_code), "\n")

ccsr_prevalence <- ccsr_long |>
  distinct(subject_id, ccsr_code) |>
  count(ccsr_code, name = "n_patients") |>
  mutate(prevalence = n_patients / n_distinct(ccsr_long$subject_id)) |>
  arrange(desc(prevalence))

psm_ccsr_codes <- ccsr_prevalence |>
  filter(prevalence >= 0.10, ccsr_code != "RSP012") |>
  pull(ccsr_code)

cat("CCSR flags (>= 10% prevalence, excl. RSP012):", length(psm_ccsr_codes), "\n")

ccsr_flags <- ccsr_long |>
  distinct(subject_id, ccsr_code) |>
  filter(ccsr_code %in% psm_ccsr_codes) |>
  mutate(present = 1) |>
  pivot_wider(
    id_cols     = subject_id,
    names_from  = ccsr_code,
    values_from = present,
    values_fill = 0
  )

unit_groups <- list(
  CVICU      = "Cardiac Vascular Intensive Care Unit (CVICU)",
  Medical    = c("Medical Intensive Care Unit (MICU)",
                 "Medical/Surgical Intensive Care Unit (MICU/SICU)"),
  SurgTrauma = c("Surgical Intensive Care Unit (SICU)",
                 "Trauma SICU (TSICU)")
)

analysis_base <- explicit_extubations |>
  filter(caregiver_unit_n >= 20) |>
  left_join(ccsr_flags, by = "subject_id") |>
  mutate(across(all_of(psm_ccsr_codes), ~ replace_na(.x, 0))) |>
  mutate(
    unit_group = case_when(
      first_careunit %in% unit_groups$CVICU      ~ "CVICU",
      first_careunit %in% unit_groups$Medical    ~ "Medical",
      first_careunit %in% unit_groups$SurgTrauma ~ "Surgical/Trauma",
      TRUE                                       ~ "Other"
    )
  )

cat("\nanalysis_base N (caregiver_unit_n >= 20):", nrow(analysis_base), "\n")
cat("Unit group sizes:\n")
print(table(analysis_base$unit_group))
cat("12-month mortality rate:", round(mean(analysis_base$mortality_12mo, na.rm = TRUE), 3), "\n")


## ── 4. TABLE 1: PATIENT CHARACTERISTICS ──────────────────────────────────────
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

## table 1c, stratified by first_careunit

tab1c <- CreateTableOne(
  vars       = table1_vars,
  strata     = "first_careunit",
  data       = patient_cohort,
  factorVars = table1_cat,
  addOverall = TRUE
)
print(tab1c, showAllLevels = FALSE, smd = TRUE)

patient_cohort %>%
  group_by(first_careunit) %>%
  summarise(
    median_vent = median(vent_hours, na.rm = TRUE),
    q25_vent    = quantile(vent_hours, 0.25, na.rm = TRUE),
    q75_vent    = quantile(vent_hours, 0.75, na.rm = TRUE),
    .groups     = "drop"
  ) %>%
  arrange(median_vent) %>%
  print()

unit_summary <- patient_cohort %>%
  group_by(first_careunit) %>%
  summarise(
    n              = n(),
    mean_sofa      = mean(sofa, na.rm = TRUE),
    mortality_12mo = mean(mortality_12mo, na.rm = TRUE),
    .groups        = "drop"
  ) %>%
  filter(first_careunit %in% c(
    "Cardiac Vascular Intensive Care Unit (CVICU)",
    "Coronary Care Unit (CCU)",
    "Medical/Surgical Intensive Care Unit (MICU/SICU)",
    "Surgical Intensive Care Unit (SICU)",
    "Neuro Surgical Intensive Care Unit (Neuro SICU)",
    "Medical Intensive Care Unit (MICU)",
    "Trauma SICU (TSICU)",
    "Neuro Intermediate",
    "Neuro Stepdown"
  )) %>%
  mutate(
    label = case_when(
      str_detect(first_careunit, "CVICU")              ~ "CVICU",
      str_detect(first_careunit, "CCU")                ~ "CCU",
      str_detect(first_careunit, "Medical/Surgical")   ~ "MICU/SICU",
      str_detect(first_careunit, "Trauma")             ~ "TSICU",
      str_detect(first_careunit, "Neuro Surgical")     ~ "Neuro SICU",
      str_detect(first_careunit, "Surgical Intensive") ~ "SICU",
      str_detect(first_careunit, "Medical Intensive")  ~ "MICU",
      str_detect(first_careunit, "Intermediate")       ~ "Neuro Intermed.",
      str_detect(first_careunit, "Stepdown")           ~ "Neuro Stepdown",
      TRUE ~ first_careunit
    )
  )

ggplot(unit_summary,
       aes(mean_sofa, mortality_12mo, size = n, label = label)) +
  geom_point(color = "#2E75B6", alpha = 0.7) +
  geom_text(vjust = -0.5, size = 3.2, color = "gray30") +
  scale_size_continuous(
    range  = c(3, 14),
    name   = "N patients",
    breaks = c(500, 1000, 2000, 3000)
  ) +
  scale_y_continuous(
    labels = scales::percent_format(accuracy = 1)
  ) +
  scale_x_continuous(limits = c(2.5, 9.5)) +
  labs(
    title    = "ICU Unit Heterogeneity: Illness Severity vs. 12-Month Mortality",
    x        = "Mean SOFA Score",
    y        = "12-Month Mortality"
  ) +
  theme_minimal() +
  theme(legend.position = "right")

ggsave("../figures/unit_sofa_mortality.png",
       width = 8, height = 4, dpi = 150)

## ── 5. CAREGIVER-LEVEL DESCRIPTIVE ───────────────────────────────────────────
##
##  One row per caregiver. primary_careunit: modal first_careunit among the
##  caregiver's analysis_base patients (most common unit among those meeting
##  the caregiver_unit_n >= 20 threshold). caregiver_summary is built from
##  analysis_base so that the caregiver-level population (NLS, partial
##  correlations in Script 03) is consistent with the MSM in Script 05.

# ── Caregiver cohort flow counts ─────────────────────────────────────────────

# 1. All caregivers with any documented extubation in patient_cohort
n_any <- patient_cohort |>
  filter(!is.na(caregiver_id)) |>
  distinct(caregiver_id) |> nrow()

# 2. Caregivers with at least one EXPLICIT event
n_explicit <- patient_cohort |>
  filter(tube_event_source == "explicit", !is.na(caregiver_id)) |>
  distinct(caregiver_id) |> nrow()

# 3. After caregiver_n >= 20 (already applied in explicit_extubations)
n_vol20 <- explicit_extubations |> distinct(caregiver_id) |> nrow()

# 4. After caregiver_unit_n >= 20 (in analysis_base)
n_unit20 <- analysis_base |> distinct(caregiver_id) |> nrow()

# 5. In primary analysis groups only (excluding Other)
n_primary <- analysis_base |>
  filter(unit_group != "Other") |>
  distinct(caregiver_id) |> nrow()

# 6. Per group
group_n <- analysis_base |>
  filter(unit_group != "Other") |>
  group_by(unit_group) |>
  summarise(n_caregivers = n_distinct(caregiver_id), .groups = "drop")

cat(sprintf("Any caregiver in patient_cohort:   %d\n", n_any))
cat(sprintf("With explicit event:               %d\n", n_explicit))
cat(sprintf("caregiver_n >= 20:                 %d\n", n_vol20))
cat(sprintf("caregiver_unit_n >= 20 (any unit): %d\n", n_unit20))
cat(sprintf("In primary analysis groups:        %d\n", n_primary))
cat("\nBy group:\n")
print(group_n)


primary_careunit_map <- analysis_base %>%
  count(caregiver_id, first_careunit) %>%
  group_by(caregiver_id) %>%
  slice_max(n, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  dplyr::select(caregiver_id, primary_careunit = first_careunit)

caregiver_summary <- analysis_base %>%
  group_by(caregiver_id) %>%
  summarise(
    caregiver_fe_rate  = first(caregiver_fe_rate),
    caregiver_n        = first(caregiver_n),
    survival_rate      = mean(survival_12mo,       na.rm = TRUE),
    mortality_rate     = mean(mortality_12mo,      na.rm = TRUE),
    hospital_expire    = mean(hospital_expire_flag, na.rm = TRUE),
    died_30d_rate      = mean(died_30d,            na.rm = TRUE),
    n_patients_sample  = n(),
    .groups            = "drop"
  ) %>%
  left_join(primary_careunit_map, by = "caregiver_id")

cat("\nCaregivers in caregiver_summary (caregiver_unit_n >= 20):", nrow(caregiver_summary), "\n")
cat("Caregiver FE rate range:",
    round(range(caregiver_summary$caregiver_fe_rate), 3), "\n")
cat("Caregiver volume range:",
    round(range(caregiver_summary$caregiver_n)), "\n")


## ── 5b. UNIT-LEVEL DESCRIPTIVE ────────────────────────────────────────────────
##
##  Per-unit breakdown of patient volume, caregiver FE rate, and survival.
##  Used to assess unit heterogeneity and inform grouping decisions for
##  the stratified analysis in Script 03.

unit_summary <- explicit_extubations %>%
  group_by(first_careunit) %>%
  summarise(
    n_patients        = n(),
    n_caregivers      = n_distinct(caregiver_id),
    mean_fe_rate      = mean(caregiver_fe_rate, na.rm = TRUE),
    median_fe_rate    = median(caregiver_fe_rate, na.rm = TRUE),
    mortality_rate    = mean(mortality_12mo, na.rm = TRUE),
    mean_sofa         = mean(sofa, na.rm = TRUE),
    mean_charlson     = mean(charlson, na.rm = TRUE),
    .groups           = "drop"
  ) %>%
  arrange(mean_fe_rate)

cat("\n── Unit-level summary (explicit cohort) ─────────────────────────────────\n")
print(unit_summary, n = Inf)
cat("────────────────────────────────────────────────────────────────────────\n")


## Caregiver summary

caregiver_summary %>%
  group_by(primary_careunit) %>%
  summarise(
    n_caregivers     = n(),
    fe_median        = round(median(caregiver_fe_rate) * 100, 1),
    fe_q25           = round(quantile(caregiver_fe_rate, 0.25) * 100, 1),
    fe_q75           = round(quantile(caregiver_fe_rate, 0.75) * 100, 1),
    fe_max           = round(max(caregiver_fe_rate) * 100, 1),
    vol_median       = round(median(caregiver_n)),
    vol_q25          = round(quantile(caregiver_n, 0.25)),
    vol_q75          = round(quantile(caregiver_n, 0.75)),
    vol_max          = max(caregiver_n),
    sample_n_median  = round(median(n_patients_sample)),
    mort_mean        = round(mean(mortality_rate) * 100, 1),
    .groups          = "drop"
  ) %>%
  print()

## ── 6. POPULATION VISUALIZATIONS ─────────────────────────────────────────────

## 6a. Caregiver FE rate distribution
ggplot(caregiver_summary, aes(caregiver_fe_rate)) +
  geom_histogram(bins = 30, fill = "#2E75B6", color = "white") +
  labs(
    title    = "Distribution of Caregiver Failed Extubation Rates",
    subtitle = "Caregivers with caregiver_unit_n >= 20",
    x        = "Caregiver FE Rate",
    y        = "Number of Caregivers"
  ) +
  theme_minimal()

## 6b. Caregiver volume vs FE rate
ggplot(caregiver_summary, aes(caregiver_n, caregiver_fe_rate, color = primary_careunit)) +
  geom_point(aes(size = n_patients_sample), alpha = 0.8) +
  labs(
    title    = "Caregiver Volume vs. Failed Extubation Rate",
    subtitle = "Bubble size = patients in analytical sample",
    x        = "Caregiver Total Extubation Volume",
    y        = "Caregiver FE Rate",
    size     = "N patients",
    color    = "Primary Care Unit"
  ) +
  theme_minimal()

ggsave("../figures/caregiver_vol_fe.png",
       width = 8, height = 3, dpi = 150)

## 6c. FE rate vs 12-month mortality (caregiver level, FE rate > 0)
caregivers_nonzero <- caregiver_summary %>%
  filter(caregiver_fe_rate > 0)

# Fit competing functional forms (used in Script 03 NLS AIC comparison)
lm_linear    <- lm(mortality_rate ~ caregiver_fe_rate,
                   data = caregivers_nonzero, weights = n_patients_sample)
lm_quadratic <- lm(mortality_rate ~ caregiver_fe_rate + I(caregiver_fe_rate^2),
                   data = caregivers_nonzero, weights = n_patients_sample)
lm_log       <- lm(log(mortality_rate) ~ caregiver_fe_rate,
                   data    = caregivers_nonzero |> filter(mortality_rate > 0),
                   weights = n_patients_sample)

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
    aes(caregiver_fe_rate, mortality_rate, size = n_patients_sample),
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
    title    = "Caregiver FE Rate vs. 12-Month Patient Mortality",
    subtitle = "Caregivers with FE rate > 0. Linear (blue), quadratic (orange), exponential (teal).",
    x        = "Caregiver Failed Extubation Rate",
    y        = "12-Month Patient Mortality Rate",
    size     = "N patients"
  ) +
  theme_minimal()

## 6d. Age distribution
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

## 6e. ICD chapter breakdown
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


## ── 7. COHORT FLOW SUMMARY ───────────────────────────────────────────────────

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


## ── 8. SAVE ───────────────────────────────────────────────────────────────────

save(patient_cohort, explicit_extubations, caregiver_summary,
     caregivers_nonzero, lm_linear, lm_quadratic, lm_log,
     unit_summary, analysis_base, psm_ccsr_codes, unit_groups,
     file = "../data/cohorts.RData")
cat("Cohorts saved to ../data/cohorts.RData\n")

