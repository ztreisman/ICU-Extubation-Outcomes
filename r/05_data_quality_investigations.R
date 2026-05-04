################################################################################
##
##  ICU Extubation Outcomes — Script 05: Data Quality Investigations
##
##  Depends on: 01_cohort_and_descriptive.R (cohorts.RData)
##
##  This script documents the investigation that established the correct
##  interpretation of the `failed_extubations` variable in last_extubations.csv.
##
##  Key finding:
##    `failed_extubations` in the last_extubations output is a patient lifetime
##    count of reintubation pairs across ALL extubation events during the
##    hospital stay — not a flag for whether the LAST extubation specifically
##    failed. Evidence:
##
##    1. Patients with failed_extubations > 1 have mean vent_hours of ~167h
##       (vs ~57h for failed_extubations = 0), consistent with multiple complete
##       intubation-extubation cycles rather than multiple failures after a
##       single extubation.
##
##    2. Among patients with failed_extubations > 1, 61% have caregiver_fe_rate
##       < 5% — if failed_extubations were a last-extubation outcome, you would
##       expect these patients to have caregivers with HIGH FE rates, not low.
##
##    3. The internal consistency check (patients with failed_extubations > 0
##       have slightly higher caregiver_fe_rate: 4.68% vs 3.97%) is in the right
##       direction but the difference is small, consistent with failed_extubations
##       being a complexity marker rather than a direct outcome.
##
##  Implication:
##    `failed_extubations` should be treated as a patient complexity covariate
##    (number of prior intubation cycles) rather than the extubation outcome.
##    The correct per-event outcome is `failed_extubation_flag` from
##    all_extubations.sql, which flags whether each specific extubation was
##    followed by reintubation within 72h.
##
##  This script is diagnostic/exploratory — it does not produce output objects
##  used by other scripts.
##
################################################################################


## ── 0. SETUP ─────────────────────────────────────────────────────────────────

library(dplyr)

load("../data/cohorts.RData")
# Loads: patient_cohort, explicit_extubations


## ── 1. FAILED_EXTUBATIONS DISTRIBUTION ───────────────────────────────────────

cat("failed_extubations distribution:\n")
print(table(patient_cohort$failed_extubations))

patient_cohort %>%
  summarise(
    mean_failed_flag    = mean(failed_extubations > 0),
    mean_fe_rate        = mean(caregiver_fe_rate, na.rm = TRUE),
    pct_zero_fe         = mean(caregiver_fe_rate == 0, na.rm = TRUE),
    n_reintubations     = sum(failed_extubations),
    n_patients          = n()
  ) %>%
  print()


## ── 2. FAILED_EXTUBATIONS BY TUBE EVENT SOURCE ───────────────────────────────
##
##  Tests whether failed_extubations rate differs between explicit and inferred
##  events — if it did, it would suggest differential detection of failures
##  rather than a real outcome difference.

patient_cohort %>%
  group_by(tube_event_source, failed_extubations > 0) %>%
  summarise(n = n(), .groups = "drop") %>%
  print()


## ── 3. MULTIPLE FAILURE PATIENTS ─────────────────────────────────────────────

cat("\nPatients with failed_extubations > 1:", 
    sum(patient_cohort$failed_extubations > 1), "\n")

patient_cohort %>%
  summarise(
    pct_with_any_failure = mean(failed_extubations > 0),
    mean_failures        = mean(failed_extubations),
    pct_multiple_failure = mean(failed_extubations > 1)
  ) %>%
  print()


## ── 4. KEY DIAGNOSTIC: VENT HOURS BY FAILED_EXTUBATIONS ─────────────────────
##
##  If failed_extubations is a lifetime count of intubation cycles, patients
##  with higher counts should have longer cumulative vent hours.
##  If it were a per-last-extubation count, vent hours would not be expected
##  to increase monotonically.

cat("\nVent hours and SOFA by failed_extubations count:\n")
explicit_extubations %>%
  group_by(failed_extubations) %>%
  summarise(
    mean_vent_hours = mean(vent_hours),
    mean_sofa       = mean(sofa),
    n               = n(),
    .groups         = "drop"
  ) %>%
  print()

# Expected pattern if lifetime count: mean_vent_hours increases with count
# (longer total time on ventilator due to multiple cycles)


## ── 5. CAREGIVER FE RATE AMONG MULTI-FAILURE PATIENTS ────────────────────────
##
##  If failed_extubations > 1 were a last-extubation outcome, patients with
##  multiple failures should be concentrated among high-FE-rate caregivers.
##  If it is a complexity marker, caregiver FE rate should be near the
##  population average.

cat("\nCaregiver FE rate for patients with failed_extubations > 1:\n")
explicit_extubations %>%
  filter(failed_extubations > 1) %>%
  summarise(
    mean_caregiver_fe_rate = mean(caregiver_fe_rate),
    pct_low_fe             = mean(caregiver_fe_rate < 0.05),
    n                      = n()
  ) %>%
  print()

# 61% of these patients have caregivers with FE rate < 5% — inconsistent
# with failed_extubations measuring a last-extubation outcome.


## ── 6. INTERNAL CONSISTENCY CHECK ────────────────────────────────────────────
##
##  Do patients with any failed_extubation have higher caregiver_fe_rate?
##  Direction should be positive if the measure has validity.

cat("\nCaregiver FE rate by any-failure status:\n")
explicit_extubations %>%
  mutate(had_failure = failed_extubations > 0) %>%
  group_by(had_failure) %>%
  summarise(
    mean_caregiver_fe_rate = mean(caregiver_fe_rate),
    n                      = n(),
    .groups                = "drop"
  ) %>%
  print()

# Small positive difference (4.68% vs 3.97%) — right direction but modest,
# consistent with failed_extubations being a complexity marker that is
# correlated with caregiver FE rate but not the primary outcome variable.

cat("\nConclusion: failed_extubations is a patient complexity marker\n")
cat("(lifetime intubation cycle count), not a last-extubation outcome.\n")
cat("Use failed_extubation_flag from all_extubations.sql for outcome analyses.\n")
