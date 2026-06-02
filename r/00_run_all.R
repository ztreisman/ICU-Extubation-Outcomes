################################################################################
##
##  ICU Extubation Outcomes — Master Run Script
##
##  Runs the complete R analysis pipeline in order.
##  Set working directory to the r/ folder before running.
##
##  Script summary:
##    01_cohort.R            Cohort construction, caregiver summary, basic figs
##    02_characterization.R  Patient/caregiver strata descriptives (CVICU focus)
##    03_nls_partial_cor.R   NLS exponential model + caregiver partial correlations
##    04_mixed_effects.R     glmmTMB + rstanarm mixed effects models
##    05_msm.R               Continuous treatment MSM (primary causal analysis)
##    06_trajectory.R        Pre-extubation RSBI/P-F trajectory (BigQuery required)
##    07_data_quality.R      Data quality investigations
##    08_build_json.R        Build analysis_outputs.json
##
##  Approximate runtime:
##    01: ~2 min    02: ~1 min    03: ~2 min
##    04: ~35 min (rstanarm)      05: ~5 min
##    06: ~5 min (BigQuery)       07: ~1 min    08: <1 min
##
##  Data inputs required (in ../data/):
##    last_extubations.csv      — from ICU_Last_Extubations.sql on BigQuery
##    patient_ccsr_long.csv     — from CCSR long query on BigQuery (see 05_msm.R)
##
##  Intermediate outputs (written to ../data/):
##    cohorts.RData             — from 01_cohort.R
##    mixed_effects_models.RData — from 04_mixed_effects.R
##    msm_results.RData         — from 05_msm.R
##    nls_results.RData         — from 03_nls_partial_cor.R
##
##  Figures written to ../figures/
##
################################################################################

# setwd("~/Projects/Extubate/r")   # adjust path as needed

cat("── Script 01: Cohort Construction ───────────────────────────────────────\n")
source("01_cohort.R")
# → cohorts.RData (includes analysis_base, psm_ccsr_codes, unit_groups)

cat("\n── Script 02: Patient/Caregiver Characterization ────────────────────────\n")
source("02_characterization.R")
# ← cohorts.RData

cat("\n── Script 03: NLS Model and Partial Correlations ────────────────────────\n")
source("03_nls_partial_cor.R")
# ← cohorts.RData  → nls_results.RData

cat("\n── Script 04: Mixed Effects Models ──────────────────────────────────────\n")
source("04_mixed_effects.R")
# ← cohorts.RData  → mixed_effects_models.RData

cat("\n── Script 05: Continuous Treatment MSM ──────────────────────────────────\n")
source("05_msm.R")
# ← cohorts.RData  → msm_results.RData

cat("\n── Script 06: Trajectory Analysis (requires BigQuery) ───────────────────\n")
source("06_trajectory.R")

cat("\n── Script 07: Data Quality Investigations ───────────────────────────────\n")
source("07_data_quality.R")

cat("\n── Script 08: Build JSON Output ─────────────────────────────────────────\n")
source("08_build_json.R")

cat("\n── Pipeline complete ────────────────────────────────────────────────────\n")
