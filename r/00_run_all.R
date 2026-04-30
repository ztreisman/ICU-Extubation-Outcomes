################################################################################
##
##  ICU Extubation Outcomes — Master Run Script
##
##  Runs the complete analysis pipeline in order.
##  Set working directory to the r/ folder before running.
##
##  Runtime:
##    Script 01: ~2 minutes
##    Script 02: ~30-40 minutes (rstanarm sampling)
##    Script 03: ~5 minutes
##    Total:     ~40 minutes
##
##  Data inputs required (in ../data/):
##    last_extubations.csv      — from ICU_Last_Extubations.sql on BigQuery
##    patient_ccsr_long.csv     — from CCSR long query on BigQuery (see Script 03)
##
##  Intermediate outputs (written to ../data/):
##    cohorts.RData             — cohort objects from Script 01
##    mixed_effects_models.RData — fitted glmmTMB and rstanarm models
##    psm_results.RData         — PSM and sensitivity analysis objects
##
##  Figures written to ../figures/
##
################################################################################

# setwd("~/Projects/Extubate/r")   # adjust as needed

cat("── Script 01: Cohort Construction and Descriptive Analysis ──────────────\n")
source("01_cohort_and_descriptive.R")

cat("\n── Script 02: Mixed Effects Models ──────────────────────────────────────\n")
source("02_logistic_and_mixed_effects.R")

cat("\n── Script 03: PSM and Sensitivity Analysis ──────────────────────────────\n")
source("03_psm_and_sensitivity.R")

cat("\n── Pipeline complete ────────────────────────────────────────────────────\n")
