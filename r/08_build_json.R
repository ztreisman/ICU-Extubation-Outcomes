# 08_build_json.R
# Build structured analysis outputs for document verification.
#
# Depends on (load order):
#   01_cohort.R        → cohorts.RData
#   04_mixed_effects.R → mixed_effects_models.RData
#   05_msm.R           → msm_results.RData
#   03_nls_partial_cor.R → nls_results.RData
#   06_trajectory.R    → trajectory data (event_slopes in environment)
#
# All RData files are loaded below; trajectory section requires event_slopes
# from a prior source("06_trajectory.R") call or BigQuery session.

library(jsonlite)
library(dplyr)

load("../data/cohorts.RData")
# patient_cohort, explicit_extubations, caregiver_summary, caregivers_nonzero,
# lm_linear, lm_quadratic, lm_log, unit_summary

load("../data/msm_results.RData")
# msm_results_list, msm_summary_tbl, group_quartile_survival,
# cvicu_noc, analysis_base, psm_ccsr_codes

load("../data/nls_results.RData")
# pcor_fe, pcor_vol, unit_pcor_results,
# nls_fit, m_max_hat, k_hat, m_max_ci, k_ci, nls_results

load("../data/mixed_effects_models.RData")
# memod_freq, memod_bayes, freq_fe_coef, freq_fe_or, freq_fe_p,
# bayes_summary, bayes_post

# ── Parameters ────────────────────────────────────────────────────────────────
comfort_threshold_days  <- 3
comfort_threshold_hours <- comfort_threshold_days * 24

# ── Derived ───────────────────────────────────────────────────────────────────
explicit_extubations <- explicit_extubations %>%
  mutate(
    comfort_extubation = hospital_expire_flag == 1 &
      !is.na(days_to_death) &
      days_to_death <= comfort_threshold_days
  )

n_comfort_total   <- sum(explicit_extubations$comfort_extubation, na.rm = TRUE)
pct_comfort       <- n_comfort_total / nrow(explicit_extubations)

n_cvicu           <- sum(explicit_extubations$first_careunit ==
                           "Cardiac Vascular Intensive Care Unit (CVICU)")
n_comfort_cvicu   <- sum(
  explicit_extubations$comfort_extubation &
    explicit_extubations$first_careunit ==
      "Cardiac Vascular Intensive Care Unit (CVICU)",
  na.rm = TRUE
)
pct_comfort_cvicu <- n_comfort_cvicu / n_cvicu

icu_units <- c(
  "Cardiac Vascular Intensive Care Unit (CVICU)",
  "Coronary Care Unit (CCU)",
  "Medical/Surgical Intensive Care Unit (MICU/SICU)",
  "Surgical Intensive Care Unit (SICU)",
  "Neuro Surgical Intensive Care Unit (Neuro SICU)",
  "Medical Intensive Care Unit (MICU)",
  "Trauma SICU (TSICU)",
  "Neuro Intermediate",
  "Neuro Stepdown"
)

# ── Trajectory model (requires event_slopes from 06_trajectory.R) ────────────
trajectory_available <- exists("event_slopes")
if (trajectory_available) {
  rsbi_glm   <- glm(
    failed ~ rsbi_at_extub + rsbi_slope,
    data   = event_slopes %>% filter(!is.na(rsbi_slope), !is.na(rsbi_at_extub)),
    family = binomial
  )
  rsbi_coefs <- coef(summary(rsbi_glm))
}

# ── Build output ──────────────────────────────────────────────────────────────
analysis_outputs <- list(

  parameters = list(
    comfort_threshold_days  = comfort_threshold_days,
    comfort_threshold_hours = comfort_threshold_hours,
    min_caregiver_unit_n    = 20L,
    min_caregivers_per_unit = 20L
  ),

  cohort = list(
    patient_cohort_n             = nrow(patient_cohort),
    explicit_extubations_n       = nrow(explicit_extubations),
    n_caregivers                 = n_distinct(explicit_extubations$caregiver_id),
    n_icu_units                  = n_distinct(explicit_extubations$first_careunit),
    mortality_12mo_explicit      = round(mean(explicit_extubations$mortality_12mo,      na.rm = TRUE), 3),
    mortality_12mo_patient       = round(mean(patient_cohort$mortality_12mo,            na.rm = TRUE), 3),
    hospital_mortality_explicit  = round(mean(explicit_extubations$hospital_expire_flag, na.rm = TRUE), 3),
    hospital_mortality_patient   = round(mean(patient_cohort$hospital_expire_flag,      na.rm = TRUE), 3),
    comfort_extubations_n        = n_comfort_total,
    comfort_extubations_pct      = round(pct_comfort, 3),
    fe_rate_median               = round(median(  explicit_extubations$caregiver_fe_rate,       na.rm = TRUE), 3),
    fe_rate_q25                  = round(quantile(explicit_extubations$caregiver_fe_rate, 0.25, na.rm = TRUE), 3),
    fe_rate_q75                  = round(quantile(explicit_extubations$caregiver_fe_rate, 0.75, na.rm = TRUE), 3),
    fe_rate_max                  = round(max(     explicit_extubations$caregiver_fe_rate,        na.rm = TRUE), 3),
    mortality_30day_explicit     = round(mean(explicit_extubations$died_30d,            na.rm = TRUE), 3)
  ),

  unit_profiles = patient_cohort %>%
    filter(first_careunit %in% icu_units) %>%
    group_by(unit = as.character(first_careunit)) %>%
    summarise(
      n                  = n(),
      age_mean           = round(mean(anchor_age,           na.rm = TRUE), 1),
      age_sd             = round(sd(anchor_age,             na.rm = TRUE), 1),
      sofa_mean          = round(mean(sofa,                 na.rm = TRUE), 1),
      sofa_sd            = round(sd(sofa,                   na.rm = TRUE), 1),
      charlson_mean      = round(mean(charlson,             na.rm = TRUE), 1),
      charlson_sd        = round(sd(charlson,               na.rm = TRUE), 1),
      vent_median        = round(median(vent_hours,          na.rm = TRUE), 1),
      vent_q25           = round(quantile(vent_hours, 0.25,  na.rm = TRUE), 1),
      vent_q75           = round(quantile(vent_hours, 0.75,  na.rm = TRUE), 1),
      mortality_12mo     = round(mean(mortality_12mo,        na.rm = TRUE), 3),
      hospital_mortality = round(mean(hospital_expire_flag,  na.rm = TRUE), 3),
      .groups            = "drop"
    ) %>%
    split(.$unit),

  msm_results = lapply(msm_results_list, function(x) {
    if (is.null(x)) return(NULL)
    list(
      group         = x$group,
      n_patients    = x$n_patients,
      n_caregivers  = x$n_caregivers,
      eff_n         = round(x$eff_n, 0),
      fe_mean       = round(x$fe_mean, 4),
      fe_median     = round(x$fe_median, 4),
      fe_iqr        = round(x$fe_iqr, 4),
      r2_denominator = round(x$r2_denominator, 3),
      or_per_5pp    = round(x$or_per_5pp, 3),
      or_lo_cl      = round(x$or_lo_cl, 3),
      or_hi_cl      = round(x$or_hi_cl, 3),
      p_cl          = round(x$p_cl, 4),
      or_lo_hc3     = round(x$or_lo_hc3, 3),
      or_hi_hc3     = round(x$or_hi_hc3, 3),
      p_hc3         = round(x$p_hc3, 4),
      secondary_outcomes = x$secondary_outcomes
    )
  }),

  comfort_sensitivity = list(
    comfort_threshold_days  = comfort_threshold_days,
    n_comfort_cvicu         = n_comfort_cvicu,
    pct_comfort_cvicu       = round(pct_comfort_cvicu, 3),
    msm_or_with_comfort     = round(msm_results_list$CVICU$or_per_5pp, 3),
    msm_ci_with_comfort     = round(c(msm_results_list$CVICU$or_lo_cl,
                                      msm_results_list$CVICU$or_hi_cl), 3),
    msm_p_with_comfort      = round(msm_results_list$CVICU$p_cl, 4),
    msm_or_no_comfort       = round(cvicu_noc$or_per_5pp, 3),
    msm_ci_no_comfort       = round(c(cvicu_noc$or_lo_cl,
                                      cvicu_noc$or_hi_cl), 3),
    msm_p_no_comfort        = round(cvicu_noc$p_cl, 4)
  ),

  quartile_mortality = group_quartile_survival %>%
    mutate(across(where(is.numeric), ~ round(., 3))) %>%
    split(.$unit_group),

  nls_fit = list(
    m_max_estimate       = round(m_max_hat, 3),
    m_max_ci_lower       = round(m_max_ci[1], 3),
    m_max_ci_upper       = round(m_max_ci[2], 3),
    k_estimate           = round(k_hat, 1),
    k_ci_lower           = round(k_ci[1], 1),
    k_ci_upper           = round(k_ci[2], 1),
    n_caregivers_nonzero = nrow(caregivers_nonzero)
  ),

  mixed_effects = list(
    glmmTMB_fe_coef           = round(freq_fe_coef, 3),
    glmmTMB_fe_or             = round(freq_fe_or,   3),
    glmmTMB_fe_p              = round(freq_fe_p,    4),
    rstanarm_fe_mean          = round(bayes_summary["caregiver_fe_scaled", "mean"],  3),
    rstanarm_fe_sd            = round(bayes_summary["caregiver_fe_scaled", "sd"],    3),
    rstanarm_fe_ci_lo         = round(bayes_summary["caregiver_fe_scaled", "2.5%"],  3),
    rstanarm_fe_ci_hi         = round(bayes_summary["caregiver_fe_scaled", "97.5%"], 3),
    rstanarm_p_lt0            = round(mean(bayes_post < 0), 3),
    re_var_caregiver_stan     = round(bayes_summary[
      "Sigma[caregiver_id:(Intercept),(Intercept)]", "mean"], 4),
    re_var_careunit_stan      = round(bayes_summary[
      "Sigma[first_careunit:(Intercept),(Intercept)]", "mean"], 3)
  ),

  partial_correlation = list(
    global_r_fe_controlling_vol  = round(pcor_fe$estimate,  3),
    global_p_fe_controlling_vol  = round(pcor_fe$p.value,   4),
    global_n                     = pcor_fe$n,
    global_r_vol_controlling_fe  = round(pcor_vol$estimate, 3),
    global_p_vol_controlling_fe  = round(pcor_vol$p.value,  4),
    by_unit = unit_pcor_results %>%
      mutate(
        primary_careunit = as.character(primary_careunit),
        across(where(is.numeric), ~ round(., 3))
      ) %>%
      split(.$primary_careunit)
  ),

  trajectory = if (trajectory_available) list(
    n_events_rsbi     = nrow(event_slopes %>% filter(!is.na(rsbi_slope))),
    n_failed_rsbi     = sum(event_slopes$failed,                                       na.rm = TRUE),
    rsbi_mean_failed  = round(mean(event_slopes$rsbi_at_extub[ event_slopes$failed],   na.rm = TRUE), 1),
    rsbi_mean_success = round(mean(event_slopes$rsbi_at_extub[!event_slopes$failed],   na.rm = TRUE), 1),
    rsbi_level_coef   = round(rsbi_coefs["rsbi_at_extub", "Estimate"],   4),
    rsbi_level_se     = round(rsbi_coefs["rsbi_at_extub", "Std. Error"], 4),
    rsbi_level_p      = round(rsbi_coefs["rsbi_at_extub", "Pr(>|z|)"],   4),
    rsbi_slope_coef   = round(rsbi_coefs["rsbi_slope",    "Estimate"],   4),
    rsbi_slope_se     = round(rsbi_coefs["rsbi_slope",    "Std. Error"], 4),
    rsbi_slope_p      = round(rsbi_coefs["rsbi_slope",    "Pr(>|z|)"],   4)
  ) else list(note = "event_slopes not available — run 06_trajectory.R first")
)

# ── Save ──────────────────────────────────────────────────────────────────────
write_json(
  analysis_outputs,
  "../data/analysis_outputs.json",
  pretty     = TRUE,
  auto_unbox = TRUE
)

cat("Saved: ../data/analysis_outputs.json\n")
cat("Top-level keys:", paste(names(analysis_outputs), collapse = ", "), "\n")
cat("Comfort threshold:", comfort_threshold_days, "days\n")
cat("Comfort extubations (total):", n_comfort_total,
    "(", round(pct_comfort * 100, 1), "%)\n")
cat("Comfort extubations (CVICU):", n_comfort_cvicu,
    "(", round(pct_comfort_cvicu * 100, 1), "%)\n")
