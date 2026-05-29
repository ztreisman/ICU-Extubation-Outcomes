# 06_build_json.R
# Build structured analysis outputs for document verification by Claude Code

library(jsonlite)
library(dplyr)

# ── Parameters ────────────────────────────────────────────────────────────────
comfort_threshold_days  <- 3   # death within N days of extubation = comfort case

# ── Derived ───────────────────────────────────────────────────────────────────
comfort_threshold_days  <- 3

explicit_extubations <- explicit_extubations %>%
  mutate(
    comfort_extubation = hospital_expire_flag == 1 &
      !is.na(days_to_death) &
      days_to_death <= comfort_threshold_days
  )

n_comfort_total <- sum(explicit_extubations$comfort_extubation, na.rm = TRUE)
pct_comfort     <- n_comfort_total / nrow(explicit_extubations)

n_cvicu         <- sum(explicit_extubations$first_careunit ==
                         "Cardiac Vascular Intensive Care Unit (CVICU)")
pct_comfort_cvicu <- n_comfort_cvicu / n_cvicu

# ── Trajectory model ──────────────────────────────────────────────────────────
rsbi_glm       <- glm(
  failed ~ rsbi_at_extub + rsbi_slope,
  data   = event_slopes %>% filter(!is.na(rsbi_slope), !is.na(rsbi_at_extub)),
  family = binomial
)
rsbi_coefs     <- coef(summary(rsbi_glm))

# ── Build output ──────────────────────────────────────────────────────────────
analysis_outputs <- list(
  
  parameters = list(
    comfort_threshold_days  = comfort_threshold_days,
    comfort_threshold_hours = comfort_threshold_hours
  ),
  
  cohort = list(
    patient_cohort_n             = nrow(patient_cohort),
    explicit_extubations_n       = nrow(explicit_extubations),
    n_caregivers                 = n_distinct(explicit_extubations$caregiver_id),
    n_icu_units                  = n_distinct(explicit_extubations$first_careunit),
    survival_12mo_explicit       = round(mean(explicit_extubations$survival_12mo,       na.rm = TRUE), 3),
    survival_12mo_patient        = round(mean(patient_cohort$survival_12mo,             na.rm = TRUE), 3),
    hospital_mortality_explicit  = round(mean(explicit_extubations$hospital_expire_flag, na.rm = TRUE), 3),
    hospital_mortality_patient   = round(mean(patient_cohort$hospital_expire_flag,      na.rm = TRUE), 3),
    comfort_extubations_n        = n_comfort_total,
    comfort_extubations_pct      = round(pct_comfort, 3),
    fe_rate_median               = round(median(  explicit_extubations$caregiver_fe_rate,        na.rm = TRUE), 3),
    fe_rate_q25                  = round(quantile(explicit_extubations$caregiver_fe_rate, 0.25,  na.rm = TRUE), 3),
    fe_rate_q75                  = round(quantile(explicit_extubations$caregiver_fe_rate, 0.75,  na.rm = TRUE), 3),
    fe_rate_max                  = round(max(     explicit_extubations$caregiver_fe_rate,         na.rm = TRUE), 3),
    mortality_30day_explicit     = round(mean(explicit_extubations$died_30d,             na.rm = TRUE), 3)
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
      survival_12mo      = round(mean(survival_12mo,         na.rm = TRUE), 3),
      hospital_mortality = round(mean(hospital_expire_flag,  na.rm = TRUE), 3),
      .groups            = "drop"
    ) %>%
    split(.$unit),
  
  psm_results = lapply(psm_results_list, function(x) list(
    group      = x$group,
    n_matched  = x$n_psm,
    fe_q25     = round(as.numeric(x$fe_q25), 3),
    fe_q75     = round(as.numeric(x$fe_q75), 3),
    surv_high  = round(x$surv_high, 3),
    surv_low   = round(x$surv_low,  3),
    or_unadj   = round(as.numeric(x$or_unadj), 3),
    ci_unadj   = round(as.numeric(x$ci_unadj), 3),
    p_unadj    = round(x$p_unadj, 4),
    or_adj     = round(as.numeric(x$or_adj),  3),
    ci_adj     = round(as.numeric(x$ci_adj),  3),
    p_adj      = round(x$p_adj,  4),
    or_ipw     = round(as.numeric(x$or_ipw),  3),
    ci_ipw     = round(as.numeric(x$ci_ipw),  3),
    p_ipw      = round(x$p_ipw,  4)
  )),
  
  comfort_sensitivity = list(
    comfort_threshold_days  = comfort_threshold_days,
    n_comfort_cvicu         = n_comfort_cvicu,
    pct_comfort_cvicu       = round(pct_comfort_cvicu, 3),
    psm_or_with_comfort     = round(as.numeric(psm_results_list$CVICU$or_adj),  3),
    psm_ci_with_comfort     = round(as.numeric(psm_results_list$CVICU$ci_adj),  3),
    psm_p_with_comfort      = round(psm_results_list$CVICU$p_adj,  4),
    psm_or_no_comfort       = round(as.numeric(cvicu_noc$or_adj),  3),
    psm_ci_no_comfort       = round(as.numeric(cvicu_noc$ci_adj),  3),
    psm_p_no_comfort        = round(cvicu_noc$p_adj,  4),
    ipw_or_with_comfort     = round(as.numeric(psm_results_list$CVICU$or_ipw), 3),
    ipw_p_with_comfort      = round(psm_results_list$CVICU$p_ipw, 4),
    ipw_or_no_comfort       = round(as.numeric(cvicu_noc$or_ipw), 3),
    ipw_p_no_comfort        = round(cvicu_noc$p_ipw, 4)
  ),
  
  quartile_survival = group_quartile_survival %>%
    mutate(across(where(is.numeric), ~ round(., 3))) %>%
    split(.$unit_group),
  
  nls_fit = list(
    b_estimate           = round(coef(nls_fit)["b"],  3),
    b_ci_lower           = round(b_ci[1], 3),
    b_ci_upper           = round(b_ci[2], 3),
    k_estimate           = round(coef(nls_fit)["k"],  1),
    k_ci_lower           = round(k_ci[1], 1),
    k_ci_upper           = round(k_ci[2], 1),
    n_caregivers_nonzero = nrow(caregivers_nonzero)
  ),
  
  mixed_effects = list(
    glmmTMB_fe_coef      = round(freq_fe_coef, 3),
    glmmTMB_fe_or        = round(freq_fe_or,   3),
    glmmTMB_fe_p         = round(freq_fe_p,    4),
    re_var_caregiver_glmmTMB  = 0,
    re_var_careunit_glmmTMB   = 0.271,
    rstanarm_fe_mean     = round(bayes_summary["caregiver_fe_scaled", "mean"],  3),
    rstanarm_fe_sd       = round(bayes_summary["caregiver_fe_scaled", "sd"],    3),
    rstanarm_fe_ci_lo    = round(bayes_summary["caregiver_fe_scaled", "2.5%"],  3),
    rstanarm_fe_ci_hi    = round(bayes_summary["caregiver_fe_scaled", "97.5%"], 3),
    rstanarm_p_lt0       = round(mean(bayes_post < 0), 3),
    re_var_caregiver_stan = round(bayes_summary[
      "Sigma[caregiver_id:(Intercept),(Intercept)]", "mean"], 4),
    re_var_careunit_stan  = round(bayes_summary[
      "Sigma[first_careunit:(Intercept),(Intercept)]", "mean"], 3)
  ),
  
  partial_correlation = list(
    
    # Global partial correlations (n >= 5, 134 caregivers)
    global_r_fe_controlling_vol  = round(pcor_fe$estimate,  3),
    global_p_fe_controlling_vol  = round(pcor_fe$p.value,   4),
    global_n                     = pcor_fe$n,
    global_r_vol_controlling_fe  = round(pcor_vol$estimate, 3),
    global_p_vol_controlling_fe  = round(pcor_vol$p.value,  4),
    
    # Unit-stratified partial correlations
    by_unit = unit_pcor_results %>%
      mutate(
        primary_careunit = as.character(primary_careunit),
        across(where(is.numeric), ~ round(., 3))
      ) %>%
      split(.$primary_careunit)
  ),
  
  trajectory = list(
    n_events_rsbi     = nrow(event_slopes %>% filter(!is.na(rsbi_slope))),
    n_failed_rsbi     = sum(event_slopes$failed,                                         na.rm = TRUE),
    rsbi_mean_failed  = round(mean(event_slopes$rsbi_at_extub[ event_slopes$failed],     na.rm = TRUE), 1),
    rsbi_mean_success = round(mean(event_slopes$rsbi_at_extub[!event_slopes$failed],     na.rm = TRUE), 1),
    rsbi_level_coef   = round(rsbi_coefs["rsbi_at_extub", "Estimate"],   4),
    rsbi_level_se     = round(rsbi_coefs["rsbi_at_extub", "Std. Error"], 4),
    rsbi_level_p      = round(rsbi_coefs["rsbi_at_extub", "Pr(>|z|)"],   4),
    rsbi_slope_coef   = round(rsbi_coefs["rsbi_slope",    "Estimate"],   4),
    rsbi_slope_se     = round(rsbi_coefs["rsbi_slope",    "Std. Error"], 4),
    rsbi_slope_p      = round(rsbi_coefs["rsbi_slope",    "Pr(>|z|)"],   4)
  )
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
cat("Comfort threshold:", comfort_threshold_days, "days (",
    comfort_threshold_hours, "h)\n")
cat("Comfort extubations (total):", n_comfort_total,
    "(", round(pct_comfort * 100, 1), "%)\n")
cat("Comfort extubations (CVICU):", n_comfort_cvicu,
    "(", round(pct_comfort_cvicu * 100, 1), "%)\n")