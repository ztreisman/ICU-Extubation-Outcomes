################################################################################
##
##  ICU Extubation Outcomes — Script 02: Mixed Effects Models
##
##  Depends on: 01_cohort_and_descriptive.R (run first, or load cohorts.RData)
##
##  Research question:
##    Does caregiver failed extubation rate predict patient 12-month survival
##    after accounting for patient-level illness severity and the three-level
##    clustering structure of the data (patients within caregivers within
##    ICU units)?
##
##  Data: explicit_extubations
##    One row per patient. Explicitly documented extubation events only.
##    Caregivers with > 10 total extubations. Genuine caregiver FE rate
##    (not imputed from overall mean).
##
##  Models:
##    memod_freq  — frequentist mixed effects (glmmTMB)
##                  Random intercepts for caregiver_id and first_careunit
##                  Fixed effects: intubation_type, illness severity covariates,
##                  caregiver_fe_rate
##
##    memod_bayes — Bayesian mixed effects (rstanarm::stan_glmer)
##                  Same structure, weakly informative priors
##                  Scaled predictors to resolve Stan initialization issues
##
##  Key findings:
##    - first_careunit RE variance: ~0.265 (glmmTMB) — ICU unit explains
##      meaningful survival variation independently of patient factors
##    - caregiver_id RE variance: ~0 — no residual caregiver-level variation
##      once unit and patient covariates are controlled
##    - caregiver_fe_rate fixed effect: coef -2.80, p = 0.097 (glmmTMB);
##      posterior mean -0.1, 95% CI (-0.2, 0.0) (rstanarm)
##    - Finding is marginal in the mixed effects framework; the PSM and
##      partial correlation analyses in Script 03 provide stronger evidence
##
##  Runtime note:
##    memod_bayes takes approximately 20-40 minutes (4 chains × 2000 iter,
##    271 caregiver groups). Run 00_run_all.R overnight or run interactively
##    and save the fitted object.
##
################################################################################


## ── 0. LIBRARIES ─────────────────────────────────────────────────────────────

library(MASS)       # load before dplyr
library(dplyr)
library(tidyr)
library(readr)
library(tibble)
library(ggplot2)
library(glmmTMB)
library(rstanarm)

set.seed(237)


## ── 1. LOAD COHORTS ───────────────────────────────────────────────────────────

load("../data/cohorts.RData")
# Loads: patient_cohort, explicit_extubations, caregiver_summary,
#        caregivers_nonzero, lm_linear, lm_quadratic, lm_log

cat("explicit_extubations N:", nrow(explicit_extubations), "\n")
cat("Unique caregivers:", n_distinct(explicit_extubations$caregiver_id), "\n")
cat("Unique ICU units:", n_distinct(explicit_extubations$first_careunit), "\n")


## ── 2. FREQUENTIST MIXED EFFECTS MODEL (glmmTMB) ─────────────────────────────
##
##  Three-level structure:
##    Level 1: patient (outcome: survival_12mo)
##    Level 2: caregiver_id (271 groups)
##    Level 3: first_careunit (12 groups)
##
##  Centered predictors used for interpretability of fixed effects.
##  Centering constants derived from explicit_extubations means.

memod_freq <- glmmTMB(
  survival_12mo ~ intubation_type +
    vent_hours_centered +
    charlson_centered +
    norepinephrine_centered +
    sofa_centered +
    caregiver_fe_rate +
    (1 | caregiver_id) +
    (1 | first_careunit),
  family = binomial,
  data   = explicit_extubations
)

summary(memod_freq)

## Extract key results
freq_fe_coef <- fixef(memod_freq)$cond["caregiver_fe_rate"]
freq_fe_se   <- sqrt(diag(vcov(memod_freq)$cond))["caregiver_fe_rate"]
freq_fe_or   <- exp(freq_fe_coef)
freq_fe_p    <- summary(memod_freq)$coefficients$cond["caregiver_fe_rate", "Pr(>|z|)"]

cat("\n── glmmTMB Results ─────────────────────────────────────────────────────\n")
cat("caregiver_fe_rate coefficient:", round(freq_fe_coef, 3), "\n")
cat("caregiver_fe_rate OR:         ", round(freq_fe_or,   3), "\n")
cat("caregiver_fe_rate p-value:    ", round(freq_fe_p,    4), "\n")
cat("\nRandom effect variances:\n")
cat("  caregiver_id:   ", round(VarCorr(memod_freq)$cond$caregiver_id[1], 6), "\n")
cat("  first_careunit: ", round(VarCorr(memod_freq)$cond$first_careunit[1], 3), "\n")
cat("────────────────────────────────────────────────────────────────────────\n")


## ── 3. BAYESIAN MIXED EFFECTS MODEL (rstanarm) ───────────────────────────────
##
##  Same structure as memod_freq but with:
##    - Scaled predictors (mean 0, SD 1) to resolve Stan initialization failures
##      caused by predictors on very different scales
##    - Weakly informative priors: normal(0, 2.5) on fixed effects,
##      decov(scale = 1) on random effect covariance
##
##  The Bayesian approach is particularly informative here because:
##    - With 271 caregivers but small per-caregiver sample sizes, the posterior
##      on caregiver RE variance is appropriately regularized toward zero rather
##      than producing a glmmTMB boundary estimate of ~0
##    - Posterior credible intervals on caregiver_fe_rate provide uncertainty
##      quantification that complements the frequentist p-value
##
##  Convergence diagnostics: all Rhat = 1.0, n_eff > 1600

memod_bayes <- stan_glmer(
  survival_12mo ~ intubation_type +
    vent_hours_scaled +
    charlson_scaled +
    norepinephrine_scaled +
    sofa_scaled +
    caregiver_fe_scaled +
    (1 | caregiver_id) +
    (1 | first_careunit),
  family           = binomial,
  data             = explicit_extubations,
  prior            = normal(0, 2.5),
  prior_covariance = decov(scale = 1),
  chains           = 4,
  iter             = 2000,
  seed             = 237
)

## Summary of key parameters
bayes_summary <- summary(
  memod_bayes,
  pars  = c("caregiver_fe_scaled",
            "Sigma[caregiver_id:(Intercept),(Intercept)]",
            "Sigma[first_careunit:(Intercept),(Intercept)]"),
  probs = c(0.025, 0.975)
)
print(bayes_summary)

## Extract posterior for caregiver_fe_scaled
bayes_post <- as.data.frame(memod_bayes)
fe_post    <- bayes_post[["caregiver_fe_scaled"]]

cat("\n── rstanarm Results ────────────────────────────────────────────────────\n")
cat("caregiver_fe_scaled posterior mean:", round(mean(fe_post), 3), "\n")
cat("caregiver_fe_scaled 95% CI:",
    round(quantile(fe_post, 0.025), 3), "to",
    round(quantile(fe_post, 0.975), 3), "\n")
cat("P(caregiver_fe_scaled < 0):", round(mean(fe_post < 0), 3), "\n")
cat("────────────────────────────────────────────────────────────────────────\n")

## Posterior density plot for caregiver_fe_scaled
ggplot(data.frame(fe_post), aes(fe_post)) +
  geom_density(fill = "#2E75B6", alpha = 0.4) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray40") +
  geom_vline(xintercept = mean(fe_post), color = "#2E75B6", linewidth = 0.8) +
  labs(
    title    = "Posterior Distribution: caregiver_fe_scaled",
    subtitle = paste0("Mean = ", round(mean(fe_post), 3),
                      ", 95% CI (", round(quantile(fe_post, 0.025), 3),
                      ", ", round(quantile(fe_post, 0.975), 3), ")"),
    x        = "Posterior estimate (standardized scale)",
    y        = "Density"
  ) +
  theme_minimal()


## ── 4. MODEL COMPARISON TABLE ─────────────────────────────────────────────────

cat("\n── Mixed Effects Model Comparison ──────────────────────────────────────\n")
cat(sprintf("%-45s %-10s %-15s %-8s\n",
            "Model", "Estimate", "95% CI", "p / P(<0)"))
cat(strrep("─", 80), "\n")
cat(sprintf("%-45s %-10s %-15s %-8s\n",
            "glmmTMB (caregiver_fe_rate, log-odds)",
            round(freq_fe_coef, 3),
            "—",
            round(freq_fe_p, 4)))
cat(sprintf("%-45s %-10s %-15s %-8s\n",
            "rstanarm (caregiver_fe_scaled, posterior)",
            round(mean(fe_post), 3),
            paste0("(", round(quantile(fe_post, 0.025), 3), ", ",
                   round(quantile(fe_post, 0.975), 3), ")"),
            paste0(round(mean(fe_post < 0) * 100, 1), "%")))
cat(strrep("─", 80), "\n")
cat("\nRandom effect variances:\n")
cat(sprintf("  %-20s glmmTMB: %-10s rstanarm posterior mean: %-10s\n",
            "caregiver_id",
            round(VarCorr(memod_freq)$cond$caregiver_id[1], 6),
            round(mean(bayes_post[["Sigma[caregiver_id:(Intercept),(Intercept)]"]]), 4)))
cat(sprintf("  %-20s glmmTMB: %-10s rstanarm posterior mean: %-10s\n",
            "first_careunit",
            round(VarCorr(memod_freq)$cond$first_careunit[1], 3),
            round(mean(bayes_post[["Sigma[first_careunit:(Intercept),(Intercept)]"]]), 3)))
cat("────────────────────────────────────────────────────────────────────────\n")


## ── 5. SAVE FITTED MODELS ─────────────────────────────────────────────────────
##
##  Save both models so they can be loaded in 00_run_all.R without re-fitting.
##  memod_bayes is particularly expensive (~20-40 min).

save(memod_freq, memod_bayes,
     file = "../data/mixed_effects_models.RData")
cat("Models saved to ../data/mixed_effects_models.RData\n")
