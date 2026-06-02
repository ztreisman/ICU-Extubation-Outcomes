################################################################################
##
##  ICU Extubation Outcomes — Script 05: Continuous Treatment MSM
##
##  Depends on: 01_cohort.R (cohorts.RData), patient_ccsr_long.csv
##  Also see:   02_characterization.R, 03_nls_partial_cor.R
##
##  Marginal structural model with stabilized Hirano-Imbens density weights:
##    sw_i = f(T_i) / f(T_i | X_i)
##  where T = caregiver_fe_rate (continuous), X = patient covariates.
##  Densities from OLS-implied normal distributions; weights trimmed at 99th
##  percentile; cluster-robust SEs by caregiver_id.
##
##  Primary outcome: 12-month mortality. Restricts to caregivers with
##  caregiver_unit_n >= 20 and unit groups with >= 20 eligible caregivers.
##
##  Limitation: zero-inflation in caregiver_fe_rate (see 02_characterization.R
##  for the volume confounding analysis). The Python DML (06_dml_cvicu.py)
##  provides a non-parametric complement without the normal density assumption.
##
################################################################################


## ── 0. LIBRARIES ─────────────────────────────────────────────────────────────

library(dplyr)
library(purrr)
library(tibble)
library(ggplot2)
library(sandwich)
library(lmtest)

set.seed(237)


## ── 1. LOAD DATA ─────────────────────────────────────────────────────────────
##
##  analysis_base, psm_ccsr_codes, and unit_groups are built in 01_cohort.R
##  and saved into cohorts.RData, so this script needs no CSV reads.

load("../data/cohorts.RData")
# Loads: patient_cohort, explicit_extubations, caregiver_summary,
#        caregivers_nonzero, lm_linear, lm_quadratic, lm_log,
#        unit_summary, analysis_base, psm_ccsr_codes, unit_groups

cat("analysis_base N:", nrow(analysis_base), "\n")
cat("Unit group sizes:\n")
print(table(analysis_base$unit_group))
cat("12-month mortality rate:",
    round(mean(analysis_base$mortality_12mo, na.rm = TRUE), 3), "\n")


## ── 2. WITHIN-GROUP QUARTILE SUMMARY (DESCRIPTIVE) ───────────────────────────

group_quartile_survival <- analysis_base |>
  filter(unit_group != "Other") |>
  group_by(unit_group) |>
  mutate(fe_quartile = ntile(caregiver_fe_rate, 4)) |>
  group_by(unit_group, fe_quartile) |>
  summarise(
    n              = n(),
    n_caregivers   = n_distinct(caregiver_id),
    fe_range       = paste0(round(min(caregiver_fe_rate) * 100, 1), "–",
                            round(max(caregiver_fe_rate) * 100, 1), "%"),
    mortality_12mo = round(mean(mortality_12mo, na.rm = TRUE), 3),
    .groups        = "drop"
  ) |>
  arrange(unit_group, fe_quartile)

cat("\n── Within-group FE rate quartile summary ",
    "───────────\n")
print(group_quartile_survival, n = Inf)
cat(paste0(rep("─", 72), collapse = ""), "\n")


## ── 3. CONTINUOUS TREATMENT MSM ──────────────────────────────────────────────

run_group_msm <- function(group_name, units, base_data) {

  group_data <- base_data |>
    filter(first_careunit %in% units) |>
    mutate(
      gender_m       = as.integer(gender == "M"),
      intub_surgical = as.integer(intubation_type == "surgical"),
      intub_med_resp = as.integer(intubation_type == "medical-respiratory")
    ) |>
    filter(
      !is.na(charlson), !is.na(sofa),
      !is.na(norepinephrine), !is.na(vent_hours)
    )

  n_pts        <- nrow(group_data)
  n_caregivers <- n_distinct(group_data$caregiver_id)
  fe_mean      <- mean(group_data$caregiver_fe_rate)
  fe_median    <- median(group_data$caregiver_fe_rate)
  fe_iqr       <- IQR(group_data$caregiver_fe_rate)
  fe_pct_zero  <- mean(group_data$caregiver_fe_rate == 0) * 100
  mort_mean    <- mean(group_data$mortality_12mo, na.rm = TRUE)

  if (n_caregivers < 20) {
    cat(sprintf("\n── %s: SKIPPED (%d caregivers < 20 minimum)\n",
                group_name, n_caregivers))
    return(NULL)
  }

  cat(sprintf(
    "\n── %s %s\n", group_name,
    paste0(rep("─", max(0, 46 - nchar(group_name))), collapse = "")
  ))
  cat(sprintf("N patients: %d  |  Caregivers: %d\n", n_pts, n_caregivers))
  cat(sprintf("FE rate — mean: %.4f  median: %.4f  IQR: %.4f\n",
              fe_mean, fe_median, fe_iqr))
  cat(sprintf("Zero FE rate: %.1f%%\n", fe_pct_zero))
  cat(sprintf("12-month mortality: %.3f\n", mort_mean))

  nonzero_ccsr <- psm_ccsr_codes[
    sapply(psm_ccsr_codes,
           function(v) var(group_data[[v]], na.rm = TRUE) > 0)
  ]
  cat(sprintf("CCSR flags (non-zero variance): %d\n", length(nonzero_ccsr)))

  cov_vars <- c(
    "anchor_age", "gender_m", "charlson", "sofa",
    "norepinephrine", "vent_hours", "intub_surgical", "intub_med_resp",
    nonzero_ccsr
  )
  if (length(units) > 1) cov_vars <- c(cov_vars, "first_careunit")

  ## Hirano-Imbens stabilized density weights
  den_lm  <- lm(
    as.formula(paste("caregiver_fe_rate ~",
                     paste(cov_vars, collapse = " + "))),
    data = group_data
  )
  num_lm  <- lm(caregiver_fe_rate ~ 1, data = group_data)

  f_den <- dnorm(group_data$caregiver_fe_rate,
                 predict(den_lm), summary(den_lm)$sigma)
  f_num <- dnorm(group_data$caregiver_fe_rate,
                 predict(num_lm), summary(num_lm)$sigma)

  sw <- f_num / pmax(f_den, 1e-10)
  sw <- pmin(sw, quantile(sw, 0.99))
  sw <- sw / mean(sw)

  eff_n  <- sum(sw)^2 / sum(sw^2)
  r2_den <- summary(den_lm)$r.squared

  cat(sprintf("Denominator R² (FE rate ~ covariates): %.3f\n", r2_den))
  cat(sprintf(
    "Effective N after weighting: %.0f  |  weight range %.2f–%.2f\n",
    eff_n, min(sw), max(sw)
  ))

  ## Covariate balance: weighted correlation with FE rate
  cont_covs <- c("anchor_age", "charlson", "sofa",
                 "norepinephrine", "vent_hours")
  cat("\nCovariate balance (|r| with FE rate, unweighted → weighted):\n")

  bal_rows <- lapply(cont_covs, function(v) {
    x  <- group_data[[v]]
    tr <- group_data$caregiver_fe_rate
    ok <- complete.cases(x, tr)
    r_unw <- cor(x[ok], tr[ok])
    r_wtd <- cov.wt(cbind(x[ok], tr[ok]),
                    wt = sw[ok], cor = TRUE)$cor[1, 2]
    cat(sprintf("  %-20s  %+.3f  →  %+.3f\n", v, r_unw, r_wtd))
    tibble(covariate    = v,
           r_unweighted = round(r_unw, 4),
           r_weighted   = round(r_wtd, 4))
  })
  bal_tbl <- bind_rows(bal_rows)

  ## Weighted logistic regression: 12-month mortality
  out_formula <- if (length(units) > 1) {
    mortality_12mo ~ caregiver_fe_rate + first_careunit
  } else {
    mortality_12mo ~ caregiver_fe_rate
  }

  msm_fit   <- glm(out_formula, data = group_data,
                   weights = sw, family = binomial)
  vcov_hc3  <- vcovHC(msm_fit, type = "HC3")
  vcov_cl   <- vcovCL(msm_fit, cluster = group_data$caregiver_id)

  fe_coef   <- coef(msm_fit)["caregiver_fe_rate"]
  fe_se_hc3 <- sqrt(vcov_hc3["caregiver_fe_rate", "caregiver_fe_rate"])
  fe_se_cl  <- sqrt(vcov_cl[ "caregiver_fe_rate", "caregiver_fe_rate"])

  delta  <- 0.05                        # effect per 5pp increase
  or_5pp <- exp(fe_coef * delta)
  lo_hc3 <- exp((fe_coef - 1.96 * fe_se_hc3) * delta)
  hi_hc3 <- exp((fe_coef + 1.96 * fe_se_hc3) * delta)
  p_hc3  <- 2 * pnorm(-abs(fe_coef / fe_se_hc3))
  lo_cl  <- exp((fe_coef - 1.96 * fe_se_cl)  * delta)
  hi_cl  <- exp((fe_coef + 1.96 * fe_se_cl)  * delta)
  p_cl   <- 2 * pnorm(-abs(fe_coef / fe_se_cl))

  cat("\nMSM: 12-month mortality ~ FE rate (cluster-robust SEs):\n")
  print(coeftest(msm_fit, vcov = vcov_cl))
  cat(sprintf(
    "\nOR per 5pp [cluster]: %.3f (%.3f–%.3f), p = %.4f\n",
    or_5pp, lo_cl, hi_cl, p_cl
  ))
  cat(sprintf(
    "OR per 5pp [HC3]:     %.3f (%.3f–%.3f), p = %.4f\n",
    or_5pp, lo_hc3, hi_hc3, p_hc3
  ))

  ## Secondary outcomes: hospital mortality, 30-day mortality
  sec_outcomes <- list()
  for (sec_outcome in c("hospital_expire_flag", "died_30d")) {
    if (!sec_outcome %in% names(group_data)) next
    sec_form <- if (length(units) > 1) {
      reformulate(c("caregiver_fe_rate", "first_careunit"),
                  response = sec_outcome)
    } else {
      reformulate("caregiver_fe_rate", response = sec_outcome)
    }
    fit_s <- tryCatch(
      glm(sec_form, data = group_data, weights = sw, family = binomial),
      error = function(e) NULL
    )
    if (is.null(fit_s)) next
    vcov_s <- vcovCL(fit_s, cluster = group_data$caregiver_id)
    fe_s   <- coef(fit_s)["caregiver_fe_rate"]
    se_s   <- sqrt(vcov_s["caregiver_fe_rate", "caregiver_fe_rate"])
    or_s   <- exp(fe_s * delta)
    lo_s   <- exp((fe_s - 1.96 * se_s) * delta)
    hi_s   <- exp((fe_s + 1.96 * se_s) * delta)
    p_s    <- 2 * pnorm(-abs(fe_s / se_s))
    lbl    <- if (sec_outcome == "hospital_expire_flag") "Hospital mortality"
              else "30-day mortality  "
    cat(sprintf("%s OR per 5pp: %.3f (%.3f–%.3f), p = %.4f\n",
                lbl, or_s, lo_s, hi_s, p_s))
    sec_outcomes[[length(sec_outcomes) + 1]] <- list(
      outcome = if (sec_outcome == "hospital_expire_flag") "hospital_mortality"
                else "died_30d",
      or  = round(or_s, 3),
      ci  = round(c(lo_s, hi_s), 3),
      p   = round(p_s, 4)
    )
  }

  list(
    group             = group_name,
    n_patients        = n_pts,
    n_caregivers      = n_caregivers,
    eff_n             = eff_n,
    fe_mean           = fe_mean,
    fe_median         = fe_median,
    fe_iqr            = fe_iqr,
    fe_pct_zero       = fe_pct_zero,
    r2_denominator    = r2_den,
    fe_coef           = fe_coef,
    or_per_5pp        = or_5pp,
    or_lo_cl          = lo_cl, or_hi_cl  = hi_cl, p_cl  = p_cl,
    or_lo_hc3         = lo_hc3, or_hi_hc3 = hi_hc3, p_hc3 = p_hc3,
    balance           = bal_tbl,
    msm_fit           = msm_fit,
    sw                = sw,
    group_data        = group_data,
    secondary_outcomes = sec_outcomes
  )
}

msm_results_list <- list(
  CVICU = run_group_msm(
    "CVICU", unit_groups$CVICU, analysis_base
  ),
  Medical = run_group_msm(
    "Medical (MICU + MICU/SICU)", unit_groups$Medical, analysis_base
  ),
  SurgTrauma = run_group_msm(
    "Surgical/Trauma (SICU + TSICU)", unit_groups$SurgTrauma, analysis_base
  )
)

cat("\n── MSM Summary Across Groups ",
    paste0(rep("─", 44), collapse = ""), "\n")
msm_summary_tbl <- bind_rows(lapply(msm_results_list, function(r) {
  if (is.null(r)) return(NULL)
  tibble(
    group         = r$group,
    n_patients    = r$n_patients,
    n_caregivers  = r$n_caregivers,
    fe_mean_pct   = round(r$fe_mean   * 100, 2),
    fe_median_pct = round(r$fe_median * 100, 2),
    r2_denom      = round(r$r2_denominator, 3),
    or_per_5pp    = round(r$or_per_5pp, 3),
    ci_lo         = round(r$or_lo_cl,   3),
    ci_hi         = round(r$or_hi_cl,   3),
    p_cl          = round(r$p_cl,        4)
  )
}))
print(msm_summary_tbl, n = Inf)
cat(paste0(rep("─", 72), collapse = ""), "\n")


## ── 4. COMFORT EXTUBATION SENSITIVITY — CVICU ────────────────────────────────

cat("\n── CVICU MSM: excluding likely comfort extubations ",
    paste0(rep("─", 23), collapse = ""), "\n")
n_comfort_cvicu <- sum(
  analysis_base$likely_comfort & analysis_base$unit_group == "CVICU",
  na.rm = TRUE
)
cat("Comfort extubations excluded from CVICU:", n_comfort_cvicu, "\n")

cvicu_noc <- run_group_msm(
  "CVICU (comfort-excluded)",
  unit_groups$CVICU,
  filter(analysis_base, !likely_comfort)
)


## (Partial correlations and NLS are in 03_nls_partial_cor.R)


## ── 5. FIGURES ────────────────────────────────────────────────────────────────

## 9a. Covariate balance plots (one per group)
for (r in msm_results_list) {
  if (is.null(r)) next

  bal <- r$balance |>
    pivot_longer(c(r_unweighted, r_weighted),
                 names_to = "type", values_to = "r") |>
    mutate(
      type = factor(
        recode(type, r_unweighted = "Unweighted", r_weighted = "Weighted"),
        levels = c("Unweighted", "Weighted")
      )
    )

  ggplot(bal, aes(y = covariate, x = r, color = type, shape = type)) +
    geom_vline(xintercept = c(-0.1, 0.1), linetype = "dashed",
               color = "gray60", linewidth = 0.5) +
    geom_vline(xintercept = 0, color = "gray40", linewidth = 0.4) +
    geom_point(size = 3.5) +
    scale_color_manual(
      values = c("Unweighted" = "#D95F02", "Weighted" = "#1D9E8E")
    ) +
    scale_shape_manual(values = c("Unweighted" = 16, "Weighted" = 17)) +
    scale_x_continuous(limits = c(-0.5, 0.5),
                       breaks = seq(-0.5, 0.5, 0.1)) +
    labs(
      title    = paste("MSM Covariate Balance —", r$group),
      subtitle = paste0(
        "Pearson r with caregiver FE rate. Dashed = ±0.10 threshold.  ",
        sprintf("N = %d patients, %d caregivers",
                r$n_patients, r$n_caregivers)
      ),
      x = "Correlation with caregiver FE rate",
      y = NULL, color = NULL, shape = NULL
    ) +
    theme_minimal() +
    theme(
      legend.position  = "top",
      panel.grid.minor = element_blank(),
      axis.text.y      = element_text(size = 10),
      plot.title       = element_text(size = 11, face = "bold"),
      plot.subtitle    = element_text(size = 8.5, color = "gray40")
    )

  fname <- paste0("../figures/msm_balance_",
                  tolower(gsub("[^a-zA-Z0-9]", "_", r$group)), ".png")
  ggsave(fname, width = 7, height = 4.5, dpi = 150)
  cat("Balance plot saved:", basename(fname), "\n")
}


## 9b. Primary forest plot (12-month mortality, cluster-robust SEs)
msm_forest_data <- bind_rows(lapply(msm_results_list, function(r) {
  if (is.null(r)) return(NULL)
  tibble(
    group  = r$group,
    or     = r$or_per_5pp,
    ci_lo  = r$or_lo_cl,
    ci_hi  = r$or_hi_cl,
    p      = r$p_cl,
    n      = r$n_patients
  )
})) |>
  mutate(
    group    = factor(group, levels = rev(unique(group))),
    sig      = p < 0.05,
    label_or = sprintf("%.3f (%.3f–%.3f)", or, ci_lo, ci_hi),
    label_p  = sprintf("p = %.3f", p),
    label_n  = sprintf("N = %d", n)
  )

ggplot(msm_forest_data,
       aes(y = group, x = or, xmin = ci_lo, xmax = ci_hi, color = sig)) +
  geom_vline(xintercept = 1, linetype = "dashed",
             color = "gray50", linewidth = 0.6) +
  geom_errorbarh(height = 0.15, linewidth = 0.8) +
  geom_point(size = 3.5) +
  geom_text(aes(x = 2.1,  label = label_or),
            hjust = 0, size = 3.2, color = "gray20") +
  geom_text(aes(x = 2.1, y = as.numeric(group) - 0.3, label = label_p),
            hjust = 0, size = 2.8, color = "gray40") +
  geom_text(aes(x = 0.22, label = label_n),
            hjust = 1, size = 2.8, color = "gray40") +
  scale_color_manual(
    values = c("FALSE" = "gray50", "TRUE" = "#2E75B6")
  ) +
  scale_x_continuous(
    limits = c(0.2, 3.2),
    breaks = c(0.25, 0.5, 1.0, 1.5, 2.0),
    trans  = "log10"
  ) +
  labs(
    title    = "MSM Results by ICU Group: 12-Month Mortality",
    subtitle = paste0(
      "OR per 5pp increase in caregiver FE rate. OR > 1 = more deaths.\n",
      "All caregivers included. Cluster-robust SEs by caregiver. Log scale."
    ),
    x = "Odds Ratio per 5pp FE rate (log scale)",
    y = NULL
  ) +
  theme_minimal() +
  theme(
    legend.position  = "none",
    panel.grid.minor = element_blank(),
    axis.text.y      = element_text(size = 10),
    plot.title       = element_text(size = 11, face = "bold"),
    plot.subtitle    = element_text(size = 8.5, color = "gray40")
  )

ggsave("../figures/msm_forest.png", width = 8, height = 3, dpi = 150)
cat("Forest plot saved: msm_forest.png\n")


## 9c. Secondary outcomes forest plot
sec_data <- bind_rows(lapply(msm_results_list, function(r) {
  if (is.null(r) || length(r$secondary_outcomes) == 0) return(NULL)
  bind_rows(lapply(r$secondary_outcomes, function(s) {
    tibble(group = r$group, outcome = s$outcome,
           or = s$or, ci_lo = s$ci[1], ci_hi = s$ci[2], p = s$p)
  }))
}))

if (nrow(sec_data) > 0) {
  sec_data <- sec_data |>
    mutate(
      group = factor(group,
                     levels = rev(levels(msm_forest_data$group))),
      outcome = recode(outcome,
                       hospital_mortality = "Hospital Mortality",
                       died_30d           = "30-Day Mortality") |>
        factor(levels = c("Hospital Mortality", "30-Day Mortality")),
      sig      = p < 0.05,
      label_or = sprintf("%.3f (%.3f–%.3f)\np = %.3f",
                         or, ci_lo, ci_hi, p)
    )

  ggplot(sec_data,
         aes(y = group, x = or, xmin = ci_lo, xmax = ci_hi,
             color = sig)) +
    facet_wrap(~outcome) +
    geom_vline(xintercept = 1, linetype = "dashed",
               color = "gray50", linewidth = 0.6) +
    geom_errorbarh(height = 0.15, linewidth = 0.8) +
    geom_point(size = 3.5) +
    geom_text(aes(x = 5.5, label = label_or),
              hjust = 0, size = 2.6, color = "gray20", lineheight = 0.9) +
    scale_color_manual(
      values = c("FALSE" = "gray50", "TRUE" = "#D95F02")
    ) +
    scale_x_continuous(
      limits = c(0.2, 9.0),
      breaks = c(0.25, 0.5, 1.0, 2.0, 4.0),
      trans  = "log10"
    ) +
    labs(
      title    = "Secondary Outcomes: Mortality by ICU Group",
      subtitle = paste0(
        "OR per 5pp increase in caregiver FE rate. OR > 1 = more deaths.\n",
        "Cluster-robust SEs by caregiver. Log scale. Orange = p < 0.05."
      ),
      x = "Odds Ratio per 5pp FE rate (log scale)",
      y = NULL
    ) +
    theme_minimal() +
    theme(
      legend.position  = "none",
      panel.grid.minor = element_blank(),
      strip.text       = element_text(size = 10, face = "bold"),
      axis.text.y      = element_text(size = 9),
      plot.title       = element_text(size = 11, face = "bold"),
      plot.subtitle    = element_text(size = 8.5, color = "gray40")
    )

  ggsave("../figures/msm_secondary_forest.png",
         width = 10, height = 3.5, dpi = 150)
  cat("Secondary forest plot saved: msm_secondary_forest.png\n")
}


## (NLS figures are in 03_nls_partial_cor.R)


## ── 6. RESULTS SUMMARY ───────────────────────────────────────────────────────

cat("\n── FINAL RESULTS SUMMARY ",
    paste0(rep("─", 49), collapse = ""), "\n")
cat("Continuous MSM (all caregivers, OR per 5pp FE rate,\n")
cat("cluster-robust SEs by caregiver_id, outcome = 12-mo mortality):\n\n")

for (r in msm_results_list) {
  if (is.null(r)) next
  cat(sprintf("  %-38s  N=%d (%d caregivers)\n",
              r$group, r$n_patients, r$n_caregivers))
  cat(sprintf("    Denom R²=%.3f  Eff.N=%.0f\n",
              r$r2_denominator, r$eff_n))
  cat(sprintf("    OR=%.3f (%.3f–%.3f), p=%.4f\n\n",
              r$or_per_5pp, r$or_lo_cl, r$or_hi_cl, r$p_cl))
}

cat("(Partial correlations and NLS results: see 03_nls_partial_cor.R)\n")
cat(paste0(rep("─", 72), collapse = ""), "\n")


## ── 7. SAVE ───────────────────────────────────────────────────────────────────

save(
  msm_results_list, msm_summary_tbl, group_quartile_survival,
  cvicu_noc, analysis_base, psm_ccsr_codes,
  file = "../data/msm_results.RData"
)
cat("\nSaved: ../data/msm_results.RData\n")
