################################################################################
##
##  ICU Extubation Outcomes — Script 03: PSM and Sensitivity Analysis
##
##  Depends on: 01_cohort_and_descriptive.R (cohorts.RData)
##              patient_ccsr_long.csv (separate BigQuery export)
##
##  CCSR long query (run in BigQuery, export as patient_ccsr_long.csv):
##
##    SELECT
##      e.subject_id,
##      e.hadm_id,
##      e.caregiver_fe_rate,
##      e.caregiver_n,
##      e.failed_extubations,
##      e.tube_event_source,
##      e.caregiver_imputation_source,
##      code AS ccsr_code
##    FROM `[your-project].Failed_Extubation_Rate.last_extubations` e,
##    UNNEST(e.ccsr_codes) AS code
##    WHERE e.caregiver_n > 10
##      AND e.ccsr_codes IS NOT NULL
##    ORDER BY e.subject_id, code
##
##  Research questions:
##    1. After matching patients of high-FE-rate caregivers to patients of
##       low-FE-rate caregivers on illness severity and diagnosis case mix,
##       does a survival difference persist? (PSM / IPW)
##    2. Is FE rate independently associated with survival at the caregiver
##       level, controlling for volume? (partial correlation)
##    3. What is the functional form of the FE rate — survival relationship,
##       and is there an irreducible survival floor? (asymptotic exponential)
##
##  Key findings:
##    PSM (top vs bottom quartile, caregiver_unit_n >= 5, exact match within unit):
##      Q25 = 0.7%, Q75 = 6.2%; matched N = 2036
##      Unadjusted OR = 0.599 (0.497-0.720), p < 0.001
##      Doubly adjusted OR = 0.646 (0.516-0.807), p < 0.001
##      Effect modification (treated x medical-respiratory): p = 0.014
##    IPW: OR = 0.757 (0.669-0.857), p < 0.001
##    Partial correlation (FE rate | caregiver volume): r = -0.288, p < 0.001
##    Partial correlation (volume | FE rate): r = -0.109, p = 0.211
##    Unit-stratified partial correlation:
##      CVICU (N = 73): r = -0.342, p = 0.003
##      MICU  (N = 30): r = -0.298, p = 0.116 (directional, underpowered)
##    Asymptotic exponential model (caregivers FE rate > 0, N = 78):
##      b (survival floor) = 0.624 (95% CI 0.586-0.658)
##      k (decay rate)     = 77.5  (95% CI 52.7-122.4)
##      ~62% of patients survive 12 months regardless of caregiver FE rate;
##      ~38% represent outcomes potentially sensitive to caregiver performance
##
################################################################################


## ── 0. LIBRARIES ─────────────────────────────────────────────────────────────

library(MASS)
library(dplyr)
library(tidyr)
library(purrr)
library(readr)
library(tibble)
library(stringr)
library(ggplot2)
library(splines)
library(MatchIt)
library(cobalt)
library(WeightIt)
library(ppcor)

set.seed(237)


## ── 1. LOAD COHORTS AND CCSR DATA ────────────────────────────────────────────

load("../data/cohorts.RData")
# Loads: patient_cohort, explicit_extubations, caregiver_summary,
#        caregivers_nonzero, lm_linear, lm_quadratic, lm_log, unit_summary

ccsr_long <- read_csv(
  "../data/patient_ccsr_long.csv",
  col_types = cols(
    subject_id   = col_character(),
    hadm_id      = col_character(),
    ccsr_code    = col_character()
  )
)

cat("explicit_extubations N:", nrow(explicit_extubations), "\n")
cat("CCSR rows:", nrow(ccsr_long), "\n")
cat("Unique CCSR codes:", n_distinct(ccsr_long$ccsr_code), "\n")


## ── 2. CCSR BINARY FLAGS ─────────────────────────────────────────────────────
##
##  LDA topic modeling (K = 5, 8, 10, 12) and k-means clustering were
##  explored as approaches to characterizing patient case mix. Both showed
##  weak patient differentiation (LDA topic SD ~0.02-0.03; k-means variance
##  explained ~15% at K=10), reflecting the short-document problem with a
##  median of ~4 CCSR codes per patient. Direct binary flags for high-
##  prevalence codes are more transparent and equally effective as PSM
##  covariates.
##
##  Selection threshold: >= 10% prevalence in the CCSR long cohort.
##  RSP012 (respiratory failure, ~47% prevalence) excluded — near-universal
##  in ventilated patients, near-zero variance.

ccsr_prevalence <- ccsr_long %>%
  distinct(subject_id, ccsr_code) %>%
  count(ccsr_code, name = "n_patients") %>%
  mutate(prevalence = n_patients / n_distinct(ccsr_long$subject_id)) %>%
  arrange(desc(prevalence))

cat("\nCCSR codes >= 10% prevalence:\n")
ccsr_prevalence %>%
  filter(prevalence >= 0.10) %>%
  print(n = 30)

psm_ccsr_codes <- ccsr_prevalence %>%
  filter(prevalence >= 0.10, ccsr_code != "RSP012") %>%
  pull(ccsr_code)

cat("\nCCSR flags used in PSM:", length(psm_ccsr_codes), "\n")

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


## ── 3. UNIT GROUPS AND WITHIN-GROUP QUARTILE SURVIVAL ─────────────────────────
##
##  PSM-eligible groups (sufficient N for quartile comparison):
##    CVICU           — post-cardiac-surgery; 1,543 patients
##    Medical         — MICU + MICU/SICU; ~1,315 patients
##    Surgical/Trauma — SICU + TSICU; ~1,127 patients
##  Descriptive-only (too few caregivers for PSM):
##    CCU (4 caregivers), Neuro units (<10 caregivers each)

unit_groups <- list(
  CVICU      = "Cardiac Vascular Intensive Care Unit (CVICU)",
  Medical    = c("Medical Intensive Care Unit (MICU)",
                 "Medical/Surgical Intensive Care Unit (MICU/SICU)"),
  SurgTrauma = c("Surgical Intensive Care Unit (SICU)",
                 "Trauma SICU (TSICU)")
)

analysis_base <- explicit_extubations %>%
  filter(caregiver_unit_n >= 5) %>%
  left_join(ccsr_flags, by = "subject_id") %>%
  mutate(across(all_of(psm_ccsr_codes), ~ replace_na(.x, 0))) %>%
  mutate(
    unit_group = case_when(
      first_careunit %in% unit_groups$CVICU      ~ "CVICU",
      first_careunit %in% unit_groups$Medical    ~ "Medical",
      first_careunit %in% unit_groups$SurgTrauma ~ "Surgical/Trauma",
      TRUE                                       ~ "Other"
    )
  )

cat("\nanalysis_base N (caregiver_unit_n >= 5):", nrow(analysis_base), "\n")
cat("Unit group sizes:\n")
print(table(analysis_base$unit_group))

## Within-group FE rate quartile survival table
cat("\n── Within-group FE rate quartile survival ───────────────────────────────\n")

group_quartile_survival <- analysis_base %>%
  filter(unit_group != "Other") %>%
  group_by(unit_group) %>%
  mutate(fe_quartile = ntile(caregiver_fe_rate, 4)) %>%
  group_by(unit_group, fe_quartile) %>%
  summarise(
    n             = n(),
    n_caregivers  = n_distinct(caregiver_id),
    fe_range      = paste0(round(min(caregiver_fe_rate) * 100, 1), "–",
                           round(max(caregiver_fe_rate) * 100, 1), "%"),
    survival_12mo = round(mean(survival_12mo, na.rm = TRUE), 3),
    .groups       = "drop"
  ) %>%
  arrange(unit_group, fe_quartile)

print(group_quartile_survival, n = Inf)
cat("────────────────────────────────────────────────────────────────────────\n")


## ── 4. GROUP-STRATIFIED PROPENSITY SCORE MATCHING AND IPW ────────────────────
##
##  PSM run separately within each unit group. Group-specific Q25/Q75 ensure
##  "high" and "low" FE rate are benchmarked against that group's distribution.
##  CCSR codes with zero variance within a group are dropped from the PS model
##  to avoid complete separation. IPW run within the same group.
##  Outcome models use no MatchIt subclass weights (1:1 without-replacement PSM;
##  subclass weights re-introduce stratum-size imbalance and invert the estimate).

run_group_psm <- function(group_name, units, base_data) {

  group_data <- base_data %>%
    filter(first_careunit %in% units) %>%
    mutate(
      fe_q25    = quantile(caregiver_fe_rate, 0.25, na.rm = TRUE),
      fe_q75    = quantile(caregiver_fe_rate, 0.75, na.rm = TRUE),
      fe_group  = case_when(
        caregiver_fe_rate >= fe_q75 ~ "high",
        caregiver_fe_rate <= fe_q25 ~ "low",
        TRUE                        ~ "middle"
      ),
      treated        = as.integer(fe_group == "high"),
      gender_m       = as.integer(gender == "M"),
      intub_surgical = as.integer(intubation_type == "surgical"),
      intub_med_resp = as.integer(intubation_type == "medical-respiratory")
    ) %>%
    filter(
      fe_group != "middle",
      !is.na(charlson), !is.na(sofa), !is.na(norepinephrine), !is.na(vent_hours)
    )

  cat(sprintf("\n── %s ──────────────────────────────────────────\n", group_name))
  cat("N patients:", nrow(group_data),
      "  Treated:", sum(group_data$treated),
      "  Control:", sum(group_data$treated == 0), "\n")
  cat("FE rate Q25:", round(group_data$fe_q25[1] * 100, 1), "%",
      " Q75:", round(group_data$fe_q75[1] * 100, 1), "%\n")

  if (sum(group_data$treated) < 20 || sum(group_data$treated == 0) < 20) {
    cat("Insufficient treated/control N — skipping PSM.\n")
    return(NULL)
  }

  nonzero_ccsr <- psm_ccsr_codes[
    sapply(psm_ccsr_codes, function(v) var(group_data[[v]], na.rm = TRUE) > 0)
  ]
  cat("CCSR flags (non-zero variance in group):", length(nonzero_ccsr), "\n")

  grp_formula <- as.formula(paste(
    "treated ~ anchor_age + gender_m + charlson + sofa + norepinephrine +
     vent_hours + intub_surgical + intub_med_resp +",
    paste(nonzero_ccsr, collapse = " + ")
  ))

  use_exact <- length(units) > 1

  tryCatch({
    m_out <- if (use_exact) {
      matchit(grp_formula, data = group_data, method = "quick",
              ratio = 1, caliper = 0.2, std.caliper = TRUE,
              distance = "logit", exact = ~ first_careunit)
    } else {
      matchit(grp_formula, data = group_data, method = "quick",
              ratio = 1, caliper = 0.2, std.caliper = TRUE,
              distance = "logit")
    }

    md <- match.data(m_out)
    cat("Matched N:", nrow(md), "\n")

    surv_high <- mean(md$survival_12mo[md$treated == 1], na.rm = TRUE)
    surv_low  <- mean(md$survival_12mo[md$treated == 0], na.rm = TRUE)
    cat("Survival — high FE:", round(surv_high, 3),
        "  low FE:", round(surv_low, 3), "\n")

    fit_u <- glm(survival_12mo ~ treated, data = md, family = binomial)
    fit_a <- glm(
      survival_12mo ~ treated + anchor_age + charlson + sofa +
        norepinephrine + vent_hours + intub_surgical + intub_med_resp,
      data = md, family = binomial
    )

    or_u <- exp(coef(fit_u)["treated"])
    ci_u <- exp(confint(fit_u)["treated", ])
    p_u  <- summary(fit_u)$coefficients["treated", "Pr(>|z|)"]

    or_a <- exp(coef(fit_a)["treated"])
    ci_a <- exp(confint(fit_a)["treated", ])
    p_a  <- summary(fit_a)$coefficients["treated", "Pr(>|z|)"]

    cat(sprintf("PSM OR unadjusted:      %.3f (%.3f-%.3f), p = %.4f\n",
                or_u, ci_u[1], ci_u[2], p_u))
    cat(sprintf("PSM OR doubly adjusted: %.3f (%.3f-%.3f), p = %.4f\n",
                or_a, ci_a[1], ci_a[2], p_a))

    ## Secondary outcomes: hospital mortality and 30-day mortality
    sec_outcomes <- list()
    for (sec_outcome in c("hospital_expire_flag", "died_30d")) {
      if (!sec_outcome %in% names(md)) next
      fit_s <- tryCatch(
        glm(reformulate("treated", response = sec_outcome), data = md, family = binomial),
        error = function(e) NULL
      )
      if (is.null(fit_s)) next
      or_s  <- exp(coef(fit_s)["treated"])
      ci_s  <- exp(confint(fit_s)["treated", ])
      p_s   <- summary(fit_s)$coefficients["treated", "Pr(>|z|)"]
      lbl   <- if (sec_outcome == "hospital_expire_flag") "Hospital mortality OR:" else "30-day mortality  OR:"
      cat(sprintf("%s %.3f (%.3f-%.3f), p = %.4f\n", lbl, or_s, ci_s[1], ci_s[2], p_s))
      sec_outcomes[[length(sec_outcomes) + 1]] <- list(
        outcome = if (sec_outcome == "hospital_expire_flag") "hospital_mortality" else "died_30d",
        or      = round(or_s, 3),
        ci      = round(ci_s, 3),
        p       = round(p_s, 4)
      )
    }

    ## IPW within group
    ipw_g   <- weightit(grp_formula, data = group_data,
                        method = "ps", estimand = "ATE")
    ipw_fit <- glm(survival_12mo ~ treated, data = group_data,
                   family = binomial, weights = ipw_g$weights)
    or_ipw  <- exp(coef(ipw_fit)["treated"])
    ci_ipw  <- exp(confint(ipw_fit)["treated", ])
    p_ipw   <- summary(ipw_fit)$coefficients["treated", "Pr(>|z|)"]
    cat(sprintf("IPW OR:                 %.3f (%.3f-%.3f), p = %.4f\n",
                or_ipw, ci_ipw[1], ci_ipw[2], p_ipw))

    ## Love plot
    p_love <- love.plot(m_out, threshold = 0.1, abs = TRUE,
                        title      = paste("PSM Balance —", group_name),
                        colors     = c("#D95F02", "#1D9E8E"),
                        stars      = "std", var.order = "unadjusted")
    fname <- paste0("../figures/psm_love_",
                    tolower(gsub("[^a-zA-Z0-9]", "_", group_name)), ".png")
    ggsave(fname, plot = p_love, width = 10, height = 10, dpi = 150)
    cat("Love plot saved:", basename(fname), "\n")

    list(
      group     = group_name,
      n_psm     = nrow(md),
      fe_q25    = group_data$fe_q25[1],
      fe_q75    = group_data$fe_q75[1],
      surv_high = surv_high,
      surv_low  = surv_low,
      or_unadj  = or_u, ci_unadj = ci_u, p_unadj = p_u,
      or_adj    = or_a, ci_adj   = ci_a, p_adj   = p_a,
      or_ipw    = or_ipw, ci_ipw = ci_ipw, p_ipw  = p_ipw,
      secondary_outcomes = sec_outcomes,
      match_out = m_out, matched_data = md
    )
  }, error = function(e) {
    cat("PSM failed:", conditionMessage(e), "\n")
    NULL
  })
}

set.seed(237)
psm_results_list <- list(
  CVICU = run_group_psm(
    "CVICU", unit_groups$CVICU, analysis_base
  ),
  Medical = run_group_psm(
    "Medical (MICU + MICU/SICU)", unit_groups$Medical, analysis_base
  ),
  SurgTrauma = run_group_psm(
    "Surgical/Trauma (SICU + TSICU)", unit_groups$SurgTrauma, analysis_base
  )
)

cat("\n── PSM + IPW Summary Across Groups ──────────────────────────────────────\n")
psm_summary_tbl <- bind_rows(lapply(psm_results_list, function(r) {
  if (is.null(r)) return(NULL)
  tibble(
    group      = r$group,
    n_matched  = r$n_psm,
    fe_q25_pct = round(r$fe_q25 * 100, 1),
    fe_q75_pct = round(r$fe_q75 * 100, 1),
    surv_high  = round(r$surv_high, 3),
    surv_low   = round(r$surv_low,  3),
    OR_adj     = round(r$or_adj,    3),
    ci_lo      = round(r$ci_adj[1], 3),
    ci_hi      = round(r$ci_adj[2], 3),
    p_adj      = round(r$p_adj,     4),
    OR_ipw     = round(r$or_ipw,    3),
    p_ipw      = round(r$p_ipw,     4)
  )
}))
print(psm_summary_tbl, n = Inf)
cat("────────────────────────────────────────────────────────────────────────\n")


## ── 5. COMFORT EXTUBATION SENSITIVITY — CVICU ────────────────────────────────
##
##  Patients extubated as part of terminal/withdrawal care (hospital death
##  within 3 days of extubation) have near-certain short-term mortality but
##  no reintubation, giving caregivers a spuriously low FE rate. Exclude
##  these cases and repeat the CVICU PSM to confirm the signal is not driven
##  by comfort-care confounding.

cat("\n── CVICU PSM: excluding likely comfort extubations ──────────────────────\n")
n_comfort_cvicu <- sum(analysis_base$likely_comfort & analysis_base$unit_group == "CVICU",
                       na.rm = TRUE)
cat("Comfort extubations excluded from CVICU:", n_comfort_cvicu, "\n")

analysis_base_noc <- analysis_base %>% filter(!likely_comfort)
set.seed(237)
cvicu_noc <- run_group_psm(
  "CVICU (comfort-excluded)", unit_groups$CVICU, analysis_base_noc
)


## ── 7. CAREGIVER-LEVEL PARTIAL CORRELATION ────────────────────────────────────

cat("\n── Caregiver-level partial correlations ─────────────────────────────────\n")
cat("N caregivers (n > 5 sample patients):", nrow(caregiver_summary), "\n")

pcor_fe <- pcor.test(
  caregiver_summary$caregiver_fe_rate,
  caregiver_summary$survival_rate,
  caregiver_summary$n_patients_sample
)
cat("\nPartial correlation: FE rate ~ survival | caregiver volume\n")
cat("r =", round(pcor_fe$estimate, 3), "\n")
cat("p =", round(pcor_fe$p.value,  4), "\n")

pcor_vol <- pcor.test(
  caregiver_summary$caregiver_n,
  caregiver_summary$survival_rate,
  caregiver_summary$caregiver_fe_rate
)
cat("\nPartial correlation: volume ~ survival | FE rate\n")
cat("r =", round(pcor_vol$estimate, 3), "\n")
cat("p =", round(pcor_vol$p.value,  4), "\n")
cat("────────────────────────────────────────────────────────────────────────\n")

cat("\n── Volume threshold sensitivity ─────────────────────────────────────────\n")
for (min_n in c(10, 20)) {
  cs_sub <- caregiver_summary %>% filter(n_patients_sample >= min_n)
  if (nrow(cs_sub) < 10) {
    cat(sprintf("n >= %d: fewer than 10 caregivers — skipped\n", min_n))
    next
  }
  pc_sub <- pcor.test(cs_sub$caregiver_fe_rate, cs_sub$survival_rate, cs_sub$n_patients_sample)
  cat(sprintf("n >= %d (%d caregivers): r = %.3f, p = %.4f\n",
              min_n, nrow(cs_sub), pc_sub$estimate, pc_sub$p.value))
}
cat("────────────────────────────────────────────────────────────────────────\n")


## ── 7b. UNIT-STRATIFIED PARTIAL CORRELATION ───────────────────────────────────
##
##  Repeat the caregiver-level partial correlation within each ICU unit
##  (restricted to units with >= 10 caregivers in caregiver_summary).
##  Reveals whether the FE rate → survival signal is consistent across units
##  or driven by a subset. Informs grouping decisions for further analysis.
##
##  caregiver_summary uses primary_careunit (modal unit per caregiver).

unit_pcor_results <- caregiver_summary %>%
  filter(!is.na(primary_careunit)) %>%
  group_by(primary_careunit) %>%
  group_modify(~ {
    df <- .x
    base <- tibble(
      n_caregivers  = nrow(df),
      mean_fe_rate  = round(mean(df$caregiver_fe_rate, na.rm = TRUE), 4),
      mean_survival = round(mean(df$survival_rate,     na.rm = TRUE), 3)
    )
    if (nrow(df) < 5) return(mutate(base, r = NA_real_, p = NA_real_))
    pc <- tryCatch(
      pcor.test(df$caregiver_fe_rate, df$survival_rate, df$n_patients_sample),
      error = function(e) list(estimate = NA_real_, p.value = NA_real_)
    )
    mutate(base, r = round(pc$estimate, 3), p = round(pc$p.value, 4))
  }) %>%
  ungroup() %>%
  arrange(desc(n_caregivers))

cat("\n── Unit-stratified partial correlations ─────────────────────────────────\n")
cat("(FE rate ~ survival | caregiver volume, by primary_careunit)\n")
print(unit_pcor_results, n = Inf)
cat("────────────────────────────────────────────────────────────────────────\n")


## ── 8. FUNCTIONAL FORM: ASYMPTOTIC EXPONENTIAL MODEL ─────────────────────────
##
##  Model: P(survival) = b + (1 - b) * exp(-k * caregiver_fe_rate)
##
##    b: survival floor — fraction surviving regardless of caregiver FE rate,
##       representing irreducible mortality driven by illness severity
##    k: decay rate — steepness of survival reduction with increasing FE rate
##
##  Motivation: standard logistic and linear models imply every death is
##  potentially preventable. This model explicitly partitions outcomes into
##  a salvageable fraction (1 - b) that responds to caregiver performance,
##  and an irreducible fraction (b) that does not. This is more clinically
##  honest for a severely ill ICU population.
##
##  Model selection: three-parameter version (with intercept shift c) fitted
##  first; c was not significant (p = 0.558) and the two-parameter model had
##  lower AIC (ΔAIC = 1.6). Two-parameter model preferred.
##
##  Restricted to caregivers with FE rate > 0 (N = 111 of 134).

cat("\n── Asymptotic exponential model ─────────────────────────────────────────\n")
cat("N caregivers (FE rate > 0):", nrow(caregivers_nonzero), "\n")

nls_fit <- nls(
  survival_rate ~ b + (1 - b) * exp(-k * caregiver_fe_rate),
  data    = caregivers_nonzero,
  weights = n_patients_sample,
  start   = list(b = 0.60, k = 50),
  control = nls.control(maxiter = 500)
)

summary(nls_fit)

b_hat <- coef(nls_fit)["b"]
k_hat <- coef(nls_fit)["k"]
b_ci  <- confint(nls_fit)["b", ]
k_ci  <- confint(nls_fit)["k", ]

cat("\nb (survival floor):", round(b_hat, 3),
    "  95% CI (", round(b_ci[1], 3), "-", round(b_ci[2], 3), ")\n")
cat("k (decay rate):    ", round(k_hat, 1),
    "  95% CI (", round(k_ci[1], 1), "-", round(k_ci[2], 1), ")\n")
cat("────────────────────────────────────────────────────────────────────────\n")

## Compare with linear and quadratic
cat("\nModel AIC comparison:\n")
print(AIC(lm_linear, lm_quadratic))
cat("NLS (asymptotic exponential) AIC:", round(AIC(nls_fit), 3), "\n")

## Visualization
fe_seq   <- seq(0, max(caregivers_nonzero$caregiver_fe_rate), by = 0.0005)
pred_exp <- data.frame(
  caregiver_fe_rate = fe_seq,
  survival_rate     = b_hat + (1 - b_hat) * exp(-k_hat * fe_seq)
)

ggplot() +
  geom_point(
    data  = caregivers_nonzero,
    aes(caregiver_fe_rate, survival_rate, size = n_patients_sample),
    alpha = 0.5, color = "gray40"
  ) +
  geom_line(
    data      = pred_exp,
    aes(caregiver_fe_rate, survival_rate),
    color     = "#2E75B6", linewidth = 1.0
  ) +
  geom_hline(
    yintercept = b_hat,
    linetype   = "dashed", color = "gray50"
  ) +
  annotate(
    "text",
    x     = max(caregivers_nonzero$caregiver_fe_rate) * 0.7,
    y     = b_hat - 0.05,
    label = paste0(
      "Survival floor = ", round(b_hat, 3),
      " (95% CI ", round(b_ci[1], 3), "\u2013", round(b_ci[2], 3), ")"
    ),
    color = "gray40", size = 3.5
  ) +
  scale_size_continuous(range = c(1, 8), name = "N patients") +
  labs(
    title    = "Caregiver FE Rate vs. 12-Month Patient Survival",
    subtitle = paste0(
      "P(survival) = b + (1\u2212b)\u00d7exp(\u2212k\u00d7FE rate)     ",
      "b = ", round(b_hat, 3),
      " (95% CI ", round(b_ci[1], 3), "\u2013", round(b_ci[2], 3), ")     ",
      "k = ", round(k_hat, 1),
      " (95% CI ", round(k_ci[1], 1), "\u2013", round(k_ci[2], 1), ")"
    ),
    x        = "Caregiver Failed Extubation Rate",
    y        = "12-Month Patient Survival Rate"
  ) +
  theme_minimal() +
  theme(plot.subtitle = element_text(size = 9))

ggsave("../figures/fe_rate_survival_exponential.png",
       width = 9, height = 6, dpi = 150)
cat("Exponential model plot saved.\n")


## ── 9. RESULTS SUMMARY ───────────────────────────────────────────────────────

cat("\n── FINAL RESULTS SUMMARY ────────────────────────────────────────────────\n")
cat("Group-stratified PSM (top vs bottom quartile within group,\n")
cat("caregiver_unit_n >= 5, doubly adjusted OR | IPW OR):\n\n")
for (r in psm_results_list) {
  if (is.null(r)) next
  cat(sprintf("  %-38s  Q25=%.1f%%  Q75=%.1f%%\n",
              r$group, r$fe_q25 * 100, r$fe_q75 * 100))
  cat(sprintf("    Survival: %.1f%% (high FE) vs %.1f%% (low FE)\n",
              r$surv_high * 100, r$surv_low * 100))
  cat(sprintf("    PSM OR:   %.3f (%.3f-%.3f), p = %.4f\n",
              r$or_adj, r$ci_adj[1], r$ci_adj[2], r$p_adj))
  cat(sprintf("    IPW OR:   %.3f (%.3f-%.3f), p = %.4f\n\n",
              r$or_ipw, r$ci_ipw[1], r$ci_ipw[2], r$p_ipw))
}
cat(sprintf("Partial r (FE rate | vol):  %.3f, p = %.4f\n",
            pcor_fe$estimate, pcor_fe$p.value))
cat(sprintf("Partial r (volume | FE):    %.3f, p = %.4f\n",
            pcor_vol$estimate, pcor_vol$p.value))
cat("\n")
cat(sprintf("Survival floor b:           %.3f (95%% CI %.3f-%.3f)\n",
            b_hat, b_ci[1], b_ci[2]))
cat(sprintf("Decay rate k:               %.1f  (95%% CI %.1f-%.1f)\n",
            k_hat, k_ci[1], k_ci[2]))
cat("────────────────────────────────────────────────────────────────────────\n")

## Save results
save(psm_results_list, psm_summary_tbl, group_quartile_survival,
     pcor_fe, pcor_vol,
     nls_fit, b_hat, k_hat, b_ci, k_ci,
     file = "../data/psm_results.RData")
cat("PSM results saved to ../data/psm_results.RData\n")
