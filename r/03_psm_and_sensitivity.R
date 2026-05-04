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
##    PSM (top vs bottom quartile FE rate):
##      Unadjusted OR = 0.660 (0.553-0.787), p < 0.001
##      Doubly adjusted OR = 0.571 (0.465-0.700), p < 0.001
##    IPW: OR = 0.699 (0.616-0.792), p < 0.001
##    Partial correlation (FE rate | caregiver volume): r = -0.363, p < 0.001
##    Partial correlation (volume | FE rate): r = -0.121, p = 0.167
##    Asymptotic exponential model:
##      b (survival floor) = 0.590 (95% CI 0.521-0.633)
##      k (decay rate)     = 47.2  (95% CI 31.2-71.0)
##      ~59% of patients survive 12 months regardless of caregiver FE rate;
##      ~41% represent outcomes potentially sensitive to caregiver performance
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
#        caregivers_nonzero, lm_linear, lm_quadratic, lm_log

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


## ── 3. PSM DATASET CONSTRUCTION ──────────────────────────────────────────────
##
##  Treatment: caregiver_fe_rate above 75th percentile ("high")
##             caregiver_fe_rate below 25th percentile ("low")
##  Middle 50% excluded to sharpen the contrast.

analysis_data <- explicit_extubations %>%
  left_join(ccsr_flags, by = "subject_id") %>%
  mutate(across(all_of(psm_ccsr_codes), ~ replace_na(.x, 0))) %>%
  mutate(
    fe_q25 = quantile(caregiver_fe_rate, 0.25, na.rm = TRUE),
    fe_q75 = quantile(caregiver_fe_rate, 0.75, na.rm = TRUE),
    fe_group = case_when(
      caregiver_fe_rate >= fe_q75 ~ "high",
      caregiver_fe_rate <= fe_q25 ~ "low",
      TRUE                        ~ "middle"
    )
  )

cat("\nFE rate Q25:", round(analysis_data$fe_q25[1], 3), "\n")
cat("FE rate Q75:", round(analysis_data$fe_q75[1], 3), "\n")
cat("Group sizes:\n")
print(table(analysis_data$fe_group))

cat("\nSurvival by FE category:\n")
analysis_data %>%
  mutate(fe_cat = case_when(
    caregiver_fe_rate < 0.05  ~ "1. <5%",
    caregiver_fe_rate < 0.15  ~ "2. 5-15%",
    caregiver_fe_rate < 0.25  ~ "3. 15-25%",
    TRUE                      ~ "4. >25%"
  )) %>%
  group_by(fe_cat) %>%
  summarise(
    n_patients    = n(),
    n_caregivers  = n_distinct(caregiver_id),
    mean_fe_rate  = mean(caregiver_fe_rate),
    survival_12mo = mean(survival_12mo, na.rm = TRUE),
    .groups       = "drop"
  ) %>%
  print()

psm_data <- analysis_data %>%
  filter(fe_group != "middle") %>%
  mutate(
    treated        = as.integer(fe_group == "high"),
    gender_m       = as.integer(gender == "M"),
    intub_surgical = as.integer(intubation_type == "surgical"),
    intub_med_resp = as.integer(intubation_type == "medical-respiratory")
  ) %>%
  filter(
    !is.na(charlson),
    !is.na(sofa),
    !is.na(norepinephrine),
    !is.na(vent_hours)
  )

cat("\nPSM dataset N:", nrow(psm_data), "\n")
cat("Treated (high FE >=", round(analysis_data$fe_q75[1], 3), "):",
    sum(psm_data$treated), "\n")
cat("Control (low FE <=", round(analysis_data$fe_q25[1], 3), "):",
    sum(psm_data$treated == 0), "\n")

ccsr_formula_terms <- paste(psm_ccsr_codes, collapse = " + ")

psm_formula <- as.formula(paste(
  "treated ~ anchor_age + gender_m + charlson + sofa + norepinephrine +
   vent_hours + intub_surgical + intub_med_resp +",
  ccsr_formula_terms
))


## ── 4. PROPENSITY SCORE MATCHING ─────────────────────────────────────────────

set.seed(237)
match_out <- matchit(
  psm_formula,
  data        = psm_data,
  method      = "quick",
  ratio       = 1,
  caliper     = 0.2,
  std.caliper = TRUE,
  distance    = "logit"
)

summary(match_out)

p_love <- love.plot(
  match_out,
  threshold  = 0.1,
  title      = "Covariate Balance Before and After PSM",
  colors     = c("#D95F02", "#1D9E8E"),
  stars      = "std",
  var.order  = "unadjusted",
  abs        = TRUE
)
print(p_love)
ggsave("../figures/psm_love_plot.png", plot = p_love,
       width = 10, height = 12, dpi = 150)
cat("Love plot saved.\n")


## ── 5. TREATMENT EFFECT ESTIMATION ───────────────────────────────────────────

matched_data <- match.data(match_out)

cat("\nMatched N:", nrow(matched_data), "\n")
cat("Survival — high FE:",
    round(mean(matched_data$survival_12mo[matched_data$treated == 1],
               na.rm = TRUE), 3), "\n")
cat("Survival — low FE:",
    round(mean(matched_data$survival_12mo[matched_data$treated == 0],
               na.rm = TRUE), 3), "\n")

outcome_unadj <- glm(
  survival_12mo ~ treated,
  data    = matched_data,
  family  = binomial,
  weights = weights
)
or_unadj <- exp(coef(outcome_unadj)["treated"])
ci_unadj <- exp(confint(outcome_unadj)["treated", ])
p_unadj  <- summary(outcome_unadj)$coefficients["treated", "Pr(>|z|)"]

cat("\nOR (unadjusted):", round(or_unadj, 3), "\n")
cat("95% CI:", round(ci_unadj, 3), "\n")
cat("p-value:", round(p_unadj, 4), "\n")

outcome_adj <- glm(
  survival_12mo ~ treated + anchor_age + charlson + sofa +
    norepinephrine + vent_hours + intub_surgical + intub_med_resp,
  data    = matched_data,
  family  = binomial,
  weights = weights
)
or_adj <- exp(coef(outcome_adj)["treated"])
ci_adj <- exp(confint(outcome_adj)["treated", ])
p_adj  <- summary(outcome_adj)$coefficients["treated", "Pr(>|z|)"]

cat("\nOR (doubly adjusted):", round(or_adj, 3), "\n")
cat("95% CI:", round(ci_adj, 3), "\n")
cat("p-value:", round(p_adj, 4), "\n")

outcome_int <- glm(
  survival_12mo ~ treated * intub_med_resp,
  data    = matched_data,
  family  = binomial,
  weights = weights
)
p_interaction <- summary(outcome_int)$coefficients[
  "treated:intub_med_resp", "Pr(>|z|)"]
cat("\nEffect modification (treated × medical-respiratory) p:",
    round(p_interaction, 4), "\n")


## ── 6. IPW SENSITIVITY ANALYSIS ──────────────────────────────────────────────

set.seed(237)
ipw_out <- weightit(
  psm_formula,
  data     = psm_data,
  method   = "ps",
  estimand = "ATE"
)
summary(ipw_out)

ipw_model <- glm(
  survival_12mo ~ treated,
  data    = psm_data,
  family  = binomial,
  weights = ipw_out$weights
)
or_ipw <- exp(coef(ipw_model)["treated"])
ci_ipw <- exp(confint(ipw_model)["treated", ])
p_ipw  <- summary(ipw_model)$coefficients["treated", "Pr(>|z|)"]

cat("\nIPW OR:", round(or_ipw, 3), "\n")
cat("IPW 95% CI:", round(ci_ipw, 3), "\n")
cat("IPW p-value:", round(p_ipw, 4), "\n")


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
    x     = max(caregivers_nonzero$caregiver_fe_rate) * 0.62,
    y     = b_hat + 0.025,
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
cat(sprintf("Survival — high FE caregiver (>= %.1f%%):  %.1f%%\n",
            analysis_data$fe_q75[1] * 100,
            mean(matched_data$survival_12mo[matched_data$treated == 1],
                 na.rm = TRUE) * 100))
cat(sprintf("Survival — low FE caregiver  (<= %.1f%%):  %.1f%%\n",
            analysis_data$fe_q25[1] * 100,
            mean(matched_data$survival_12mo[matched_data$treated == 0],
                 na.rm = TRUE) * 100))
cat("\n")
cat(sprintf("PSM OR (unadjusted):        %.3f (%.3f-%.3f), p = %.4f\n",
            or_unadj, ci_unadj[1], ci_unadj[2], p_unadj))
cat(sprintf("PSM OR (doubly adjusted):   %.3f (%.3f-%.3f), p = %.4f\n",
            or_adj, ci_adj[1], ci_adj[2], p_adj))
cat(sprintf("IPW OR:                     %.3f (%.3f-%.3f), p = %.4f\n",
            or_ipw, ci_ipw[1], ci_ipw[2], p_ipw))
cat(sprintf("Effect modification p:      %.4f\n", p_interaction))
cat("\n")
cat(sprintf("Partial r (FE rate | vol):  %.3f, p = %.4f\n",
            pcor_fe$estimate, pcor_fe$p.value))
cat(sprintf("Partial r (volume | FE):    %.3f, p = %.4f\n",
            pcor_vol$estimate, pcor_vol$p.value))
cat("\n")
cat(sprintf("Survival floor b:           %.3f (95%% CI %.3f-%.3f)\n",
            b_hat, b_ci[1], b_ci[2]))
cat(sprintf("Decay rate k:               %.1f  (95%% CI %.1f-%.1f)\n",
            k_hat, k_ci[1], k_ci[2]))
cat("Interpretation: ~59% survive regardless of caregiver FE rate;\n")
cat("~41% represent outcomes potentially sensitive to caregiver performance.\n")
cat("────────────────────────────────────────────────────────────────────────\n")

## Save results
save(match_out, matched_data, ipw_out,
     or_unadj, ci_unadj, p_unadj,
     or_adj, ci_adj, p_adj,
     or_ipw, ci_ipw, p_ipw,
     pcor_fe, pcor_vol,
     nls_fit, b_hat, k_hat, b_ci, k_ci,
     file = "../data/psm_results.RData")
cat("PSM results saved to ../data/psm_results.RData\n")
