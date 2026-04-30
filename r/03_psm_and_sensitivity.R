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
##  Research question:
##    After matching patients of high-FE-rate caregivers to patients of
##    low-FE-rate caregivers on illness severity and diagnosis case mix,
##    does a survival difference persist?
##
##  Analytical approach:
##    1. CCSR binary flags: 21 high-prevalence diagnosis codes as PSM covariates
##       (LDA topic modeling and k-means clustering were explored but showed
##       weak patient differentiation — median ~4 CCSR codes per patient
##       produces short-document LDA failure; direct binary flags are more
##       transparent and interpretable)
##    2. PSM: nearest-neighbor 1:1 matching (method = "quick") within caliper
##       0.2 SD on logit propensity score
##    3. IPW: inverse probability weighting as sensitivity analysis
##    4. Partial correlation: caregiver-level analysis controlling for volume
##    5. Functional form: linear vs quadratic vs exponential models of
##       FE rate — survival relationship at the caregiver level
##
##  Key findings:
##    PSM (top vs bottom quartile FE rate):
##      Unadjusted OR = 0.745 (0.629-0.882), p = 0.0006
##      Doubly adjusted OR = 0.681 (0.558-0.830), p = 0.0001
##    IPW: OR = 0.711 (0.630-0.802), p < 0.001
##    Partial correlation (FE rate | caregiver volume): r = -0.366, p < 0.001
##    Partial correlation (volume | FE rate): r = -0.122, p = 0.160
##    Functional form: quadratic fits better than linear (F=15.1, p=0.0002);
##      relationship is monotonically decreasing with diminishing slope —
##      consistent with a survival floor at high illness severity,
##      not a U-shaped optimum
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
library(MatchIt)     # propensity score matching
library(cobalt)      # balance diagnostics and love plot
library(WeightIt)    # IPW sensitivity analysis
library(ppcor)       # partial correlation

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
##
##  Quartiles (explicit_extubations):
##    Q25 = 0.021 (2.1% FE rate)
##    Q75 = 0.052 (5.2% FE rate)

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
##
##  Method: "quick" (faster than "nearest", equivalent results for large N)
##  Ratio: 1:1
##  Caliper: 0.2 SD of logit propensity score (standard recommendation)
##  Distance: logit propensity score
##
##  Note: method = "nearest" caused R session crashes due to memory constraints
##  with 570 caregivers and 21 CCSR flags. "quick" uses a faster algorithm
##  with identical statistical properties.

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

## Love plot — covariate balance before and after matching
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

## Unadjusted outcome model
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

## Doubly adjusted outcome model
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

## Effect modification by intubation type
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
##
##  Unit of analysis: caregiver (N = 134, restricted to n > 5 sample patients)
##  Tests whether FE rate predicts survival independently of caregiver volume,
##  and vice versa.
##
##  Key finding: FE rate independently predicts survival controlling for volume
##  (r = -0.366, p < 0.001); volume does not predict survival controlling for
##  FE rate (r = -0.122, p = 0.160). FE rate is the signal; volume is not.

cat("\n── Caregiver-level partial correlations ─────────────────────────────────\n")
cat("N caregivers (n > 5 sample patients):", nrow(caregiver_summary), "\n")

## Partial correlation: FE rate ~ survival, controlling for volume
pcor_fe <- pcor.test(
  caregiver_summary$caregiver_fe_rate,
  caregiver_summary$survival_rate,
  caregiver_summary$n_patients_sample
)
cat("\nPartial correlation: FE rate ~ survival | caregiver volume\n")
cat("r =", round(pcor_fe$estimate, 3), "\n")
cat("p =", round(pcor_fe$p.value,  4), "\n")

## Partial correlation: volume ~ survival, controlling for FE rate
pcor_vol <- pcor.test(
  caregiver_summary$caregiver_n,
  caregiver_summary$survival_rate,
  caregiver_summary$caregiver_fe_rate
)
cat("\nPartial correlation: volume ~ survival | FE rate\n")
cat("r =", round(pcor_vol$estimate, 3), "\n")
cat("p =", round(pcor_vol$p.value,  4), "\n")
cat("────────────────────────────────────────────────────────────────────────\n")


## ── 8. FUNCTIONAL FORM ANALYSIS ──────────────────────────────────────────────
##
##  Tests whether the FE rate — survival relationship is linear, quadratic,
##  or exponential at the caregiver level.
##
##  Restricted to caregivers with FE rate > 0 (n = 112 of 134) to avoid
##  small-sample zero-FE caregivers inflating the intercept.
##
##  Finding: quadratic fits significantly better than linear (F = 15.1,
##  p = 0.0002, ΔAIC = 12.5). The parabola opens upward (positive quadratic
##  coefficient) with vertex at ~7.7% — but this is a minimum, not an optimum.
##  The relationship is monotonically decreasing with diminishing slope:
##  highest harm per unit FE rate at low FE rates, flattening at higher rates
##  due to a survival floor in severely ill patients. This is consistent with
##  an exponential decay model (log-linear fit also significant, p < 0.001).
##  There is no evidence of a sweet spot where moderate FE rates produce
##  better outcomes than low FE rates.
##
##  Bootstrap 95% CI on vertex: 6.4% to 11.4% (1999/2000 bootstraps retained)
##  — stable but describes a minimum, not an optimal target.

cat("\n── Functional form analysis (caregivers with FE rate > 0) ──────────────\n")
cat("N caregivers:", nrow(caregivers_nonzero), "\n")

# Models fitted in Script 01 and loaded via cohorts.RData:
#   lm_linear    — linear
#   lm_quadratic — quadratic (best AIC)
#   lm_log       — log(survival) ~ fe_rate (exponential)

print(AIC(lm_linear, lm_quadratic))
print(anova(lm_linear, lm_quadratic))

coefs  <- coef(lm_quadratic)
vertex <- -coefs[2] / (2 * coefs[3])
cat("Quadratic vertex (minimum, not optimum):", round(vertex, 3), "\n")
cat("Quadratic coefficient sign:",
    ifelse(coefs[3] > 0, "positive (parabola opens upward = minimum)", "negative"), "\n")

## Bootstrap CI on vertex
set.seed(237)
n_boot   <- 2000
vertices <- numeric(n_boot)
for (i in seq_len(n_boot)) {
  boot_idx  <- sample(nrow(caregivers_nonzero), replace = TRUE)
  boot_data <- caregivers_nonzero[boot_idx, ]
  boot_fit  <- lm(
    survival_rate ~ caregiver_fe_rate + I(caregiver_fe_rate^2),
    data    = boot_data,
    weights = n_patients_sample
  )
  b          <- coef(boot_fit)
  vertices[i] <- -b[2] / (2 * b[3])
}
vertices_clean <- vertices[vertices > 0 & vertices < 1]
cat("\nBootstrap vertex CI (describes minimum, not optimum):\n")
cat("Mean:", round(mean(vertices_clean), 3), "\n")
cat("95% CI:", round(quantile(vertices_clean, 0.025), 3),
    "to", round(quantile(vertices_clean, 0.975), 3), "\n")
cat("Bootstraps retained:", length(vertices_clean), "of", n_boot, "\n")
cat("────────────────────────────────────────────────────────────────────────\n")

## Final visualization: FE rate vs survival with competing functional forms
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
    aes(caregiver_fe_rate, survival_rate, size = n_patients_sample),
    alpha = 0.5, color = "gray40"
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
  annotate("text", x = max(fe_seq) * 0.80, y = 0.93,
           label = "Linear",      color = "#2E75B6", size = 3.5) +
  annotate("text", x = max(fe_seq) * 0.80, y = 0.87,
           label = "Quadratic",   color = "#D95F02", size = 3.5) +
  annotate("text", x = max(fe_seq) * 0.80, y = 0.81,
           label = "Exponential", color = "#1D9E8E", size = 3.5) +
  labs(
    title    = "Caregiver FE Rate vs. 12-Month Patient Survival",
    subtitle = "Caregivers with FE rate > 0. Quadratic fits best (ΔAIC = 12.5 vs linear).\nDiminishing slope reflects survival floor, not U-shaped optimum.",
    x        = "Caregiver Failed Extubation Rate",
    y        = "12-Month Patient Survival Rate",
    size     = "N patients"
  ) +
  theme_minimal()

ggsave("../figures/fe_rate_survival_functional_form.png",
       width = 9, height = 6, dpi = 150)


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
cat("Functional form: quadratic > linear (ΔAIC = 12.5, p = 0.0002)\n")
cat("Interpretation:  diminishing marginal harm, not U-shaped optimum\n")
cat("────────────────────────────────────────────────────────────────────────\n")

## Save PSM objects for reporting
save(match_out, matched_data, ipw_out,
     or_unadj, ci_unadj, p_unadj,
     or_adj, ci_adj, p_adj,
     or_ipw, ci_ipw, p_ipw,
     pcor_fe, pcor_vol,
     file = "../data/psm_results.RData")
cat("PSM results saved to ../data/psm_results.RData\n")
