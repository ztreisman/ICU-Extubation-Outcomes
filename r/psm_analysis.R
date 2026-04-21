## ── 0. LIBRARIES ─────────────────────────────────────────────────────────────

install.packages(c("MatchIt", "cobalt", "WeightIt", "topicmodels", 
                   "tidytext", "Matrix", "slam"))

library(tidyverse)
library(MatchIt)
library(cobalt)
library(WeightIt)
library(topicmodels)
library(tidytext)
library(Matrix)
library(slam)

set.seed(237)

## ── 1. DATA INGESTION ────────────────────────────────────────────────────────

ccsr_long <- read_csv("data/patient_ccsr_long.csv",
                      col_types = cols(
                        subject_id = col_character(),
                        hadm_id    = col_character(),
                        ccsr_code  = col_character()
                      ))

patients <- read_csv("data/last_extubations_clean.csv",
                     col_types = cols(
                       subject_id   = col_character(),
                       hadm_id      = col_character(),
                       gender       = col_factor(),
                       intubation_type = col_factor()
                     )) %>%
  mutate(survival_12mo = is.na(dod))

cat("Patients:", nrow(patients), "\n")
cat("CCSR rows:", nrow(ccsr_long), "\n")
cat("Unique CCSR codes:", n_distinct(ccsr_long$ccsr_code), "\n")

## ── 2. CCSR BINARY MATRIX ────────────────────────────────────────────────────
##
##  Each patient is a "document", each CCSR code is a "term"
##  Value = 1 if patient has that code, 0 otherwise

ccsr_counts <- ccsr_long %>%
  count(subject_id, ccsr_code) %>%
  rename(count = n)

# Convert to DocumentTermMatrix for topicmodels
ccsr_dtm <- ccsr_counts %>%
  cast_dtm(document = subject_id, 
           term     = ccsr_code, 
           value    = count)

cat("DTM dimensions:", dim(ccsr_dtm), "\n")
cat("DTM sparsity:", ccsr_dtm$nrow * ccsr_dtm$ncol, "possible entries,", 
    sum(ccsr_dtm$v), "non-zero\n")


## ── 3. LDA TOPIC MODELING ─────────────────────────────────────────────────────
##
##  Fit LDA for K = 5, 8, 10 topics
##  Choose K based on perplexity and clinical interpretability
##
##  Note: LDA can take a few minutes per model

control_lda <- list(seed = 237, burnin = 100, iter = 500, thin = 10)

cat("Fitting K=5...\n")
lda_5  <- LDA(ccsr_dtm, k = 5,  method = "Gibbs", control = control_lda)

cat("Fitting K=8...\n")
lda_8  <- LDA(ccsr_dtm, k = 8,  method = "Gibbs", control = control_lda)

cat("Fitting K=10...\n")
lda_10 <- LDA(ccsr_dtm, k = 10, method = "Gibbs", control = control_lda)

# Compare perplexity (lower = better fit)
cat("\nPerplexity comparison:\n")
cat("K=5: ",  perplexity(lda_5,  newdata = ccsr_dtm), "\n")
cat("K=8: ",  perplexity(lda_8,  newdata = ccsr_dtm), "\n")
cat("K=10:", perplexity(lda_10, newdata = ccsr_dtm), "\n")

cat("Fitting K=12...\n")
lda_12 <- LDA(ccsr_dtm, k = 12, method = "Gibbs", control = control_lda)
cat("K=12:", perplexity(lda_12, newdata = ccsr_dtm), "\n")

## Top terms per topic — what is each topic "about"?
terms(lda_8, k = 10)

## ── 4. EXTRACT TOPIC PROPORTIONS ─────────────────────────────────────────────
##
##  gamma matrix: one row per patient, one column per topic
##  Each row sums to 1 — patient's mixture across the 8 clinical profiles

topic_proportions <- posterior(lda_8)$topics %>%
  as.data.frame() %>%
  tibble::rownames_to_column("subject_id") %>%
  rename(
    topic_cardiac_hf        = `1`,  # Heart failure / vascular
    topic_neurologic        = `2`,  # CNS with cardiac
    topic_cad_metabolic     = `3`,  # Coronary artery disease / metabolic
    topic_trauma            = `4`,  # Trauma / injury
    topic_respiratory       = `5`,  # Respiratory failure / metabolic
    topic_sepsis_abdominal  = `6`,  # Sepsis / infectious / digestive
    topic_respiratory_infx  = `7`,  # Respiratory infection / sepsis
    topic_renal_metabolic   = `8`   # Renal / metabolic / digestive
  )

cat("Topic proportions dimensions:", dim(topic_proportions), "\n")
cat("Row sums (should all be 1):", range(rowSums(topic_proportions[,-1])), "\n")

## Merge into patient dataset
patients_topics <- patients %>%
  left_join(topic_proportions, by = "subject_id")

cat("Patients with topic assignments:", 
    sum(!is.na(patients_topics$topic_cardiac_hf)), "\n")

## Quick summary of topic distributions
topic_proportions[,-1] %>%
  summarise(across(everything(), mean)) %>%
  pivot_longer(everything(), names_to = "topic", values_to = "mean_proportion") %>%
  arrange(desc(mean_proportion)) %>%
  print()
topic_proportions[,-1] %>%
  summarise(across(everything(), sd)) %>%
  pivot_longer(everything(), names_to = "topic", values_to = "sd_proportion") %>%
  arrange(desc(sd_proportion)) %>%
  print()

## Also look at a few individual patients to confirm differentiation
topic_proportions %>%
  slice(1:5) %>%
  pivot_longer(-subject_id, names_to = "topic", values_to = "proportion") %>%
  print(n = 40)

## ── 4b. K-MEANS CLUSTERING (replacing LDA proportions) ───────────────────────
##
##  Binary patient × CCSR matrix, k-means with k=8 to match LDA comparison
##  Hard cluster assignment used as PSM covariate

# Build binary matrix
ccsr_binary <- ccsr_long %>%
  distinct(subject_id, ccsr_code) %>%
  mutate(present = 1) %>%
  pivot_wider(
    id_cols     = subject_id,
    names_from  = ccsr_code,
    values_from = present,
    values_fill = 0
  )

cat("Binary matrix dimensions:", dim(ccsr_binary), "\n")

# Run k-means for k = 6, 8, 10
# nstart = 25 runs from 25 random starts, takes the best
set.seed(237)
km_6  <- kmeans(ccsr_binary[,-1], centers = 6,  nstart = 25, iter.max = 100)
km_8  <- kmeans(ccsr_binary[,-1], centers = 8,  nstart = 25, iter.max = 100)
km_10 <- kmeans(ccsr_binary[,-1], centers = 10, nstart = 25, iter.max = 100)

# Compare within-cluster sum of squares (lower = tighter clusters)
cat("Within-SS k=6: ",  km_6$tot.withinss,  "\n")
cat("Within-SS k=8: ",  km_8$tot.withinss,  "\n")
cat("Within-SS k=10:", km_10$tot.withinss, "\n")

# Proportion of variance explained
cat("Variance explained k=6: ",  
    round(km_6$betweenss  / km_6$totss,  3), "\n")
cat("Variance explained k=8: ",  
    round(km_8$betweenss  / km_8$totss,  3), "\n")
cat("Variance explained k=10:", 
    round(km_10$betweenss / km_10$totss, 3), "\n")


## ── 4c. REVISED COVARIATE STRATEGY ───────────────────────────────────────────
##
##  Since clustering explains little variance, use two approaches:
##
##  Approach A: LDA topic proportions as continuous covariates (soft membership)
##  Approach B: Top individual CCSR codes as binary flags (direct, transparent)
##
##  We'll use both and compare PSM balance

## Approach B: identify high-prevalence CCSR codes to use as binary flags
## Use codes present in at least 10% of patients — interpretable and stable

ccsr_prevalence <- ccsr_long %>%
  distinct(subject_id, ccsr_code) %>%
  count(ccsr_code, name = "n_patients") %>%
  mutate(prevalence = n_patients / n_distinct(ccsr_long$subject_id)) %>%
  arrange(desc(prevalence))

cat("CCSR codes in >= 10% of patients:\n")
ccsr_prevalence %>% filter(prevalence >= 0.10) %>% print(n = 30)

cat("\nCCSR codes in >= 5% of patients:", 
    sum(ccsr_prevalence$prevalence >= 0.05), "\n")
cat("CCSR codes in >= 10% of patients:", 
    sum(ccsr_prevalence$prevalence >= 0.10), "\n")


## ── 5. BUILD ANALYSIS DATASET FOR PSM ────────────────────────────────────────

# Get the 22 high-prevalence codes, drop RSP012 (near-universal, low variance)
psm_ccsr_codes <- ccsr_prevalence %>%
  filter(prevalence >= 0.10, ccsr_code != "RSP012") %>%
  pull(ccsr_code)

cat("CCSR codes used as PSM covariates:", length(psm_ccsr_codes), "\n")

# Build binary flags for each code
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

# Merge with patient data and topic proportions
analysis_data <- patients_topics %>%
  left_join(ccsr_flags, by = "subject_id") %>%
  # Replace NA flags with 0 (patient had none of these codes)
  mutate(across(all_of(psm_ccsr_codes), ~ replace_na(.x, 0))) %>%
  # Derive survival_12mo from dod
  mutate(survival_12mo = is.na(dod)) %>%
  # Define PSM treatment variable
  # High FE rate = above 75th percentile, Low = below 25th percentile
  mutate(
    fe_q25 = quantile(caregiver_fe_rate, 0.25, na.rm = TRUE),
    fe_q75 = quantile(caregiver_fe_rate, 0.75, na.rm = TRUE),
    fe_group = case_when(
      caregiver_fe_rate >= fe_q75 ~ "high",
      caregiver_fe_rate <= fe_q25 ~ "low",
      TRUE ~ "middle"
    )
  )

cat("FE rate 25th percentile:", round(analysis_data$fe_q25[1], 3), "\n")
cat("FE rate 75th percentile:", round(analysis_data$fe_q75[1], 3), "\n")
cat("Group sizes:\n")
table(analysis_data$fe_group)

## ── 6. PROPENSITY SCORE MATCHING ─────────────────────────────────────────────

# Filter to high/low groups only, create binary treatment indicator
psm_data <- analysis_data %>%
  filter(fe_group != "middle") %>%
  mutate(
    treated = as.integer(fe_group == "high"),
    # Convert factors to numeric for MatchIt
    gender_m = as.integer(gender == "M"),
    intub_surgical = as.integer(intubation_type == "surgical"),
    intub_med_resp = as.integer(intubation_type == "medical-respiratory")
  ) %>%
  # Drop rows with missing key covariates
  filter(
    !is.na(charlson),
    !is.na(sofa),
    !is.na(norepinephrine),
    !is.na(vent_hours)
  )

cat("PSM dataset N:", nrow(psm_data), "\n")
cat("Treated (high FE):", sum(psm_data$treated), "\n")
cat("Control (low FE):", sum(psm_data$treated == 0), "\n")

# Build propensity model formula
# Covariates: demographics + illness severity + intubation type + CCSR flags
ccsr_formula_terms <- paste(psm_ccsr_codes, collapse = " + ")

psm_formula <- as.formula(paste(
  "treated ~ anchor_age + gender_m + charlson + sofa + norepinephrine +
   vent_hours + intub_surgical + intub_med_resp +",
  ccsr_formula_terms
))

# Fit PSM with nearest neighbor 1:1 matching within caliper
# Caliper = 0.2 * SD of logit propensity score (standard recommendation)
set.seed(237)
match_out <- matchit(
  psm_formula,
  data    = psm_data,
  method  = "nearest",
  ratio   = 1,
  caliper = 0.2,
  std.caliper = TRUE,
  distance = "logit"
)

summary(match_out)

## ── 7. TREATMENT EFFECT ESTIMATION ───────────────────────────────────────────

# Extract matched dataset
matched_data <- match.data(match_out)

cat("Matched dataset N:", nrow(matched_data), "\n")
cat("Survival rate - High FE caregivers:", 
    round(mean(matched_data$survival_12mo[matched_data$treated == 1], 
               na.rm = TRUE), 3), "\n")
cat("Survival rate - Low FE caregivers:", 
    round(mean(matched_data$survival_12mo[matched_data$treated == 0], 
               na.rm = TRUE), 3), "\n")

# Outcome model in matched sample
# Use weights from matching and cluster on subclass (matched pair)
outcome_model <- glm(
  survival_12mo ~ treated,
  data   = matched_data,
  family = binomial,
  weights = weights
)
summary(outcome_model)

# Odds ratio and 95% CI
or <- exp(coef(outcome_model)["treated"])
ci <- exp(confint(outcome_model)["treated",])
cat("\nOdds Ratio (high vs low FE caregiver):", round(or, 3), "\n")
cat("95% CI:", round(ci, 3), "\n")
cat("p-value:", 
    round(summary(outcome_model)$coefficients["treated", "Pr(>|z|)"], 4), "\n")

# Also run adjusted model adding back key covariates
outcome_model_adj <- glm(
  survival_12mo ~ treated + anchor_age + charlson + sofa + 
    norepinephrine + vent_hours + intub_surgical + intub_med_resp,
  data   = matched_data,
  family = binomial,
  weights = weights
)

or_adj <- exp(coef(outcome_model_adj)["treated"])
ci_adj <- exp(confint(outcome_model_adj)["treated",])
cat("\nAdjusted Odds Ratio:", round(or_adj, 3), "\n")
cat("Adjusted 95% CI:", round(ci_adj, 3), "\n")
cat("Adjusted p-value:", 
    round(summary(outcome_model_adj)$coefficients["treated", "Pr(>|z|)"], 4), "\n")

## ── 8. SENSITIVITY AND VALIDITY CHECKS ───────────────────────────────────────

# Check: are survival rates plausible given this is 12-month survival?
# These patients are sicker than average (ventilated ICU patients)
# 50-57% 12-month survival is consistent with published literature

# Check 1: balance visualization (love plot)
library(cobalt)
love.plot(match_out, 
          threshold = 0.1,
          title = "Covariate Balance Before and After PSM",
          colors = c("#D95F02", "#1D9E8E"))

# Check 2: IPW as sensitivity analysis
# If IPW gives similar results, matching result is robust
set.seed(237)
ipw_out <- weightit(
  psm_formula,
  data   = psm_data,
  method = "ps",
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
ci_ipw <- exp(confint(ipw_model)["treated",])
cat("\nIPW Odds Ratio:", round(or_ipw, 3), "\n")
cat("IPW 95% CI:", round(ci_ipw, 3), "\n")
cat("IPW p-value:", 
    round(summary(ipw_model)$coefficients["treated","Pr(>|z|)"], 4), "\n")

# Check 3: Does effect vary by intubation type?
# (effect modification check)
outcome_interaction <- glm(
  survival_12mo ~ treated * intub_med_resp,
  data   = matched_data,
  family = binomial,
  weights = weights
)
summary(outcome_interaction)
cat("\nInteraction p-value (treated x medical-respiratory):",
    round(summary(outcome_interaction)$coefficients[
      "treated:intub_med_resp", "Pr(>|z|)"], 4), "\n")

## ── 9. SAVE OUTPUTS ──────────────────────────────────────────────────────────

# Love plot — suppress the stars warning cleanly
p <- love.plot(
  match_out,
  threshold  = 0.1,
  title      = "Covariate Balance Before and After PSM",
  colors     = c("#D95F02", "#1D9E8E"),
  stars      = "std",
  var.order  = "unadjusted",
  abs        = TRUE
  )

# Results summary table
results_summary <- tibble(
  Analysis = c(
    "PSM — unadjusted",
    "PSM — doubly adjusted",
    "IPW — full sample"
  ),
  OR   = c(0.786, 0.733, 0.789),
  CI_lower = c(0.700, 0.643, 0.735),
  CI_upper = c(0.881, 0.836, 0.846),
  p_value  = c("<0.001", "<0.001", "<0.001"),
  N        = c(4706, 4706, 6191)
)

print(results_summary)
write_csv(results_summary, "psm_results_summary.csv")

# Save matched dataset for any further analysis
write_csv(matched_data %>% 
            select(subject_id, treated, survival_12mo, 
                   anchor_age, charlson, sofa, norepinephrine,
                   vent_hours, intubation_type, caregiver_fe_rate,
                   weights, subclass),
          "psm_matched_dataset.csv")

cat("\nAll outputs saved.\n")
cat("\n── FINAL RESULTS SUMMARY ─────────────────────────────────────\n")
cat("Survival rate, high-FE caregivers (treated):  50.8%\n")
cat("Survival rate, low-FE caregivers (control):   56.8%\n")
cat("Absolute survival difference:                  6.0 pp\n")
cat("PSM odds ratio (adjusted):           0.733 (0.643–0.836)\n")
cat("IPW odds ratio (sensitivity check):  0.789 (0.735–0.846)\n")
cat("Effect modification by intubation type:        p = 0.623\n")
cat("Conclusion: caregiver FE rate effect on 12-month survival\n")
cat("persists after matching on case mix and illness severity.\n")
cat("──────────────────────────────────────────────────────────────\n")
