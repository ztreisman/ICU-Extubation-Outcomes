################################################################################
##
##  ICU Extubation Outcomes — Script 02: Patient and Caregiver
##  Characterization by FE Rate Strata
##
##  Depends on: 05_msm.R (msm_results.RData — provides analysis_base and
##              psm_ccsr_codes)
##
##  Purpose: Describe how patient case mix and caregiver profiles differ across
##  FE rate strata within each ICU group. Characterizes the volume-confounding
##  of the 0%-FE caregiver stratum identified in the MSM analysis.
##
##  Key question: Do patients of 0%-FE caregivers have better outcomes because
##  of systematic case mix differences (lower illness severity, simpler
##  procedures) or because of genuine caregiver performance differences?
##
################################################################################


## ── 0. SETUP ─────────────────────────────────────────────────────────────────

library(dplyr)
library(tidyr)
library(tibble)
library(stringr)
library(ggplot2)

## ppcor (loaded by script 03) pulls in MASS, which masks dplyr::select.
## Pin select to dplyr explicitly so it survives the session.
select <- dplyr::select

load("../data/cohorts.RData")
# Loads: analysis_base, psm_ccsr_codes, unit_groups, explicit_extubations, ...


## ── 1. PATIENT CHARACTERISTICS BY FE RATE QUARTILE ───────────────────────────
##
##  For each ICU group, summarise patient demographics and illness severity
##  across FE rate quartiles. This replicates the "Table 2" split by exposure
##  group and allows direct inspection of the imbalance the MSM addresses.

cat("\n── Patient characteristics by FE rate quartile ──────────────────────────\n")

quartile_char <- analysis_base |>
  filter(unit_group != "Other") |>
  group_by(unit_group) |>
  mutate(fe_quartile = ntile(caregiver_fe_rate, 4)) |>
  ungroup() |>
  group_by(unit_group, fe_quartile) |>
  summarise(
    n              = n(),
    n_caregivers   = n_distinct(caregiver_id),
    fe_range       = paste0(round(min(caregiver_fe_rate) * 100, 1), "–",
                            round(max(caregiver_fe_rate) * 100, 1), "%"),
    age_mean_sd    = sprintf("%.1f (%.1f)",
                             mean(anchor_age,    na.rm = TRUE),
                             sd(anchor_age,      na.rm = TRUE)),
    pct_male       = round(mean(gender == "M",   na.rm = TRUE) * 100, 1),
    sofa_mean_sd   = sprintf("%.1f (%.1f)",
                             mean(sofa,          na.rm = TRUE),
                             sd(sofa,            na.rm = TRUE)),
    charlson_mean_sd = sprintf("%.1f (%.1f)",
                               mean(charlson,    na.rm = TRUE),
                               sd(charlson,      na.rm = TRUE)),
    vent_h_med_iqr = sprintf("%d [%d–%d]",
                             as.integer(median(vent_hours, na.rm = TRUE)),
                             as.integer(quantile(vent_hours, 0.25, na.rm = TRUE)),
                             as.integer(quantile(vent_hours, 0.75, na.rm = TRUE))),
    pct_surgical   = round(
      mean(intubation_type == "surgical", na.rm = TRUE) * 100, 1
    ),
    mortality_12mo = round(mean(mortality_12mo, na.rm = TRUE) * 100, 1),
    .groups        = "drop"
  ) |>
  arrange(unit_group, fe_quartile)

print(quartile_char, n = Inf, width = Inf)
cat(paste0(rep("─", 72), collapse = ""), "\n")


## ── 2. CVICU DETAILED CHARACTERIZATION ───────────────────────────────────────
##
##  CVICU receives detailed treatment because:
##    - The MSM signal (if any) is CVICU-specific
##    - 64.7% of CVICU patients are with 0%-FE caregivers
##    - Q25 = 0%, so the "Q1" group is really a categorical 0%-FE stratum
##
##  Strata: exactly 0% vs 0–median vs median–Q75 vs >Q75 of the
##  non-zero FE rate distribution.

cvicu_base <- analysis_base |>
  filter(unit_group == "CVICU") |>
  mutate(
    fe_nonzero_q50 = quantile(caregiver_fe_rate[caregiver_fe_rate > 0],
                               0.50, na.rm = TRUE),
    fe_nonzero_q75 = quantile(caregiver_fe_rate[caregiver_fe_rate > 0],
                               0.75, na.rm = TRUE),
    fe_stratum = case_when(
      caregiver_fe_rate == 0                       ~ "0%",
      caregiver_fe_rate <= fe_nonzero_q50          ~ ">0–median",
      caregiver_fe_rate <= fe_nonzero_q75          ~ "median–Q75",
      TRUE                                         ~ ">Q75"
    ),
    fe_stratum = factor(fe_stratum,
                        levels = c("0%", ">0–median", "median–Q75", ">Q75"))
  )

cat("\n── CVICU patient characterization by FE rate stratum ───────────────────\n")
cat("FE rate: nonzero median =",
    round(unique(cvicu_base$fe_nonzero_q50) * 100, 2), "%,",
    "nonzero Q75 =",
    round(unique(cvicu_base$fe_nonzero_q75) * 100, 2), "%\n\n")

cvicu_char <- cvicu_base |>
  group_by(fe_stratum) |>
  summarise(
    n              = n(),
    n_caregivers   = n_distinct(caregiver_id),
    fe_mean_pct    = round(mean(caregiver_fe_rate) * 100, 2),
    age_mean_sd    = sprintf("%.1f (%.1f)",
                             mean(anchor_age,  na.rm = TRUE),
                             sd(anchor_age,    na.rm = TRUE)),
    pct_male       = round(mean(gender == "M", na.rm = TRUE) * 100, 1),
    sofa_mean_sd   = sprintf("%.1f (%.1f)",
                             mean(sofa,        na.rm = TRUE),
                             sd(sofa,          na.rm = TRUE)),
    charlson_mean_sd = sprintf("%.1f (%.1f)",
                               mean(charlson,  na.rm = TRUE),
                               sd(charlson,    na.rm = TRUE)),
    vent_h_med_iqr = sprintf("%d [%d–%d]",
                             as.integer(median(vent_hours, na.rm = TRUE)),
                             as.integer(quantile(vent_hours, 0.25, na.rm = TRUE)),
                             as.integer(quantile(vent_hours, 0.75, na.rm = TRUE))),
    pct_surgical   = round(
      mean(intubation_type == "surgical", na.rm = TRUE) * 100, 1
    ),
    mortality_12mo_pct = round(mean(mortality_12mo, na.rm = TRUE) * 100, 1),
    .groups        = "drop"
  )

print(cvicu_char, n = Inf, width = Inf)
cat(paste0(rep("─", 72), collapse = ""), "\n")


## ── 3. TOP CCSR DIAGNOSES BY CVICU FE STRATUM ───────────────────────────────
##
##  Prevalence of each CCSR diagnostic flag by FE rate stratum.
##  Differences here explain the imbalance the MSM is adjusting for.

cat("\n── CCSR diagnosis prevalence by CVICU FE stratum (>= 15%) ─────────────\n")

ccsr_by_stratum <- cvicu_base |>
  select(subject_id, fe_stratum, all_of(psm_ccsr_codes)) |>
  pivot_longer(all_of(psm_ccsr_codes),
               names_to = "ccsr_code", values_to = "present") |>
  group_by(fe_stratum, ccsr_code) |>
  summarise(
    n_patients  = sum(present, na.rm = TRUE),
    prevalence  = round(mean(present, na.rm = TRUE), 3),
    .groups     = "drop"
  ) |>
  filter(prevalence >= 0.15) |>
  arrange(fe_stratum, desc(prevalence))

print(ccsr_by_stratum, n = Inf)
cat(paste0(rep("─", 72), collapse = ""), "\n")

## CCSR codes that vary most across strata (largest range of prevalence)
cat("\n── CCSR codes with largest prevalence variation across strata ──────────\n")

ccsr_variation <- cvicu_base |>
  select(fe_stratum, all_of(psm_ccsr_codes)) |>
  pivot_longer(all_of(psm_ccsr_codes),
               names_to = "ccsr_code", values_to = "present") |>
  group_by(ccsr_code, fe_stratum) |>
  summarise(prevalence = mean(present, na.rm = TRUE), .groups = "drop") |>
  group_by(ccsr_code) |>
  summarise(
    prev_range  = max(prevalence) - min(prevalence),
    prev_min    = round(min(prevalence), 3),
    prev_max    = round(max(prevalence), 3),
    .groups     = "drop"
  ) |>
  arrange(desc(prev_range)) |>
  head(15)

print(ccsr_variation, n = Inf)
cat(paste0(rep("─", 72), collapse = ""), "\n")


## ── 4. CAREGIVER PROFILES BY FE RATE STRATUM (CVICU) ─────────────────────────
##
##  Caregiver-level summary: volume, mean patient characteristics, outcome.
##  Separates whether the 0%-FE stratum represents high-volume specialists
##  (low FE because they only do straightforward cases) or simply providers
##  who happen to have had no failures yet.

cat("\n── CVICU caregiver profiles by FE rate stratum ─────────────────────────\n")

caregiver_fe_strata <- cvicu_base |>
  distinct(caregiver_id, caregiver_fe_rate, fe_stratum)

caregiver_profiles <- cvicu_base |>
  group_by(caregiver_id) |>
  summarise(
    fe_stratum     = first(fe_stratum),
    n_patients     = n(),
    fe_rate        = first(caregiver_fe_rate),
    age_mean       = round(mean(anchor_age,  na.rm = TRUE), 1),
    sofa_mean      = round(mean(sofa,        na.rm = TRUE), 2),
    pct_surgical   = round(
      mean(intubation_type == "surgical", na.rm = TRUE) * 100, 1
    ),
    mortality_12mo = round(mean(mortality_12mo, na.rm = TRUE), 3),
    .groups        = "drop"
  )

caregiver_summary_by_stratum <- caregiver_profiles |>
  group_by(fe_stratum) |>
  summarise(
    n_caregivers   = n(),
    vol_med_iqr    = sprintf("%d [%d–%d]",
                             as.integer(median(n_patients)),
                             as.integer(quantile(n_patients, 0.25)),
                             as.integer(quantile(n_patients, 0.75))),
    vol_max        = max(n_patients),
    age_mean       = round(mean(age_mean),   1),
    sofa_mean      = round(mean(sofa_mean),  2),
    pct_surgical   = round(mean(pct_surgical), 1),
    mortality_mean = round(mean(mortality_12mo), 3),
    .groups        = "drop"
  )

print(caregiver_summary_by_stratum, n = Inf, width = Inf)
cat(paste0(rep("─", 72), collapse = ""), "\n")


## ── 5. FIGURES ────────────────────────────────────────────────────────────────

## 5a. SOFA and age distributions by CVICU FE stratum (boxplot)
##  vent_hours is log10-transformed before pivoting; its distribution is
##  heavily right-skewed (range ~1–500h) and collapses on a linear scale.
cvicu_plot_long <- cvicu_base |>
  mutate(log10_vent_hours = log10(vent_hours)) |>
  select(fe_stratum, sofa, anchor_age, log10_vent_hours, charlson) |>
  pivot_longer(-fe_stratum, names_to = "variable", values_to = "value") |>
  mutate(
    variable = recode(variable,
                      sofa             = "SOFA score",
                      anchor_age       = "Age (years)",
                      log10_vent_hours = "Vent hours (log₁₀)",
                      charlson         = "Charlson index")
  )

ggplot(cvicu_plot_long,
       aes(x = fe_stratum, y = value, fill = fe_stratum)) +
  facet_wrap(~variable, scales = "free_y", nrow = 2) +
  geom_boxplot(outlier.size = 0.5, outlier.alpha = 0.4, linewidth = 0.5) +
  scale_fill_brewer(palette = "Blues", direction = 1) +
  labs(
    title    = "CVICU Patient Characteristics by Caregiver FE Rate Stratum",
    subtitle = paste0(
      "0% = exactly zero FE rate caregiver. Other strata from non-zero ",
      "distribution.\nBoxplot: median, IQR, 1.5×IQR whiskers."
    ),
    x = "Caregiver FE rate stratum",
    y = NULL
  ) +
  theme_minimal() +
  theme(
    legend.position  = "none",
    strip.text       = element_text(size = 10, face = "bold"),
    axis.text.x      = element_text(size = 9),
    plot.title       = element_text(size = 11, face = "bold"),
    plot.subtitle    = element_text(size = 8.5, color = "gray40")
  )

ggsave("../figures/cvicu_strata_covariates.png",
       width = 9, height = 6, dpi = 150)
cat("CVICU strata covariate plot saved.\n")


## 5b. 12-month mortality by FE rate stratum, all ICU groups

mort_by_quartile <- analysis_base |>
  filter(unit_group != "Other") |>
  group_by(unit_group) |>
  mutate(fe_quartile = ntile(caregiver_fe_rate, 4),
         fe_label = sprintf("Q%d", fe_quartile)) |>
  group_by(unit_group, fe_quartile, fe_label) |>
  summarise(
    n              = n(),
    mortality_mean = mean(mortality_12mo, na.rm = TRUE),
    mortality_se   = sd(mortality_12mo,  na.rm = TRUE) / sqrt(n()),
    .groups        = "drop"
  )

ggplot(mort_by_quartile,
       aes(x = fe_label, y = mortality_mean,
           ymin = mortality_mean - 1.96 * mortality_se,
           ymax = mortality_mean + 1.96 * mortality_se)) +
  facet_wrap(~unit_group) +
  geom_col(fill = "#4472C4", alpha = 0.75, width = 0.6) +
  geom_errorbar(width = 0.2, linewidth = 0.6, color = "gray30") +
  geom_text(aes(label = sprintf("n=%d", n), y = 0.01),
            size = 2.8, color = "white", fontface = "bold", vjust = 0) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1),
                     expand = expansion(mult = c(0, 0.05))) +
  labs(
    title    = "12-Month Mortality by Caregiver FE Rate Quartile",
    subtitle = paste0(
      "Q1 = lowest FE rate caregivers. Error bars = 95% CI (normal approx).\n",
      "Unadjusted; see MSM results in Script 03 for covariate-adjusted estimates."
    ),
    x = "Caregiver FE rate quartile",
    y = "12-Month Mortality"
  ) +
  theme_minimal() +
  theme(
    strip.text       = element_text(size = 10, face = "bold"),
    plot.title       = element_text(size = 11, face = "bold"),
    plot.subtitle    = element_text(size = 8.5, color = "gray40"),
    panel.grid.major.x = element_blank()
  )

ggsave("../figures/mortality_by_quartile.png",
       width = 10, height = 4, dpi = 150)
cat("Mortality by quartile plot saved.\n")
