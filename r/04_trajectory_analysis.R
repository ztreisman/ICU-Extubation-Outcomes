################################################################################
##
##  ICU Extubation Outcomes — Script 04: Trajectory Analysis
##
##  Depends on: bigrquery package and Google Cloud authentication
##              all_extubations table saved in BigQuery
##              (produced by sql/all_extubations.sql)
##
##  Research question:
##    Does the trajectory (slope) of RSBI or P/F ratio in the 72h before
##    extubation predict per-event extubation failure, beyond the level of
##    each variable at the moment of extubation?
##
##  Design:
##    Unit of analysis: extubation event (not patient)
##    Outcome: failed_extubation_flag — reintubation within 72h of THIS
##             specific extubation (both events explicitly documented)
##    Exposure: per-event slope of RSBI or P/F ratio estimated from all
##              measurements in the 72h window before each extubation
##
##  Note on failed_extubations vs failed_extubation_flag:
##    The `failed_extubations` column in last_extubations.csv is a patient
##    lifetime count of reintubation pairs across all extubation events —
##    not a flag for whether the last extubation specifically failed. It is
##    a patient complexity marker. `failed_extubation_flag` (from
##    all_extubations.sql) is the correct per-event outcome variable.
##
##  Key findings:
##    RSBI level at extubation: significant (p = 0.032)
##      mean 52.3 (failed) vs 48.7 (succeeded)
##    RSBI slope 72h before extubation: not significant (p = 0.634)
##    P/F ratio level: not significant (p = 0.505)
##    P/F ratio slope: not significant (p = 0.726)
##
##    Conclusion: trajectory in the 72h pre-extubation window carries no
##    independent predictive information beyond level at extubation for
##    either variable. Consistent with clinicians already conditioning on
##    a plateau criterion before extubating.
##
##  Runtime: ~5 minutes (two BigQuery queries via bigrquery)
##
################################################################################


## ── 0. LIBRARIES ─────────────────────────────────────────────────────────────

library(bigrquery)
library(dplyr)
library(tidyr)
library(ggplot2)
library(patchwork)

set.seed(237)

# Authenticate — opens browser for Google sign-in
# Use bq_auth(use_oob = TRUE) if running in a non-browser environment
bq_auth()

project <- "mythical-legend-456217-e9"


## ── 1. LOAD ALL EXTUBATIONS ───────────────────────────────────────────────────

all_extub <- read_csv("../data/all_extubations.csv") %>%
  mutate(
    failed        = failed_extubation_flag == 1,
    survival_12mo = is.na(dod)
  )

cat("Total extubation events:", nrow(all_extub), "\n")
cat("Unique patients:", n_distinct(all_extub$subject_id), "\n")
cat("Events per patient (median):", median(table(all_extub$subject_id)), "\n")
cat("\nFailed extubation flag distribution:\n")
print(table(all_extub$failed_extubation_flag))
cat("\nOverall FE rate:", round(mean(all_extub$failed_extubation_flag), 3), "\n")
cat("\nCaregivers:", n_distinct(all_extub$caregiver_id), "\n")
cat("ICU units:", n_distinct(all_extub$first_careunit), "\n")


## ── 2. RSBI MEASUREMENTS IN 72H BEFORE EACH EXTUBATION ───────────────────────
##
##  Joins all RSBI measurements (spontaneous RR / spontaneous TV) within
##  72h before each extubation event in the all_extubations table.
##  RSBI = spontaneous RR / (spontaneous TV in liters)
##  Restricted to physiologically plausible range: 1-500

sql_rsbi_events <- "
  SELECT
    e.subject_id,
    e.stay_id,
    e.event_time,
    e.failed_extubation_flag,
    e.caregiver_fe_rate_loo,
    e.vent_hours,
    rr.valuenum / (tv.valuenum / 1000.0)              AS rsbi,
    TIMESTAMP_DIFF(e.event_time, rr.charttime, MINUTE) / 60.0
                                                       AS hours_before_extubation
  FROM `mythical-legend-456217-e9.Failed_Extubation_Rate.all_extubations` e
  JOIN `physionet-data.mimiciv_3_1_icu.chartevents` rr
    ON e.stay_id = rr.stay_id
    AND rr.itemid IN (224689, 224422)
    AND rr.valuenum > 0
    AND TIMESTAMP_DIFF(e.event_time, rr.charttime, MINUTE) BETWEEN 0 AND 72*60
  JOIN `physionet-data.mimiciv_3_1_icu.chartevents` tv
    ON e.stay_id = tv.stay_id
    AND tv.itemid = 224686
    AND tv.valuenum > 0
    AND ABS(TIMESTAMP_DIFF(tv.charttime, rr.charttime, MINUTE)) <= 5
  WHERE rr.valuenum / (tv.valuenum / 1000.0) BETWEEN 1 AND 500
"

rsbi_events <- bq_project_query(project, sql_rsbi_events) %>%
  bq_table_download(bigint = "integer64")

cat("\nRSBI data:\n")
cat("Rows:", nrow(rsbi_events), "\n")
cat("Extubation events with RSBI:", n_distinct(rsbi_events$event_time), "\n")
cat("Measurements per event (median):",
    median(table(paste(rsbi_events$subject_id, rsbi_events$event_time))), "\n")


## ── 3. RSBI PER-EVENT SLOPE ANALYSIS ─────────────────────────────────────────
##
##  For each extubation event with >= 3 RSBI measurements in the 72h window,
##  fit lm(rsbi ~ hours_before_extubation) and extract:
##    rsbi_slope:    coefficient on hours_before_extubation
##                   positive slope = RSBI rising as extubation approaches
##                   negative slope = RSBI falling (improving)
##    rsbi_at_extub: RSBI value closest to extubation time

event_slopes <- rsbi_events %>%
  mutate(event_id = paste(subject_id, event_time)) %>%
  group_by(event_id, subject_id, failed_extubation_flag,
           caregiver_fe_rate_loo, vent_hours) %>%
  filter(n() >= 3) %>%
  summarise(
    rsbi_slope    = coef(lm(rsbi ~ hours_before_extubation))[2],
    rsbi_at_extub = rsbi[which.min(hours_before_extubation)],
    rsbi_mean     = mean(rsbi),
    n_obs         = n(),
    .groups       = "drop"
  ) %>%
  mutate(failed = failed_extubation_flag == 1)

cat("\nRSBI slope analysis:\n")
cat("Events with >= 3 measurements:", nrow(event_slopes), "\n")
cat("Failed:", sum(event_slopes$failed), "\n")
cat("Succeeded:", sum(!event_slopes$failed), "\n")

cat("\nSlope distribution by outcome:\n")
event_slopes %>%
  group_by(failed) %>%
  summarise(
    mean_slope    = mean(rsbi_slope, na.rm = TRUE),
    median_slope  = median(rsbi_slope, na.rm = TRUE),
    mean_at_extub = mean(rsbi_at_extub, na.rm = TRUE),
    n             = n(),
    .groups       = "drop"
  ) %>%
  print()

cat("\nLogistic regression: failed ~ rsbi_at_extub + rsbi_slope\n")
glm(failed ~ rsbi_at_extub + rsbi_slope,
    data   = event_slopes %>% filter(!is.na(rsbi_slope), !is.na(rsbi_at_extub)),
    family = binomial) %>%
  summary()


## ── 4. P/F RATIO MEASUREMENTS IN 72H BEFORE EACH EXTUBATION ─────────────────
##
##  Joins P/F ratio (pao2fio2ratio) from the mimiciv derived blood gas table
##  within 72h before each extubation event. Joins on hadm_id since the bg
##  table does not have stay_id.

sql_pf_events <- "
  SELECT
    e.subject_id,
    e.stay_id,
    e.hadm_id,
    e.event_time,
    e.failed_extubation_flag,
    e.caregiver_fe_rate_loo,
    e.vent_hours,
    bg.pao2fio2ratio                                   AS pf_ratio,
    TIMESTAMP_DIFF(e.event_time, bg.charttime, MINUTE) / 60.0
                                                       AS hours_before_extubation
  FROM `mythical-legend-456217-e9.Failed_Extubation_Rate.all_extubations` e
  JOIN `physionet-data.mimiciv_3_1_derived.bg` bg
    ON e.hadm_id = bg.hadm_id
    AND bg.pao2fio2ratio BETWEEN 50 AND 600
    AND TIMESTAMP_DIFF(e.event_time, bg.charttime, MINUTE) BETWEEN 0 AND 72*60
"

pf_events <- bq_project_query(project, sql_pf_events) %>%
  bq_table_download(bigint = "integer64")

cat("\nP/F ratio data:\n")
cat("Rows:", nrow(pf_events), "\n")
cat("Events with P/F:", n_distinct(paste(pf_events$subject_id, pf_events$event_time)), "\n")
cat("Measurements per event (median):",
    median(table(paste(pf_events$subject_id, pf_events$event_time))), "\n")


## ── 5. P/F RATIO PER-EVENT SLOPE ANALYSIS ────────────────────────────────────

pf_event_slopes <- pf_events %>%
  mutate(event_id = paste(subject_id, event_time)) %>%
  group_by(event_id, subject_id, failed_extubation_flag,
           caregiver_fe_rate_loo, vent_hours) %>%
  filter(n() >= 3) %>%
  summarise(
    pf_slope    = coef(lm(pf_ratio ~ hours_before_extubation))[2],
    pf_at_extub = pf_ratio[which.min(hours_before_extubation)],
    pf_mean     = mean(pf_ratio),
    n_obs       = n(),
    .groups     = "drop"
  ) %>%
  mutate(failed = failed_extubation_flag == 1)

cat("\nP/F slope analysis:\n")
cat("Events with >= 3 measurements:", nrow(pf_event_slopes), "\n")
cat("Failed:", sum(pf_event_slopes$failed), "\n")
cat("Succeeded:", sum(!pf_event_slopes$failed), "\n")

cat("\nSlope distribution by outcome:\n")
pf_event_slopes %>%
  group_by(failed) %>%
  summarise(
    mean_slope    = mean(pf_slope, na.rm = TRUE),
    median_slope  = median(pf_slope, na.rm = TRUE),
    mean_at_extub = mean(pf_at_extub, na.rm = TRUE),
    n             = n(),
    .groups       = "drop"
  ) %>%
  print()

cat("\nLogistic regression: failed ~ pf_at_extub + pf_slope\n")
glm(failed ~ pf_at_extub + pf_slope,
    data   = pf_event_slopes %>% filter(!is.na(pf_slope), !is.na(pf_at_extub)),
    family = binomial) %>%
  summary()


## ── 6. TRAJECTORY VISUALIZATION ──────────────────────────────────────────────
##
##  Population mean ± 95% CI for RSBI and P/F ratio at hourly bins in the
##  72h before extubation, split by extubation outcome.
##
##  Note: each hourly bin draws from a different subset of events (patients
##  do not all have measurements at every hour). The lines show population
##  cross-sections, not individual trajectories. The regression analysis in
##  sections 3 and 5 uses proper per-event slopes.

rsbi_summary <- rsbi_events %>%
  mutate(
    failed   = failed_extubation_flag == 1,
    hour_bin = floor(hours_before_extubation)
  ) %>%
  filter(hour_bin <= 72) %>%
  group_by(hour_bin, failed) %>%
  summarise(
    mean_rsbi = mean(rsbi),
    se_rsbi   = sd(rsbi) / sqrt(n()),
    n         = n(),
    .groups   = "drop"
  ) %>%
  filter(n >= 30)

pf_summary <- pf_events %>%
  mutate(
    failed   = failed_extubation_flag == 1,
    hour_bin = floor(hours_before_extubation)
  ) %>%
  filter(hour_bin <= 72) %>%
  group_by(hour_bin, failed) %>%
  summarise(
    mean_pf = mean(pf_ratio),
    se_pf   = sd(pf_ratio) / sqrt(n()),
    n       = n(),
    .groups = "drop"
  ) %>%
  filter(n >= 30)

p_rsbi <- ggplot(rsbi_summary,
                 aes(hour_bin, mean_rsbi,
                     color = failed, fill = failed, group = failed)) +
  geom_ribbon(aes(ymin = mean_rsbi - 1.96 * se_rsbi,
                  ymax = mean_rsbi + 1.96 * se_rsbi),
              alpha = 0.2, color = NA) +
  geom_line(linewidth = 0.8) +
  scale_x_reverse(
    breaks = c(72, 48, 24, 0),
    labels = c("72h", "48h", "24h", "Extubation")
  ) +
  scale_color_manual(
    values = c("FALSE" = "#1D9E8E", "TRUE" = "#D95F02"),
    labels = c("FALSE" = "Succeeded", "TRUE" = "Failed")
  ) +
  scale_fill_manual(
    values = c("FALSE" = "#1D9E8E", "TRUE" = "#D95F02"),
    labels = c("FALSE" = "Succeeded", "TRUE" = "Failed")
  ) +
  labs(
    title = "RSBI",
    x     = NULL, y = "Mean RSBI",
    color = "Extubation", fill = "Extubation"
  ) +
  theme_minimal() +
  theme(legend.position = "none")

p_pf <- ggplot(pf_summary,
               aes(hour_bin, mean_pf,
                   color = failed, fill = failed, group = failed)) +
  geom_ribbon(aes(ymin = mean_pf - 1.96 * se_pf,
                  ymax = mean_pf + 1.96 * se_pf),
              alpha = 0.2, color = NA) +
  geom_line(linewidth = 0.8) +
  scale_x_reverse(
    breaks = c(72, 48, 24, 0),
    labels = c("72h", "48h", "24h", "Extubation")
  ) +
  scale_color_manual(
    values = c("FALSE" = "#1D9E8E", "TRUE" = "#D95F02"),
    labels = c("FALSE" = "Succeeded", "TRUE" = "Failed")
  ) +
  scale_fill_manual(
    values = c("FALSE" = "#1D9E8E", "TRUE" = "#D95F02"),
    labels = c("FALSE" = "Succeeded", "TRUE" = "Failed")
  ) +
  labs(
    title = "P/F Ratio",
    x     = NULL, y = "Mean P/F Ratio",
    color = "Extubation", fill = "Extubation"
  ) +
  theme_minimal() +
  theme(legend.position = "right")

p_rsbi + p_pf +
  plot_layout(guides = "collect") +
  plot_annotation(
    title    = "Physiological Trajectories Before Extubation by Outcome",
    subtitle = paste0(
      "Per-extubation-event data (RSBI: ",
      n_distinct(rsbi_events$event_time),
      " events; P/F: ",
      n_distinct(paste(pf_events$subject_id, pf_events$event_time)),
      " events).\n",
      "Hourly bins, n \u2265 30 per bin. Mean \u00b1 95% CI. ",
      "Failed = reintubated within 72h of this specific extubation."
    ),
    theme = theme(
      plot.title    = element_text(size = 13, face = "bold"),
      plot.subtitle = element_text(size = 8.5, color = "gray40")
    )
  ) &
  theme(legend.position = "right")

ggsave("../figures/preextubation_trajectories_by_event_outcome.png",
       width = 11, height = 5.5, dpi = 150)
cat("Trajectory plot saved.\n")
