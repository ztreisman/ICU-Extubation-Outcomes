################################################################################
##
##  ICU Extubation Outcomes — Script 03: NLS Exponential Model and
##  Caregiver-Level Partial Correlations
##
##  Depends on: 01_cohort.R (cohorts.RData)
##
##  Research questions:
##    1. Is caregiver FE rate independently associated with patient mortality
##       at the caregiver level, after controlling for caregiver volume?
##       (partial correlation)
##    2. What is the functional form of the FE rate — mortality relationship
##       at the caregiver level, and is there an identifiable mortality
##       ceiling? (asymptotic exponential NLS)
##
##  NLS model: P(mortality) = m_max * (1 - exp(-k * caregiver_fe_rate))
##    m_max: mortality ceiling — maximum achievable mortality as FE rate → ∞;
##           the fraction of outcomes potentially attributable to caregiver
##           performance. Equivalent to (1 - survival floor b).
##    k:     rate of increase — steepness of mortality rise with FE rate.
##
##  Both analyses operate at the caregiver level (one row per caregiver,
##  using caregiver_summary from cohorts.RData). Restricted to caregivers
##  with FE rate > 0 for the NLS fit.
##
##  Saves: ../data/nls_results.RData
##
################################################################################


## ── 0. LIBRARIES ─────────────────────────────────────────────────────────────

library(dplyr)
library(ggplot2)
library(ppcor)
library(stringr)


## ── 1. LOAD DATA ─────────────────────────────────────────────────────────────

load("../data/cohorts.RData")
# Loads: caregiver_summary, caregivers_nonzero,
#        lm_linear, lm_quadratic, lm_log, unit_summary


## ── 2. CAREGIVER-LEVEL PARTIAL CORRELATIONS ──────────────────────────────────
##
##  Partial correlation between caregiver FE rate and patient mortality rate,
##  controlling for caregiver volume (n_patients_sample). Then the reverse:
##  volume ~ mortality controlling for FE rate. If FE rate carries independent
##  information beyond volume, the first partial r should be significant while
##  the second should not.

cat("\n── Caregiver-level partial correlations ",
    paste0(rep("─", 33), collapse = ""), "\n")
cat("N caregivers (caregiver_unit_n >= 20):", nrow(caregiver_summary), "\n")

pcor_fe <- pcor.test(
  caregiver_summary$caregiver_fe_rate,
  caregiver_summary$mortality_rate,
  caregiver_summary$n_patients_sample
)
cat("\nPartial r: FE rate ~ mortality | caregiver volume\n")
cat("r =", round(pcor_fe$estimate, 3), "  p =", round(pcor_fe$p.value, 4), "\n")

pcor_vol <- pcor.test(
  caregiver_summary$caregiver_n,
  caregiver_summary$mortality_rate,
  caregiver_summary$caregiver_fe_rate
)
cat("\nPartial r: volume ~ mortality | FE rate\n")
cat("r =", round(pcor_vol$estimate, 3), "  p =", round(pcor_vol$p.value, 4), "\n")

cat("\nVolume threshold sensitivity:\n")
for (min_n in c(10, 20)) {
  cs_sub <- caregiver_summary |> filter(n_patients_sample >= min_n)
  if (nrow(cs_sub) < 10) {
    cat(sprintf("n >= %d: fewer than 10 caregivers — skipped\n", min_n))
    next
  }
  pc_sub <- pcor.test(cs_sub$caregiver_fe_rate, cs_sub$mortality_rate,
                      cs_sub$n_patients_sample)
  cat(sprintf("n >= %d (%d caregivers): r = %.3f, p = %.4f\n",
              min_n, nrow(cs_sub), pc_sub$estimate, pc_sub$p.value))
}

unit_pcor_results <- caregiver_summary |>
  filter(!is.na(primary_careunit)) |>
  group_by(primary_careunit) |>
  group_modify(~ {
    df <- .x
    base <- tibble(
      n_caregivers   = nrow(df),
      mean_fe_rate   = round(mean(df$caregiver_fe_rate, na.rm = TRUE), 4),
      mean_mortality = round(mean(df$mortality_rate,    na.rm = TRUE), 3)
    )
    if (nrow(df) < 5) return(mutate(base, r = NA_real_, p = NA_real_))
    pc <- tryCatch(
      pcor.test(df$caregiver_fe_rate, df$mortality_rate,
                df$n_patients_sample),
      error = function(e) list(estimate = NA_real_, p.value = NA_real_)
    )
    mutate(base, r = round(pc$estimate, 3), p = round(pc$p.value, 4))
  }) |>
  ungroup() |>
  arrange(desc(n_caregivers))

cat("\nUnit-stratified partial correlations (FE rate ~ mortality | volume):\n")
print(unit_pcor_results, n = Inf)
cat(paste0(rep("─", 72), collapse = ""), "\n")


## ── 3. ASYMPTOTIC EXPONENTIAL NLS MODEL ──────────────────────────────────────
##
##  P(mortality) = m_max * (1 - exp(-k * caregiver_fe_rate))
##
##  Restricted to caregivers with FE rate > 0 (zero-FE caregivers cannot
##  inform the shape of the dose-response curve). Model is weighted by
##  n_patients_sample so high-volume caregivers exert proportionally more
##  influence on the fit.
##
##  AIC compared against linear, quadratic, and log-linear forms fit in
##  01_cohort.R (stored in cohorts.RData as lm_linear, lm_quadratic, lm_log).

cat("\n── Asymptotic exponential model (mortality ceiling) ",
    paste0(rep("─", 22), collapse = ""), "\n")
cat("N caregivers (FE rate > 0):", nrow(caregivers_nonzero), "\n")

nls_fit <- nls(
  mortality_rate ~ m_max * (1 - exp(-k * caregiver_fe_rate)),
  data    = caregivers_nonzero,
  weights = n_patients_sample,
  start   = list(m_max = 0.40, k = 50),
  control = nls.control(maxiter = 500)
)

summary(nls_fit)

m_max_hat <- coef(nls_fit)["m_max"]
k_hat     <- coef(nls_fit)["k"]
m_max_ci  <- confint(nls_fit)["m_max", ]
k_ci      <- confint(nls_fit)["k", ]

cat(sprintf("\nm_max (mortality ceiling): %.3f  95%% CI (%.3f–%.3f)\n",
            m_max_hat, m_max_ci[1], m_max_ci[2]))
cat(sprintf("k     (rate of increase):  %.1f   95%% CI (%.1f–%.1f)\n",
            k_hat, k_ci[1], k_ci[2]))

cat("\nModel AIC comparison:\n")
print(AIC(lm_linear, lm_quadratic))
cat("NLS (asymptotic exponential) AIC:", round(AIC(nls_fit), 3), "\n")

## Unit-stratified NLS
nls_by_unit <- caregiver_summary |>
  filter(caregiver_fe_rate > 0) |>
  group_by(primary_careunit) |>
  summarise(n = n(), fe_range = max(caregiver_fe_rate) - min(caregiver_fe_rate),
            .groups = "drop") |>
  filter(n >= 10, fe_range > 0.02) |>
  pull(primary_careunit)

nls_results <- caregiver_summary |>
  filter(caregiver_fe_rate > 0, primary_careunit %in% nls_by_unit) |>
  group_by(primary_careunit) |>
  group_map(function(data, key) {
    unit <- as.character(key$primary_careunit)
    tryCatch({
      fit <- nls(
        mortality_rate ~ m_max * (1 - exp(-k * caregiver_fe_rate)),
        data = data, weights = n_patients_sample,
        start = list(m_max = 0.40, k = 50),
        control = nls.control(maxiter = 500)
      )
      ci <- tryCatch(confint(fit),
                     error = function(e) matrix(
                       NA, 2, 2, dimnames = list(c("m_max","k"),
                                                 c("2.5%","97.5%"))))
      tibble(unit = unit, n = nrow(data),
             m_max       = round(coef(fit)["m_max"],     3),
             m_max_ci_lo = round(ci["m_max","2.5%"],     3),
             m_max_ci_hi = round(ci["m_max","97.5%"],    3),
             k           = round(coef(fit)["k"],         1),
             k_ci_lo     = round(ci["k","2.5%"],         1),
             k_ci_hi     = round(ci["k","97.5%"],        1),
             converged   = TRUE)
    }, error = function(e) tibble(
      unit = unit, n = nrow(data),
      m_max = NA, m_max_ci_lo = NA, m_max_ci_hi = NA,
      k = NA, k_ci_lo = NA, k_ci_hi = NA,
      converged = FALSE
    ))
  }) |>
  bind_rows()

cat("\nUnit-stratified NLS results:\n")
print(nls_results)
cat(paste0(rep("─", 72), collapse = ""), "\n")


## ── 4. FIGURES ────────────────────────────────────────────────────────────────

## 4a. Global NLS fit
fe_seq   <- seq(0, max(caregivers_nonzero$caregiver_fe_rate), by = 0.0005)
pred_exp <- data.frame(
  caregiver_fe_rate = fe_seq,
  mortality_rate    = m_max_hat * (1 - exp(-k_hat * fe_seq))
)

ggplot() +
  geom_point(data = caregivers_nonzero,
             aes(caregiver_fe_rate, mortality_rate,
                 size = n_patients_sample, color = primary_careunit),
             alpha = 0.8) +
  geom_line(data = pred_exp,
            aes(caregiver_fe_rate, mortality_rate),
            color = "#2E75B6", linewidth = 1.0) +
  geom_hline(yintercept = m_max_hat, linetype = "dashed", color = "gray50") +
  annotate("text",
           x     = max(caregivers_nonzero$caregiver_fe_rate) * 0.7,
           y     = m_max_hat + 0.03,
           label = paste0(
             "Mortality ceiling = ", round(m_max_hat, 3),
             " (95% CI ", round(m_max_ci[1], 3), "–",
             round(m_max_ci[2], 3), ")"
           ),
           color = "gray40", size = 3.5) +
  scale_size_continuous(range = c(1, 8), name = "N patients") +
  labs(
    title    = "Caregiver FE Rate vs. 12-Month Patient Mortality",
    subtitle = paste0(
      "P(mortality) = m_max×(1−exp(−k×FE rate))     ",
      "m_max = ", round(m_max_hat, 3),
      " (95% CI ", round(m_max_ci[1], 3), "–", round(m_max_ci[2], 3), ")     ",
      "k = ", round(k_hat, 1),
      " (95% CI ", round(k_ci[1], 1), "–", round(k_ci[2], 1), ")"
    ),
    x = "Caregiver Failed Extubation Rate",
    y = "12-Month Patient Mortality Rate"
  ) +
  theme_minimal() +
  theme(plot.subtitle = element_text(size = 9))

ggsave("../figures/fe_rate_mortality_exponential.png",
       width = 9, height = 6, dpi = 150)
cat("Global NLS figure saved.\n")


## 4b. Stratified NLS figure
fe_seq      <- seq(0, 0.25, length.out = 200)
m_max_cvicu <- nls_results$m_max[
  nls_results$unit == "Cardiac Vascular Intensive Care Unit (CVICU)"
]
k_cvicu     <- nls_results$k[
  nls_results$unit == "Cardiac Vascular Intensive Care Unit (CVICU)"
]

caregiver_plot <- caregiver_summary |>
  filter(caregiver_fe_rate > 0) |>
  mutate(unit_label = case_when(
    str_detect(primary_careunit, "CVICU")           ~ "CVICU",
    str_detect(primary_careunit, "Medical Intens")  ~ "MICU",
    str_detect(primary_careunit, "Trauma")          ~ "TSICU",
    TRUE                                            ~ "Other"
  ))

ggplot() +
  geom_point(data = caregiver_plot,
             aes(caregiver_fe_rate, mortality_rate,
                 size = n_patients_sample, color = unit_label),
             alpha = 0.6) +
  geom_line(data = bind_rows(
    tibble(caregiver_fe_rate = fe_seq,
           mortality = coef(nls_fit)["m_max"] *
             (1 - exp(-coef(nls_fit)["k"] * fe_seq)),
           model = sprintf("Cross-unit (m_max=%.3f)",
                           coef(nls_fit)["m_max"])),
    tibble(caregiver_fe_rate = fe_seq,
           mortality = m_max_cvicu * (1 - exp(-k_cvicu * fe_seq)),
           model = sprintf("CVICU only (m_max=%.3f)", m_max_cvicu))
  ), aes(caregiver_fe_rate, mortality, linetype = model),
  linewidth = 1.0, color = "gray20") +
  geom_hline(yintercept = m_max_cvicu, linetype = "dotted",
             color = "#2E75B6", linewidth = 0.7) +
  geom_hline(yintercept = coef(nls_fit)["m_max"], linetype = "dotted",
             color = "gray40", linewidth = 0.7) +
  annotate("text", x = 0.22, y = m_max_cvicu - 0.03,
           label = sprintf("CVICU ceiling = %.3f", m_max_cvicu),
           size = 3, color = "#2E75B6", hjust = 1) +
  annotate("text", x = 0.22, y = coef(nls_fit)["m_max"] + 0.03,
           label = sprintf("Global ceiling = %.3f", coef(nls_fit)["m_max"]),
           size = 3, color = "gray40", hjust = 1) +
  scale_color_manual(values = c("CVICU" = "#2E75B6", "MICU" = "#D95F02",
                                "TSICU" = "#1D9E8E", "Other" = "gray60")) +
  scale_size_continuous(range = c(1, 5), guide = "none") +
  scale_x_continuous(labels = scales::percent_format(accuracy = 1)) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1),
                     limits = c(0, 1)) +
  labs(
    title    = "Caregiver FE Rate vs. 12-Month Patient Mortality",
    subtitle = paste0(
      "Asymptotic exponential fit: global and CVICU-only. ",
      "Bubble size = N patients."
    ),
    x     = "Caregiver Failed Extubation Rate",
    y     = "12-Month Patient Mortality Rate",
    color = "Unit", linetype = "Model"
  ) +
  theme_minimal() +
  theme(legend.position = "right",
        plot.title      = element_text(size = 11, face = "bold"),
        plot.subtitle   = element_text(size = 8.5, color = "gray40"))

ggsave("../figures/nls_stratified.png", width = 8, height = 5, dpi = 150)
cat("Stratified NLS figure saved.\n")


## ── 5. SAVE ───────────────────────────────────────────────────────────────────

save(
  pcor_fe, pcor_vol, unit_pcor_results,
  nls_fit, m_max_hat, k_hat, m_max_ci, k_ci,
  nls_results,
  file = "../data/nls_results.RData"
)
cat("\nSaved: ../data/nls_results.RData\n")

