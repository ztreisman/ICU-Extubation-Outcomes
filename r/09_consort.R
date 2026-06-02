################################################################################
##
##  ICU Extubation Outcomes — Script 09: CONSORT Diagram
##
##  Depends on: 01_cohort.R (cohorts.RData)
##
##  All counts from the SQL level down are computed from loaded R objects.
##  Counts that predate last_extubations.csv (SQL pipeline internals) are
##  defined as named parameters at the top — update them if the query changes
##  or when reusing for a different dataset (e.g., eICU).
##
##  Output: ../figures/consort_diagram.png
##
################################################################################


## ── 0. SETUP ─────────────────────────────────────────────────────────────────

library(dplyr)
library(ggplot2)

load("../data/cohorts.RData")


## ── 1. UPSTREAM PARAMETERS (SQL-level, not in R objects) ─────────────────────
##
##  These come from the SQL pipeline and are displayed in the script 01 output.
##  Update here if the source query or database version changes.

N_VENT_EVENTS   <- 27729L  # InvasiveVent >= 4h events in MIMIC-IV v3.1
                            # (all-extubation total before dedup)
N_RAW_SQL       <- 12668L  # last_extubations.csv row count (1 per patient)


## ── 2. COMPUTE COUNTS FROM R OBJECTS ─────────────────────────────────────────

n_cohort        <- nrow(patient_cohort)
n_explicit      <- sum(patient_cohort$tube_event_source == "explicit")
n_inferred      <- sum(patient_cohort$tube_event_source == "inferred")
n_analytic      <- nrow(explicit_extubations)
n_base          <- nrow(analysis_base)

n_cvicu         <- sum(analysis_base$unit_group == "CVICU")
n_medical       <- sum(analysis_base$unit_group == "Medical")
n_surg          <- sum(analysis_base$unit_group == "Surgical/Trauma")
n_other         <- sum(analysis_base$unit_group == "Other")

## Derived exclusion counts
n_excl_covar    <- N_RAW_SQL - n_cohort          # incomplete covariates
n_excl_inferred <- n_inferred                    # no caregiver link
n_excl_vol      <- n_explicit - n_analytic       # caregiver_n < 20 or non-std ICU
n_excl_unit     <- n_analytic - n_base           # caregiver_unit_n < 20

## Thresholds (for labels — keep in sync with 01_cohort.R and 05_msm.R)
THRESH_CAREGIVER_N    <- 20L
THRESH_CAREGIVER_UNIT <- 20L


## ── 3. LAYOUT PARAMETERS ─────────────────────────────────────────────────────
##
##  All positions in arbitrary units. Adjust to taste.
##  x_main = horizontal center of main flow column
##  x_excl = horizontal center of exclusion column

x_main  <- 3.5
x_excl  <- 7.8
canvas_x <- c(0, 11)

# Y centers for main-flow boxes (decreasing)
y <- list(
  sql    = 13.0,
  cohort = 11.0,
  expl   = 8.8,
  analyt = 6.6,
  base   = 4.4,
  groups = 2.2
)

# Box dimensions
w_main  <- 3.4   # main flow box half-width → full width = 2 * w_main
w_excl  <- 2.8
w_group <- 2.0
h <- list(
  sql    = 0.55,
  cohort = 0.75,
  expl   = 0.55,
  analyt = 0.75,
  base   = 0.75,
  group  = 0.65
)


## ── 4. LABEL STRINGS ─────────────────────────────────────────────────────────

lbl <- list(
  sql    = sprintf("MIMIC-IV v3.1\nLast extubation per patient\n(N = %s)",
                   format(N_RAW_SQL, big.mark = ",")),

  cohort = sprintf("patient_cohort\n(N = %s)\nExplicit: %s  |  Inferred: %s",
                   format(n_cohort,   big.mark = ","),
                   format(n_explicit, big.mark = ","),
                   format(n_inferred, big.mark = ",")),

  expl   = sprintf("Explicit extubations\n(N = %s)",
                   format(n_explicit, big.mark = ",")),

  analyt = sprintf("explicit_extubations\ncaregiver_n ≥ %d, standard ICU\n(N = %s)",
                   THRESH_CAREGIVER_N,
                   format(n_analytic, big.mark = ",")),

  base   = sprintf("analysis_base\ncaregiver_unit_n ≥ %d\n(N = %s)",
                   THRESH_CAREGIVER_UNIT,
                   format(n_base, big.mark = ",")),

  ## Exclusion labels
  e_covar  = sprintf("Excluded  N = %s\nIncomplete covariates",
                     format(n_excl_covar, big.mark = ",")),

  e_inf    = sprintf("Excluded  N = %s\nInferred events\n(no caregiver link)",
                     format(n_excl_inferred, big.mark = ",")),

  e_vol    = sprintf("Excluded  N = %s\ncaregiver_n < %d or\nnon-standard ICU unit",
                     format(n_excl_vol, big.mark = ","),
                     THRESH_CAREGIVER_N),

  e_unit   = sprintf("Excluded  N = %s\ncaregiver_unit_n < %d",
                     format(n_excl_unit, big.mark = ","),
                     THRESH_CAREGIVER_UNIT),

  ## Unit group labels
  cvicu    = sprintf("CVICU\nN = %s", format(n_cvicu,   big.mark = ",")),
  medical  = sprintf("Medical\nMICU + MICU/SICU\nN = %s",
                     format(n_medical, big.mark = ",")),
  surg     = sprintf("Surgical/Trauma\nSICU + TSICU\nN = %s",
                     format(n_surg,   big.mark = ",")),
  other    = sprintf("Other ICUs\nN = %s\n(descriptive only)",
                     format(n_other,  big.mark = ","))
)


## ── 5. BUILD DIAGRAM DATA ────────────────────────────────────────────────────

## Boxes: (xc, yc, half-width, half-height, label, type)
make_box <- function(xc, yc, hw, hh, label, type = "main") {
  tibble(xc = xc, yc = yc, hw = hw, hh = hh, label = label, type = type)
}

# Y of exclusion boxes = midpoint between adjacent main flow boxes
ye <- list(
  covar = (y$sql    + y$cohort) / 2,
  inf   = (y$cohort + y$expl)   / 2,
  vol   = (y$expl   + y$analyt) / 2,
  unit  = (y$analyt + y$base)   / 2
)

## Half-heights for exclusion boxes (3 lines → taller)
hh_excl <- 0.55

boxes <- bind_rows(
  ## Main flow
  make_box(x_main, y$sql,    w_main, h$sql,    lbl$sql,    "main"),
  make_box(x_main, y$cohort, w_main, h$cohort, lbl$cohort, "main"),
  make_box(x_main, y$expl,   w_main, h$expl,   lbl$expl,   "main"),
  make_box(x_main, y$analyt, w_main, h$analyt, lbl$analyt, "main"),
  make_box(x_main, y$base,   w_main, h$base,   lbl$base,   "main"),
  ## Exclusion boxes
  make_box(x_excl, ye$covar, w_excl, hh_excl,  lbl$e_covar, "excl"),
  make_box(x_excl, ye$inf,   w_excl, hh_excl,  lbl$e_inf,   "excl"),
  make_box(x_excl, ye$vol,   w_excl, hh_excl,  lbl$e_vol,   "excl"),
  make_box(x_excl, ye$unit,  w_excl, hh_excl,  lbl$e_unit,  "excl"),
  ## Unit groups (evenly spread)
  make_box(1.2,  y$groups, w_group, h$group, lbl$cvicu,   "group"),
  make_box(3.5,  y$groups, w_group, h$group, lbl$medical, "group"),
  make_box(5.8,  y$groups, w_group, h$group, lbl$surg,    "group"),
  make_box(8.1,  y$groups, w_group, h$group, lbl$other,   "other")
)

## Arrows: (x0, y0, x1, y1, type)
make_arr <- function(x0, y0, x1, y1, type = "main") {
  tibble(x0 = x0, y0 = y0, x1 = x1, y1 = y1, type = type)
}

arrows_df <- bind_rows(
  ## Vertical main flow
  make_arr(x_main, y$sql    - h$sql,    x_main, y$cohort + h$cohort),
  make_arr(x_main, y$cohort - h$cohort, x_main, y$expl   + h$expl),
  make_arr(x_main, y$expl   - h$expl,   x_main, y$analyt + h$analyt),
  make_arr(x_main, y$analyt - h$analyt, x_main, y$base   + h$base),
  ## Base → unit groups (horizontal spread then down)
  make_arr(1.2,  y$base - h$base, 1.2,  y$groups + h$group, "group"),
  make_arr(3.5,  y$base - h$base, 3.5,  y$groups + h$group, "group"),
  make_arr(5.8,  y$base - h$base, 5.8,  y$groups + h$group, "group"),
  make_arr(8.1,  y$base - h$base, 8.1,  y$groups + h$group, "other"),
  ## Horizontal connectors to exclusion boxes (right-angle elbows)
  ## These go: main flow right edge → horizontal to excl box left edge
  make_arr(x_main + w_main, ye$covar, x_excl - w_excl, ye$covar, "excl"),
  make_arr(x_main + w_main, ye$inf,   x_excl - w_excl, ye$inf,   "excl"),
  make_arr(x_main + w_main, ye$vol,   x_excl - w_excl, ye$vol,   "excl"),
  make_arr(x_main + w_main, ye$unit,  x_excl - w_excl, ye$unit,  "excl")
)

## Horizontal cap segments (short tick from main flow box to elbow y)
caps_df <- bind_rows(
  make_arr(x_main + w_main, y$cohort, x_main + w_main, ye$covar, "cap"),
  make_arr(x_main + w_main, y$cohort, x_main + w_main, ye$inf,   "cap"),
  make_arr(x_main + w_main, y$expl,   x_main + w_main, ye$vol,   "cap"),
  make_arr(x_main + w_main, y$analyt, x_main + w_main, ye$unit,  "cap")
)

## Horizontal spreader at bottom of analysis_base → unit groups
spreader_y <- y$base - h$base
spreader_df <- tibble(
  x0 = 1.2, y0 = spreader_y, x1 = 8.1, y1 = spreader_y
)


## ── 6. BUILD THE PLOT ────────────────────────────────────────────────────────

p <- ggplot() +

  ## Horizontal spreader
  geom_segment(data = spreader_df,
               aes(x = x0, y = y0, xend = x1, yend = y1),
               colour = "grey40", linewidth = 0.4) +

  ## Cap segments (no arrow)
  geom_segment(data = caps_df,
               aes(x = x0, y = y0, xend = x1, yend = y1),
               colour = "grey40", linewidth = 0.4) +

  ## Exclusion arrows
  geom_segment(data = filter(arrows_df, type == "excl"),
               aes(x = x0, y = y0, xend = x1, yend = y1),
               colour = "grey40", linewidth = 0.4,
               arrow = arrow(length = unit(0.15, "cm"),
                             type = "closed")) +

  ## Main-flow arrows
  geom_segment(data = filter(arrows_df, type == "main"),
               aes(x = x0, y = y0, xend = x1, yend = y1),
               colour = "grey30", linewidth = 0.5,
               arrow = arrow(length = unit(0.18, "cm"),
                             type = "closed")) +

  ## Group arrows
  geom_segment(data = filter(arrows_df, type %in% c("group", "other")),
               aes(x = x0, y = y0, xend = x1, yend = y1),
               colour = "grey30", linewidth = 0.5,
               arrow = arrow(length = unit(0.15, "cm"),
                             type = "closed")) +

  ## Box outlines
  geom_rect(data = boxes,
            aes(xmin = xc - hw, xmax = xc + hw,
                ymin = yc - hh, ymax = yc + hh,
                fill = type, colour = type),
            linewidth = 0.4) +

  ## Box text
  geom_text(data = boxes,
            aes(x = xc, y = yc, label = label),
            size = 2.7, lineheight = 1.2, vjust = 0.5) +

  ## Styling
  scale_fill_manual(values = c(
    main  = "white",
    excl  = "#FFF3E0",   # pale amber for exclusions
    group = "#E8F4FD",   # pale blue for analysis groups
    other = "#F5F5F5"    # light grey for "other" (not primary analysis)
  )) +
  scale_colour_manual(values = c(
    main  = "grey30",
    excl  = "#E6A020",
    group = "#2E75B6",
    other = "grey60"
  )) +
  coord_cartesian(xlim = canvas_x,
                  ylim = c(y$groups - h$group - 0.3, y$sql + h$sql + 0.5)) +
  theme_void() +
  theme(legend.position = "none",
        plot.margin = margin(8, 8, 8, 8))


## ── 7. SAVE ───────────────────────────────────────────────────────────────────

ggsave("../figures/consort_diagram.png",
       plot = p, width = 9, height = 12, dpi = 200)
cat("Saved: ../figures/consort_diagram.png\n")
cat(sprintf("Cohort flow: %s → %s explicit → %s analytic → %s analysis_base\n",
            format(n_cohort, big.mark = ","),
            format(n_explicit, big.mark = ","),
            format(n_analytic, big.mark = ","),
            format(n_base, big.mark = ",")))
