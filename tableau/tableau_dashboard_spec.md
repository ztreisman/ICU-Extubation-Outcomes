# Tableau Dashboard Design Specification
## ICU Extubation Outcomes: Does Provider Volume Predict Patient Survival?

**Platform:** Tableau Desktop (publish to Tableau Public)
**Data sources:** `synthetic_last_extubations.csv`, `synthetic_caregiver_rate.csv`
**Dashboard count:** 3
**Estimated build time:** 6–10 hours total

---

## Data Source Setup

### Connection
Connect both CSVs as separate data sources. Name them:
- `Patient Records` → `synthetic_last_extubations.csv`
- `Caregiver Summary` → `synthetic_caregiver_rate.csv`

### Field type overrides (Patient Records)
After connecting, set these manually in the Data Source tab:

| Field | Set type to |
|---|---|
| `subject_id` | String |
| `caregiver_id` | String |
| `survival_12mo` | Boolean |
| `hospital_expire_flag` | Number (treat as Dimension) |
| `failed_extubations` | Number (whole) |
| `intubation_type` | Dimension |
| `primary_diagnosis` | Dimension |
| `tube_event_source` | Dimension |
| `gender` | Dimension |

### Calculated fields to create upfront (Patient Records)

Create these before building any sheets via **Analysis > Create Calculated Field**.

**`Failed Extubation Flag`**
```
[failed_extubations] > 0
```

**`Survival Label`**
```
IF [survival_12mo] = TRUE THEN "Survived 12 mo" ELSE "Did Not Survive" END
```

**`Intubation Type Label`**
```
CASE [intubation_type]
  WHEN "surgical" THEN "Surgical"
  WHEN "medical-respiratory" THEN "Medical – Respiratory"
  WHEN "medical-non-respiratory" THEN "Medical – Non-Respiratory"
  ELSE [intubation_type]
END
```

**`Caregiver Volume Tier`**
```
IF [caregiver_n] < 20 THEN "1 – Low (<20)"
ELSEIF [caregiver_n] < 50 THEN "2 – Medium (20–49)"
ELSEIF [caregiver_n] < 100 THEN "3 – High (50–99)"
ELSE "4 – Very High (100+)"
END
```
*(Prefix numbers force correct sort order)*

---

## DASHBOARD 1: Patient Cohort Overview

**Purpose:** Establish clinical context — who are these patients and what are baseline outcome rates?
**Story beat:** Data storytelling starts with cohort description, not modeling.
**Size:** 1200 × 850 px

### Layout (Tiled)

```
┌──────────────────────────────────────────────────────────────┐
│  HEADER: Title + subtitle + 4 KPI tiles                      │  ~120px
├────────────────────┬─────────────────────────────────────────┤
│  Sheet 1a          │  Sheet 1b                               │
│  Age Distribution  │  Intubation Type × Outcome              │  ~320px
│  (histogram)       │  (stacked bar)                          │
├────────────────────┴──────────────┬──────────────────────────┤
│  Sheet 1c                         │  Sheet 1d                │
│  ICD Chapter Breakdown            │  Acuity Scatter          │
│  (horizontal bar)                 │  (SOFA vs Charlson)      │  ~320px
└───────────────────────────────────┴──────────────────────────┘
```

### KPI Tiles

Build four single-number sheets (Mark type: Text, no axes, no headers):

| Tile | Formula | Format |
|---|---|---|
| Total Patients | `COUNTD([subject_id])` | "800 Patients" |
| Failed Extubation Rate | `AVG(INT([Failed Extubation Flag])) * 100` | "14.5% Failed Extubation" |
| 12-Mo Survival | `AVG(INT([survival_12mo])) * 100` | "73.9% Survived 12 Mo" |
| Unique Caregivers | `COUNTD([caregiver_id])` | "80 Caregivers" |

Format KPI numbers at 28pt bold. Place in a horizontal container in the header area.

---

### Sheet 1a: Age Distribution

- **Mark type:** Bar
- **Columns:** `FLOOR([anchor_age] / 5) * 5` (5-year bins)
- **Rows:** `COUNT([subject_id])`
- **Color:** `[Survival Label]` — teal (#1D9E8E) = survived, coral (#D95F02) = did not survive
- **Sort:** Columns ascending
- **Title:** "Age Distribution by 12-Month Survival"
- **X-axis label:** "Age (5-year bins)" | **Y-axis label:** "Patient Count"
- Remove all gridlines (Format > Lines > None for rows and columns)

---

### Sheet 1b: Intubation Type × Outcome

- **Mark type:** Bar (stacked)
- **Columns:** `[Intubation Type Label]`
- **Rows:** `COUNTD([subject_id])`
- **Color:** `[Survival Label]` (same palette as 1a)
- **Sort columns:** Descending by count
- **Title:** "Outcomes by Intubation Type"
- Add a reference line on Y axis: `AVG(INT([Failed Extubation Flag]))` across the table, dashed gray, label "Overall FE Rate"
- **Tooltip:** Intubation type, N, FE rate %, 12-mo survival %

---

### Sheet 1c: ICD Chapter Breakdown

- **Mark type:** Bar (horizontal)
- **Rows:** `[primary_diagnosis]`
- **Columns:** `COUNTD([subject_id])`
- **Color:** `AVG(INT([Failed Extubation Flag])) * 100` — continuous orange-blue diverging (low FE = blue, high = orange)
- **Sort:** Rows descending by count
- **Title:** "Patient Volume & FE Rate by Primary ICD Chapter"
- **Tooltip:** Chapter, N patients, FE rate %, 12-mo survival %
- Color legend labeled "FE Rate (%)"

---

### Sheet 1d: Acuity Scatter

- **Mark type:** Circle
- **Columns:** `[charlson]`
- **Rows:** `[sofa]`
- **Color:** `[Survival Label]`
- **Opacity:** 50% | **Size:** Small-fixed
- **Title:** "Illness Severity: SOFA vs. Charlson Comorbidity Index"
- Add reference lines at `MEDIAN([charlson])` (vertical) and `MEDIAN([sofa])` (horizontal) to create quadrants. Style: light gray, dashed.
- **Tooltip:** Age, SOFA, Charlson, intubation type, outcome

---

### Dashboard 1 Filter (visible to viewer)

- Field: `[Intubation Type Label]`
- Control type: Checkbox list, "All" selected by default
- Apply to: All worksheets using Patient Records

---

## DASHBOARD 2: Provider Volume & Outcomes

**Purpose:** The central research finding. Show the volume-outcome relationship, its scatter, and its consistency across patient tiers.
**Data source:** Caregiver Summary (primary)
**Size:** 1200 × 900 px

### Layout

```
┌──────────────────────────────────────────────────────────────┐
│  HEADER: Title + 2 KPI tiles (Low-vol vs High-vol FE rates)  │  ~100px
├─────────────────────────────┬────────────────────────────────┤
│  Sheet 2a (MAIN)            │  Sheet 2b                      │
│  Volume vs. FE Rate scatter │  Volume vs. 12-Mo Survival     │
│  with polynomial trend      │  scatter with linear trend     │  ~380px
├─────────────────────────────┴────────────────────────────────┤
│  Sheet 2c                                                     │
│  Volume tier box plots — FE rate spread per tier             │  ~300px
└──────────────────────────────────────────────────────────────┘
```

### KPI Tiles (header — these are the headline finding)

Build in Caregiver Summary source. Place side-by-side with an arrow or "vs." separator:

| Tile | Formula | Note |
|---|---|---|
| Low-Volume FE Rate | `AVG([observed_fe_rate])` filtered to `[caregiver_n] < 20` | Use a fixed filter on the sheet |
| High-Volume FE Rate | `AVG([observed_fe_rate])` filtered to `[caregiver_n] >= 50` | Use a fixed filter on the sheet |

Format at 32pt bold. The gap between these two numbers is the story.

---

### Sheet 2a: Volume vs. FE Rate

- **Data source:** Caregiver Summary
- **Mark type:** Circle
- **Columns:** `[caregiver_n]`
- **Rows:** `[observed_fe_rate]`
- **Size:** `[n_patients_in_sample]` (bubble size = sample representation)
- **Color:** `[observed_fe_rate]` — continuous diverging (low = blue, high = orange)
- **Opacity:** 75%
- **Title:** "Caregiver Volume vs. Observed Failed Extubation Rate"
- **X-axis:** "Caregiver Extubation Volume" | **Y-axis:** "Observed FE Rate"

**Trend line:**
- Analysis > Trend Lines > Show Trend Lines
- Model: Polynomial, degree 2
- Show Confidence Bands: yes
- Uncheck "Show Tooltips for Trend Lines"

**Reference line:** `AVG([observed_fe_rate])`, horizontal dashed, labeled "Overall Mean"

**Manual annotation** (Worksheet > Annotate > Area, upper-left region of chart):
> "Low-volume caregivers show higher and more variable FE rates. Variation narrows as volume increases."

**Tooltip:** Caregiver ID, Volume, FE Rate (%), 12-Mo Survival Rate (%), Patients in Sample

---

### Sheet 2b: Volume vs. 12-Mo Survival Rate

- **Data source:** Caregiver Summary
- **Mark type:** Circle
- **Columns:** `[caregiver_n]`
- **Rows:** `[survival_12mo_rate]`
- **Color:** `[survival_12mo_rate]` — sequential blue (darker = higher survival)
- **Size:** `[n_patients_in_sample]`
- **Opacity:** 70%
- **Title:** "Caregiver Volume vs. Patient 12-Month Survival Rate"
- **Trend line:** Linear
- **Reference line:** `AVG([survival_12mo_rate])`, dashed

---

### Sheet 2c: Volume Tier Box Plots

- **Data source:** Caregiver Summary
- **Mark type:** Circle (individual caregiver dots)
- **Columns:** `[Caregiver Volume Tier]`
  *(Create this field in Caregiver Summary with the same formula as in Patient Records)*
- **Rows:** `[observed_fe_rate]`
- **Color:** `[Caregiver Volume Tier]` — 4 sequential blues
- **Sort columns:** Alphabetical (the number prefixes enforce Low→Very High order)
- **Title:** "FE Rate Distribution by Volume Tier"

**Add box plot overlay:**
- Analytics pane (right side) > drag "Box Plot" onto the view
- This overlays median, IQR, and whiskers on the individual dots

**Tooltip:** Volume tier, N caregivers in tier, Median FE Rate, Min/Max

---

### Dashboard 2 Filter

- Field: `[Caregiver Volume Tier]` (from Caregiver Summary)
- Control: Radio button (single select + "All")
- Apply to: Sheets 2a and 2b only (NOT 2c, which is the tier breakdown itself)

---

## DASHBOARD 3: Risk Score Simulator

**Purpose:** Interactive logistic regression simulator. Sliders adjust patient parameters; dashboard shows estimated 12-month survival probability and factor contributions.
**Story beat:** "This is the kind of tool MDCalc builds. I built the statistical model; this shows I can translate it into a product."
**Size:** 1200 × 900 px

### How the simulator works

Tableau implements the logistic equation directly as calculated fields driven by parameters. Coefficients are stylized to match published literature effect directions — **not literal MIMIC-derived values** — keeping the simulator fully independent of restricted data.

---

### Step 1: Create Parameters

**Data pane > right-click > Create Parameter** for each:

| Parameter Name | Data Type | Min | Max | Default | Step | Control |
|---|---|---|---|---|---|---|
| `P: Charlson Index` | Integer | 0 | 15 | 5 | 1 | Slider |
| `P: SOFA Score` | Integer | 0 | 20 | 8 | 1 | Slider |
| `P: Vasopressor Dose` | Float | 0.00 | 0.99 | 0.00 | 0.05 | Slider |
| `P: Vent Hours` | Integer | 1 | 500 | 96 | 12 | Slider |
| `P: Caregiver FE Rate` | Float | 0.00 | 0.50 | 0.13 | 0.01 | Slider |
| `P: Intubation Type` | String | — | — | "Surgical" | — | Dropdown |

For `P: Intubation Type`, set allowable values to: `Surgical`, `Medical – Respiratory`, `Medical – Non-Respiratory`

Right-click each parameter > **Show Parameter** to make it visible on the dashboard.

---

### Step 2: Simulator Calculated Fields

Create all of these in the **Patient Records** data source.

**`Sim: Intubation Effect`**
```
CASE [P: Intubation Type]
  WHEN "Surgical" THEN -0.30
  WHEN "Medical – Respiratory" THEN 0.10
  WHEN "Medical – Non-Respiratory" THEN 0.20
  ELSE 0
END
```

**`Sim: Linear Predictor`**
```
1.70
+ [Sim: Intubation Effect]
+ (-0.18) * ([P: Charlson Index] - 5)
+ (-0.20) * ([P: SOFA Score] - 8)
+ (-1.80) * ([P: Vasopressor Dose] - 0.05)
+ (-0.004) * ([P: Vent Hours] - 138)
+ (-3.0) * [P: Caregiver FE Rate]
+ (5.0) * [P: Caregiver FE Rate] * [P: Caregiver FE Rate]
```
*(Centering constants match synthetic data means: Charlson=5, SOFA=8, Norepi=0.05, Vent=138h)*

**`Sim: Predicted Survival`**
```
1 / (1 + EXP(-[Sim: Linear Predictor]))
```

**`Sim: Survival Pct`**
```
[Sim: Predicted Survival] * 100
```

**`Sim: Risk Category`**
```
IF [Sim: Predicted Survival] >= 0.80 THEN "Lower Risk"
ELSEIF [Sim: Predicted Survival] >= 0.60 THEN "Moderate Risk"
ELSE "Higher Risk"
END
```

**Factor contribution fields** (deviation from population mean — used in waterfall chart):

**`Contrib: Charlson`**
```
(-0.18) * ([P: Charlson Index] - 5)
```

**`Contrib: SOFA`**
```
(-0.20) * ([P: SOFA Score] - 8)
```

**`Contrib: Vasopressor`**
```
(-1.80) * ([P: Vasopressor Dose] - 0.05)
```

**`Contrib: Vent Hours`**
```
(-0.004) * ([P: Vent Hours] - 138)
```

**`Contrib: Caregiver FE Rate`**
```
((-3.0) * [P: Caregiver FE Rate] + (5.0) * [P: Caregiver FE Rate] * [P: Caregiver FE Rate])
- ((-3.0) * 0.13 + (5.0) * 0.13 * 0.13)
```
*(Subtracts baseline at mean FE rate 0.13, so a mean-rate caregiver contributes 0)*

**`Contrib: Intubation Type`**
```
[Sim: Intubation Effect] - (-0.30)
```
*(Deviation from Surgical baseline)*

---

### Step 3: Simulator Sheets

#### Sheet 3a: Predicted Probability (headline number)

- **Mark type:** Text
- **Rows:** `MIN([Sim: Survival Pct])`
- Format the text mark: 40pt bold, centered
- **Color:** `MIN([Sim: Risk Category])` — green / amber / red matching style guide
- No axes, no headers, no gridlines
- Add a second text mark row with `MIN([Sim: Risk Category])` at 14pt below the number

This sheet shows a single updating number like **"76.3%"**.

#### Sheet 3b: Probability Gauge Bar

- **Mark type:** Bar
- **Rows:** `MIN([Sim: Predicted Survival])`
- **Columns:** `MIN(1)` (constant anchor)
- **Color:** `MIN([Sim: Risk Category])` — same palette
- Fix Y-axis range: 0 to 1.0
- **Reference line:** At 0.739 (population survival mean), dashed gray, label "Population avg."
- Make the bar wide — this is purely decorative reinforcement of the number in 3a

#### Sheet 3c: Factor Contributions (horizontal bar chart)

This is the most technically interesting piece — it shows *why* the estimate is where it is.

**Setup approach:** Create a small third data source `simulator_factors.csv` with exactly 6 rows:

```csv
factor_label,factor_order,factor_field
Charlson Index,1,charlson
SOFA Score,2,sofa
Vasopressor Dose,3,vasopressor
Vent Hours,4,vent
Intubation Type,5,intubation
Caregiver FE Rate,6,caregiver
```

Connect as `Factor Labels`. Then build:

- **Mark type:** Bar (horizontal)
- **Rows:** `[factor_label]` from Factor Labels (sort by `factor_order`)
- **Columns:** Calculated field using `CASE [factor_field]` to route each label to its contribution field
- **Color:** Diverging — positive contributions (help survival) in teal, negative in coral
- **Reference line:** At 0 — vertical line marking "no effect"
- **Title:** "What's Driving This Estimate?"
- **Subtitle / annotation:** "Bars show deviation from population-average patient. Teal = better than average; coral = worse."
- **Tooltip:** Factor name, contribution value (formatted as log-odds deviation)

**Simpler alternative** if blending feels complex: Build 6 tiny horizontal bar sheets (one per factor), stack them vertically in a layout container. Less elegant but no data blending required.

---

### Dashboard 3 Layout

```
┌───────────────────────────────────────────────────────────────┐
│  TITLE: "ICU Extubation Risk Estimator"     [disclaimer text] │  ~60px
├──────────────────────┬────────────────────────────────────────┤
│  PARAMETERS          │  Sheet 3a: "76.3%"  (big number)      │
│  (parameter          │  [Lower Risk]                          │
│   controls shown     │                                        │
│   here via           │  Sheet 3b: gauge bar                   │
│   Show Parameter)    │                                        │  ~420px
│                      │  "Estimated 12-month survival"         │
│  P: Charlson Index   │  "Population average: 73.9%"          │
│  P: SOFA Score       │                                        │
│  P: Vasopressor      │                                        │
│  P: Vent Hours       │                                        │
│  P: Caregiver FE     │                                        │
│  P: Intubation Type  │                                        │
├──────────────────────┴────────────────────────────────────────┤
│  Sheet 3c: Factor contributions bar chart                     │  ~300px
│  "What's driving this estimate?"                              │
└───────────────────────────────────────────────────────────────┘
```

**Required disclaimer text object:**
> *"This simulator uses a logistic regression model calibrated to a synthetic dataset generated from published ICU extubation literature. It is for educational and portfolio demonstration purposes only and must not be used for clinical decision-making."*

---

## Navigation Buttons

Add across all three dashboards using **Objects > Navigation** (Tableau Desktop 2020.1+):

- `📋 Cohort Overview` → Dashboard 1
- `👩‍⚕️ Provider Volume` → Dashboard 2
- `🎛 Risk Simulator` → Dashboard 3

Place in a horizontal container pinned to the top of each dashboard. Style with accent color (#2E75B6), white text, consistent 40px height.

---

## Visual Style Guide

| Element | Value |
|---|---|
| Background | White (#FFFFFF) |
| Header band | Light gray (#F5F5F5) |
| Accent / nav | Steel blue (#2E75B6) |
| Survived / lower risk | Teal (#1D9E8E) |
| Not survived / higher risk | Coral (#D95F02) |
| Moderate risk | Amber (#F39C12) |
| Gridlines | None (remove everywhere) |
| Chart borders | 1px light gray (#E0E0E0) |
| Title font | Tableau Book or Arial, 14pt bold |
| Label font | Tableau Book or Arial, 10pt |
| Tooltips | Customize all; remove sheet name |

---

## Tableau Public Publishing Checklist

- [ ] Hide all worksheet tabs — only dashboards visible
- [ ] Customize URL slug: `icu-extubation-provider-volume-outcomes`
- [ ] Verify Dashboard 1 is first (used as thumbnail)
- [ ] Write Tableau Public description (2–3 sentences; note synthetic data)
- [ ] Add tags: `healthcare`, `ICU`, `logistic-regression`, `clinical`, `data-science`
- [ ] Test all parameter sliders after publishing
- [ ] Test navigation buttons after publishing

---

## Interview Talking Points

**"Walk me through this."**
> "I built this on a synthetic cohort modeled after MIMIC-IV ICU extubation data — I couldn't publish the real records due to the PhysioNet user agreement, so I generated synthetic data preserving the distributional structure from published literature. The question is whether provider volume predicts patient outcomes after extubation — a volume-outcome relationship well-documented in surgery, less studied in ventilator weaning."

**"Why synthetic data?"**
> "Knowing what you can't share is as important as knowing how to analyze it. I parameterized the cohort from Fernandez et al. 2024 for failed extubation rates and MIMIC-III published demographics for age and comorbidity distributions — so the numbers are grounded even though no patient records were used."

**"What did the analysis find?"**
> "Low-volume caregivers showed higher and more variable failed extubation rates, consistent with a learning curve effect. We looked for case mix as an explanation using LASSO across ICD chapters — it didn't account for the variation. The volume effect persists after diagnosis adjustment."

**"How does the risk simulator work?"**
> "It implements the logistic regression equation directly as Tableau calculated fields driven by parameter sliders. The coefficients are stylized to match published effect directions — not the literal MIMIC-fitted values, which keeps the simulator fully independent of restricted data. The contribution chart shows how each factor deviates from a population-average patient, which is exactly the kind of explainability framing clinical decision-support tools need."
