"""
Central configuration for ICU Extubation Outcomes — Deep Learning Pipeline.

Mirrors the R pipeline's analytical choices (set.seed(237), same cohort
filters, same CCSR threshold) to keep results directly comparable.
"""

from pathlib import Path

# ── Paths ─────────────────────────────────────────────────────────────────────
DATA_DIR    = Path("../data")
FIGURES_DIR = Path("figures")
MODELS_DIR  = Path("models")

LAST_EXTUBATIONS_CSV = DATA_DIR / "last_extubations.csv"
CCSR_LONG_CSV        = DATA_DIR / "patient_ccsr_long.csv"  # unnested fallback

# ── Reproducibility ────────────────────────────────────────────────────────────
SEED = 237  # matches set.seed(237) throughout the R pipeline

# ── Cohort filters (matching explicit_extubations in 01_cohort_and_descriptive.R)
MIN_CAREGIVER_N = 10   # caregiver_n > 10

# ── Features ──────────────────────────────────────────────────────────────────
CONTINUOUS_FEATURES = [
    "anchor_age",         # patient age at admission
    "charlson",           # Charlson comorbidity index
    "sofa",               # first-day SOFA score
    "norepinephrine",     # max norepinephrine-equivalent dose (vasopressor load)
    "vent_hours",         # total ventilation hours (log-transformed in preprocessing)
    "caregiver_fe_rate",  # caregiver's historical failed extubation rate
    "caregiver_n",        # caregiver total extubation volume (log-transformed)
    "failed_extubations", # patient lifetime reintubation count (complexity marker)
]

# Note: vent_hours and caregiver_n are right-skewed; log1p transform applied.

CATEGORICAL_FEATURES = [
    "gender",            # M / F  → embedding dim 2
    "intubation_type",   # surgical / medical-respiratory / medical-non-respiratory
    "primary_diagnosis", # ICD chapter (~15 levels)
    "first_careunit",    # 12 ICU units  (analogous to random intercept in glmmTMB)
]

TARGET = "survival_12mo"

# ── CCSR binary flags ─────────────────────────────────────────────────────────
# Threshold and exclusions match 03_psm_and_sensitivity.R
CCSR_PREVALENCE_THRESHOLD = 0.10
CCSR_EXCLUDE = ["RSP012"]  # near-universal (47%) in ventilated patients; near-zero variance

# ── Train / val / test split ───────────────────────────────────────────────────
TRAIN_FRAC = 0.70
VAL_FRAC   = 0.15
TEST_FRAC  = 0.15

# ── MLP hyperparameters ───────────────────────────────────────────────────────
HIDDEN_DIMS   = [256, 128, 64]
DROPOUT_RATE  = 0.30
LEARNING_RATE = 1e-3
WEIGHT_DECAY  = 1e-4
BATCH_SIZE    = 128
MAX_EPOCHS    = 300
PATIENCE      = 20    # early stopping — stop if val loss doesn't improve

# Embedding dimension rule: min(50, (n_categories + 1) // 2)
# Applied per categorical feature in data_utils.py
