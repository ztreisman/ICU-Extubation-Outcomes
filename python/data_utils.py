"""
Data loading and feature engineering for the ICU Extubation DL pipeline.

Produces the same analytical cohort as 01_cohort_and_descriptive.R:
  - explicit_extubations: tube_event_source == 'explicit',
    caregiver_fe_rate not NULL, caregiver_n > 10
  - survival_12mo: derived from dod (matching R's is.na(dod))
  - CCSR binary flags: >= 10% prevalence, excluding RSP012

Outputs a PreparedData namedtuple with:
  - X_train, X_val, X_test  : float32 numpy arrays (model-ready)
  - y_train, y_val, y_test  : int32 numpy arrays (0/1)
  - feature_names           : list of feature names in column order
  - cat_dims                : list of (n_categories,) per categorical feature
  - cat_emb_dims            : list of embedding dims per categorical feature
  - cat_start_idx           : index in X where categorical columns start
  - ccsr_cols               : list of CCSR column names used
  - encoders                : dict of fitted LabelEncoders for categoricals
  - scaler                  : fitted StandardScaler for continuous features
  - df_explicit             : the unencoded explicit_extubations DataFrame
"""

import ast
import re
from collections import namedtuple
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.model_selection import train_test_split
from sklearn.preprocessing import LabelEncoder, StandardScaler

import config

PreparedData = namedtuple("PreparedData", [
    "X_train", "X_val", "X_test",
    "y_train", "y_val", "y_test",
    "feature_names",
    "cat_dims", "cat_emb_dims", "cat_start_idx",
    "ccsr_cols",
    "encoders", "scaler",
    "df_explicit",
])


# ── 1. Load raw CSV ────────────────────────────────────────────────────────────

def load_raw(path: Path = config.LAST_EXTUBATIONS_CSV) -> pd.DataFrame:
    dtype_map = {
        "subject_id":                  str,
        "hadm_id":                     str,
        "caregiver_id":                str,
        "gender":                      str,
        "intubation_type":             str,
        "tube_event_source":           str,
        "caregiver_imputation_source": str,
        "first_careunit":              str,
        "primary_diagnosis":           str,
        "recorded_intubation":         str,
        "hospital_expire_flag":        "Int64",
    }
    df = pd.read_csv(
        path,
        dtype=dtype_map,
        parse_dates=["event_time", "dod"],
        low_memory=False,
    )
    # Derive 12-month survival: dead if dod is within 12 months of extubation
    # Matches R's is.na(dod) for patients with no recorded death;
    # patients who died > 12 months after extubation are also survivors.
    df["survival_12mo"] = (
        df["dod"].isna() |
        ((df["dod"] - df["event_time"].dt.normalize()).dt.days > 365)
    ).astype(int)
    return df


# ── 2. Parse ccsr_codes array column ──────────────────────────────────────────
# BigQuery CSV export represents ARRAY<STRING> in several possible formats;
# this parser handles the most common ones.

def _parse_ccsr_cell(val) -> list[str]:
    """Return a list of CCSR code strings from a single cell value."""
    if pd.isna(val) or val == "" or val == "[]":
        return []
    s = str(val).strip()
    # BigQuery-style: ['RSP011', 'CAR002', ...]
    if s.startswith("["):
        try:
            codes = ast.literal_eval(s)
            return [c.strip() for c in codes if isinstance(c, str)]
        except Exception:
            pass
    # Bare comma-separated: RSP011,CAR002
    return [c.strip() for c in re.split(r"[,\s]+", s) if c.strip()]


def parse_ccsr_codes(df: pd.DataFrame) -> pd.DataFrame:
    """Parse the ccsr_codes array column in-place; returns df."""
    if "ccsr_codes" not in df.columns:
        df["ccsr_codes_list"] = [[] for _ in range(len(df))]
    else:
        df["ccsr_codes_list"] = df["ccsr_codes"].apply(_parse_ccsr_cell)
    return df


# ── 3. Build CCSR binary flags (mirrors 03_psm_and_sensitivity.R) ─────────────

def build_ccsr_flags(
    df: pd.DataFrame,
    threshold: float = config.CCSR_PREVALENCE_THRESHOLD,
    exclude: list[str] = config.CCSR_EXCLUDE,
    ccsr_long_path: Path = config.CCSR_LONG_CSV,
) -> tuple[pd.DataFrame, list[str]]:
    """
    Add binary CCSR flag columns to df.

    If patient_ccsr_long.csv exists, use it (more reliable unnested source).
    Otherwise parse ccsr_codes_list from the array column.

    Returns (df_with_flags, selected_ccsr_cols).
    """
    n_patients = len(df)

    # Prefer the separate unnested CSV if available
    if ccsr_long_path.exists():
        ccsr_long = pd.read_csv(
            ccsr_long_path,
            dtype={"subject_id": str, "hadm_id": str, "ccsr_code": str},
            usecols=["subject_id", "ccsr_code"],
        )
        wide = (
            ccsr_long[ccsr_long["subject_id"].isin(df["subject_id"])]
            .drop_duplicates(["subject_id", "ccsr_code"])
            .assign(present=1)
            .pivot(index="subject_id", columns="ccsr_code", values="present")
            .fillna(0)
            .astype(int)
        )
    else:
        # Fall back to parsing the ccsr_codes_list column
        if "ccsr_codes_list" not in df.columns:
            df = parse_ccsr_codes(df)
        all_codes = {code for codes in df["ccsr_codes_list"] for code in codes}
        wide = pd.DataFrame(
            {
                code: df["ccsr_codes_list"].apply(lambda lst: int(code in lst))
                for code in all_codes
            },
            index=df.index,
        )
        wide.index = df["subject_id"].values

    # Compute prevalence in this cohort
    prevalence = wide.mean()
    selected = prevalence[
        (prevalence >= threshold) & (~prevalence.index.isin(exclude))
    ].index.tolist()
    selected.sort()

    print(f"\nCCSR codes >= {threshold*100:.0f}% prevalence: {len(selected)} (excluding {exclude})")

    # Merge flags back to df on subject_id
    flag_df = wide[selected].reset_index().rename(columns={"index": "subject_id"})
    flag_df.columns.name = None
    df = df.merge(flag_df, on="subject_id", how="left")
    for col in selected:
        df[col] = df[col].fillna(0).astype(int)

    ccsr_cols = selected
    return df, ccsr_cols


# ── 4. Apply cohort filters (matching explicit_extubations in R) ───────────────

def filter_explicit_cohort(df: pd.DataFrame) -> pd.DataFrame:
    before = len(df)
    df = df[
        (df["tube_event_source"] == "explicit") &
        (df["caregiver_fe_rate"].notna()) &
        (df["caregiver_n"] > config.MIN_CAREGIVER_N) &
        (df["sofa"].notna()) &
        (df["charlson"].notna()) &
        (df["primary_diagnosis"].notna()) &
        (df["vent_hours"].notna())
    ].copy()
    print(f"\nCohort filter: {before:,} → {len(df):,} patients (explicit_extubations)")
    print(f"  12-month survival rate: {df['survival_12mo'].mean():.3f}")
    return df


# ── 5. Feature encoding ────────────────────────────────────────────────────────

def _emb_dim(n_cat: int) -> int:
    """Standard embedding dimension: min(50, (n+1)//2)."""
    return min(50, (n_cat + 1) // 2)


def encode_features(
    df_train: pd.DataFrame,
    df_val: pd.DataFrame,
    df_test: pd.DataFrame,
    ccsr_cols: list[str],
) -> tuple:
    """
    Fit encoders on training set; apply to val and test.

    Continuous features are log1p-transformed (vent_hours, caregiver_n)
    before scaling to address right-skew.

    Returns:
        X_train, X_val, X_test  : float32 numpy arrays
        feature_names           : column names in order
        cat_dims                : number of categories per categorical feature
        cat_emb_dims            : embedding dim per categorical feature
        cat_start_idx           : column index where categoricals begin
        encoders                : dict of LabelEncoder per categorical feature
        scaler                  : fitted StandardScaler
    """
    # ── Continuous ─────────────────────────────────────────────────────────────
    log_transform_cols = {"vent_hours", "caregiver_n"}
    cont_cols = config.CONTINUOUS_FEATURES

    def _prep_continuous(df: pd.DataFrame) -> np.ndarray:
        out = []
        for col in cont_cols:
            vals = df[col].fillna(df[col].median()).values.astype(float)
            if col in log_transform_cols:
                vals = np.log1p(vals)
            out.append(vals)
        return np.column_stack(out)

    X_cont_train = _prep_continuous(df_train)
    X_cont_val   = _prep_continuous(df_val)
    X_cont_test  = _prep_continuous(df_test)

    scaler = StandardScaler()
    X_cont_train = scaler.fit_transform(X_cont_train)
    X_cont_val   = scaler.transform(X_cont_val)
    X_cont_test  = scaler.transform(X_cont_test)

    # ── Categorical (label-encoded integers for embedding lookup) ───────────────
    encoders = {}
    X_cat_parts_train, X_cat_parts_val, X_cat_parts_test = [], [], []
    cat_dims, cat_emb_dims = [], []

    for col in config.CATEGORICAL_FEATURES:
        le = LabelEncoder()
        train_vals = df_train[col].fillna("Unknown").astype(str)
        le.fit(train_vals)

        # Unknown categories in val/test get mapped to a reserved 'Unknown' class
        all_classes = list(le.classes_) + ["Unknown"]
        le.classes_ = np.array(all_classes)
        n_cat = len(le.classes_)

        def _safe_transform(series, le=le):
            vals = series.fillna("Unknown").astype(str)
            return np.array([
                le.transform([v])[0] if v in le.classes_ else le.transform(["Unknown"])[0]
                for v in vals
            ])

        X_cat_parts_train.append(_safe_transform(df_train[col]))
        X_cat_parts_val.append(_safe_transform(df_val[col]))
        X_cat_parts_test.append(_safe_transform(df_test[col]))

        cat_dims.append(n_cat)
        cat_emb_dims.append(_emb_dim(n_cat))
        encoders[col] = le

    X_cat_train = np.column_stack(X_cat_parts_train).astype(np.int64)
    X_cat_val   = np.column_stack(X_cat_parts_val).astype(np.int64)
    X_cat_test  = np.column_stack(X_cat_parts_test).astype(np.int64)

    # ── CCSR multi-hot ─────────────────────────────────────────────────────────
    X_ccsr_train = df_train[ccsr_cols].values.astype(np.float32)
    X_ccsr_val   = df_val[ccsr_cols].values.astype(np.float32)
    X_ccsr_test  = df_test[ccsr_cols].values.astype(np.float32)

    # ── Concatenate: [continuous | CCSR] + [categorical] ───────────────────────
    # Categorical kept separate because they need embedding lookup, not dense ops.
    # Convention: X[:, :cat_start_idx] = continuous+CCSR, X[:, cat_start_idx:] = cats
    cat_start_idx = X_cont_train.shape[1] + X_ccsr_train.shape[1]

    X_train = np.concatenate([X_cont_train, X_ccsr_train, X_cat_train], axis=1).astype(np.float32)
    X_val   = np.concatenate([X_cont_val,   X_ccsr_val,   X_cat_val],   axis=1).astype(np.float32)
    X_test  = np.concatenate([X_cont_test,  X_ccsr_test,  X_cat_test],  axis=1).astype(np.float32)

    feature_names = (
        [f"{c}_log1p_scaled" if c in log_transform_cols else f"{c}_scaled" for c in cont_cols]
        + [f"ccsr_{c}" for c in ccsr_cols]
        + [f"cat_{c}" for c in config.CATEGORICAL_FEATURES]
    )

    return (
        X_train, X_val, X_test,
        feature_names, cat_dims, cat_emb_dims, cat_start_idx,
        encoders, scaler,
    )


# ── 6. Main entry point ────────────────────────────────────────────────────────

def prepare_data() -> PreparedData:
    """
    Full data preparation pipeline. Returns a PreparedData namedtuple.
    Call this from any training or evaluation script.
    """
    np.random.seed(config.SEED)

    print("Loading last_extubations.csv ...")
    df = load_raw()
    print(f"  Raw rows: {len(df):,}")

    df = parse_ccsr_codes(df)
    df = filter_explicit_cohort(df)
    df, ccsr_cols = build_ccsr_flags(df)

    # ── Train / val / test split (stratified on survival_12mo) ─────────────────
    df_train, df_temp = train_test_split(
        df,
        test_size=(1 - config.TRAIN_FRAC),
        stratify=df[config.TARGET],
        random_state=config.SEED,
    )
    val_size_relative = config.VAL_FRAC / (config.VAL_FRAC + config.TEST_FRAC)
    df_val, df_test = train_test_split(
        df_temp,
        test_size=(1 - val_size_relative),
        stratify=df_temp[config.TARGET],
        random_state=config.SEED,
    )

    print(f"\nSplit sizes — train: {len(df_train)}, val: {len(df_val)}, test: {len(df_test)}")
    for name, split in [("train", df_train), ("val", df_val), ("test", df_test)]:
        rate = split[config.TARGET].mean()
        print(f"  {name} survival rate: {rate:.3f}")

    (
        X_train, X_val, X_test,
        feature_names, cat_dims, cat_emb_dims, cat_start_idx,
        encoders, scaler,
    ) = encode_features(df_train, df_val, df_test, ccsr_cols)

    y_train = df_train[config.TARGET].values.astype(np.int32)
    y_val   = df_val[config.TARGET].values.astype(np.int32)
    y_test  = df_test[config.TARGET].values.astype(np.int32)

    print(f"\nFeature matrix shape: {X_train.shape[1]} columns")
    print(f"  Continuous + CCSR: {cat_start_idx}")
    print(f"  Categorical (embedded): {len(config.CATEGORICAL_FEATURES)}")
    print(f"  Cat dims: {list(zip(config.CATEGORICAL_FEATURES, cat_dims))}")
    print(f"  Emb dims: {list(zip(config.CATEGORICAL_FEATURES, cat_emb_dims))}")

    return PreparedData(
        X_train=X_train, X_val=X_val, X_test=X_test,
        y_train=y_train, y_val=y_val, y_test=y_test,
        feature_names=feature_names,
        cat_dims=cat_dims, cat_emb_dims=cat_emb_dims,
        cat_start_idx=cat_start_idx,
        ccsr_cols=ccsr_cols,
        encoders=encoders, scaler=scaler,
        df_explicit=df,
    )
