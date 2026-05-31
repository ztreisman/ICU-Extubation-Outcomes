"""
Script 06 — CVICU-only Double/Debiased ML

Same cross-fit DML estimator as Script 05, restricted to
Cardiac Vascular Intensive Care Unit (CVICU) patients.

Motivation: the global DML (Script 05) pools all ICU units and finds a null
result (theta = -0.172, p = 0.581), diluted by the null Medical and
Surgical/Trauma groups. PSM stratified by unit finds a significant signal in
CVICU (OR = 0.581, p = 0.027; IPW OR = 0.656, p = 0.005). This script
isolates the causal estimate using only within-CVICU variation in caregiver
FE rate, after partialling out patient complexity within that unit.

Key differences from Script 05:
  - Filtered to CVICU patients only (~1,536 in the analytical cohort)
  - first_careunit excluded from the confounder matrix (zero variance in a
    single-unit analysis; including it would add nothing)
  - CCSR prevalence threshold reapplied to the CVICU subset (CVICU is 84%
    circulatory, so the prevalent codes differ from the full cohort)
  - min_samples_leaf reduced to 10 (fold training sets ~1,200 vs ~3,000
    globally; keeps GBM from becoming too shallow)

Outputs:
  figures/dml_cvicu.png  — partialled-out scatter + residual distribution
"""

import os
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats
from sklearn.ensemble import GradientBoostingClassifier, GradientBoostingRegressor
from sklearn.metrics import roc_auc_score
from sklearn.model_selection import KFold
from sklearn.preprocessing import StandardScaler

sys.path.insert(0, os.path.dirname(__file__))
import config
from data_utils import (
    build_ccsr_flags,
    filter_explicit_cohort,
    load_raw,
    parse_ccsr_codes,
)

CVICU_LABEL = "Cardiac Vascular Intensive Care Unit (CVICU)"
N_FOLDS = 5
SEED = config.SEED

# Exclude first_careunit — constant within CVICU, adds no information
CATEGORICAL_CONFOUNDERS = [c for c in config.CATEGORICAL_FEATURES if c != "first_careunit"]


# ── Feature preparation ────────────────────────────────────────────────────────

def build_Wmat(df: pd.DataFrame, ccsr_cols: list[str]) -> tuple:
    log_cols = {"vent_hours", "caregiver_n"}
    cont_conf = [f for f in config.CONTINUOUS_FEATURES if f != "caregiver_fe_rate"]

    cont_parts, cont_names = [], []
    for col in cont_conf:
        vals = df[col].fillna(df[col].median()).values.astype(float)
        if col in log_cols:
            vals = np.log1p(vals)
        cont_parts.append(vals)
        cont_names.append(col + ("_log1p" if col in log_cols else ""))

    W_cont = StandardScaler().fit_transform(np.column_stack(cont_parts))

    cat_parts, cat_names = [], []
    for col in CATEGORICAL_CONFOUNDERS:
        dummies = pd.get_dummies(
            df[col].fillna("Unknown").astype(str), prefix=col, drop_first=False
        )
        cat_parts.append(dummies.values.astype(float))
        cat_names.extend(dummies.columns.tolist())

    W_ccsr = df[ccsr_cols].values.astype(float)
    ccsr_names = [f"ccsr_{c}" for c in ccsr_cols]

    W = np.concatenate([W_cont, W_ccsr] + cat_parts, axis=1)
    feature_names = cont_names + ccsr_names + cat_names

    D = df["caregiver_fe_rate"].values.astype(float)
    return W, D, feature_names


# ── Cross-fit DML ──────────────────────────────────────────────────────────────

def _gbm_outcome(**kw):
    return GradientBoostingClassifier(
        n_estimators=300, max_depth=4, learning_rate=0.05,
        subsample=0.8, min_samples_leaf=10, random_state=SEED, **kw
    )


def _gbm_treatment(**kw):
    return GradientBoostingRegressor(
        n_estimators=300, max_depth=4, learning_rate=0.05,
        subsample=0.8, min_samples_leaf=10, random_state=SEED, **kw
    )


def run_dml(df: pd.DataFrame, ccsr_cols: list[str]) -> dict:
    W, D, feature_names = build_Wmat(df, ccsr_cols)
    Y = df[config.TARGET].values.astype(float)
    n = len(Y)

    Y_res = np.empty(n)
    D_res = np.empty(n)

    kf = KFold(n_splits=N_FOLDS, shuffle=True, random_state=SEED)

    print(f"\nCross-fitting {N_FOLDS} folds (N = {n:,}) ...")
    for fold, (tr_idx, te_idx) in enumerate(kf.split(W), 1):
        W_tr, W_te = W[tr_idx], W[te_idx]
        Y_tr, D_tr = Y[tr_idx], D[tr_idx]
        Y_te, D_te = Y[te_idx], D[te_idx]

        out_m = _gbm_outcome()
        out_m.fit(W_tr, Y_tr)
        Y_hat = out_m.predict_proba(W_te)[:, 1]

        trt_m = _gbm_treatment()
        trt_m.fit(W_tr, D_tr)
        D_hat = trt_m.predict(W_te)

        Y_res[te_idx] = Y_te - Y_hat
        D_res[te_idx] = D_te - D_hat

        auc = roc_auc_score(Y_te, Y_hat)
        d_r2 = 1 - np.var(D_te - D_hat) / np.var(D_te)
        print(f"  Fold {fold}: outcome AUC = {auc:.4f} | treatment R² = {d_r2:.4f}")

    theta = np.dot(D_res, Y_res) / np.dot(D_res, D_res)

    psi = D_res * (Y_res - theta * D_res)
    denom = np.sum(D_res ** 2)
    se = np.sqrt(np.sum(psi ** 2) / denom ** 2)

    z = theta / se
    p = 2 * (1 - stats.norm.cdf(abs(z)))

    return {
        "theta": theta, "se": se, "z": z, "p": p,
        "ci_lo": theta - 1.96 * se,
        "ci_hi": theta + 1.96 * se,
        "D": D, "Y": Y,
        "D_res": D_res, "Y_res": Y_res,
        "feature_names": feature_names,
        "n": n,
    }


# ── Results reporting ──────────────────────────────────────────────────────────

def print_results(res: dict) -> None:
    theta = res["theta"]
    print("\n" + "═" * 60)
    print("DML — CVICU only — caregiver_fe_rate → 12-month survival")
    print("═" * 60)
    print(f"  N patients : {res['n']:,}")
    print(f"  θ̂         = {theta:+.4f}  (ΔP(survival) per unit of FE rate)")
    print(f"  SE         = {res['se']:.4f}")
    print(f"  95% CI     : [{res['ci_lo']:+.4f}, {res['ci_hi']:+.4f}]")
    print(f"  z = {res['z']:.3f},  p = {res['p']:.4f}")

    iqr_lo, iqr_hi = np.percentile(res["D"], [25, 75])
    iqr_range = iqr_hi - iqr_lo
    iqr_effect = theta * iqr_range
    print(f"\n  CVICU FE rate IQR : {iqr_lo:.4f} → {iqr_hi:.4f}  (range {iqr_range:.4f})")
    print(f"  IQR survival diff : {iqr_effect:+.3f} ({iqr_effect*100:+.1f} pp)")

    p_bar = res["Y"].mean()
    logodds_approx = theta / (p_bar * (1 - p_bar))
    print(f"\n  Mean CVICU 12-mo survival : {p_bar:.3f}")
    print(f"  Approx log-odds equivalent: {logodds_approx:+.3f}")
    print(f"\n  For comparison:")
    print(f"    Global DML (all units)  : θ̂ = -0.172, p = 0.581")
    print(f"    CVICU PSM (OR)          : 0.581, p = 0.027")
    print(f"    CVICU IPW (OR)          : 0.656, p = 0.005")


# ── Figure ─────────────────────────────────────────────────────────────────────

def plot_results(res: dict) -> None:
    config.FIGURES_DIR.mkdir(exist_ok=True)
    D_res, Y_res, theta = res["D_res"], res["Y_res"], res["theta"]

    fig, axes = plt.subplots(1, 2, figsize=(12, 5))
    fig.suptitle(f"DML — CVICU only (N = {res['n']:,})", fontsize=12)

    # Panel 1: binned partialled-out scatter
    ax = axes[0]
    quantiles = np.linspace(0, 1, 22)
    bin_edges = np.quantile(D_res, quantiles)
    bin_idx = np.digitize(D_res, bin_edges[1:-1])
    n_bins = len(bin_edges) - 1

    centers, means, counts = [], [], []
    for i in range(n_bins):
        mask = bin_idx == i
        if mask.sum() > 0:
            centers.append((bin_edges[i] + bin_edges[i + 1]) / 2)
            means.append(Y_res[mask].mean())
            counts.append(mask.sum())

    ax.scatter(centers, means, s=[c / 2 for c in counts],
               alpha=0.75, color="steelblue", zorder=3)
    x_line = np.linspace(D_res.min(), D_res.max(), 200)
    ax.plot(x_line, theta * x_line, color="tomato", lw=2,
            label=f"θ̂ = {theta:+.4f}  (p = {res['p']:.3f})")
    ax.axhline(0, color="gray", lw=0.8, ls="--")
    ax.axvline(0, color="gray", lw=0.8, ls="--")
    ax.set_xlabel("D̃  (FE rate residual after partialling out confounders)")
    ax.set_ylabel("Ỹ  (survival residual)")
    ax.set_title("Partialled-out relationship\n(binned means, size ∝ n)")
    ax.legend(fontsize=9)

    # Panel 2: treatment residual distribution
    ax = axes[1]
    ax.hist(D_res, bins=50, color="steelblue", alpha=0.7, edgecolor="white")
    ax.axvline(0, color="tomato", lw=1.5, ls="--")
    ax.set_xlabel("D̃  (FE rate residual)")
    ax.set_ylabel("Count")
    ax.set_title("Treatment residual\n(within-CVICU variation unexplained by patient factors)")

    plt.tight_layout()
    outpath = config.FIGURES_DIR / "dml_cvicu.png"
    plt.savefig(outpath, dpi=150, bbox_inches="tight")
    print(f"\nFigure saved to {outpath}")


# ── Entry point ────────────────────────────────────────────────────────────────

def main() -> dict:
    print("=" * 60)
    print("Script 06 — CVICU-only DML")
    print("=" * 60)

    df = load_raw()
    print(f"Loading last_extubations.csv ...")
    print(f"  Raw rows: {len(df):,}")

    df = parse_ccsr_codes(df)
    df = filter_explicit_cohort(df)

    df_cvicu = df[df["first_careunit"] == CVICU_LABEL].copy()
    print(f"\nCVICU filter: {len(df):,} → {len(df_cvicu):,} patients")
    print(f"  12-mo survival : {df_cvicu['survival_12mo'].mean():.3f}")
    print(f"  FE rate mean   : {df_cvicu['caregiver_fe_rate'].mean():.4f}")
    print(f"  FE rate median : {df_cvicu['caregiver_fe_rate'].median():.4f}")
    print(f"  Caregivers     : {df_cvicu['caregiver_id'].nunique()}")

    df_cvicu, ccsr_cols = build_ccsr_flags(df_cvicu)

    res = run_dml(df_cvicu, ccsr_cols)
    print_results(res)
    plot_results(res)

    return res


if __name__ == "__main__":
    main()
