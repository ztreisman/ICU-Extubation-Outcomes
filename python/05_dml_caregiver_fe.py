"""
Script 05 — Double/Debiased ML (DML) for Caregiver FE Rate Effect

Estimates the causal effect of caregiver_fe_rate on 12-month mortality using
the partially linear DML estimator (Robinson 1988, Chernozhukov et al. 2018).

Structural model:
    Y = θ·D + g(W) + ε
    D = m(W) + v

    Y = mortality_12mo (binary, treated as continuous in the linear projection)
    D = caregiver_fe_rate (continuous treatment)
    W = all other patient / caregiver / unit confounders
    θ = average partial effect of FE rate on P(mortality)

Procedure:
    1. K-fold cross-fitting: for each held-out fold, train nuisance models on
       the remaining folds and produce out-of-sample residuals
         Ỹ = Y − Ê[Y|W]    (outcome residual)
         D̃ = D − Ê[D|W]    (treatment residual)
    2. Pool across folds: θ̂ = Σ(D̃·Ỹ) / Σ(D̃²)  (Frisch-Waugh)
    3. Sandwich (HC) standard error from the influence-function score

Nuisance models: GradientBoosting — classifier for E[Y|W], regressor for E[D|W].

θ is on the probability scale (change in P(mortality) per unit of caregiver_fe_rate).
The R glmmTMB coefficient is on the log-odds scale; rough comparison:
    log-odds coef × p̄(1−p̄) ≈ probability-scale effect at the mean mortality rate.

Outputs:
    figures/dml_caregiver_fe.png  — partialled-out scatter + residual distribution
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

N_FOLDS = 5
SEED = config.SEED


# ── Feature preparation ────────────────────────────────────────────────────────

def build_Wmat(df: pd.DataFrame, ccsr_cols: list[str]):
    """
    Build confounder matrix W (all features except caregiver_fe_rate).

    Continuous features are log1p-transformed where skewed, then standardised.
    Categoricals are one-hot encoded (GBM does not use embeddings).
    Returns W (ndarray), D (ndarray), feature_names (list).
    """
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
    for col in config.CATEGORICAL_FEATURES:
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
        subsample=0.8, min_samples_leaf=20, random_state=SEED, **kw
    )


def _gbm_treatment(**kw):
    return GradientBoostingRegressor(
        n_estimators=300, max_depth=4, learning_rate=0.05,
        subsample=0.8, min_samples_leaf=20, random_state=SEED, **kw
    )


def run_dml(df: pd.DataFrame, ccsr_cols: list[str]) -> dict:
    W, D, feature_names = build_Wmat(df, ccsr_cols)
    Y = df[config.TARGET].values.astype(float)
    n = len(Y)

    Y_res = np.empty(n)
    D_res = np.empty(n)
    Y_hat_all = np.empty(n)

    kf = KFold(n_splits=N_FOLDS, shuffle=True, random_state=SEED)

    print(f"\nCross-fitting {N_FOLDS} folds ...")
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
        Y_hat_all[te_idx] = Y_hat

        auc = roc_auc_score(Y_te, Y_hat)
        d_var_explained = 1 - np.var(D_te - D_hat) / np.var(D_te)
        print(f"  Fold {fold}: outcome AUC = {auc:.4f} | treatment R² = {d_var_explained:.4f}")

    # Frisch-Waugh point estimate
    theta = np.dot(D_res, Y_res) / np.dot(D_res, D_res)

    # Influence-function sandwich SE
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
    }


# ── Results reporting ──────────────────────────────────────────────────────────

def print_results(res: dict) -> None:
    theta = res["theta"]
    print("\n" + "═" * 60)
    print("DML — caregiver_fe_rate → 12-month mortality")
    print("═" * 60)
    print(f"  θ̂   = {theta:+.4f}  (ΔP(mortality) per unit of FE rate)")
    print(f"  SE   = {res['se']:.4f}")
    print(f"  95% CI: [{res['ci_lo']:+.4f}, {res['ci_hi']:+.4f}]")
    print(f"  z = {res['z']:.3f},  p = {res['p']:.4f}")

    iqr_lo, iqr_hi = np.percentile(res["D"], [25, 75])
    iqr_range = iqr_hi - iqr_lo
    iqr_effect = theta * iqr_range
    print(f"\n  FE rate IQR: {iqr_lo:.4f} → {iqr_hi:.4f}  (range {iqr_range:.4f})")
    print(f"  IQR mortality difference: {iqr_effect:+.3f} ({iqr_effect*100:+.1f} pp)")

    # Rough log-odds comparison: coef_logodds ≈ θ / (p̄·(1−p̄))
    p_bar = res["Y"].mean()
    logodds_approx = theta / (p_bar * (1 - p_bar))
    print(f"\n  Approximate log-odds equivalent: {logodds_approx:+.3f}")


# ── Figure ─────────────────────────────────────────────────────────────────────

def plot_results(res: dict) -> None:
    config.FIGURES_DIR.mkdir(exist_ok=True)
    D_res, Y_res, theta = res["D_res"], res["Y_res"], res["theta"]

    fig, axes = plt.subplots(1, 2, figsize=(12, 5))

    # Panel 1: binned partialled-out scatter
    ax = axes[0]
    quantiles = np.linspace(0, 1, 22)
    bin_edges = np.quantile(D_res, quantiles)
    bin_idx = np.digitize(D_res, bin_edges[1:-1])
    n_bins = len(bin_edges) - 1

    centers = [(bin_edges[i] + bin_edges[i + 1]) / 2 for i in range(n_bins)]
    means   = [Y_res[bin_idx == i].mean() for i in range(n_bins)]
    counts  = [(bin_idx == i).sum() for i in range(n_bins)]

    ax.scatter(centers, means, s=[c / 2 for c in counts],
               alpha=0.75, color="steelblue", zorder=3)
    x_line = np.linspace(D_res.min(), D_res.max(), 200)
    ax.plot(x_line, theta * x_line, color="tomato", lw=2,
            label=f"θ̂ = {theta:+.4f}  (p = {res['p']:.3f})")
    ax.axhline(0, color="gray", lw=0.8, ls="--")
    ax.axvline(0, color="gray", lw=0.8, ls="--")
    ax.set_xlabel("D̃  (FE rate residual after partialling out confounders)")
    ax.set_ylabel("Ỹ  (mortality residual)")
    ax.set_title("Partialled-out relationship\n(binned means, size ∝ n)")
    ax.legend(fontsize=9)

    # Panel 2: treatment residual distribution
    ax = axes[1]
    ax.hist(D_res, bins=60, color="steelblue", alpha=0.7, edgecolor="white")
    ax.axvline(0, color="tomato", lw=1.5, ls="--")
    ax.set_xlabel("D̃  (FE rate residual)")
    ax.set_ylabel("Count")
    ax.set_title("Treatment residual distribution\n"
                 "(variation in FE rate not explained by confounders)")

    plt.tight_layout()
    outpath = config.FIGURES_DIR / "dml_caregiver_fe.png"
    plt.savefig(outpath, dpi=150, bbox_inches="tight")
    print(f"\nFigure saved to {outpath}")


# ── Entry point ────────────────────────────────────────────────────────────────

def main() -> dict:
    print("=" * 60)
    print("Script 05 — Double/Debiased ML (DML)")
    print("=" * 60)

    df = load_raw()
    print(f"Loading last_extubations.csv ...")
    print(f"  Raw rows: {len(df):,}")

    df = parse_ccsr_codes(df)
    df = filter_explicit_cohort(df)
    df, ccsr_cols = build_ccsr_flags(df)

    res = run_dml(df, ccsr_cols)
    print_results(res)
    plot_results(res)

    return res


if __name__ == "__main__":
    main()
