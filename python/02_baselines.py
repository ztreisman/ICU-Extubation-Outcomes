"""
Script 02 — Baseline Models

Logistic regression and gradient-boosted tree baselines.
These serve as performance floors for the MLP:
  - If MLP can't beat logistic regression, the DL architecture adds no value.
  - If gradient boosting beats MLP, the data likely doesn't reward non-linearity
    at this sample size.

Mirrors the frequentist modeling in 02_logistic_and_mixed_effects.R but
without the caregiver-level random effects (those are captured by the
first_careunit embedding in the MLP).

Outputs:
  - Baseline AUC-ROC and AUC-PR printed to console
  - Logistic regression coefficients (for comparison to glmmTMB estimates)
  - ROC curve plot saved to figures/
"""

import os
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from sklearn.ensemble import GradientBoostingClassifier
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import (
    roc_auc_score, average_precision_score,
    roc_curve, precision_recall_curve,
    brier_score_loss,
)
import warnings

sys.path.insert(0, os.path.dirname(__file__))
import config
from data_utils import prepare_data

warnings.filterwarnings("ignore", category=FutureWarning)
np.random.seed(config.SEED)


def evaluate(name: str, y_true, y_prob) -> dict:
    auc    = roc_auc_score(y_true, y_prob)
    ap     = average_precision_score(y_true, y_prob)
    brier  = brier_score_loss(y_true, y_prob)
    print(f"\n{name}")
    print(f"  AUC-ROC : {auc:.4f}")
    print(f"  AUC-PR  : {ap:.4f}")
    print(f"  Brier   : {brier:.4f}")
    return {"auc": auc, "ap": ap, "brier": brier, "name": name, "y_prob": y_prob}


def main():
    print("=" * 60)
    print("Script 02 — Baseline Models")
    print("=" * 60)

    data = prepare_data()

    # Use only continuous + CCSR columns for LR (no label-encoded cats)
    # to match logistic regression assumptions.
    # Categorical columns are integers but aren't meaningful as ordinal.
    X_tr  = data.X_train[:, :data.cat_start_idx]
    X_va  = data.X_val[:, :data.cat_start_idx]
    X_te  = data.X_test[:, :data.cat_start_idx]
    y_tr  = data.y_train
    y_va  = data.y_val
    y_te  = data.y_test

    results = []

    # ── Logistic regression (L2-regularized) ──────────────────────────────────
    print("\nFitting logistic regression ...")
    lr = LogisticRegression(
        C=1.0, max_iter=1000, random_state=config.SEED, solver="lbfgs"
    )
    lr.fit(X_tr, y_tr)
    lr_prob_val  = lr.predict_proba(X_va)[:, 1]
    lr_prob_test = lr.predict_proba(X_te)[:, 1]

    results.append(evaluate("Logistic Regression (val)",  y_va, lr_prob_val))
    results.append(evaluate("Logistic Regression (test)", y_te, lr_prob_test))

    # Print coefficients for the continuous features (comparable to glmmTMB)
    cont_feature_names = [n for n in data.feature_names if not n.startswith("cat_")][:len(config.CONTINUOUS_FEATURES)]
    print("\nLogistic regression coefficients (continuous features):")
    n_cont = len(config.CONTINUOUS_FEATURES)
    for name, coef in zip(cont_feature_names[:n_cont], lr.coef_[0][:n_cont]):
        print(f"  {name:<35} {coef:+.4f}")

    # ── Gradient boosting ──────────────────────────────────────────────────────
    # Include categorical columns (label-encoded as integers) — GBM is invariant
    # to monotone transforms of input features, so label encoding is fine.
    X_tr_full = data.X_train
    X_va_full = data.X_val
    X_te_full = data.X_test

    print("\nFitting gradient boosting classifier ...")
    gb = GradientBoostingClassifier(
        n_estimators=300,
        max_depth=4,
        learning_rate=0.05,
        subsample=0.8,
        min_samples_leaf=20,
        random_state=config.SEED,
        validation_fraction=0.1,
        n_iter_no_change=20,
        tol=1e-4,
        verbose=0,
    )
    gb.fit(X_tr_full, y_tr)
    gb_prob_val  = gb.predict_proba(X_va_full)[:, 1]
    gb_prob_test = gb.predict_proba(X_te_full)[:, 1]

    results.append(evaluate("Gradient Boosting (val)",  y_va, gb_prob_val))
    results.append(evaluate("Gradient Boosting (test)", y_te, gb_prob_test))

    print(f"\nGradient boosting stopped at {gb.n_estimators_} estimators")

    # ── Top GB feature importances ─────────────────────────────────────────────
    importances = gb.feature_importances_
    idx = np.argsort(importances)[::-1][:15]
    print("\nTop 15 Gradient Boosting feature importances:")
    for rank, i in enumerate(idx):
        fname = data.feature_names[i] if i < len(data.feature_names) else f"col_{i}"
        print(f"  {rank+1:2d}. {fname:<40} {importances[i]:.4f}")

    # ── ROC curve plot ─────────────────────────────────────────────────────────
    config.FIGURES_DIR.mkdir(exist_ok=True)
    fig, axes = plt.subplots(1, 2, figsize=(12, 5))

    ax = axes[0]
    for label, prob in [
        ("Logistic Reg", lr_prob_test),
        ("Grad Boosting", gb_prob_test),
    ]:
        fpr, tpr, _ = roc_curve(y_te, prob)
        auc = roc_auc_score(y_te, prob)
        ax.plot(fpr, tpr, label=f"{label} (AUC={auc:.3f})")
    ax.plot([0, 1], [0, 1], "k--", alpha=0.4)
    ax.set_xlabel("False Positive Rate")
    ax.set_ylabel("True Positive Rate")
    ax.set_title("ROC Curves — Baseline Models (test set)")
    ax.legend()

    ax = axes[1]
    baseline_rate = y_te.mean()
    for label, prob in [
        ("Logistic Reg", lr_prob_test),
        ("Grad Boosting", gb_prob_test),
    ]:
        prec, rec, _ = precision_recall_curve(y_te, prob)
        ap = average_precision_score(y_te, prob)
        ax.plot(rec, prec, label=f"{label} (AP={ap:.3f})")
    ax.axhline(baseline_rate, linestyle="--", color="k", alpha=0.4, label=f"Baseline ({baseline_rate:.2f})")
    ax.set_xlabel("Recall")
    ax.set_ylabel("Precision")
    ax.set_title("Precision-Recall Curves — Baseline Models (test set)")
    ax.legend()

    plt.tight_layout()
    out_path = config.FIGURES_DIR / "baseline_roc_pr_curves.png"
    plt.savefig(out_path, dpi=150, bbox_inches="tight")
    print(f"\nFigure saved to {out_path}")

    return results, {"lr": lr, "gb": gb}


if __name__ == "__main__":
    main()
