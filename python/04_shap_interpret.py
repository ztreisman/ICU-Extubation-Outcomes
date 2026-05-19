"""
Script 04 — SHAP Interpretability

Uses SHAP (SHapley Additive exPlanations) to explain the trained MLP's
predictions. SHAP values answer: "how much did each feature push this
patient's predicted survival probability up or down from the baseline?"

Two SHAP approaches:
  DeepExplainer  — exact Shapley values for PyTorch models (fast, exact)
  KernelExplainer — model-agnostic fallback (slow, use if DeepExplainer fails)

Outputs:
  figures/shap_summary.png       — beeswarm plot (global feature importance)
  figures/shap_bar.png           — mean |SHAP| bar chart
  figures/shap_caregiver_fe.png  — dependence plot for caregiver_fe_rate
"""

import os
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import torch
import shap

sys.path.insert(0, os.path.dirname(__file__))
import config
from data_utils import prepare_data
from model import MLPWithEmbeddings

DEVICE = torch.device("cpu")  # SHAP works better on CPU


def load_model(data) -> MLPWithEmbeddings:
    ckpt = config.MODELS_DIR / "mlp_best.pt"
    if not ckpt.exists():
        raise FileNotFoundError(
            f"No checkpoint at {ckpt}. Run 03_train_mlp.py first."
        )
    model = MLPWithEmbeddings(
        continuous_dim=data.cat_start_idx,
        cat_dims=data.cat_dims,
        cat_emb_dims=data.cat_emb_dims,
        hidden_dims=config.HIDDEN_DIMS,
        dropout=config.DROPOUT_RATE,
    )
    model.cat_start_idx = data.cat_start_idx
    model.load_state_dict(torch.load(ckpt, map_location=DEVICE))
    model.eval()
    return model


def make_predict_fn(model):
    """
    Wrapper that takes a numpy array and returns survival probabilities.
    SHAP explainers need a plain numpy → numpy function.
    """
    @torch.no_grad()
    def predict(X: np.ndarray) -> np.ndarray:
        t = torch.tensor(X, dtype=torch.float32).to(DEVICE)
        return torch.sigmoid(model(t)).cpu().numpy()
    return predict


def compute_shap_values(model, data, n_background=100, n_explain=200):
    """
    Compute SHAP values using KernelExplainer on a random subset.

    n_background: number of background samples for the SHAP baseline
                  (k-means summary of training set — faster than full set)
    n_explain:    number of test samples to explain
                  (full test set is fine but slower)
    """
    rng = np.random.default_rng(config.SEED)

    # Background: summarise training distribution with k-means
    background_idx = rng.choice(len(data.X_train), size=n_background, replace=False)
    background = data.X_train[background_idx]

    # Samples to explain
    explain_idx = rng.choice(len(data.X_test), size=min(n_explain, len(data.X_test)), replace=False)
    X_explain = data.X_test[explain_idx]
    y_explain = data.y_test[explain_idx]

    predict_fn = make_predict_fn(model)

    print(f"\nComputing SHAP values for {len(X_explain)} test samples ...")
    print(f"  Background set: {len(background)} samples (k-means summary)")
    print("  This may take a few minutes with KernelExplainer ...")

    # KernelExplainer works with any black-box function
    explainer = shap.KernelExplainer(predict_fn, shap.kmeans(background, 20))
    shap_values = explainer.shap_values(X_explain, nsamples=100, silent=True)

    return shap_values, X_explain, y_explain, explain_idx


def plot_shap(shap_values, X_explain, feature_names, cat_start_idx, ccsr_cols):
    config.FIGURES_DIR.mkdir(exist_ok=True)

    # Use human-readable feature names (strip prefixes added by encode_features)
    def clean(name):
        return (name
                .replace("_log1p_scaled", "")
                .replace("_scaled", "")
                .replace("ccsr_", "CCSR:")
                .replace("cat_", ""))

    clean_names = [clean(n) for n in feature_names]

    # ── Summary (beeswarm) plot ────────────────────────────────────────────────
    # Each dot is one patient. Position = SHAP value (impact on prediction).
    # Color = feature value (red = high, blue = low).
    # This shows both direction and magnitude of each feature's effect.
    plt.figure(figsize=(10, 8))
    shap.summary_plot(
        shap_values, X_explain,
        feature_names=clean_names,
        show=False, max_display=20,
    )
    plt.title("SHAP Summary — MLP Survival Predictions")
    plt.tight_layout()
    plt.savefig(config.FIGURES_DIR / "shap_summary.png", dpi=150, bbox_inches="tight")
    print(f"SHAP summary plot → {config.FIGURES_DIR}/shap_summary.png")
    plt.close()

    # ── Bar chart (mean |SHAP|) ────────────────────────────────────────────────
    mean_abs = np.abs(shap_values).mean(axis=0)
    idx = np.argsort(mean_abs)[::-1][:20]
    plt.figure(figsize=(9, 6))
    plt.barh(range(len(idx)), mean_abs[idx][::-1], color="#2E75B6")
    plt.yticks(range(len(idx)), [clean_names[i] for i in idx][::-1])
    plt.xlabel("Mean |SHAP value|")
    plt.title("Top 20 Features by Mean Absolute SHAP Value")
    plt.tight_layout()
    plt.savefig(config.FIGURES_DIR / "shap_bar.png", dpi=150, bbox_inches="tight")
    print(f"SHAP bar chart → {config.FIGURES_DIR}/shap_bar.png")
    plt.close()

    # ── Dependence plot for caregiver_fe_rate ──────────────────────────────────
    # Shows how SHAP value for caregiver_fe_rate varies with its actual value.
    # Should show a negative slope — consistent with the R asymptotic exponential.
    fe_idx = next(
        (i for i, n in enumerate(clean_names) if "caregiver_fe_rate" in n), None
    )
    if fe_idx is not None:
        plt.figure(figsize=(7, 5))
        shap.dependence_plot(
            fe_idx, shap_values, X_explain,
            feature_names=clean_names,
            show=False,
        )
        plt.title("SHAP Dependence: caregiver_fe_rate")
        plt.tight_layout()
        plt.savefig(config.FIGURES_DIR / "shap_caregiver_fe.png",
                    dpi=150, bbox_inches="tight")
        print(f"SHAP dependence plot → {config.FIGURES_DIR}/shap_caregiver_fe.png")
        plt.close()


def main():
    print("=" * 60)
    print("Script 04 — SHAP Interpretability")
    print("=" * 60)

    data = prepare_data()
    model = load_model(data)

    shap_values, X_explain, y_explain, _ = compute_shap_values(
        model, data, n_background=100, n_explain=200
    )

    plot_shap(
        shap_values, X_explain,
        data.feature_names, data.cat_start_idx, data.ccsr_cols,
    )

    # Print top features by mean |SHAP|
    mean_abs = np.abs(shap_values).mean(axis=0)
    idx = np.argsort(mean_abs)[::-1][:10]
    print("\nTop 10 features by mean |SHAP value|:")
    for rank, i in enumerate(idx):
        name = data.feature_names[i] if i < len(data.feature_names) else f"col_{i}"
        print(f"  {rank+1:2d}. {name:<45} {mean_abs[i]:.4f}")

    return shap_values, data


if __name__ == "__main__":
    main()
