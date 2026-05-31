"""
Script 00 — Run Full Pipeline

Runs scripts 02 → 03 → 04 → 05 → 06 in order, mirroring 00_run_all.R.

Usage (from the python/ directory):
    python 00_run_all.py                        # full pipeline
    python 00_run_all.py --skip-shap            # skip SHAP (faster)
    python 00_run_all.py --skip-dml             # skip both DML scripts
    python 00_run_all.py --skip-shap --skip-dml # baselines + MLP only
"""

import argparse
import sys
import os

sys.path.insert(0, os.path.dirname(__file__))


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--skip-shap", action="store_true",
                        help="Skip SHAP interpretability (saves ~5-10 min)")
    parser.add_argument("--skip-dml", action="store_true",
                        help="Skip DML causal estimation (saves ~2-3 min)")
    args = parser.parse_args()

    print("\n" + "=" * 60)
    print("ICU Extubation Outcomes — Python DL Pipeline")
    print("=" * 60)

    # Shared data preparation (done once, reused across scripts)
    from data_utils import prepare_data
    data = prepare_data()

    # ── Script 02: Baselines ───────────────────────────────────────────────────
    print("\n\n" + "─" * 60)
    print("STEP 1/3 — Baseline models (logistic regression, gradient boosting)")
    print("─" * 60)
    from importlib import import_module
    baselines = import_module("02_baselines")
    baseline_results, baseline_models = baselines.main()

    # ── Script 03: MLP ────────────────────────────────────────────────────────
    print("\n\n" + "─" * 60)
    print("STEP 2/3 — MLP with entity embeddings (PyTorch)")
    print("─" * 60)
    mlp_mod = import_module("03_train_mlp")
    model, _, mlp_results = mlp_mod.main()

    # ── Script 04: SHAP ───────────────────────────────────────────────────────
    if not args.skip_shap:
        print("\n\n" + "─" * 60)
        print("STEP 3/3 — SHAP interpretability")
        print("─" * 60)
        shap_mod = import_module("04_shap_interpret")
        shap_mod.main()
    else:
        print("\nSkipping SHAP (--skip-shap flag set).")

    # ── Script 05: DML (global) ───────────────────────────────────────────────
    if not args.skip_dml:
        print("\n\n" + "─" * 60)
        print("STEP 4/5 — Double/Debiased ML (global, all units)")
        print("─" * 60)
        dml_mod = import_module("05_dml_caregiver_fe")
        dml_mod.main()

        # ── Script 06: DML (CVICU only) ───────────────────────────────────────
        print("\n\n" + "─" * 60)
        print("STEP 5/5 — Double/Debiased ML (CVICU only)")
        print("─" * 60)
        dml_cvicu_mod = import_module("06_dml_cvicu")
        dml_cvicu_mod.main()
    else:
        print("\nSkipping DML (--skip-dml flag set).")

    print("\n\n" + "=" * 60)
    print("Pipeline complete.")
    print(f"  Figures: python/figures/")
    print(f"  Model:   python/models/mlp_best.pt")
    print("=" * 60)


if __name__ == "__main__":
    main()
