"""
Script 03 — Train MLP with Entity Embeddings (PyTorch)

Key concepts illustrated:
  - Entity embeddings for categorical variables
  - Mini-batch gradient descent with Adam optimizer
  - Learning rate scheduling (ReduceLROnPlateau)
  - Early stopping to prevent overfitting
  - Weighted loss for class imbalance
  - Model checkpointing (save best val-loss weights)
  - AUC-ROC, AUC-PR, Brier score evaluation

Outputs:
  - models/mlp_best.pt          — best checkpoint
  - models/mlp_history.npz      — training history arrays
  - figures/mlp_training_curves.png
"""

import os
import sys
import time

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import torch
import torch.nn as nn
from torch.utils.data import DataLoader, TensorDataset
from sklearn.metrics import roc_auc_score, average_precision_score, brier_score_loss

sys.path.insert(0, os.path.dirname(__file__))
import config
from data_utils import prepare_data
from model import MLPWithEmbeddings

torch.manual_seed(config.SEED)
np.random.seed(config.SEED)
DEVICE = torch.device("cuda" if torch.cuda.is_available() else "cpu")


def make_loader(X, y, shuffle, batch_size=config.BATCH_SIZE):
    ds = TensorDataset(torch.tensor(X, dtype=torch.float32),
                       torch.tensor(y, dtype=torch.float32))
    return DataLoader(ds, batch_size=batch_size, shuffle=shuffle, num_workers=0)


def train_epoch(model, loader, optimizer, criterion):
    model.train()
    total, n = 0.0, 0
    for xb, yb in loader:
        xb, yb = xb.to(DEVICE), yb.to(DEVICE)
        optimizer.zero_grad()
        loss = criterion(model(xb), yb)
        loss.backward()
        nn.utils.clip_grad_norm_(model.parameters(), 1.0)
        optimizer.step()
        total += loss.item() * len(yb); n += len(yb)
    return total / n


@torch.no_grad()
def eval_epoch(model, loader, criterion):
    model.eval()
    total, n = 0.0, 0
    logits_all, labels_all = [], []
    for xb, yb in loader:
        xb, yb = xb.to(DEVICE), yb.to(DEVICE)
        out = model(xb)
        total += criterion(out, yb).item() * len(yb); n += len(yb)
        logits_all.append(out.cpu()); labels_all.append(yb.cpu())
    probs  = torch.sigmoid(torch.cat(logits_all)).numpy()
    labels = torch.cat(labels_all).numpy()
    return total / n, roc_auc_score(labels, probs)


def train(data):
    config.MODELS_DIR.mkdir(exist_ok=True)
    config.FIGURES_DIR.mkdir(exist_ok=True)

    train_loader = make_loader(data.X_train, data.y_train, shuffle=True)
    val_loader   = make_loader(data.X_val,   data.y_val,   shuffle=False)

    pos_rate   = data.y_train.mean()
    pos_weight = torch.tensor([(1 - pos_rate) / pos_rate], dtype=torch.float32).to(DEVICE)
    criterion  = nn.BCEWithLogitsLoss(pos_weight=pos_weight)
    print(f"\npos_weight: {pos_weight.item():.3f}  (mortality rate: {pos_rate:.3f})")

    model = MLPWithEmbeddings(
        continuous_dim=data.cat_start_idx,
        cat_dims=data.cat_dims,
        cat_emb_dims=data.cat_emb_dims,
        hidden_dims=config.HIDDEN_DIMS,
        dropout=config.DROPOUT_RATE,
    )
    model.cat_start_idx = data.cat_start_idx
    model.to(DEVICE)
    print(f"Parameters: {sum(p.numel() for p in model.parameters() if p.requires_grad):,}")

    optimizer = torch.optim.Adam(model.parameters(),
                                  lr=config.LEARNING_RATE,
                                  weight_decay=config.WEIGHT_DECAY)
    scheduler = torch.optim.lr_scheduler.ReduceLROnPlateau(
        optimizer, mode="min", factor=0.5, patience=10)

    history = {"train_loss": [], "val_loss": [], "val_auc": [], "lr": []}
    best_val_loss, best_epoch, patience_counter = float("inf"), 0, 0
    ckpt = config.MODELS_DIR / "mlp_best.pt"
    t0 = time.time()

    print(f"\nTraining up to {config.MAX_EPOCHS} epochs (early stopping patience={config.PATIENCE})")
    for epoch in range(1, config.MAX_EPOCHS + 1):
        tl = train_epoch(model, train_loader, optimizer, criterion)
        vl, va = eval_epoch(model, val_loader, criterion)
        lr = optimizer.param_groups[0]["lr"]
        history["train_loss"].append(tl)
        history["val_loss"].append(vl)
        history["val_auc"].append(va)
        history["lr"].append(lr)
        scheduler.step(vl)

        if vl < best_val_loss:
            best_val_loss, best_epoch, patience_counter = vl, epoch, 0
            torch.save(model.state_dict(), ckpt)
        else:
            patience_counter += 1

        if epoch % 20 == 0 or epoch <= 5:
            print(f"  ep {epoch:4d} | train={tl:.4f} | val={vl:.4f} | auc={va:.4f} "
                  f"| lr={lr:.1e} | {time.time()-t0:.0f}s")

        if patience_counter >= config.PATIENCE:
            print(f"\nEarly stop at epoch {epoch}  (best epoch {best_epoch})")
            break

    model.load_state_dict(torch.load(ckpt, map_location=DEVICE))
    np.savez(config.MODELS_DIR / "mlp_history.npz",
             **{k: np.array(v) for k, v in history.items()})

    # Training curves
    fig, axes = plt.subplots(1, 3, figsize=(14, 4))
    axes[0].plot(history["train_loss"], label="train")
    axes[0].plot(history["val_loss"],   label="val")
    axes[0].axvline(best_epoch - 1, color="red", linestyle="--", alpha=0.6,
                    label=f"best ep {best_epoch}")
    axes[0].set(title="Loss", xlabel="Epoch", ylabel="BCE")
    axes[0].legend()
    axes[1].plot(history["val_auc"])
    axes[1].axvline(best_epoch - 1, color="red", linestyle="--", alpha=0.6)
    axes[1].set(title="Val AUC-ROC", xlabel="Epoch")
    axes[2].semilogy(history["lr"])
    axes[2].set(title="Learning Rate", xlabel="Epoch")
    plt.tight_layout()
    plt.savefig(config.FIGURES_DIR / "mlp_training_curves.png", dpi=150, bbox_inches="tight")
    print(f"Training curves → {config.FIGURES_DIR}/mlp_training_curves.png")

    return model, history


@torch.no_grad()
def evaluate_test(model, data):
    model.eval()
    probs = torch.sigmoid(
        model(torch.tensor(data.X_test, dtype=torch.float32).to(DEVICE))
    ).cpu().numpy()
    y = data.y_test
    auc   = roc_auc_score(y, probs)
    ap    = average_precision_score(y, probs)
    brier = brier_score_loss(y, probs)
    print(f"\n── MLP test-set results ───────────────────────────────")
    print(f"  AUC-ROC : {auc:.4f}")
    print(f"  AUC-PR  : {ap:.4f}")
    print(f"  Brier   : {brier:.4f}")
    return {"auc": auc, "ap": ap, "brier": brier, "probs": probs}


def main():
    print("=" * 60)
    print("Script 03 — MLP with Entity Embeddings")
    print("=" * 60)
    data = prepare_data()
    model, history = train(data)
    results = evaluate_test(model, data)
    return model, data, results


if __name__ == "__main__":
    main()
