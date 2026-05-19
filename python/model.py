"""
PyTorch model architectures for ICU extubation survival prediction.

Two models:
  MLPWithEmbeddings  — the main model: entity embeddings for categorical
                       features + dense layers for continuous/CCSR features.
  TabNetWrapper      — thin wrapper around pytorch-tabnet for comparison.

Entity embedding approach:
  Guo & Berkhahn (2016) "Entity Embeddings of Categorical Variables"
  https://arxiv.org/abs/1604.06737

  Each categorical variable gets its own Embedding layer. The embedding
  vectors are learned end-to-end with the rest of the network. This is
  the deep learning analogue of a random intercept in mixed-effects models:
  first_careunit embeddings ≈ ICU-unit random effects in glmmTMB.
"""

import torch
import torch.nn as nn
import torch.nn.functional as F
from torch import Tensor


# ── Helper ─────────────────────────────────────────────────────────────────────

class ResidualBlock(nn.Module):
    """
    Two-layer residual block with BatchNorm and Dropout.

    The skip connection (identity shortcut) helps gradient flow in deeper
    networks. For tabular data with ~4K rows this is optional, but it's
    good practice to know.
    """
    def __init__(self, dim: int, dropout: float):
        super().__init__()
        self.block = nn.Sequential(
            nn.Linear(dim, dim),
            nn.BatchNorm1d(dim),
            nn.GELU(),
            nn.Dropout(dropout),
            nn.Linear(dim, dim),
            nn.BatchNorm1d(dim),
        )
        self.act = nn.GELU()
        self.dropout = nn.Dropout(dropout)

    def forward(self, x: Tensor) -> Tensor:
        return self.act(self.dropout(self.block(x)) + x)


# ── Main model ─────────────────────────────────────────────────────────────────

class MLPWithEmbeddings(nn.Module):
    """
    Feedforward network with entity embeddings for categorical variables.

    Input structure (from data_utils.encode_features):
      x[:, :cat_start_idx]  — continuous + CCSR features (float32)
      x[:, cat_start_idx:]  — label-encoded categoricals (int, cast inside forward)

    Architecture:
      Embedding → flatten → concat with continuous → hidden layers → output

    Args:
        continuous_dim   : number of continuous + CCSR features
        cat_dims         : list of n_categories per categorical feature
        cat_emb_dims     : list of embedding dimensions per categorical feature
        hidden_dims      : list of hidden layer sizes, e.g. [256, 128, 64]
        dropout          : dropout probability (same across all layers)
        use_residual     : if True, use ResidualBlocks instead of plain Linear
    """
    def __init__(
        self,
        continuous_dim: int,
        cat_dims: list[int],
        cat_emb_dims: list[int],
        hidden_dims: list[int] = (256, 128, 64),
        dropout: float = 0.30,
        use_residual: bool = True,
    ):
        super().__init__()

        self.cat_start_idx = continuous_dim  # set externally after construction
        self.n_cats = len(cat_dims)

        # Entity embedding layers — one per categorical feature
        self.embeddings = nn.ModuleList([
            nn.Embedding(n_cat + 1, emb_dim)  # +1 for unknown/padding category
            for n_cat, emb_dim in zip(cat_dims, cat_emb_dims)
        ])

        total_emb_dim = sum(cat_emb_dims)
        input_dim = continuous_dim + total_emb_dim

        # Input projection → first hidden layer
        layers = [
            nn.Linear(input_dim, hidden_dims[0]),
            nn.BatchNorm1d(hidden_dims[0]),
            nn.GELU(),
            nn.Dropout(dropout),
        ]

        # Hidden layers — either plain or residual
        for i in range(len(hidden_dims) - 1):
            in_dim, out_dim = hidden_dims[i], hidden_dims[i + 1]
            if use_residual and in_dim == out_dim:
                layers.append(ResidualBlock(in_dim, dropout))
            else:
                layers += [
                    nn.Linear(in_dim, out_dim),
                    nn.BatchNorm1d(out_dim),
                    nn.GELU(),
                    nn.Dropout(dropout),
                ]

        self.network = nn.Sequential(*layers)
        self.head = nn.Linear(hidden_dims[-1], 1)  # single logit output

        self._init_weights()

    def _init_weights(self):
        """Kaiming initialization for Linear layers; uniform for embeddings."""
        for module in self.modules():
            if isinstance(module, nn.Linear):
                nn.init.kaiming_normal_(module.weight, nonlinearity="relu")
                if module.bias is not None:
                    nn.init.zeros_(module.bias)
            elif isinstance(module, nn.Embedding):
                nn.init.uniform_(module.weight, -0.05, 0.05)

    def forward(self, x: Tensor) -> Tensor:
        # Split continuous and categorical parts
        x_cont = x[:, :self.cat_start_idx]                  # (B, continuous_dim)
        x_cat  = x[:, self.cat_start_idx:].long()           # (B, n_cats)

        # Look up embeddings and flatten
        emb_parts = [
            self.embeddings[i](x_cat[:, i])
            for i in range(self.n_cats)
        ]
        x_emb = torch.cat(emb_parts, dim=1)                 # (B, total_emb_dim)

        # Concatenate and pass through MLP
        x_combined = torch.cat([x_cont, x_emb], dim=1)      # (B, input_dim)
        h = self.network(x_combined)
        return self.head(h).squeeze(1)                       # (B,) — raw logits

    @torch.no_grad()
    def predict_proba(self, x: Tensor) -> Tensor:
        """Returns survival probabilities in [0, 1]."""
        self.eval()
        return torch.sigmoid(self.forward(x))


# ── TabNet wrapper ─────────────────────────────────────────────────────────────

class TabNetSurvivalWrapper:
    """
    Thin wrapper around pytorch_tabnet.TabNetClassifier for comparison.

    TabNet uses sequential attention to select features at each step,
    providing built-in interpretability (feature importance masks).
    It treats all inputs as a flat continuous matrix, so categoricals
    are label-encoded (not embedded), which is less expressive but
    requires no special architecture.

    Usage:
        wrapper = TabNetSurvivalWrapper()
        wrapper.fit(X_train, y_train, X_val, y_val)
        preds   = wrapper.predict_proba(X_test)[:, 1]
    """
    def __init__(self, cat_idxs: list[int] = None, cat_dims: list[int] = None):
        from pytorch_tabnet.tab_model import TabNetClassifier
        self.model = TabNetClassifier(
            cat_idxs=cat_idxs or [],
            cat_dims=cat_dims or [],
            cat_emb_dim=3,
            n_d=16, n_a=16,       # width of decision step output
            n_steps=5,            # number of sequential attention steps
            gamma=1.5,            # coefficient for feature reusage
            n_independent=2,
            n_shared=2,
            momentum=0.02,
            mask_type="sparsemax",
            verbose=1,
            seed=237,
        )

    def fit(self, X_train, y_train, X_val, y_val, max_epochs=200, patience=20):
        self.model.fit(
            X_train=X_train.astype(float),
            y_train=y_train,
            eval_set=[(X_val.astype(float), y_val)],
            eval_name=["val"],
            eval_metric=["auc"],
            max_epochs=max_epochs,
            patience=patience,
            batch_size=128,
            virtual_batch_size=32,
        )

    def predict_proba(self, X):
        return self.model.predict_proba(X.astype(float))

    @property
    def feature_importances_(self):
        return self.model.feature_importances_
