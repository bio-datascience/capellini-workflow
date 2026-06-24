"""Numerical transformations: CLR, row normalization, shrinkage correlation."""

from __future__ import annotations

import numpy as np
import pandas as pd


def clr(X, pseudocount: float = 1e-6, eps: float = 1e-12):
    """Row-wise centered log-ratio transform.

    Args:
        X: Nonnegative array-like of shape (n_samples, n_features).
        pseudocount: Constant added before the log to handle zeros.
        eps: Numerical stability constant inside the log.

    Returns:
        CLR-transformed numpy array.
    """
    X = np.asarray(X, dtype=float)
    X = np.clip(X, 0.0, None) + pseudocount
    log_x = np.log(X + eps)
    return log_x - log_x.mean(axis=1, keepdims=True)


def row_normalize(W, eps: float = 1e-12):
    """Row-normalize a matrix so each row sums to 1.

    Args:
        W: Matrix-like.
        eps: Small constant to avoid division by zero.

    Returns:
        Row-normalized numpy array.
    """
    W = np.asarray(W, dtype=float)
    rs = W.sum(axis=1, keepdims=True)
    return W / (rs + eps)


def schaefer_strimmer_corr(X: pd.DataFrame) -> tuple[pd.DataFrame, dict]:
    """Schäfer-Strimmer shrinkage correlation estimator.

    Args:
        X: Samples x features DataFrame (at least 3 samples required).

    Returns:
        Tuple (shrinkage_correlation_DataFrame, diagnostics_dict) with keys
        ``n_samples``, ``n_features``, ``lambda_var``, ``lambda_corr``.
    """
    values = X.to_numpy(dtype=float)
    n, p = values.shape
    if n < 3:
        raise ValueError(f"Need at least 3 samples for shrinkage correlation, got {n}")

    # variance shrinkage parameter (computed for diagnostics)
    w = (values - values.mean(axis=0)) ** 2
    w_bar = np.mean(w, axis=0)
    var_unb = (n / (n - 1)) * w_bar
    var_s = (n / (n - 1) ** 3) * np.sum((w - w_bar) ** 2, axis=0)
    med_var = np.median(var_unb)
    denom_var = np.sum((var_unb - med_var) ** 2)
    lambda_var = 1.0 if denom_var <= 1e-15 else min(1.0, float(np.sum(var_s) / denom_var))

    # off-diagonal shrinkage parameter
    sd = np.std(values, axis=0, ddof=1)
    X_st = np.divide(values, sd, out=np.zeros_like(values), where=sd > 1e-12)
    X_c_st = X_st - X_st.mean(axis=0)
    w_st = X_c_st.T @ X_c_st
    w_st_sq = (X_c_st ** 2).T @ (X_c_st ** 2)
    w_bar_st = w_st / n
    var_s_st = (n / (n - 1) ** 3) * (w_st_sq - 2 * w_bar_st * w_st + n * w_bar_st ** 2)
    corr_unb_st = (n / (n - 1)) * w_bar_st
    denom_corr = np.sum(corr_unb_st ** 2) - np.sum(np.diag(corr_unb_st) ** 2)
    numer_corr = np.sum(var_s_st) - np.sum(np.diag(var_s_st))
    lambda_corr = 1.0 if denom_corr <= 1e-15 else min(1.0, float(numer_corr / denom_corr))

    corr_X = np.nan_to_num(np.corrcoef(values.T), nan=0.0, posinf=0.0, neginf=0.0)
    corr_shrink = (1.0 - lambda_corr) * corr_X
    np.fill_diagonal(corr_shrink, 1.0)

    out = pd.DataFrame(corr_shrink, index=X.columns, columns=X.columns)
    info = {
        "n_samples": int(n),
        "n_features": int(p),
        "lambda_var": float(lambda_var),
        "lambda_corr": float(lambda_corr),
    }
    return out, info
