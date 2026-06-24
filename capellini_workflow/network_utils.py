"""Network-level utilities: residual message passing, CRISPR smoothing, taxonomy kernels, abundance helpers.

This is the slim core of the CAPELLINI network stage. The math follows the
paper (W̃ = (1−α) W + α K_vir W K_bac, residual propagation Z* = Z + η (cross
- Z)). Helper plumbing (sample alignment, orientation auto-detection, CRISPR
binarisation, taxonomy clean-up) lives below the math.
"""

from __future__ import annotations

import numpy as np
import pandas as pd

from capellini_workflow.transforms import clr, row_normalize
from capellini_workflow.taxonomy import (
    apply_custom_renames,
    clean_df_ids,
    clean_index_ids,
    parse_bool_series,
    rename_clostridium_sensu_stricto,
    sanitize_index,
)


# ── Validation helpers (kept verbatim) ────────────────────────────────────────

def _validate_inputs(V: np.ndarray, B: np.ndarray, W: np.ndarray, lam: float, n_steps: int) -> None:
    """Validate shapes and parameter ranges for message-passing functions."""
    n, q = V.shape
    n_b, p = B.shape
    if n != n_b:
        raise ValueError(f"V and B must have the same number of samples, got {n} and {n_b}.")
    if W.shape != (q, p):
        raise ValueError(f"W_vh_smooth must have shape ({q}, {p}), got {W.shape}.")
    if lam < 0:
        raise ValueError(f"lam must be nonnegative, got {lam}.")
    if n_steps < 1:
        raise ValueError(f"n_steps must be at least 1, got {n_steps}.")


def _build_common_inputs(
    V_df: pd.DataFrame,
    B_df: pd.DataFrame,
    W_vh_smooth_df: pd.DataFrame,
    pseudocount: float,
    eps: float,
) -> tuple:
    """Align samples, CLR-transform V and B, subset W to the common features.

    Returns:
        (V0, B0, V_clr_df, B_clr_df, X_clr_df, W_sub) — raw & CLR DataFrames + the
        W block restricted to overlapping virus rows / bacteria columns.
    """
    common_samples = V_df.index.intersection(B_df.index)
    if len(common_samples) == 0:
        raise ValueError("No overlapping sample IDs between V_df and B_df.")
    V0 = V_df.loc[common_samples].copy()
    B0 = B_df.loc[common_samples].copy()

    V_clr_df = pd.DataFrame(clr(V0, pseudocount=pseudocount, eps=eps), index=V0.index, columns=V0.columns)
    B_clr_df = pd.DataFrame(clr(B0, pseudocount=pseudocount, eps=eps), index=B0.index, columns=B0.columns)

    viruses = V_clr_df.columns.intersection(W_vh_smooth_df.index)
    bacteria = B_clr_df.columns.intersection(W_vh_smooth_df.columns)
    if len(viruses) == 0:
        raise ValueError("No overlapping viruses between V_clr_df.columns and W_vh_smooth_df.index.")
    if len(bacteria) == 0:
        raise ValueError("No overlapping bacteria between B_clr_df.columns and W_vh_smooth_df.columns.")

    V_clr_df = V_clr_df.loc[:, viruses].copy()
    B_clr_df = B_clr_df.loc[:, bacteria].copy()
    X_clr_df = pd.concat([V_clr_df, B_clr_df], axis=1)
    W_sub = W_vh_smooth_df.loc[viruses, bacteria].copy()
    return V0, B0, V_clr_df, B_clr_df, X_clr_df, W_sub


def orient_W_viruses_by_bacteria(
    W_df: pd.DataFrame,
    V_df: pd.DataFrame,
    B_df: pd.DataFrame,
    verbose: bool = True,
) -> pd.DataFrame:
    """Detect and correct the orientation of W to be viruses x bacteria."""
    n_vir_row = len(set(W_df.index.astype(str)) & set(map(str, V_df.columns)))
    n_bac_col = len(set(W_df.columns.astype(str)) & set(map(str, B_df.columns)))
    n_bac_row = len(set(W_df.index.astype(str)) & set(map(str, B_df.columns)))
    n_vir_col = len(set(W_df.columns.astype(str)) & set(map(str, V_df.columns)))
    if verbose:
        print("W rows overlapping viruses:", n_vir_row, "/", V_df.shape[1])
        print("W cols overlapping bacteria:", n_bac_col, "/", B_df.shape[1])
        print("W rows overlapping bacteria:", n_bac_row, "/", B_df.shape[1])
        print("W cols overlapping viruses:", n_vir_col, "/", V_df.shape[1])
    if n_bac_row > n_vir_row and n_vir_col > n_bac_col:
        print("Detected bacteria x viruses orientation; transposing W.")
        W_df = W_df.T
    return W_df


# ── Message passing ──────────────────────────────────────────────────────────

def _residual_message_passing(
    V: np.ndarray,
    B: np.ndarray,
    W: np.ndarray,
    lam: float,
    n_steps: int,
) -> tuple[np.ndarray, np.ndarray]:
    """Residual additive cross-domain message passing (paper Eq. 4-5).

    Z*_v = Z_v + η (Z_b P_h − Z_v),   Z*_b = Z_b + η (Z_v P_v − Z_b)
    """
    V = np.asarray(V, dtype=float)
    B = np.asarray(B, dtype=float)
    W = np.asarray(W, dtype=float)
    _validate_inputs(V, B, W, lam=lam, n_steps=n_steps)

    P_vh = row_normalize(W)        # viruses → bacteria  (q × p)
    P_hv = row_normalize(W.T)      # bacteria → viruses  (p × q)

    Vt, Bt = V.copy(), B.copy()
    for _ in range(n_steps):
        Rv = Vt - Vt.mean(axis=0, keepdims=True)
        Rb = Bt - Bt.mean(axis=0, keepdims=True)
        Vt = Vt + lam * (Rb @ P_hv)
        Bt = Bt + lam * (Rv @ P_vh)
    return Vt, Bt


def build_xstar_from_smoothed_crispr(
    V_df: pd.DataFrame,
    B_df: pd.DataFrame,
    W_vh_smooth_df: pd.DataFrame,
    *,
    pseudocount: float = 1e-6,
    lam: float = 0.5,
    n_steps: int = 1,
    eps: float = 1e-12,
) -> dict[str, pd.DataFrame]:
    """Residual X* pipeline: align samples, CLR, propagate via W̃.

    Returns a dict with V_clr, B_clr, X_clr, V_star, B_star, X_star,
    W_smooth_aligned.
    """
    V0, B0, V_clr_df, B_clr_df, X_clr_df, W_sub = _build_common_inputs(
        V_df, B_df, W_vh_smooth_df, pseudocount=pseudocount, eps=eps
    )
    V_star, B_star = _residual_message_passing(
        V_clr_df.to_numpy(), B_clr_df.to_numpy(), W_sub.to_numpy(),
        lam=lam, n_steps=n_steps,
    )
    V_star_df = pd.DataFrame(V_star, index=V_clr_df.index, columns=V_clr_df.columns)
    B_star_df = pd.DataFrame(B_star, index=B_clr_df.index, columns=B_clr_df.columns)
    X_star_df = pd.concat([V_star_df, B_star_df], axis=1)
    return {
        "V_clr": V_clr_df,
        "B_clr": B_clr_df,
        "X_clr": X_clr_df,
        "V_star": V_star_df,
        "B_star": B_star_df,
        "X_star": X_star_df,
        "W_smooth_aligned": W_sub,
    }


# ── Taxonomy kernel + CRISPR smoothing ────────────────────────────────────────

def build_taxonomy_kernel(
    ids,
    tax_df: pd.DataFrame,
    ranks,
    weights=None,
    fill_value: str = "",
    normalize_rows: bool = True,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Build a similarity kernel where K[i,j] = Σ w_k · I(rank_k[i] == rank_k[j]).

    Args:
        ids: Feature IDs to include.
        tax_df: Taxonomy table (rows = IDs, columns include ``ranks``).
        ranks: Ordered rank columns (deepest last).
        weights: Optional per-rank weights; default 1..len(ranks). Normalized
            to sum to 1.
        fill_value: String treated as missing — pairs sharing this value at
            a given rank are NOT counted as a match.
        normalize_rows: Row-normalize the resulting kernel.

    Returns:
        (K, aligned_tax) — kernel as DataFrame indexed by the overlap of
        ``ids`` and ``tax_df.index``, plus the cleaned taxonomy slice.
    """
    ids = pd.Index(ids)
    ranks = list(ranks)
    missing = [r for r in ranks if r not in tax_df.columns]
    if missing:
        raise ValueError(f"Missing taxonomy columns: {missing}")

    common = ids.intersection(tax_df.index)
    if len(common) == 0:
        raise ValueError("No overlap between ids and taxonomy index")

    weights = (np.arange(1, len(ranks) + 1, dtype=float) if weights is None
               else np.asarray(weights, dtype=float))
    if len(weights) != len(ranks):
        raise ValueError("weights must have same length as ranks")
    weights = weights / weights.sum()

    tax = tax_df.loc[common, ranks].fillna(fill_value).astype(str)
    tax = tax.replace({"nan": fill_value, "None": fill_value, "NA": fill_value,
                       "N/A": fill_value, "unknown": fill_value, "unclassified": fill_value})
    for r in ranks:
        tax[r] = tax[r].str.replace(r"^[a-z]__+", "", regex=True)

    X = tax.to_numpy()
    n = X.shape[0]
    K = np.zeros((n, n), dtype=float)
    for k in range(len(ranks)):
        col = X[:, k]
        K += weights[k] * ((col[:, None] == col[None, :]) & (col[:, None] != fill_value)).astype(float)
    np.fill_diagonal(K, 1.0)
    if normalize_rows:
        rs = K.sum(axis=1, keepdims=True)
        rs[rs == 0] = 1.0
        K = K / rs

    return pd.DataFrame(K, index=common, columns=common), tax


def smooth_crispr_bac_vir(
    crispr_df: pd.DataFrame,
    K_bac: pd.DataFrame,
    K_vir: pd.DataFrame,
    alpha: float = 0.95,
) -> pd.DataFrame:
    """W̃ = (1 − α) W + α (K_bac · W · K_vir), restricted to common rows/cols."""
    bac = crispr_df.index.intersection(K_bac.index)
    vir = crispr_df.columns.intersection(K_vir.index)
    if len(bac) == 0 or len(vir) == 0:
        raise ValueError("No overlap between crispr_df and the taxonomy kernels.")
    W = crispr_df.loc[bac, vir].to_numpy(dtype=float)
    Kb = K_bac.loc[bac, bac].to_numpy()
    Kv = K_vir.loc[vir, vir].to_numpy()
    W_smooth = (1 - alpha) * W + alpha * (Kb @ W @ Kv)
    return pd.DataFrame(W_smooth, index=bac, columns=vir)


def build_binary_crispr_matrix(
    raw_crispr_path: str,
    bacteria_features,
    virus_features,
    transpose_after_load: bool = True,
) -> pd.DataFrame:
    """Load a raw CRISPR network CSV and build a binary bacteria × viruses matrix.

    The output is reindexed onto the requested ``bacteria_features`` /
    ``virus_features`` (missing rows/cols are zero) and the values are clipped
    to {0, 1}.
    """
    df = pd.read_csv(raw_crispr_path, index_col=0)
    if transpose_after_load:
        df = df.T
    df.index = df.index.astype(str)
    df.columns = df.columns.astype(str)

    bac_idx = pd.Index([str(b) for b in bacteria_features])
    vir_idx = pd.Index([str(v) for v in virus_features])

    # We start from a CRISPR matrix in (vir × bac) after the optional transpose
    # above; transpose once more to match the (bac × vir) output convention.
    src = df.T
    out = pd.DataFrame(0, index=bac_idx, columns=vir_idx, dtype=float)
    rows = bac_idx.intersection(src.index)
    cols = vir_idx.intersection(src.columns)
    if len(rows) and len(cols):
        out.loc[rows, cols] = src.loc[rows, cols].values
    out[out != 0] = 1
    return out


def get_hierarchies(df_b: pd.DataFrame, df_v: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Prepare bacteria/virus taxonomy frames for taxonomy-kernel building.

    Bacteria are reindexed by ``progenomes_taxid_genus`` (drop NaN, dedup).
    Viruses are reindexed by ``lev0`` (dedup).
    """
    b = df_b.copy()
    v = df_v.copy()
    v = v.set_index("lev0")
    b["progenomes_taxid_genus"] = b["progenomes_taxid_genus"].astype(str).str.split(".").str[0]
    b = b.set_index("progenomes_taxid_genus")
    b = b.loc[[i for i in b.index if "nan" not in str(i)]]
    b = b.loc[~b.index.duplicated(keep="first")]
    v = v.loc[~v.index.duplicated(keep="first")]
    b["progenomes_taxid_genus"] = b.index
    v["lev0"] = v.index
    return b, v


def build_smoothed_crispr_for_study(
    raw_crispr_path: str,
    bacteria_features,
    virus_features,
    tax_bac_path: str,
    tax_vir_path: str,
    bacterial_ranks,
    viral_ranks,
    bacterial_weights,
    viral_weights,
    alpha: float = 0.95,
    transpose_after_load: bool = True,
) -> dict[str, pd.DataFrame]:
    """Orchestrate the smoothed-CRISPR build for a single study.

    Returns the artefacts (CRISPR binary, K_bac, K_vir, smoothed W, aligned
    taxonomy frames). The caller is responsible for persisting them.
    """
    tax_bac = pd.read_csv(tax_bac_path, index_col=0)
    tax_vir = pd.read_csv(tax_vir_path, index_col=0)
    tax_bac_aligned, tax_vir_aligned = get_hierarchies(tax_bac, tax_vir)

    crispr_binary = build_binary_crispr_matrix(
        raw_crispr_path=raw_crispr_path,
        bacteria_features=bacteria_features,
        virus_features=virus_features,
        transpose_after_load=transpose_after_load,
    )
    K_bac, _ = build_taxonomy_kernel(crispr_binary.index, tax_bac_aligned,
                                     ranks=bacterial_ranks, weights=bacterial_weights)
    K_vir, _ = build_taxonomy_kernel(crispr_binary.columns, tax_vir_aligned,
                                     ranks=viral_ranks, weights=viral_weights)
    crispr_smooth = smooth_crispr_bac_vir(crispr_binary, K_bac, K_vir, alpha=alpha)

    return {
        "crispr_binary": crispr_binary,
        "K_bac": K_bac,
        "K_vir": K_vir,
        "crispr_smooth": crispr_smooth,
        "bac_tax_aligned": tax_bac_aligned,
        "vir_tax_aligned": tax_vir_aligned,
    }


def crispr_matrix_aggregate_viruses(
    df_crispr: pd.DataFrame,
    vir_tax: pd.DataFrame,
    *,
    bac_col: int = 0,
    vir_col: int = 1,
    vir_rank: str = "lev0",
) -> pd.DataFrame:
    """Aggregate a SpacePHARER predictions TSV into a (bac_taxid × vOTU) matrix.

    Bacterial spacer IDs of the form ``...>TAXID...`` are parsed to the
    leading taxid; viral contigs are mapped to ``vir_rank`` via ``vir_tax``.
    """
    def _parse_bac_taxid(s):
        if pd.isna(s):
            return None
        try:
            return int(str(s).split(">")[1].split(".")[0])
        except (IndexError, ValueError):
            return None

    df = df_crispr[[bac_col, vir_col]].copy()
    df["bac_taxid"] = df[bac_col].map(_parse_bac_taxid)
    df["contig"] = df[vir_col].astype(str).str.strip()
    df = df.dropna(subset=["bac_taxid", "contig"])
    df["bac_taxid"] = df["bac_taxid"].astype(int)

    M = pd.crosstab(df["bac_taxid"], df["contig"]).astype(int)

    if vir_rank not in vir_tax.columns:
        raise KeyError(f"vir_tax missing column {vir_rank!r}")
    vt = vir_tax[[vir_rank]].copy()
    vt.index = vt.index.astype(str).str.strip()
    vt = vt.loc[~vt.index.duplicated(keep="first")]
    labels = pd.Series(M.columns.astype(str), index=M.columns).map(vt[vir_rank])

    valid = labels.notna() & (labels.astype(str).str.strip() != "")
    M = M.loc[:, valid.values]
    labels = labels.loc[valid].astype(str).str.strip()
    out = M.T.groupby(labels).sum().T
    out.index = out.index.astype(str)
    return out


# ── Abundance preprocessing (kept verbatim per user instruction) ──────────────

def aggregate_otu_columns_by_rank_skip_nan(
    otu: pd.DataFrame,
    tax: pd.DataFrame,
    rank: str,
) -> pd.DataFrame:
    """Aggregate OTU/ASV columns to a taxonomy rank, skipping NaN labels."""
    common_asvs = otu.columns.intersection(tax.index)
    if len(common_asvs) == 0:
        raise ValueError("No common ASV/OTU IDs between abundance columns and taxonomy index.")
    otu2 = otu.loc[:, common_asvs]
    labels = tax.loc[common_asvs, rank]
    keep = labels.notna() & (labels.astype(str).str.strip() != "")
    otu2 = otu2.loc[:, keep.values]
    labels = labels.loc[keep].astype(str).str.strip()
    return otu2.T.groupby(labels).sum().T


def prepare_bacteria_genus_abundance(
    otu: pd.DataFrame,
    tax: pd.DataFrame,
    rank: str = "target_taxids",
) -> pd.DataFrame:
    """Aggregate bacteria OTUs to ``rank``, sanitize the resulting column index."""
    out = aggregate_otu_columns_by_rank_skip_nan(otu, tax, rank)
    out = rename_clostridium_sensu_stricto(out)
    out.columns = sanitize_index(out.columns)
    out = apply_custom_renames(out)
    return out


def remove_disease_columns_from_virus_abundance(V: pd.DataFrame) -> pd.DataFrame:
    """Drop phenotype/metadata columns that sometimes leak into viral abundance tables."""
    bad_cols = {
        "disease", "Disease", "disease_original", "disease_binary",
        "phenotype", "Phenotype", "label", "Label",
        "class", "Class", "group", "Group",
        "sample_id", "SampleID", "subject_id", "SubjectNo",
        "Reads", "reads",
    }
    drop = [c for c in V.columns if str(c) in bad_cols]
    if drop:
        print("Dropping non-feature columns from viral abundance:", drop)
        return V.drop(columns=drop)
    return V.copy()


def prevalence_filter_df(df: pd.DataFrame, prevalence: float = 0.10, verbose: bool = True) -> pd.DataFrame:
    """Keep features present in at least ``prevalence`` × n_samples samples."""
    n_samples = df.shape[0]
    min_samples = max(int(prevalence * n_samples), 1)
    keep = (df > 0).sum(axis=0) >= min_samples
    if verbose:
        print(f"prevalence={prevalence} -> min_samples={min_samples}/{n_samples}; "
              f"kept {int(keep.sum())}/{df.shape[1]} features")
    return df.loc[:, keep]


def align_abundance_from_metadata(
    virus_abundance: pd.DataFrame,
    bacteria_abundance: pd.DataFrame,
    metadata: pd.DataFrame,
    *,
    keep_col: str = "keep_for_analysis",
    virus_id_col: str = "virus_sample_id",
    bacteria_id_col: str = "bacteria_sample_id",
    final_index_col: str = "virus_sample_id",
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """Align V/B with the standardized metadata (filter, reorder, rename index)."""
    V_raw = clean_df_ids(virus_abundance)
    B_raw = clean_df_ids(bacteria_abundance)
    meta = metadata.copy()
    for col in [virus_id_col, bacteria_id_col]:
        if col not in meta.columns:
            raise ValueError(f"metadata is missing required column: {col}")
    if keep_col is not None and keep_col in meta.columns:
        keep = parse_bool_series(meta[keep_col])
        print(f"metadata filter {keep_col}: kept {int(keep.sum())} / {len(keep)} rows")
        meta = meta.loc[keep].copy()
    else:
        print(f"metadata filter {keep_col!r} not found; keeping all {meta.shape[0]} rows")
    if "analysis_order" in meta.columns:
        meta = meta.sort_values("analysis_order").copy()
    for col in [virus_id_col, bacteria_id_col, final_index_col, "sample_id", "subject_id"]:
        if col in meta.columns:
            meta[col] = clean_index_ids(meta[col])

    virus_ids = list(meta[virus_id_col].astype(str))
    bacteria_ids = list(meta[bacteria_id_col].astype(str))
    missing_v = sorted(set(virus_ids) - set(V_raw.index.astype(str)))
    missing_b = sorted(set(bacteria_ids) - set(B_raw.index.astype(str)))
    if missing_v:
        raise ValueError(f"{len(missing_v)} viral sample IDs missing in viral abundance. Examples: {missing_v[:10]}")
    if missing_b:
        raise ValueError(f"{len(missing_b)} bacterial sample IDs missing in bacterial abundance. Examples: {missing_b[:10]}")

    V = V_raw.loc[virus_ids].copy()
    B = B_raw.loc[bacteria_ids].copy()
    final_index = list(meta[final_index_col].astype(str)) if final_index_col in meta.columns else virus_ids
    V.index = final_index
    B.index = final_index
    meta = meta.reset_index(drop=True)
    return V, B, meta
