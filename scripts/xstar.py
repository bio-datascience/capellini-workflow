#!/usr/bin/env python3
"""Residual X* propagation + Schaefer-Strimmer correlations on Z*.

Snakemake script. Migrated from capellini/stages/network.py build_xstar_one.
"""

from __future__ import annotations

import argparse
import logging
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

import pandas as pd

from capellini_workflow.io import read_table, write_df
from capellini_workflow.transforms import schaefer_strimmer_corr
from capellini_workflow.network_utils import (
    build_xstar_from_smoothed_crispr,
    orient_W_viruses_by_bacteria,
)

logger = logging.getLogger(__name__)


def build_xstar(
    v_processed: str,
    b_processed: str,
    crispr_smooth_vir_bac: str,
    output_dir: str,
    pseudocount: float = 1e-6,
    lam: float = 0.5,
    n_steps: int = 1,
) -> dict[str, Path]:
    """Residual X* propagation + Schaefer-Strimmer correlations on Z*."""
    out_dir = Path(output_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    V = read_table(v_processed)
    B = read_table(b_processed)
    W = read_table(crispr_smooth_vir_bac)
    W = orient_W_viruses_by_bacteria(W, V, B)

    res = build_xstar_from_smoothed_crispr(
        V_df=V, B_df=B, W_vh_smooth_df=W,
        pseudocount=pseudocount, lam=lam, n_steps=n_steps,
    )

    # Shrinkage on Z* = [B_star  V_star]
    V_pref = res["V_star"].add_prefix("V__")
    B_pref = res["B_star"].add_prefix("B__")
    joint_star = pd.concat([B_pref, V_pref], axis=1)
    corr_star_full, info_star = schaefer_strimmer_corr(joint_star)
    corr_star_bv = corr_star_full.loc[B_pref.columns, V_pref.columns]

    paths = {
        "V_clr": out_dir / "V_clr.csv.gz",
        "B_clr": out_dir / "B_clr.csv.gz",
        "V_star": out_dir / "V_star.csv.gz",
        "B_star": out_dir / "B_star.csv.gz",
        "X_star": out_dir / "X_star.csv.gz",
        "shrinkage_xstar_full": out_dir / "shrinkage_xstar_corr_full.csv.gz",
        "shrinkage_xstar_BV": out_dir / "shrinkage_xstar_corr_BV.csv.gz",
    }
    write_df(res["V_clr"], paths["V_clr"])
    write_df(res["B_clr"], paths["B_clr"])
    write_df(res["V_star"], paths["V_star"])
    write_df(res["B_star"], paths["B_star"])
    write_df(res["X_star"], paths["X_star"])
    write_df(corr_star_full, paths["shrinkage_xstar_full"])
    write_df(corr_star_bv, paths["shrinkage_xstar_BV"])
    logger.info("X* + Z*-shrinkage written to %s (lambda=%s)", out_dir, info_star)
    return paths


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Build X* propagation")
    parser.add_argument("--v-processed", required=True)
    parser.add_argument("--b-processed", required=True)
    parser.add_argument("--crispr-smooth-vir-bac", required=True)
    parser.add_argument("--output-dir", required=True)
    parser.add_argument("--pseudocount", type=float, default=1e-6)
    parser.add_argument("--lam", type=float, default=0.5)
    parser.add_argument("--n-steps", type=int, default=1)
    args = parser.parse_args()

    logging.basicConfig(level=logging.INFO)
    build_xstar(
        v_processed=args.v_processed,
        b_processed=args.b_processed,
        crispr_smooth_vir_bac=args.crispr_smooth_vir_bac,
        output_dir=args.output_dir,
        pseudocount=args.pseudocount,
        lam=args.lam,
        n_steps=args.n_steps,
    )
