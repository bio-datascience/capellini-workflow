#!/usr/bin/env python3
"""Schaefer-Strimmer shrinkage correlations on Z = [B_clr  V_clr].

Snakemake script. Migrated from capellini/stages/network.py build_shrinkage_one.
"""

from __future__ import annotations

import argparse
import logging
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

import pandas as pd

from capellini_workflow.io import read_table, write_df
from capellini_workflow.transforms import clr, schaefer_strimmer_corr

logger = logging.getLogger(__name__)


def build_shrinkage(
    v_processed: str,
    b_processed: str,
    output_dir: str,
    pseudocount: float = 1e-6,
) -> dict[str, Path]:
    """SS shrinkage correlation on Z = [B_clr  V_clr]."""
    out_dir = Path(output_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    V = read_table(v_processed)
    B = read_table(b_processed).loc[V.index]

    V_clr = pd.DataFrame(clr(V, pseudocount=pseudocount), index=V.index, columns=V.columns)
    B_clr = pd.DataFrame(clr(B, pseudocount=pseudocount), index=B.index, columns=B.columns)
    V_pref = V_clr.add_prefix("V__")
    B_pref = B_clr.add_prefix("B__")
    joint = pd.concat([B_pref, V_pref], axis=1)

    corr_df, info = schaefer_strimmer_corr(joint)
    bv = corr_df.loc[B_pref.columns, V_pref.columns]

    paths = {
        "corr_full": out_dir / "shrinkage_corr_full.csv.gz",
        "corr_bv": out_dir / "shrinkage_corr_BV.csv.gz",
    }
    write_df(corr_df, paths["corr_full"])
    write_df(bv, paths["corr_bv"])
    logger.info("shrinkage written (lambda=%s)", info)
    return paths


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Shrinkage correlations")
    parser.add_argument("--v-processed", required=True)
    parser.add_argument("--b-processed", required=True)
    parser.add_argument("--output-dir", required=True)
    parser.add_argument("--pseudocount", type=float, default=1e-6)
    args = parser.parse_args()

    logging.basicConfig(level=logging.INFO)
    build_shrinkage(
        v_processed=args.v_processed,
        b_processed=args.b_processed,
        output_dir=args.output_dir,
        pseudocount=args.pseudocount,
    )
