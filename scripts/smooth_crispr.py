#!/usr/bin/env python3
"""Smooth the raw CRISPR matrix via taxonomy kernels.

Snakemake script. Migrated from capellini/stages/network.py build_smooth_crispr_one.
"""

from __future__ import annotations

import argparse
import json
import logging
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

import pandas as pd

from capellini_workflow.io import read_table, write_df
from capellini_workflow.network_utils import build_smoothed_crispr_for_study

logger = logging.getLogger(__name__)


def build_smooth_crispr(
    v_processed: str,
    b_processed: str,
    crispr_raw: str,
    tax_bac_for_smoothing: str,
    tax_vir: str,
    output_dir: str,
    bacterial_ranks: list | None = None,
    viral_ranks: list | None = None,
    bacterial_weights: list | None = None,
    viral_weights: list | None = None,
    alpha: float = 0.95,
    transpose_after_load: bool = True,
) -> dict[str, Path]:
    """Smooth the raw CRISPR matrix via taxonomy kernels; persist all artefacts."""
    if bacterial_ranks is None:
        bacterial_ranks = ["Phylum", "Class", "Order", "Family", "progenomes_taxid_genus"]
    if viral_ranks is None:
        viral_ranks = ["lev8", "lev7", "lev6", "lev5", "lev4", "lev3", "lev2", "lev1", "lev0"]
    if bacterial_weights is None:
        bacterial_weights = [1, 2, 3, 6, 8]
    if viral_weights is None:
        viral_weights = [1, 1, 2, 3, 4, 6, 8, 10, 12]

    out_dir = Path(output_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    V = read_table(v_processed)
    B = read_table(b_processed)

    art = build_smoothed_crispr_for_study(
        raw_crispr_path=crispr_raw,
        bacteria_features=B.columns,
        virus_features=V.columns,
        tax_bac_path=tax_bac_for_smoothing,
        tax_vir_path=tax_vir,
        bacterial_ranks=bacterial_ranks,
        viral_ranks=viral_ranks,
        bacterial_weights=bacterial_weights,
        viral_weights=viral_weights,
        alpha=alpha,
        transpose_after_load=transpose_after_load,
    )

    paths: dict[str, Path] = {}
    for name in ("crispr_binary", "crispr_smooth", "K_bac", "K_vir"):
        p = out_dir / f"{name}.csv.gz"
        write_df(art[name], p)
        paths[name] = p
    p_vb = out_dir / "crispr_smooth_vir_bac.csv.gz"
    write_df(art["crispr_smooth"].T, p_vb)
    paths["crispr_smooth_vir_bac"] = p_vb
    logger.info("smoothed CRISPR written to %s", out_dir)
    return paths


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Smooth CRISPR matrix")
    parser.add_argument("--v-processed", required=True)
    parser.add_argument("--b-processed", required=True)
    parser.add_argument("--crispr-raw", required=True)
    parser.add_argument("--tax-bac-for-smoothing", required=True)
    parser.add_argument("--tax-vir", required=True)
    parser.add_argument("--output-dir", required=True)
    parser.add_argument("--bacterial-ranks", type=str, default=None,
                        help="JSON list of bacterial taxonomy ranks")
    parser.add_argument("--viral-ranks", type=str, default=None,
                        help="JSON list of viral taxonomy ranks")
    parser.add_argument("--bacterial-weights", type=str, default=None,
                        help="JSON list of bacterial weights")
    parser.add_argument("--viral-weights", type=str, default=None,
                        help="JSON list of viral weights")
    parser.add_argument("--alpha", type=float, default=0.95)
    parser.add_argument("--transpose-after-load", action="store_true", default=True)
    args = parser.parse_args()

    logging.basicConfig(level=logging.INFO)
    build_smooth_crispr(
        v_processed=args.v_processed,
        b_processed=args.b_processed,
        crispr_raw=args.crispr_raw,
        tax_bac_for_smoothing=args.tax_bac_for_smoothing,
        tax_vir=args.tax_vir,
        output_dir=args.output_dir,
        bacterial_ranks=json.loads(args.bacterial_ranks) if args.bacterial_ranks else None,
        viral_ranks=json.loads(args.viral_ranks) if args.viral_ranks else None,
        bacterial_weights=json.loads(args.bacterial_weights) if args.bacterial_weights else None,
        viral_weights=json.loads(args.viral_weights) if args.viral_weights else None,
        alpha=args.alpha,
        transpose_after_load=args.transpose_after_load,
    )
