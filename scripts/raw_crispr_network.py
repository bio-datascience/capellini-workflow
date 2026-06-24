#!/usr/bin/env python3
"""Aggregate SpacePHARER predictions into a (bac x vOTU) CRISPR matrix.

Snakemake script. Migrated from capellini/stages/network.py build_raw_crispr_one.
"""

from __future__ import annotations

import argparse
import logging
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

import pandas as pd

from capellini_workflow.io import read_table, write_df
from capellini_workflow.network_utils import crispr_matrix_aggregate_viruses

logger = logging.getLogger(__name__)


def build_raw_crispr(
    phage_host_predictions: str,
    tax_vir: str,
    output_dir: str,
    aggregate_viral_rank: str = "lev0",
) -> Path:
    """Aggregate the SpacePHARER predictions into a (bac x vOTU) CRISPR matrix."""
    out_dir = Path(output_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    df_pred = pd.read_csv(phage_host_predictions, sep="\t", header=None, comment="#")
    vir_tax = read_table(tax_vir)
    crispr = crispr_matrix_aggregate_viruses(df_pred, vir_tax, vir_rank=aggregate_viral_rank)

    out = out_dir / "crispr_net.csv.gz"
    write_df(crispr, out)
    logger.info("raw CRISPR matrix written: %s", out)
    return out


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Build raw CRISPR network")
    parser.add_argument("--phage-host-predictions", required=True)
    parser.add_argument("--tax-vir", required=True)
    parser.add_argument("--output-dir", required=True)
    parser.add_argument("--aggregate-viral-rank", default="lev0")
    args = parser.parse_args()

    logging.basicConfig(level=logging.INFO)
    build_raw_crispr(
        phage_host_predictions=args.phage_host_predictions,
        tax_vir=args.tax_vir,
        output_dir=args.output_dir,
        aggregate_viral_rank=args.aggregate_viral_rank,
    )
