#!/usr/bin/env python3
"""Build common abundance matrices (V_processed, B_processed).

Snakemake script. Migrated from capellini/stages/network.py build_common_abundance_one.
"""

from __future__ import annotations

import argparse
import logging
import sys
from pathlib import Path

# Add parent so capellini_workflow is importable when run from the workflow root
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

import pandas as pd

from capellini_workflow.io import read_table, write_df
from capellini_workflow.taxonomy import clean_df_ids, load_bacteria_taxonomy
from capellini_workflow.network_utils import (
    align_abundance_from_metadata,
    prepare_bacteria_genus_abundance,
    prevalence_filter_df,
    remove_disease_columns_from_virus_abundance,
)

logger = logging.getLogger(__name__)


def build_common_abundance(
    virus_abundance_raw: str,
    bacteria_otu: str,
    bacteria_taxonomy: str,
    metadata_path: str,
    output_dir: str,
    bacteria_taxonomy_rank: str = "target_taxids",
    keep_column: str = "keep_for_analysis",
    prevalence: float = 0.10,
) -> dict[str, Path]:
    """Load raw inputs, align, filter, and write processed abundance matrices."""
    logger.info("common abundance: loading raw inputs")
    V_raw = read_table(virus_abundance_raw)
    V_raw = remove_disease_columns_from_virus_abundance(V_raw)
    V_raw = clean_df_ids(V_raw)

    B_otu = read_table(bacteria_otu)
    B_tax = load_bacteria_taxonomy(bacteria_taxonomy)
    B_genus = prepare_bacteria_genus_abundance(B_otu, B_tax, rank=bacteria_taxonomy_rank)

    meta = read_table(metadata_path, index_col=None)
    V, B, meta_aligned = align_abundance_from_metadata(
        virus_abundance=V_raw,
        bacteria_abundance=B_genus,
        metadata=meta,
        keep_col=keep_column,
    )
    V = prevalence_filter_df(V, prevalence=prevalence)
    B = prevalence_filter_df(B, prevalence=prevalence)

    out_dir = Path(output_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    paths = {
        "V": out_dir / "V_processed.csv",
        "B": out_dir / "B_processed.csv",
        "meta": out_dir / "metadata_aligned.csv",
    }
    write_df(V, paths["V"])
    write_df(B, paths["B"])
    write_df(meta_aligned, paths["meta"])
    logger.info("common abundance written to %s", out_dir)
    return paths


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Build common abundance matrices")
    parser.add_argument("--virus-abundance-raw", required=True)
    parser.add_argument("--bacteria-otu", required=True)
    parser.add_argument("--bacteria-taxonomy", required=True)
    parser.add_argument("--metadata", required=True)
    parser.add_argument("--output-dir", required=True)
    parser.add_argument("--bacteria-taxonomy-rank", default="target_taxids")
    parser.add_argument("--keep-column", default="keep_for_analysis")
    parser.add_argument("--prevalence", type=float, default=0.10)
    args = parser.parse_args()

    logging.basicConfig(level=logging.INFO)
    build_common_abundance(
        virus_abundance_raw=args.virus_abundance_raw,
        bacteria_otu=args.bacteria_otu,
        bacteria_taxonomy=args.bacteria_taxonomy,
        metadata_path=args.metadata,
        output_dir=args.output_dir,
        bacteria_taxonomy_rank=args.bacteria_taxonomy_rank,
        keep_column=args.keep_column,
        prevalence=args.prevalence,
    )
