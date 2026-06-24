#!/usr/bin/env python3
"""Extract GCA target list from the harmonized taxonomy table.

Helper script for the protein_clusters Snakemake rule.
"""

from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd


def extract_gca_list(
    silva_fixed_path: str,
    output_path: str,
    species_level: bool = False,
) -> None:
    """Read silva_fixed and write one GCA per line to output_path."""
    silva_fixed = pd.read_csv(silva_fixed_path, index_col=0)

    if species_level:
        gca_target_set = (
            set(silva_fixed["GCA_species"].dropna())
            | set(silva_fixed["GCA_family"].dropna())
        )
    else:
        gca_target_set = (
            set(silva_fixed["GCA_genus"].dropna())
            | set(silva_fixed["GCA_family"].dropna())
        )

    out = Path(output_path)
    out.parent.mkdir(parents=True, exist_ok=True)
    with open(out, "w") as f:
        for gca in sorted(gca_target_set):
            f.write(f"{gca}\n")

    print(f"Wrote {len(gca_target_set)} GCA targets to {out}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Extract GCA target list")
    parser.add_argument("--silva-fixed", required=True)
    parser.add_argument("--output", required=True)
    parser.add_argument("--species-level", type=str, default="false",
                        help="'true' or 'false'")
    args = parser.parse_args()

    species_level = args.species_level.lower() in ("true", "1", "yes")
    extract_gca_list(args.silva_fixed, args.output, species_level)
