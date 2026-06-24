"""CAPELLINI Snakemake workflow.

CRISPR-Abundance Phage-Evidence Linkage for Leveraging Interaction Network Inference.

Usage:
    snakemake --configfile config/ibd.yaml --cores 4
    snakemake --configfile config/ibd.yaml --dry-run
"""

import os
from pathlib import Path

# Directory containing the Snakefile (workflow root)
WORKFLOW_DIR = Path(workflow.basedir)
SCRIPTS_DIR = WORKFLOW_DIR / "scripts"

# ── Validate required config keys ────────────────────────────────────────────
for key in ("base", "download_path", "virus_fasta_name", "metadata_path"):
    if not config.get(key):
        raise ValueError(f"Required config key '{key}' is missing or empty.")


# ── Derived paths (mirrors CapelliniConfig.__post_init__) ────────────────────
def _dir_suffix(direction):
    return {"forward": "F", "reverse": "R", "paired": "P"}.get(direction, "F")


_BASE = config["base"]
_SUFFIX = _dir_suffix(config.get("direction", "forward"))
_DADA2 = _BASE + "/DADA2 output"
_MMSEQ = _BASE + "/MMSeqs2 Output"
_SP = _BASE + "/SpacePHARER output"
_PROCS = _BASE + "/Procs Estimations"
_NET = _BASE + "/Enhanced Networks"
_INPUT_FASTA = _BASE + "/Inputs/Fasta Collection"


# ── Include rule modules ─────────────────────────────────────────────────────
include: "rules/dada2.smk"
include: "rules/harmonize.smk"
include: "rules/spacepharer.smk"
include: "rules/protein_clusters.smk"
include: "rules/network.smk"


# ── Build the target list based on run_* flags ───────────────────────────────
def _network_targets():
    targets = []
    if config.get("run_common_abundance", True):
        targets.append(_NET + "/common/V_processed.csv")
        targets.append(_NET + "/common/B_processed.csv")
    if config.get("run_shrinkage_correlations", True):
        targets.append(_NET + "/shrinkage/shrinkage_corr_BV.csv.gz")
    if config.get("run_raw_crispr_networks", True):
        targets.append(_NET + "/crispr_raw/crispr_net.csv.gz")
    if config.get("run_smooth_crispr", True):
        targets.append(_NET + "/crispr_smooth/crispr_smooth_vir_bac.csv.gz")
    if config.get("run_xstar", True):
        targets.append(_NET + "/xstar/X_star.csv.gz")
        targets.append(_NET + "/xstar/shrinkage_xstar_corr_BV.csv.gz")
    return targets


rule all:
    """Default target: run the complete CAPELLINI pipeline."""
    input:
        # Stage 1: DADA2 — the ASV FASTA now stays in DADA2 output/ so the
        # harmonizer's --input-path can pick up both files together.
        _DADA2 + f"/OTU_table_{_SUFFIX}.csv",
        _DADA2 + f"/taxonomy_table_{_SUFFIX}.csv",
        _DADA2 + f"/ASV_sequences_{_SUFFIX}.fasta",
        # Stage 2: Harmonize (progenomes-harmonizer)
        _DADA2 + f"/taxonomy_table_{_SUFFIX}_withIDs.csv",
        # Stage 3: SpacePHARER
        _SP + "/output/phage_host_predictions.tsv",
        # Stage 4: Protein clustering (both count + binary matrices)
        _PROCS + "/pc_matrix.csv",
        _PROCS + "/pb_matrix.csv",
        # Stage 5: Network
        *_network_targets(),
