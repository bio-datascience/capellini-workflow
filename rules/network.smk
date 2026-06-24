"""Network rules: common abundance, shrinkage, raw CRISPR, smooth CRISPR, X*."""

import json as _json


def _dir_suffix(direction):
    return {"forward": "F", "reverse": "R", "paired": "P"}.get(direction, "F")


_NET_ROOT = config["base"] + "/Enhanced Networks"
_DADA2 = config["base"] + "/DADA2 output"
_SUFFIX = _dir_suffix(config["direction"])
_SP_FOLDER = config["base"] + "/SpacePHARER output"


# --- Resolve network input paths (same logic as CapelliniConfig.__post_init__) ---
_bacteria_otu = config.get("bacteria_otu", "") or f"{_DADA2}/OTU_table_{_SUFFIX}.csv"
_bacteria_taxonomy = config.get("bacteria_taxonomy", "") or f"{_DADA2}/taxonomy_table_{_SUFFIX}_withIDs.csv"
_phage_host_predictions = config.get("phage_host_predictions", "") or f"{_SP_FOLDER}/output/phage_host_predictions.tsv"
_tax_bac_for_smoothing = config.get("tax_bac_for_smoothing", "") or f"{_DADA2}/taxonomy_table_{_SUFFIX}_withIDs.csv"

# These may contain un-interpolated Python f-strings from legacy configs — detect and fix
import re
_FSTR_RE = re.compile(r'^\s*f["\'].*["\']\s*$')
if _FSTR_RE.match(_bacteria_otu):
    _bacteria_otu = f"{_DADA2}/OTU_table_{_SUFFIX}.csv"
if _FSTR_RE.match(_bacteria_taxonomy):
    _bacteria_taxonomy = f"{_DADA2}/taxonomy_table_{_SUFFIX}_withIDs.csv"
if _FSTR_RE.match(_phage_host_predictions):
    _phage_host_predictions = f"{_SP_FOLDER}/output/phage_host_predictions.tsv"
if _FSTR_RE.match(_tax_bac_for_smoothing):
    _tax_bac_for_smoothing = f"{_DADA2}/taxonomy_table_{_SUFFIX}_withIDs.csv"


rule common_abundance:
    """Build common abundance matrices (V_processed, B_processed)."""
    input:
        virus_abundance_raw=config["virus_abundance_raw"],
        bacteria_otu=_bacteria_otu,
        bacteria_taxonomy=_bacteria_taxonomy,
        metadata=config["metadata_path"],
    output:
        v_processed=_NET_ROOT + "/common/V_processed.csv",
        b_processed=_NET_ROOT + "/common/B_processed.csv",
        meta_aligned=_NET_ROOT + "/common/metadata_aligned.csv",
    params:
        output_dir=_NET_ROOT + "/common",
        bacteria_taxonomy_rank=config.get("bacteria_taxonomy_rank", "target_taxids"),
        keep_column=config.get("keep_column", "keep_for_analysis"),
        prevalence=config.get("prevalence", 0.10),
        script=str(SCRIPTS_DIR / "common_abundance.py"),
    shell:
        """
        python "{params.script}" \
            --virus-abundance-raw "{input.virus_abundance_raw}" \
            --bacteria-otu "{input.bacteria_otu}" \
            --bacteria-taxonomy "{input.bacteria_taxonomy}" \
            --metadata "{input.metadata}" \
            --output-dir "{params.output_dir}" \
            --bacteria-taxonomy-rank "{params.bacteria_taxonomy_rank}" \
            --keep-column "{params.keep_column}" \
            --prevalence {params.prevalence}
        """


rule shrinkage_correlations:
    """Schaefer-Strimmer shrinkage correlations on Z = [B_clr V_clr]."""
    input:
        v_processed=_NET_ROOT + "/common/V_processed.csv",
        b_processed=_NET_ROOT + "/common/B_processed.csv",
    output:
        corr_full=_NET_ROOT + "/shrinkage/shrinkage_corr_full.csv.gz",
        corr_bv=_NET_ROOT + "/shrinkage/shrinkage_corr_BV.csv.gz",
    params:
        output_dir=_NET_ROOT + "/shrinkage",
        pseudocount=config.get("pseudocount", 1e-6),
        script=str(SCRIPTS_DIR / "shrinkage_correlations.py"),
    shell:
        """
        python "{params.script}" \
            --v-processed "{input.v_processed}" \
            --b-processed "{input.b_processed}" \
            --output-dir "{params.output_dir}" \
            --pseudocount {params.pseudocount}
        """


rule raw_crispr_network:
    """Aggregate SpacePHARER predictions into a (bac x vOTU) CRISPR matrix."""
    input:
        predictions=_phage_host_predictions,
        tax_vir=config["tax_vir"],
    output:
        crispr_net=_NET_ROOT + "/crispr_raw/crispr_net.csv.gz",
    params:
        output_dir=_NET_ROOT + "/crispr_raw",
        aggregate_viral_rank=config.get("aggregate_viral_rank", "lev0"),
        script=str(SCRIPTS_DIR / "raw_crispr_network.py"),
    shell:
        """
        python "{params.script}" \
            --phage-host-predictions "{input.predictions}" \
            --tax-vir "{input.tax_vir}" \
            --output-dir "{params.output_dir}" \
            --aggregate-viral-rank "{params.aggregate_viral_rank}"
        """


rule smooth_crispr:
    """Smooth the raw CRISPR matrix via taxonomy kernels."""
    input:
        v_processed=_NET_ROOT + "/common/V_processed.csv",
        b_processed=_NET_ROOT + "/common/B_processed.csv",
        crispr_raw=_NET_ROOT + "/crispr_raw/crispr_net.csv.gz",
        tax_bac=_tax_bac_for_smoothing,
        tax_vir=config["tax_vir"],
    output:
        crispr_smooth=_NET_ROOT + "/crispr_smooth/crispr_smooth.csv.gz",
        crispr_smooth_vir_bac=_NET_ROOT + "/crispr_smooth/crispr_smooth_vir_bac.csv.gz",
    params:
        output_dir=_NET_ROOT + "/crispr_smooth",
        # Serialize the list params as JSON so smooth_crispr.py can json.loads()
        # them. Snakemake otherwise stringifies a Python list by joining elements
        # with spaces, which would silently corrupt these inputs.
        bacterial_ranks=_json.dumps(config.get("bacterial_ranks", ["Phylum", "Class", "Order", "Family", "progenomes_taxid_genus"])),
        viral_ranks=_json.dumps(config.get("viral_ranks", ["lev8", "lev7", "lev6", "lev5", "lev4", "lev3", "lev2", "lev1", "lev0"])),
        bacterial_weights=_json.dumps(config.get("bacterial_weights", [1, 2, 3, 6, 8])),
        viral_weights=_json.dumps(config.get("viral_weights", [1, 1, 2, 3, 4, 6, 8, 10, 12])),
        alpha=config.get("crispr_smooth_alpha", 0.95),
        transpose="--transpose-after-load" if config.get("transpose_raw_crispr_after_load", True) else "",
        script=str(SCRIPTS_DIR / "smooth_crispr.py"),
    shell:
        """
        python "{params.script}" \
            --v-processed "{input.v_processed}" \
            --b-processed "{input.b_processed}" \
            --crispr-raw "{input.crispr_raw}" \
            --tax-bac-for-smoothing "{input.tax_bac}" \
            --tax-vir "{input.tax_vir}" \
            --output-dir "{params.output_dir}" \
            --bacterial-ranks '{params.bacterial_ranks}' \
            --viral-ranks '{params.viral_ranks}' \
            --bacterial-weights '{params.bacterial_weights}' \
            --viral-weights '{params.viral_weights}' \
            --alpha {params.alpha} \
            {params.transpose}
        """


rule xstar:
    """Residual X* propagation + Schaefer-Strimmer correlations on Z*."""
    input:
        v_processed=_NET_ROOT + "/common/V_processed.csv",
        b_processed=_NET_ROOT + "/common/B_processed.csv",
        crispr_smooth_vir_bac=_NET_ROOT + "/crispr_smooth/crispr_smooth_vir_bac.csv.gz",
    output:
        x_star=_NET_ROOT + "/xstar/X_star.csv.gz",
        shrinkage_bv=_NET_ROOT + "/xstar/shrinkage_xstar_corr_BV.csv.gz",
    params:
        output_dir=_NET_ROOT + "/xstar",
        pseudocount=config.get("pseudocount", 1e-6),
        lam=config.get("lam", 0.5),
        n_steps=config.get("n_steps", 1),
        script=str(SCRIPTS_DIR / "xstar.py"),
    shell:
        """
        python "{params.script}" \
            --v-processed "{input.v_processed}" \
            --b-processed "{input.b_processed}" \
            --crispr-smooth-vir-bac "{input.crispr_smooth_vir_bac}" \
            --output-dir "{params.output_dir}" \
            --pseudocount {params.pseudocount} \
            --lam {params.lam} \
            --n-steps {params.n_steps}
        """
