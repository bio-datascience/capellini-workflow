"""ProCs Maker rule: extract proteins, cluster, build pc/pb matrices."""


def _dir_suffix(direction):
    return {"forward": "F", "reverse": "R", "paired": "P"}.get(direction, "F")


rule protein_clusters:
    """Run the ProCs Maker pipeline (bi-modal: bacterial + viral)."""
    input:
        bac_gca_table=config["base"] + "/DADA2 output/taxonomy_table_{suffix}_withIDs.csv".format(
            suffix=_dir_suffix(config["direction"])
        ),
        viral_contigs=config["base"] + "/Inputs/Fasta Collection/" + config["virus_fasta_name"],
    output:
        # procs_maker always emits both matrices in bi-modal mode.
        pc_matrix=config["base"] + "/Procs Estimations/pc_matrix.csv",
        pb_matrix=config["base"] + "/Procs Estimations/pb_matrix.csv",
    params:
        output_dir=config["base"] + "/Procs Estimations",
        download_path=config["download_path"],
        batch_size=config.get("batch_size", 1500),
        filter_flag="--filter-1bac-1vir" if config.get("filter_1bac_1vir", False) else "",
        species_flag="--species-level" if config.get("species_level", False) else "",
    shell:
        """
        procs_maker \
            --bac-gca-table "{input.bac_gca_table}" \
            --viral-contigs "{input.viral_contigs}" \
            --output-dir "{params.output_dir}" \
            --download-path "{params.download_path}" \
            --batch-size {params.batch_size} \
            {params.species_flag} \
            {params.filter_flag}
        """
