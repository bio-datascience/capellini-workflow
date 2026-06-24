"""SpacePHARER rule: spacer extraction, filtering, DB creation, prediction."""


def _dir_suffix(direction):
    return {"forward": "F", "reverse": "R", "paired": "P"}.get(direction, "F")


rule spacepharer:
    """Run the SpacePHARER stage."""
    input:
        silva_fixed=config["base"] + "/DADA2 output/taxonomy_table_{suffix}_withIDs.csv".format(
            suffix=_dir_suffix(config["direction"])
        ),
        virus_fasta=config["base"] + "/Inputs/Fasta Collection/" + config["virus_fasta_name"],
    output:
        predictions=config["base"] + "/SpacePHARER output/output/phage_host_predictions.tsv",
    params:
        sp_folder=config["base"] + "/SpacePHARER output",
        input_fasta_folder=config["base"] + "/Inputs/Fasta Collection",
        download_path=config["download_path"],
        species_level="--species-level" if config.get("species_level", False) else "",
        min_n_spacers=config.get("min_n_spacers", 3),
        min_length=config.get("min_length", 23),
        max_length=config.get("max_length", 47),
        fdr=config.get("fdr", 0.05),
        script=str(SCRIPTS_DIR / "spacepharer_wrapper.py"),
    shell:
        """
        python "{params.script}" \
            --silva-fixed "{input.silva_fixed}" \
            --virus-fasta "{input.virus_fasta}" \
            --sp-folder "{params.sp_folder}" \
            --input-fasta-folder "{params.input_fasta_folder}" \
            --download-path "{params.download_path}" \
            --output-tsv "{output.predictions}" \
            --min-n-spacers {params.min_n_spacers} \
            --min-length {params.min_length} \
            --max-length {params.max_length} \
            --fdr {params.fdr} \
            {params.species_level}
        """
