"""DADA2 rule: run the R pipeline to produce OTU table, taxonomy table, and ASV FASTA.

The R script writes the ASV FASTA into the same DADA2 output/ folder as the
taxonomy table. We deliberately do NOT move it elsewhere — the harmonizer's
``--input-path`` mode expects both files to live side-by-side.
"""


def _dir_suffix(direction):
    return {"forward": "F", "reverse": "R", "paired": "P"}.get(direction, "F")


rule dada2:
    """Run the DADA2 R pipeline."""
    input:
        bacterial_raw=config["bacterial_raw_fasta_folder"],
        silva_ref=config["silva_ref_path"],
        silva_taxmap=config["silva_taxmap_path"],
    output:
        otu_table=config["base"] + "/DADA2 output/OTU_table_{suffix}.csv".format(
            suffix=_dir_suffix(config["direction"])
        ),
        taxonomy_table=config["base"] + "/DADA2 output/taxonomy_table_{suffix}.csv".format(
            suffix=_dir_suffix(config["direction"])
        ),
        asv_fasta=config["base"] + "/DADA2 output/ASV_sequences_{suffix}.fasta".format(
            suffix=_dir_suffix(config["direction"])
        ),
    params:
        output_dir=config["base"] + "/DADA2 output",
        direction=config["direction"],
        fasta_generation="TRUE" if config.get("fasta_generation", True) else "FALSE",
        # Toggle removeBimeraDenovo from the YAML config. Defaults to FALSE
        # to preserve the validated legacy ASV count when the key is absent
        # (e.g. when using configs that predate this knob).
        chimera_removal="TRUE" if config.get("chimera_removal", False) else "FALSE",
        # SILVA species-assignment reference for addSpecies(). Declared as a
        # *param* (not an input) because the R script downloads it on demand —
        # as an input Snakemake would abort before it could be fetched. It always
        # lives in download_path with every other reference; not configurable.
        species_path=(
            config["download_path"].rstrip("/") + "/silva_species_assignment_v138.1.fa.gz"
        ),
        species_url=config.get(
            "silva_species_url",
            "https://zenodo.org/records/4587955/files/silva_species_assignment_v138.1.fa.gz?download=1",
        ),
        r_script=str(SCRIPTS_DIR / "dada2_pipe.R"),
    shell:
        """
        mkdir -p "{params.output_dir}"
        Rscript "{params.r_script}" \
            "{input.bacterial_raw}" \
            "{params.output_dir}" \
            "{input.silva_ref}" \
            "{input.silva_taxmap}" \
            "{params.direction}" \
            "{params.fasta_generation}" \
            "{params.chimera_removal}" \
            "{params.species_path}" \
            "{params.species_url}"
        """
