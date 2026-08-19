"""Harmonize rule: run progenomes-harmonizer (NCBI mapping + 3-layer assignment).

Uses the harmonizer's ``--input-path`` mode: it auto-detects the taxonomy
table and the ASV FASTA from the same folder (DADA2 output/), so we only have
to point it at one directory.
"""


def _dir_suffix(direction):
    return {"forward": "F", "reverse": "R", "paired": "P"}.get(direction, "F")


rule harmonize:
    """Run progenomes-harmonizer."""
    # Falls back to every core when `cores` is unset (or empty — `config.get`
    # with a default would return "" there and Snakemake needs an int).
    # Snakemake caps this by the --cores given on the command line anyway.
    threads: config.get("cores") or (os.cpu_count() or 8)
    input:
        # Declared as inputs so Snakemake tracks them as dependencies — but
        # the harmonizer rediscovers them itself from --input-path.
        taxonomy_table=config["base"] + "/DADA2 output/taxonomy_table_{suffix}.csv".format(
            suffix=_dir_suffix(config["direction"])
        ),
        asv_fasta=config["base"] + "/DADA2 output/ASV_sequences_{suffix}.fasta".format(
            suffix=_dir_suffix(config["direction"])
        ),
    output:
        silva_fixed=config["base"] + "/DADA2 output/taxonomy_table_{suffix}_withIDs.csv".format(
            suffix=_dir_suffix(config["direction"])
        ),
    params:
        input_path=config["base"] + "/DADA2 output",
        download_path=config["download_path"],
        mmseq_folder=config["base"] + "/MMSeqs2 Output",
        direction=config["direction"],
        min_bitscore=config.get("min_bitscore", 50),
        max_matches=config.get("max_matches", 20),
        ncbi_taxdmp_url=config.get("ncbi_taxdmp_url", "https://ftp.ncbi.nih.gov/pub/taxonomy/taxdmp.zip"),
        species_level="--species-level" if config.get("species_level", False) else "",
        ref_removal="--ref-removal" if config.get("ref_removal", True) else "",
        # Optional pre-built proGenomes4 16S reference (else built with barrnap).
        reference_16s='--reference-16s "{}"'.format(config["reference_16s"]) if config.get("reference_16s", "") else "",
        # proGenomes4 GCA->taxid map (optional; auto-downloaded if omitted).
        pg4_ncbi_taxonomy='--pg4-ncbi-taxonomy "{}"'.format(config["pg4_ncbi_taxonomy"]) if config.get("pg4_ncbi_taxonomy", "") else "",
        # barrnap 16S-build params (parallel barrnap over split genome chunks).
        barrnap_kingdom=config.get("barrnap_kingdom", "bac"),
        barrnap_evalue=config.get("barrnap_evalue", 1e-6),
        barrnap_reject=config.get("barrnap_reject", 0.25),
    shell:
        """
        mkdir -p "{params.mmseq_folder}"
        progenomes-harmonize \
            --input-path "{params.input_path}" \
            --download-path "{params.download_path}" \
            --mmseq-folder "{params.mmseq_folder}" \
            --output "{output.silva_fixed}" \
            --direction "{params.direction}" \
            --threads {threads} \
            --min-bitscore {params.min_bitscore} \
            --max-matches {params.max_matches} \
            --ncbi-taxdmp-url "{params.ncbi_taxdmp_url}" \
            --barrnap-kingdom "{params.barrnap_kingdom}" \
            --barrnap-evalue {params.barrnap_evalue} \
            --barrnap-reject {params.barrnap_reject} \
            {params.species_level} \
            {params.ref_removal} \
            {params.reference_16s} \
            {params.pg4_ncbi_taxonomy}
        """
