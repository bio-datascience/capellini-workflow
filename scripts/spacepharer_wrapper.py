#!/usr/bin/env python3
"""SpacePHARER wrapper: spacer extraction, DB creation, prediction, and statistics.

Snakemake script entry point. Migrated from capellini/stages/spacepharer.py.
All parameters are passed explicitly (no CapelliniConfig dependency).
"""

from __future__ import annotations

import argparse
import gzip
import logging
import os
import re
import shutil
import subprocess
import sys
import tempfile
import urllib.request
from collections import Counter, defaultdict
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import networkx as nx
import numpy as np
import pandas as pd
import seaborn as sns
from Bio import SeqIO

logger = logging.getLogger(__name__)

# proGenomes4 representative genome contigs — the source MinCED extracts CRISPR
# spacers from, and the proGenomes4 GCA->taxid map for target filtering.
PG4_GENOMES_URL = (
    "https://progenomes.embl.de/data/repGenomes/pg4_genomes_representatives.fna.gz"
)
PG4_NCBI_TAXONOMY_URL = "https://progenomes.embl.de/data/pg4_ncbi_taxonomy.tsv.gz"

# Bundled spacers collection path (proGenomes4-derived).
_BUNDLED_SPACERS = (
    Path(__file__).resolve().parent.parent
    / "capellini_workflow" / "data" / "references" / "spacers"
    / "spacers_CompleteCollection.fasta"
)


def sh(cmd: str, desc: str = "") -> subprocess.CompletedProcess:
    """Run a shell command, raising on failure."""
    if desc:
        print(desc)
    print(f"Executing command: {cmd}")
    try:
        r = subprocess.run(cmd, shell=True, capture_output=True, text=True, check=True)
    except subprocess.CalledProcessError as e:
        msg = (
            f"Command failed with code {e.returncode}\n"
            f"STDOUT:\n{e.stdout}\nSTDERR:\n{e.stderr}"
        )
        raise RuntimeError(msg) from e
    if r.stdout:
        print(r.stdout)
    if r.stderr.strip():
        print("STDERR (tool messages):\n" + r.stderr)
    return r


def check_and_install_spacepharer() -> tuple[str, str]:
    """Check that spacepharer and minced are on PATH; install via conda if missing."""
    if not (shutil.which("spacepharer") and shutil.which("minced")):
        logger.info("Installing spacepharer and minced via conda ...")
        subprocess.run(
            "conda install -y -c conda-forge -c bioconda spacepharer minced",
            shell=True,
            check=True,
        )
        logger.info("SpacePHARER and MinCED installed")
    else:
        logger.info("SpacePHARER and MinCED already installed")

    spacepharer_path = shutil.which("spacepharer") or "spacepharer"
    minced_path = shutil.which("minced") or "minced"
    return spacepharer_path, minced_path


class SpacePHARERWorkflow:
    """Wrapper around the SpacePHARER + MinCED command-line tools."""

    def __init__(self, workdir: str, spacerdir: str) -> None:
        self.wd = Path(workdir)
        self.sd = Path(spacerdir)
        self.spacepharer_path = shutil.which("spacepharer") or "spacepharer"
        self.minced_path = shutil.which("minced") or "minced"
        for d in ("spacers", "databases", "output", "tmp"):
            (self.wd / d).mkdir(parents=True, exist_ok=True)
        self.sd.mkdir(parents=True, exist_ok=True)

    def extract_spacers(
        self,
        fasta_path: str | Path,
        min_n_spacers: int,
        min_length: int,
        max_length: int,
        tag: str = "spacers_CompleteCollection",
    ) -> Path:
        """Extract CRISPR spacers from a FASTA using MinCED."""
        raw_txt = self.sd / "spacers" / f"{tag}.txt"
        raw_gff = self.sd / "spacers" / f"{tag}.gff"
        fa_out1 = self.sd / "spacers" / f"{tag}_spacers.fa"
        fa_out2 = self.sd / "spacers" / f"{tag}.fasta"

        (self.sd / "spacers").mkdir(parents=True, exist_ok=True)

        cmd = (
            f'"{self.minced_path}" -spacers -minNR {min_n_spacers} '
            f'-minRL {min_length} -maxRL {max_length} '
            f'"{fasta_path}" "{raw_txt}" "{raw_gff}"'
        )
        sh(cmd, "MinCED - Extracting CRISPR spacers")

        if fa_out1.exists():
            fa_out1.rename(fa_out2)
        return fa_out2

    def make_db(
        self,
        fasta_file: str | Path,
        dbname: str,
        is_spacer: bool = False,
        rev: bool = False,
    ) -> Path:
        """Create a SpacePHARER database from a FASTA file."""
        db_path = self.wd / "databases" / dbname
        tmp_path = self.wd / "tmp"

        flag = ""
        if is_spacer:
            flag += " --extractorf-spacer 1"
        if rev:
            flag += " --reverse-fragments 1"

        cmd = (
            f'"{self.spacepharer_path}" createsetdb '
            f'"{fasta_file}" "{db_path}" "{tmp_path}"{flag}'
        )
        sh(cmd, f"Creating DB {dbname}")
        return db_path

    def predict(
        self,
        spacerDB: Path,
        viralDB: Path,
        viralctrlDB: Path,
        out: str = "phage_host_predictions.tsv",
        fdr: float = 0.05,
    ) -> Path:
        """Run SpacePHARER predictmatch."""
        out_path = self.wd / "output" / out
        tmp_path = self.wd / "tmp"
        cmd = (
            f'"{self.spacepharer_path}" predictmatch '
            f'"{spacerDB}" "{viralDB}" "{viralctrlDB}" "{out_path}" "{tmp_path}" '
            f"--fdr {fdr}"
        )
        sh(cmd, "Running SpacePHARER predictmatch")
        return out_path

    def quick_stats(self, tsv: str | Path) -> None:
        """Print quick interaction statistics from a prediction TSV."""
        with open(str(tsv)) as fh:
            lines = [ln.strip() for ln in fh if ln.strip() and not ln.startswith("#")]
        hosts = len({ln.split("\t")[0] for ln in lines})
        phages = len({ln.split("\t")[1] for ln in lines})
        print(f"Interactions: {len(lines)}  Unique Hosts: {hosts}  Unique Phages: {phages}")


def get_spacers_collection(
    download_path: str,
    sp_folder: str,
    min_n_spacers: int = 3,
    min_length: int = 23,
    max_length: int = 47,
    remove_decomp_fasta: bool = True,
    bacContigs_reference_url: str = PG4_GENOMES_URL,
    regenerate: bool = False,
) -> Path:
    """Return path to spacers_CompleteCollection.fasta, using the bundled file if present."""
    spacers_path = Path(download_path) / "spacers"
    spacers_collection = spacers_path / "spacers_CompleteCollection.fasta"

    wf = SpacePHARERWorkflow(workdir=sp_folder, spacerdir=str(spacers_path))

    # Bundled file takes priority (must be non-empty to be trusted — an empty
    # placeholder must NOT block the MinCED regeneration fallback below).
    if _BUNDLED_SPACERS.exists() and _BUNDLED_SPACERS.stat().st_size > 0 and not regenerate:
        logger.info("Bundled spacers_CompleteCollection.fasta found - using it")
        spacers_path.mkdir(parents=True, exist_ok=True)
        if not spacers_collection.exists():
            shutil.copy2(str(_BUNDLED_SPACERS), str(spacers_collection))
        return spacers_collection

    if spacers_collection.exists() and spacers_collection.stat().st_size > 0:
        logger.info("Complete spacers collection found - skipping generation")
        return spacers_collection

    # Regenerate the complete spacer collection from the proGenomes4 genome
    # contigs with MinCED, in bounded batches. MinCED needs a plain FASTA and
    # cannot stream gzip, so we stream the .gz and hand it batch-by-batch — the
    # full ~150-200 GB uncompressed archive is never materialised (peak scratch
    # is ~one batch). CRISPR arrays live within a single contig, so batching at
    # contig boundaries yields the same spacers as a single whole-file run.
    filename = os.path.basename(bacContigs_reference_url)
    contigs_ref = Path(download_path) / filename
    if not contigs_ref.exists():
        logger.info("Downloading proGenomes4 genome contigs (~large) ...")
        urllib.request.urlretrieve(bacContigs_reference_url, str(contigs_ref))

    spacers_path.mkdir(parents=True, exist_ok=True)
    minced_path = shutil.which("minced") or "minced"
    opener = gzip.open if contigs_ref.suffix == ".gz" else open
    batch_bytes = 4 * 1024 ** 3
    n_batches = 0

    from capellini_workflow.io import atomic_write

    with tempfile.TemporaryDirectory(prefix="pg4_spacers.") as tmp:
        tmp_dir = Path(tmp)
        batch_fna = tmp_dir / "batch.fna"

        def run_batch(final_out) -> None:
            nonlocal n_batches
            n_batches += 1
            out_txt = tmp_dir / "batch.txt"
            out_gff = tmp_dir / "batch.gff"
            spacers_fa = tmp_dir / "batch_spacers.fa"  # MinCED names it <txt-stem>_spacers.fa
            for p in (out_txt, out_gff, spacers_fa):
                p.unlink(missing_ok=True)
            cmd = (
                f'"{minced_path}" -spacers -minNR {min_n_spacers} '
                f'-minRL {min_length} -maxRL {max_length} '
                f'"{batch_fna}" "{out_txt}" "{out_gff}"'
            )
            sh(cmd, f"MinCED - extracting spacers (batch {n_batches})")
            if spacers_fa.exists():
                with open(spacers_fa) as src:
                    shutil.copyfileobj(src, final_out)

        with atomic_write(spacers_collection) as final_out:
            bf = open(batch_fna, "w")
            cur = 0
            with opener(str(contigs_ref), "rt") as fh:
                for line in fh:
                    if line.startswith(">") and cur >= batch_bytes:
                        bf.close()
                        run_batch(final_out)
                        bf = open(batch_fna, "w")
                        cur = 0
                    bf.write(line)
                    cur += len(line)
            bf.close()
            if cur > 0:
                run_batch(final_out)

    logger.info("Spacer collection written: %s (%d batches)", spacers_collection, n_batches)
    return spacers_collection


def filter_target_spacers(
    spacers_collection: Path,
    ncbi_id_target_set_int: set,
    input_fasta_folder: str | Path,
    gca_to_taxid: dict[str, int],
) -> Path:
    """Filter spacers_CompleteCollection.fasta to cohort-specific NCBI IDs.

    proGenomes4 spacer headers begin with the versioned GCA accession
    (``GCA_000005825.2_...``), not the taxid, so the taxid is resolved from the
    GCA -> taxid map before comparing against the target set.
    """
    from progenomes_harmonizer.reference import parse_gca
    from capellini_workflow.io import atomic_write

    target_spacers_path = Path(input_fasta_folder) / "target_spacers.fasta"

    logger.info("Filtering spacers to %s target NCBI IDs", len(ncbi_id_target_set_int))
    # Atomic write so an interrupted run can't leave a half-filtered
    # target_spacers.fasta that the next run silently reuses.
    with atomic_write(target_spacers_path) as fasta_out:
        for record in SeqIO.parse(str(spacers_collection), "fasta"):
            gca = parse_gca(record.id)
            taxid = gca_to_taxid.get(gca) if gca else None
            if taxid in ncbi_id_target_set_int:
                SeqIO.write(record, fasta_out, "fasta")

    logger.info("Filtered spacers written: %s", target_spacers_path)
    return target_spacers_path


def _target_taxid_set(silva_fixed: pd.DataFrame, species_level: bool) -> set:
    """Cohort target taxids used to filter the CRISPR spacer collection.

    ``species_level`` only ever *adds* precision: species taxids are layered on
    top of the genus ones, never replacing them. ``addSpecies()`` assigns a
    species only on an exact, unambiguous match, so most ASVs have none —
    building the set as ``species | family`` would silently drop every
    genus-level target and make spacer selection coarser, the same trap fixed in
    the harmonizer's ``target_taxids`` chain.

    Args:
        silva_fixed: Harmonized taxonomy table with the progenomes_taxid_* columns.
        species_level: Whether to include the species layer.

    Returns:
        Set of integer NCBI taxids to target.
    """
    def _col(name: str) -> set:
        if name not in silva_fixed:
            return set()
        return set(silva_fixed[name].dropna().astype(int))

    taxids = _col("progenomes_taxid_genus") | _col("progenomes_taxid_family")
    if species_level:
        taxids |= _col("progenomes_taxid_species")
    return taxids


def run_spacepharer(
    silva_fixed_path: str,
    virus_fasta_path: str,
    sp_folder: str,
    input_fasta_folder: str,
    download_path: str,
    output_tsv: str,
    species_level: bool = False,
    min_n_spacers: int = 3,
    min_length: int = 23,
    max_length: int = 47,
    fdr: float = 0.05,
    remove_decomp_fasta: bool = True,
    bacContigs_reference_url: str = PG4_GENOMES_URL,
    pg4_ncbi_taxonomy_path: str = "",
    pg4_ncbi_taxonomy_url: str = PG4_NCBI_TAXONOMY_URL,
    regenerate_spacers: bool = False,
) -> None:
    """Run the full SpacePHARER stage: spacer collection, filtering, and prediction."""
    logger.info("SpacePHARER: starting")

    virus_fasta = Path(virus_fasta_path)
    if not virus_fasta.exists():
        raise FileNotFoundError(f"Virus FASTA not found: {virus_fasta}")

    check_and_install_spacepharer()

    # proGenomes4 spacer headers carry the versioned GCA (not the taxid), so load
    # the GCA -> taxid map used to filter spacers to the cohort's target taxa.
    from progenomes_harmonizer.reference import load_gca_to_taxid
    from capellini_workflow.io import download_if_missing
    if not pg4_ncbi_taxonomy_path:
        pg4_ncbi_taxonomy_path = str(Path(download_path) / Path(pg4_ncbi_taxonomy_url).name)
    download_if_missing(pg4_ncbi_taxonomy_url, pg4_ncbi_taxonomy_path,
                        label="proGenomes4 ncbi taxonomy")
    gca_to_taxid = load_gca_to_taxid(pg4_ncbi_taxonomy_path)

    silva_fixed = pd.read_csv(silva_fixed_path, index_col=0)

    # Build target sets
    ncbi_id_target_set_int = _target_taxid_set(silva_fixed, species_level)

    print(f"Resolution: {'species' if species_level else 'genus'}-level")
    print(f"Total unique NCBI target IDs: {len(ncbi_id_target_set_int)}")

    # Get spacers collection
    spacers_collection = get_spacers_collection(
        download_path=download_path,
        sp_folder=sp_folder,
        min_n_spacers=min_n_spacers,
        min_length=min_length,
        max_length=max_length,
        remove_decomp_fasta=remove_decomp_fasta,
        bacContigs_reference_url=bacContigs_reference_url,
        regenerate=regenerate_spacers,
    )

    # Filter to target taxa
    target_spacers_path = filter_target_spacers(
        spacers_collection, ncbi_id_target_set_int, input_fasta_folder, gca_to_taxid
    )

    # Run SpacePHARER
    spacers_path = Path(download_path) / "spacers"
    wf = SpacePHARERWorkflow(workdir=sp_folder, spacerdir=str(spacers_path))

    spacerDB = wf.make_db(target_spacers_path, "spacerDB", is_spacer=True)
    viralDB = wf.make_db(virus_fasta, "viralDB")
    controlDB = wf.make_db(virus_fasta, "viralDB_control", rev=True)

    tsv = wf.predict(spacerDB, viralDB, controlDB, out=Path(output_tsv).name, fdr=fdr)
    wf.quick_stats(tsv)

    logger.info("SpacePHARER stage complete")


def compute_spacepharer_stats(
    silva_fixed_path: str,
    topBitScore_path: str,
    virus_fasta_path: str,
    sp_folder: str,
    input_fasta_folder: str,
    download_path: str,
    species_level: bool = False,
    min_bitscore: int = 50,
) -> None:
    """Print the full pipeline statistics summary (Sections 1-3 from notebook)."""
    from progenomes_harmonizer.reference import parse_gca, load_gca_to_taxid
    from capellini_workflow.io import download_if_missing

    silva_fixed = pd.read_csv(silva_fixed_path, index_col=0)
    topBitScore_df = pd.read_csv(topBitScore_path, index_col=0)

    sp_folder = Path(sp_folder)
    virus_fasta = Path(virus_fasta_path)
    spacers_collection = Path(download_path) / "spacers" / "spacers_CompleteCollection.fasta"
    reference_16S_path = Path(input_fasta_folder) / "pg4_16s.fasta"

    # proGenomes4 headers carry the versioned GCA, not the taxid — resolve taxids
    # from the GCA -> taxid map.
    pg4_tax_path = Path(download_path) / Path(PG4_NCBI_TAXONOMY_URL).name
    download_if_missing(PG4_NCBI_TAXONOMY_URL, pg4_tax_path, label="proGenomes4 ncbi taxonomy")
    gca_to_taxid = load_gca_to_taxid(pg4_tax_path)

    n_viral_contigs = sum(1 for _ in SeqIO.parse(str(virus_fasta), "fasta"))

    n_16s_records = 0
    unique_16s_ncbi_ids: set = set()
    with open(str(reference_16S_path), "r") as fh:
        for rec in SeqIO.parse(fh, "fasta"):
            n_16s_records += 1
            taxid = gca_to_taxid.get(parse_gca(rec.id))
            if taxid is not None:
                unique_16s_ncbi_ids.add(taxid)
    n_16s_ncbi_ids = len(unique_16s_ncbi_ids)

    n_complete_spacers = sum(1 for _ in SeqIO.parse(str(spacers_collection), "fasta")) if spacers_collection.exists() else 0

    ncbi_id_target_set_int = _target_taxid_set(silva_fixed, species_level)

    ncbi_id_target_set_str = set(str(i) for i in ncbi_id_target_set_int)
    target_hits = topBitScore_df[topBitScore_df["NCBI ID"].isin(ncbi_id_target_set_str)]
    n_bac_genomes = target_hits["Genome Accession ID"].dropna().nunique()
    n_bac_genomes_ncbi = target_hits["NCBI ID"].dropna().nunique()

    target_spacers_path = Path(input_fasta_folder) / "target_spacers.fasta"
    n_bac_spacers = 0
    ncbi_ids_with_spacers: set = set()
    if target_spacers_path.exists():
        for record in SeqIO.parse(str(target_spacers_path), "fasta"):
            n_bac_spacers += 1
            taxid = gca_to_taxid.get(parse_gca(record.id))
            if taxid is not None:
                ncbi_ids_with_spacers.add(taxid)
    n_ncbi_with_spacers = len(ncbi_ids_with_spacers)

    n_asvs_total = len(silva_fixed)
    n_asvs_matched = topBitScore_df["query"].nunique()
    pct_matched = 100 * n_asvs_matched / n_asvs_total
    mean_bitscore = topBitScore_df["bitscore"].mean()
    median_bitscore = topBitScore_df["bitscore"].median()

    n_assigned_species = silva_fixed["progenomes_taxid_species"].notna().sum()
    n_assigned_genus = silva_fixed["progenomes_taxid_genus"].notna().sum()
    n_assigned_family = silva_fixed["progenomes_taxid_family"].notna().sum()
    n_assigned_any = silva_fixed[
        ["progenomes_taxid_species", "progenomes_taxid_genus", "progenomes_taxid_family"]
    ].notna().any(axis=1).sum()

    # SpacePHARER results
    pred_path = sp_folder / "output" / "phage_host_predictions.tsv"
    colnames_9 = ["spacer_id", "phage_id", "pvalue", "aln_len", "mismatch",
                  "qstart", "qend", "qprot_aln", "sprot_aln"]
    sp_df = pd.read_csv(
        str(pred_path), sep="\t", comment="#", header=None, names=colnames_9,
        engine="python", dtype=str, on_bad_lines="skip"
    ).dropna(how="all").reset_index(drop=True)
    for c in ["pvalue", "aln_len", "mismatch", "qstart", "qend"]:
        sp_df[c] = pd.to_numeric(sp_df[c], errors="coerce")
    sp_df["host_id"] = sp_df["spacer_id"].str.replace(r"_spacer_\d+$", "", regex=True)

    pairs = sp_df[["host_id", "phage_id"]].drop_duplicates().reset_index(drop=True)
    pairs = pairs.merge(
        sp_df.groupby(["host_id", "phage_id"])["spacer_id"].nunique().rename("n_spacers").reset_index(),
        on=["host_id", "phage_id"]
    )
    pairs = pairs.merge(
        sp_df.groupby(["host_id", "phage_id"])["pvalue"].min().rename("best_pval").reset_index(),
        on=["host_id", "phage_id"]
    )

    host_degree_sp = pairs.groupby("host_id")["phage_id"].nunique().rename("degree")
    virus_degree_sp = pairs.groupby("phage_id")["host_id"].nunique().rename("degree")
    n_hosts = host_degree_sp.shape[0]
    n_viruses = virus_degree_sp.shape[0]
    n_interactions = pairs.shape[0]
    n_spacer_hits = len(sp_df)
    density = n_interactions / (n_hosts * n_viruses) if n_hosts * n_viruses > 0 else 0.0

    print("=" * 64)
    print("  PIPELINE STATISTICS SUMMARY")
    print("=" * 64)
    print(f"\n  [INPUT]")
    print(f"  Total 16S ASVs                          : {n_asvs_total:>8,}")
    print(f"  Viral contigs                           : {n_viral_contigs:>8,}")
    print(f"\n  [REFERENCE FILES]")
    print(f"  16S ProGenomes ref records              : {n_16s_records:>8,}")
    print(f"  16S ProGenomes ref unique NCBI IDs      : {n_16s_ncbi_ids:>8,}")
    print(f"  ProGenomes spacers (complete collection): {n_complete_spacers:>8,}")
    print(f"\n  [TARGET SET]")
    print(f"  Bacterial genomes (unique GCAs)         : {n_bac_genomes:>8,}")
    print(f"  Bacterial genomes (unique NCBI IDs)     : {n_bac_genomes_ncbi:>8,}")
    print(f"  Bacterial spacers (cohort-specific)     : {n_bac_spacers:>8,}")
    print(f"\n  [MMSeqs2 - 16S vs proGenomes4 16S]")
    print(f"  ASVs with >=1 hit (bitscore >= {min_bitscore})      : {n_asvs_matched:>8,}  ({pct_matched:.1f}%)")
    print(f"  Mean / Median bitscore                  : {mean_bitscore:>6.1f}  / {median_bitscore:.1f}")
    print(f"\n  [NCBI TAXID ASSIGNMENT - 3-layer]")
    print(f"  ASVs assigned (any layer)               : {n_assigned_any:>8,}  ({100*n_assigned_any/n_asvs_total:.1f}%)")
    print(f"    Layer 1 - species                     : {n_assigned_species:>8,}")
    print(f"    Layer 2 - genus                       : {n_assigned_genus:>8,}")
    print(f"    Layer 3 - family                      : {n_assigned_family:>8,}")
    print(f"  NCBI IDs with >=1 spacer                : {n_ncbi_with_spacers:>8,}")
    print(f"\n  [SPACEPHARER NETWORK]")
    print(f"  Unique hosts (bacteria)                 : {n_hosts:>8,}")
    print(f"  Unique viruses (phages/contigs)         : {n_viruses:>8,}")
    print(f"  Unique host-phage interactions          : {n_interactions:>8,}")
    print(f"  Raw spacer-hit rows                     : {n_spacer_hits:>8,}")
    print(f"  Network density                         : {density:>8.5f}")
    print("=" * 64)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Run SpacePHARER stage")
    parser.add_argument("--silva-fixed", required=True, help="Path to silva_fixed CSV")
    parser.add_argument("--virus-fasta", required=True, help="Path to virus FASTA")
    parser.add_argument("--sp-folder", required=True, help="SpacePHARER working directory")
    parser.add_argument("--input-fasta-folder", required=True, help="Input FASTA folder")
    parser.add_argument("--download-path", required=True, help="Downloads directory")
    parser.add_argument("--output-tsv", required=True, help="Output predictions TSV path")
    parser.add_argument("--species-level", action="store_true", default=False)
    parser.add_argument("--min-n-spacers", type=int, default=3)
    parser.add_argument("--min-length", type=int, default=23)
    parser.add_argument("--max-length", type=int, default=47)
    parser.add_argument("--fdr", type=float, default=0.05)
    parser.add_argument("--remove-decomp-fasta", action="store_true", default=True)
    parser.add_argument("--pg4-ncbi-taxonomy", default="",
                        help="Path to pg4_ncbi_taxonomy.tsv(.gz) (GCA->taxid map). "
                             "Defaults to <download-path>/pg4_ncbi_taxonomy.tsv.gz; "
                             "downloaded if missing.")
    args = parser.parse_args()

    logging.basicConfig(level=logging.INFO)
    run_spacepharer(
        silva_fixed_path=args.silva_fixed,
        virus_fasta_path=args.virus_fasta,
        sp_folder=args.sp_folder,
        input_fasta_folder=args.input_fasta_folder,
        download_path=args.download_path,
        output_tsv=args.output_tsv,
        species_level=args.species_level,
        min_n_spacers=args.min_n_spacers,
        min_length=args.min_length,
        max_length=args.max_length,
        fdr=args.fdr,
        remove_decomp_fasta=args.remove_decomp_fasta,
        pg4_ncbi_taxonomy_path=args.pg4_ncbi_taxonomy,
    )
