# capellini-workflow

**C**RISPR-**A**bundance **P**hage-**E**vidence **L**inkage for **L**everaging **I**nteraction **N**etwork **I**nference.

Snakemake orchestration of the full CAPELLINI pipeline. This package wires
three independent tools - `progenomes-harmonizer`, `procs-maker`, and the
network inference scripts - into a reproducible, file-based workflow.

The pipeline runs 9 rules in sequence:

1. **dada2** - 16S amplicon denoising and taxonomy assignment (R/DADA2)
2. **harmonize** - 3-layer NCBI/GCA taxonomy assignment (`progenomes-harmonizer`)
3. **spacepharer** - CRISPR spacer extraction and phage-host prediction (SpacePHARER + MinCED)
4. **protein_clusters** - bacterial/viral protein extraction and clustering (`procs-maker`)
5. **common_abundance** - prevalence-filtered, aligned abundance matrices
6. **shrinkage_correlations** - Schaefer-Strimmer shrinkage correlation on joint CLR-transformed [bacteria × virus] matrix
7. **raw_crispr_network** - aggregated CRISPR spacer interaction matrix
8. **smooth_crispr** - taxonomy-kernel-smoothed CRISPR matrix
9. **xstar** - residual message-passing propagation + shrinkage on Z\*

Rules 5–9 are individually togglable. Each rule tracks its outputs as Snakemake
files so the workflow can resume from any interruption point.

---

## Installation

Install in order - the workflow depends on both upstream tools:

```bash
pip install git+https://github.com/bio-datascience/progenomes-harmonizer.git
pip install git+https://github.com/bio-datascience/procs-maker.git
pip install git+https://github.com/bio-datascience/capellini-workflow.git
```

### External tools (must be on PATH)

| Tool | Used by |
|---|---|
| `mmseqs` / `mmseqs2` | harmonize, spacepharer |
| `spacepharer` | spacepharer |
| `minced` | spacepharer (CRISPR spacer extraction) |
| `prodigal` | protein_clusters (viral gene prediction) |
| `barrnap` | harmonize — only when building the proGenomes4 16S reference (skipped if `reference_16s` is bundled/supplied) |
| `Rscript` | dada2 (with packages: DADA2, phyloseq, data.table) |

---

## Usage

The primary interface is the interactive terminal UI:

```bash
capellini
```

This opens a menu-driven interface where you can load a config file, set the
number of cores, preview the DAG, and run the full pipeline or individual stages, 
all without writing Snakemake commands directly.

You can also invoke Snakemake directly if you prefer:

```bash
# Check the DAG without running anything
snakemake --configfile config/ibd.yaml --dry-run

# Run the full pipeline
snakemake --configfile config/ibd.yaml --cores 8

# Run up to a specific stage (triggers all upstream dependencies)
snakemake --configfile config/ibd.yaml --cores 8 protein_clusters

# Resume after an interruption
snakemake --configfile config/ibd.yaml --cores 8 --rerun-incomplete
```

---

## Configuration

Each study gets its own YAML config file. Copy `config/blank.yaml` as a
starting point and fill in the required paths.

```bash
cp config/blank.yaml config/mystudy.yaml
```

### Required

These five keys must be set before the pipeline can run:

| Key | Description |
|---|---|
| `base` | Root output directory. All pipeline outputs go under this path. |
| `download_path` | Directory where large reference files are cached (proGenomes4 archives, NCBI taxonomy, SpacePHARER database). |
| `virus_fasta_name` | Filename of the viral contigs FASTA, expected at `<base>/<virus_fasta_name>`. |
| `metadata_path` | Path to the sample metadata CSV. Must contain a `keep_column` column (see `keep_column` below) to mark samples included in the analysis. |
| `bacterial_raw_fasta_folder` | Folder containing the raw 16S FASTQ files (one per sample), used by the DADA2 rule. |

### Global flags

| Key | Default | Description |
|---|---|---|
| `cores` | `1` | Number of CPU cores passed to Snakemake. Also set interactively via the `capellini` UI. |
| `species_level` | `false` | Use species-level GCA resolution throughout the pipeline (harmonizer + procs-maker). Default is genus-level. |
| `ref_removal` | `true` | Delete large downloaded reference archives after they have been used (saves disk space). Individual rules may re-download on the next run if needed. |
| `fresh_start` | `false` | Force re-run of all rules, ignoring existing outputs. Equivalent to deleting all outputs and restarting. |

### DADA2

| Key | Default | Description |
|---|---|---|
| `direction` | `"forward"` | Sequencing direction. One of `"forward"`, `"reverse"`, `"paired"`. Determines which FASTQ files are read and the output filename suffix (`F`, `R`, or `P`). |
| `bacteria_fasta_name` | `"16S_DADA2_bacteria.fasta"` | Filename for the ASV FASTA written by the DADA2 rule, used by the harmonize rule. |
| `fasta_generation` | `true` | Whether to write the ASV FASTA alongside the OTU/taxonomy tables. Disable only if you are providing a pre-built FASTA. |
| `chimera_removal` | `false` | Run `removeBimeraDenovo` (DADA2 chimera filtering) after denoising. Can significantly reduce ASV count (~73% in some datasets). Off by default - enable explicitly if your protocol requires it. |
| `silva_ref_path` | `""` | Path to the SILVA classifier FASTA (e.g. `silva_nr99_v138.1_train_set.fa.gz`). Required by the DADA2 rule for `assignTaxonomy`. |
| `silva_taxmap_path` | `""` | Path to the SILVA taxonomy map (e.g. `tax_slv_ssu_138.1.txt`). Required for species-level taxonomy assignment. |

### Harmonizer (progenomes-harmonizer)

| Key | Default | Description |
|---|---|---|
| `min_bitscore` | `50` | Minimum MMseqs2 bitscore for a hit to be accepted during 16S mapping. |
| `max_matches` | `20` | Maximum MMseqs2 hits retained per ASV query. |
| `ncbi_taxdmp_url` | NCBI FTP | Override URL for the NCBI taxonomy dump archive. |
| `reference_16s` | `""` | Optional path to a pre-built proGenomes4 16S reference (`pg4_16s.fasta`). When set, the harmonize rule skips building it with barrnap (and the ~46 GB genome-contigs download). Leave empty to use the bundled reference or build on demand. |
| `pg4_ncbi_taxonomy` | `""` | Optional path to the proGenomes4 `pg4_ncbi_taxonomy.tsv(.gz)` (`GCA → taxid` map). Auto-downloaded to `download_path` if empty. |

### SpacePHARER

| Key | Default | Description |
|---|---|---|
| `min_n_spacers` | `3` | Minimum number of spacers a CRISPR array must have to be retained. |
| `min_length` | `23` | Minimum spacer length (bp). |
| `max_length` | `47` | Maximum spacer length (bp). |
| `fdr` | `0.05` | False discovery rate threshold for SpacePHARER phage-host predictions. |
| `keep_spacers_collection` | `true` | Keep the intermediate spacers FASTA after SpacePHARER completes. |
| `remove_decomp_fasta` | `true` | Delete the decompressed proGenomes4 contigs FASTA after spacer search. |

> Spacepharer uses the **bundled** proGenomes4 spacer collection and joins on
> NCBI taxids. Only if that collection is absent does it regenerate spacers from
> the proGenomes4 genome contigs (`pg4_genomes_representatives.fna.gz`) with
> MinCED.

### Protein clustering (procs-maker)

| Key | Default | Description |
|---|---|---|
| `batch_size` | `1500` | Genome batch size for streaming the proGenomes4 protein archive. Reduce if you run out of RAM. |
| `filter_1bac_1vir` | `false` | Keep only protein clusters containing at least one bacterial and one viral protein (cross-domain clusters only). |
| `save_single_bacgenome_collection` | `false` | Save individual per-genome protein FASTA files in addition to the combined collection. |
| `keep_coords` | `false` | Keep Prodigal's `coords.gbk` gene annotation file after viral protein extraction. |
| `remove_collections` | `false` | Delete intermediate protein FASTA collections after clustering. |

### Network

#### Stage toggles

| Key | Default | Description |
|---|---|---|
| `run_common_abundance` | `true` | Run the common-abundance alignment step (rule 5). |
| `run_shrinkage_correlations` | `true` | Run Schaefer-Strimmer shrinkage correlation (rule 6). |
| `run_raw_crispr_networks` | `true` | Build the raw CRISPR interaction matrix (rule 7). |
| `run_smooth_crispr` | `true` | Apply taxonomy-kernel smoothing to the CRISPR matrix (rule 8). |
| `run_xstar` | `true` | Run residual message-passing propagation (rule 9). |

#### Abundance and filtering

| Key | Default | Description |
|---|---|---|
| `prevalence` | `0.10` | Minimum prevalence threshold (fraction of samples) for retaining a taxon in the abundance matrices. |
| `keep_column` | `"keep_for_analysis"` | Column name in the metadata CSV that marks which samples to include (`1` = include, `0` = exclude). |
| `pseudocount` | `1.0e-6` | Pseudocount added before CLR transformation to avoid log(0). |

#### Bacterial taxonomy smoothing

| Key | Default | Description |
|---|---|---|
| `bacteria_taxonomy_rank` | `"target_taxids"` | Column from the harmonized taxonomy table used to map bacteria to the CRISPR network. |
| `bacterial_ranks` | `["Phylum","Class","Order","Family","progenomes_taxid_genus"]` | Taxonomic ranks used to build the bacterial taxonomy kernel for CRISPR smoothing. |
| `bacterial_weights` | `[1, 2, 3, 6, 8]` | Weights assigned to each bacterial rank for the smoothing kernel (higher = more influence). Must match the length of `bacterial_ranks`. |

#### Viral taxonomy smoothing

| Key | Default | Description |
|---|---|---|
| `viral_ranks` | `["lev8",…,"lev0"]` | Viral taxonomy levels (from the viral taxonomy table) used for the smoothing kernel. |
| `viral_weights` | `[1, 1, 2, 3, 4, 6, 8, 10, 12]` | Weights for each viral rank. Must match the length of `viral_ranks`. |
| `aggregate_viral_rank` | `"lev0"` | Viral rank at which to aggregate the CRISPR matrix before smoothing. |

#### CRISPR propagation (X\*)

| Key | Default | Description |
|---|---|---|
| `crispr_smooth_alpha` | `0.95` | Smoothing weight for the taxonomy kernel applied to the CRISPR matrix (0 = no smoothing, 1 = full kernel). |
| `lam` | `0.5` | Lambda regularization for the X\* message-passing step. |
| `n_steps` | `1` | Number of propagation steps in the X\* algorithm. |
| `preserve_scale` | `false` | Normalize the X\* output matrix to preserve the scale of the input. |
| `transpose_raw_crispr_after_load` | `true` | Transpose the raw CRISPR matrix after loading (bacteria × viruses → viruses × bacteria). Set based on your SpacePHARER output orientation. |

#### Pre-computed input overrides

These keys are all optional. When left empty, the workflow derives each path
automatically from `base` and the upstream rule outputs. Set them explicitly
only if you want to inject pre-computed matrices and skip upstream rules.

| Key | Description |
|---|---|
| `virus_abundance_raw` | Path to a pre-computed viral abundance matrix (skips abundance derivation). |
| `bacteria_otu` | Path to a pre-computed bacterial OTU table. |
| `bacteria_taxonomy` | Path to a pre-computed bacterial taxonomy table. |
| `phage_host_predictions` | Path to a pre-computed SpacePHARER output TSV. |
| `tax_bac_for_smoothing` | Path to the taxonomy table used for bacterial CRISPR kernel smoothing. |
| `tax_vir` | Path to the viral taxonomy table. |

---

## Outputs

All outputs go under the `base` directory:

```
<base>/
├── DADA2 output/
│   ├── OTU_table_F.csv               (ASV count table)
│   ├── taxonomy_table_F.csv          (DADA2 taxonomy assignments)
│   ├── ASV_sequences_F.fasta         (representative 16S sequences)
│   └── taxonomy_table_F_withIDs.csv  (harmonized table with GCA accessions)
├── SpacePHARER output/
│   └── phage_host_predictions.tsv    (CRISPR-based phage-host pairs)
├── Procs Estimations/
│   ├── pc_matrix.csv                 (protein-count matrix)
│   └── pb_matrix.csv                 (protein presence/absence matrix)
└── Enhanced Networks/
    ├── common_abundance_*.csv        (prevalence-filtered abundance matrices)
    ├── shrinkage_correlations_*.csv  (Schaefer-Strimmer correlation matrix)
    ├── raw_crispr_network_*.csv      (aggregated CRISPR interaction matrix)
    ├── smooth_crispr_*.csv           (taxonomy-smoothed CRISPR matrix)
    └── xstar_*.csv                   (X* residual propagation output)
```

