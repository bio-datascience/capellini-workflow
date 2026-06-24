#!/usr/bin/env Rscript

library(dada2)
library(phyloseq)
library(data.table)
library(NetCoMi)

# ---- reproducibility ----
# Single seed reused before every stochastic step.
# Change here if you want to verify run-to-run stability.
SEED <- 42L
set.seed(SEED)

# ---- get command line arguments ----
args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 6) {
  stop(
    "Usage: Rscript dada2_pipeline.R <input_path> <output_path> <ref_path> <taxmap_file> <direction> <fasta_generation> [chimera_removal]\n",
    "Example: Rscript dada2_pipeline.R /path/to/fastq /path/to/out /path/to/ref.fa.gz /path/to/taxmap.txt paired TRUE FALSE\n"
  )
}

input_path       <- args[1]
output_path      <- args[2]
ref_path         <- args[3]
taxmap_file      <- args[4]
direction        <- args[5]                  # "forward" | "reverse" | "paired"
fasta_generation <- as.logical(args[6])      # "TRUE"/"FALSE" -> TRUE/FALSE

# 7th arg (optional, default FALSE) — toggle removeBimeraDenovo. See
# tests/test_chimera_removal.R for the rationale: enabling chimera removal
# drops ~70% of detected ASVs on gut 16S, which is methodologically correct
# (they're PCR artefacts) but is disabled by default to match the validated
# legacy output.
chimera_removal  <- if (length(args) >= 7) as.logical(args[7]) else FALSE
if (is.na(chimera_removal)) chimera_removal <- FALSE
message("chimera_removal flag: ", chimera_removal)

if (!dir.exists(output_path)) dir.create(output_path, recursive = TRUE)

# ---- helpers ----
find_fastqs <- function(path) {
  if (dir.exists(path)) {
    files <- list.files(path, pattern="\\.(fastq|fq)(\\.gz)?$", full.names = TRUE, ignore.case = TRUE)
  } else if (file.exists(path) && grepl("\\.(fastq|fq)(\\.gz)?$", path, ignore.case = TRUE)) {
    files <- normalizePath(path)
  } else {
    stop("No FASTQ/FQ files found. Check 'input_path'.\nGiven path: ", path)
  }
  if (length(files) == 0) stop("No FASTQ/FQ files matched in: ", path)
  files
}

basename_noext <- function(x) sub("\\.(fastq|fq)(\\.gz)?$", "", basename(x), ignore.case = TRUE)

# ---- add NCBI taxid from SILVA taxonomy file (tax_slv_ssu_138.1.txt) ----
# NOTA: usa tax_slv_ssu_138.1.txt (path -> NCBI taxid), NON taxmap_slv_ssu_ref_nr_138.1.txt
# Il taxmap contiene taxid interni SILVA, non taxid NCBI.
# Il file tax_slv_ssu ha colonne: path, taxid (NCBI reale), rank
add_taxid_from_silva_taxmap <- function(tax_df, taxmap_file,
                                        ranks = c("Kingdom","Phylum","Class","Order","Family","Genus")) {
  stopifnot(file.exists(taxmap_file))
  
  ranks_present <- intersect(ranks, colnames(tax_df))
  if (length(ranks_present) == 0) stop("Nessuno dei rank attesi trovato in tax_df.")
  
  # Leggi il file tax_slv_ssu (NON taxmap): colonne path, taxid, rank (senza header)
  dt <- data.table::fread(taxmap_file, sep = "\t", header = FALSE, quote = "",
                          col.names = c("path", "taxid", "rank", "remark", "release"),
                          fill = TRUE)
  dt[, path  := trimws(path)]
  dt[, taxid := suppressWarnings(as.integer(taxid))]
  
  # Tieni solo un taxid per path (già unici in questo file, ma per sicurezza)
  dt_unique <- dt[!is.na(path) & path != "" & !is.na(taxid),
                  .(taxid = taxid[1], rank = rank[1]), by = path]
  data.table::setkey(dt_unique, path)
  
  # Costruttore di path SILVA dal formato DADA2:
  # DADA2 assegna Kingdom="Bacteria" ma SILVA usa "Bacteria" come domain.
  # Il path SILVA include il nome del gruppo finale seguito da ";"
  # es: "Bacteria;Firmicutes;Clostridia;Lachnospirales;Lachnospiraceae;Blautia;"
  make_silva_path <- function(vals) {
    vals <- as.character(vals)
    vals <- vals[!is.na(vals) & vals != "" & vals != "NA"]
    if (length(vals) == 0) return(NA_character_)
    paste0(paste(vals, collapse = ";"), ";")
  }
  
  n_rows <- nrow(tax_df)
  taxids        <- rep(NA_integer_, n_rows)
  matched_ranks <- rep(NA_character_, n_rows)
  
  for (i in seq_len(n_rows)) {
    row_vals <- as.character(tax_df[i, ranks_present, drop = TRUE])
    
    # Tenta il match dal rank più granulare verso il più alto.
    # Accetta un rank superiore solo se il rank inferiore è genuinamente mancante.
    for (k in seq_along(ranks_present)) {
      n_keep   <- length(ranks_present) - k + 1
      sub_vals <- row_vals[seq_len(n_keep)]
      path_try <- make_silva_path(sub_vals)
      if (is.na(path_try)) next
      
      hit <- dt_unique[path_try, on = "path", nomatch = NULL]
      if (nrow(hit) == 0) next
      
      # Se k > 1 (match su rank superiore al più granulare disponibile):
      # accetta solo se i rank più granulari sono tutti NA/vuoti
      if (k > 1) {
        deeper_vals <- row_vals[(n_keep + 1):length(ranks_present)]
        all_missing <- all(is.na(deeper_vals) | deeper_vals == "" | deeper_vals == "NA")
        if (!all_missing) {
          # L'ASV ha un genus (o rank più granulare) definito ma non ha matchato ->
          # non assegnare il taxid di un rank superiore (sarebbe sbagliato)
          break
        }
      }
      
      taxids[i]        <- hit$taxid[1]
      matched_ranks[i] <- hit$rank[1]
      break
    }
  }
  
  tax_df$NCBI_taxid        <- taxids
  tax_df$taxid_matched_rank <- matched_ranks  # colonna diagnostica
  tax_df
}

# ---- locate input fastqs ----
all_fastqs <- find_fastqs(input_path)

seqtab <- NULL
sample.names <- NULL
dir_suffix <- NULL

# ---- DADA2 pipeline (single-end / paired-end) ----
if (tolower(direction) %in% c("forward", "reverse")) {
  
  # ---------- SINGLE-END ----------
  dir_token <- switch(tolower(direction),
                      "forward" = "R1",
                      "reverse" = "R2")
  
  keep_idx <- grepl(paste0("_", dir_token, "(_|\\.)"), basename(all_fastqs), ignore.case = TRUE)
  fnSel <- sort(all_fastqs[keep_idx])
  
  if (length(fnSel) == 0) {
    warning("R1/R2 direction descriptors not found in filenames; using all files as single-end.")
    fnSel <- sort(all_fastqs)
  }
  
  message("Processing ", direction, " reads (", dir_token, "): ", length(fnSel), " files.")
  
  sample.names <- sub("\\..*$", "", basename(fnSel))
  dir_suffix <- ifelse(tolower(direction) == "forward", "F", "R")
  
  filt_dir <- file.path(output_path, "filtered")
  if (!dir.exists(filt_dir)) dir.create(filt_dir, recursive = TRUE)
  
  filtSel <- file.path(filt_dir, paste0(sample.names, "_", dir_suffix, "_filt.fastq.gz"))
  names(filtSel) <- sample.names
  
  out <- filterAndTrim(
    fnSel, filtSel,
    truncLen = 240,
    trimLeft = 1,
    maxN = 0,
    maxEE = 2,
    truncQ = 2,
    rm.phix = TRUE,
    compress = TRUE,
    multithread = TRUE
  )
  
  if (nrow(out) == 0) stop("filterAndTrim returned no rows; check input files.")
  if (all(out[, "reads.out"] == 0)) stop("All reads were filtered out. Adjust truncLen/maxEE or verify reads quality.")
  
  set.seed(SEED)
  err <- learnErrors(filtSel, multithread = TRUE)
  set.seed(SEED)
  dada_res <- dada(filtSel, err = err, multithread = TRUE)
  seqtab <- makeSequenceTable(dada_res)

  # Chimera removal — gated by the YAML config flag `chimera_removal`.
  # Enabling typically drops ~70% of detected ASVs on gut 16S (they're PCR
  # bimeras; see tests/test_chimera_removal.R).
  if (isTRUE(chimera_removal)) {
    message("Running removeBimeraDenovo (single-end branch)")
    seqtab <- removeBimeraDenovo(seqtab, method = "consensus", multithread = TRUE, verbose = TRUE)
  } else {
    message("Skipping removeBimeraDenovo (chimera_removal = FALSE)")
  }

} else if (tolower(direction) == "paired") {
  
  # ---------- PAIRED-END ----------
  f_idx <- grepl("_R1(_|\\.)", basename(all_fastqs), ignore.case = TRUE)
  r_idx <- grepl("_R2(_|\\.)", basename(all_fastqs), ignore.case = TRUE)
  fnFs <- sort(all_fastqs[f_idx])
  fnRs <- sort(all_fastqs[r_idx])
  
  if (length(fnFs) == 0 || length(fnRs) == 0) {
    stop("Paired-end mode: couldn't find both R1 and R2 files (need filenames with '_R1_' and '_R2_').")
  }
  
  sampF <- sub("_R1.*$", "", basename_noext(fnFs))
  sampR <- sub("_R2.*$", "", basename_noext(fnRs))
  common <- intersect(sampF, sampR)
  if (length(common) == 0) stop("Paired-end mode: no matching R1/R2 samples after name harmonization.")
  
  fnFs <- fnFs[match(common, sampF)]
  fnRs <- fnRs[match(common, sampR)]
  if (any(is.na(fnFs)) || any(is.na(fnRs))) stop("Paired-end mode: failed to align R1/R2 by sample names.")
  
  sample.names <- common
  dir_suffix <- "P"
  
  filt_dir <- file.path(output_path, "filtered")
  if (!dir.exists(filt_dir)) dir.create(filt_dir, recursive = TRUE)
  
  filtFs <- file.path(filt_dir, paste0(sample.names, "_F_filt.fastq.gz"))
  filtRs <- file.path(filt_dir, paste0(sample.names, "_R_filt.fastq.gz"))
  names(filtFs) <- sample.names
  names(filtRs) <- sample.names
  
  out <- filterAndTrim(
    fnFs, filtFs,
    fnRs, filtRs,
    truncLen = c(240, 240),
    trimLeft = c(1, 1),
    maxN = 0,
    maxEE = c(2, 2),
    truncQ = 2,
    rm.phix = TRUE,
    compress = TRUE,
    multithread = TRUE
  )
  
  if (nrow(out) == 0) stop("filterAndTrim (paired) returned no rows; check input files.")
  if (all(out[, "reads.out"] == 0)) stop("All reads were filtered out in paired mode.")
  
  set.seed(SEED)
  errF <- learnErrors(filtFs, multithread = TRUE)
  set.seed(SEED)
  errR <- learnErrors(filtRs, multithread = TRUE)
  
  set.seed(SEED)
  dadaFs <- dada(filtFs, err = errF, multithread = TRUE)
  set.seed(SEED)
  dadaRs <- dada(filtRs, err = errR, multithread = TRUE)
  
  mergers <- mergePairs(dadaFs, filtFs, dadaRs, filtRs, verbose = TRUE)
  seqtab <- makeSequenceTable(mergers)

  # Chimera removal — gated by the YAML config flag `chimera_removal`.
  # Enabling typically drops ~70% of detected ASVs on gut 16S (they're PCR
  # bimeras; see tests/test_chimera_removal.R).
  if (isTRUE(chimera_removal)) {
    message("Running removeBimeraDenovo (paired-end branch)")
    seqtab <- removeBimeraDenovo(seqtab, method = "consensus", multithread = TRUE, verbose = TRUE)
  } else {
    message("Skipping removeBimeraDenovo (chimera_removal = FALSE)")
  }

} else {
  stop("direction must be one of: 'forward', 'reverse', or 'paired'")
}

# ---------- Taxonomy assignment + phyloseq (ASV-level) ----------
# minBoot is the bootstrap confidence threshold (out of 100) for rank assignment.
# dada2's default is 50; we use 60 to match the original capellini package and
# the notebook output that has been validated to run end-to-end through the
# pipeline (including the SpacePHARER findpam step). Raising this value to 80
# was tried but produces ~4x fewer ASVs (5.5k vs 20k on the IBD test set),
# which in turn yields a much smaller target_spacers.fasta and triggers an
# intermittent SpacePHARER findpam Bus error on macOS for this corpus.
set.seed(SEED)
taxa <- assignTaxonomy(seqtab, ref_path, minBoot = 60, multithread = TRUE)

sample_data_df <- data.frame(SampleID = sample.names, SampleType = "Type1", row.names = sample.names)
physeq <- phyloseq(
  otu_table(seqtab, taxa_are_rows = FALSE),  # samples x ASV
  tax_table(taxa),                           # ASV x ranks
  sample_data(sample_data_df)
)

# ---------- taxonomy_table: 1 ASV per riga + NCBI taxid ----------
tax_table_df <- as.data.frame(tax_table(physeq))
tax_table_df <- add_taxid_from_silva_taxmap(tax_table_df, taxmap_file)

# aggiorna phyloseq con taxid dentro
tax_table(physeq) <- tax_table(as.matrix(tax_table_df))

cat("TaxID matched (ASV-level): ",
    sum(!is.na(tax_table_df$NCBI_taxid)), "/", nrow(tax_table_df),
    " (", round(mean(!is.na(tax_table_df$NCBI_taxid))*100, 1), "%)\n")

# ---------- OTU_table (ASV-level) ----------
# campioni x ASV (colonne = sequenze ASV)
OTU_table <- as.data.frame(otu_table(physeq))
write.csv(
  OTU_table,
  file.path(output_path, paste0("OTU_table_", dir_suffix, ".csv")),
  row.names = TRUE
)

# ---------- taxonomy_table (ASV-level, con Genus e taxid) ----------
# 1 riga = 1 ASV; ASV nello stesso Genus avranno gli stessi valori Kingdom..Genus (e spesso stesso NCBI_taxid)
write.csv(
  tax_table_df,
  file.path(output_path, paste0("taxonomy_table_", dir_suffix, ".csv")),
  row.names = TRUE
)

# ---------- Genus_Table (abbondanze aggregate per Genus) ----------
# Qui integriamo renameTaxa() SOLO per gestire bene gli unclassified a Genus,
# senza toccare i nomi ASV (che restano le sequenze)

physeq_clean <- renameTaxa(
  physeq,
  pat = "<name>",
  unclass = c("unclassified", "Unclassified", "Incertae Sedis")
)

tax_clean_df <- as.data.frame(tax_table(physeq_clean))

# genus labels “pulite” per aggregazione
genus_labels <- as.character(tax_clean_df$Genus)
genus_labels[is.na(genus_labels) | genus_labels == ""] <- "unclassified"

# Mappa ASV -> Genus_label (usando rownames = ASV seq)
asv_ids <- colnames(OTU_table)
genus_map <- setNames(genus_labels[match(asv_ids, rownames(tax_clean_df))], asv_ids)
genus_map[is.na(genus_map) | genus_map == ""] <- "unclassified"

# Somma ASV per Genus (senza perdere OTU_table ASV-level)
Genus_Table <- t(rowsum(t(as.matrix(OTU_table)), group = genus_map))

# Rendi i nomi Genus univoci se necessario
colnames(Genus_Table) <- make.unique(colnames(Genus_Table))

write.csv(
  Genus_Table,
  file.path(output_path, paste0("Genus_Table_", dir_suffix, ".csv")),
  row.names = TRUE
)

# ---------- salva phyloseq ASV-level con taxid dentro ----------
saveRDS(physeq, file = file.path(output_path, paste0("phyloseq_ASV_", dir_suffix, ".rds")))

# ---------- FASTA export (ASV sequences + NCBI taxid/genus in header) ----------
if (isTRUE(fasta_generation)) {
  seqs <- taxa_names(physeq)  # sequenze ASV reali (taxa_names = ASV sequences)
  if (length(seqs) == 0L) {
    warning("FASTA export skipped: no taxa_names found in physeq.")
  } else {
    fasta_file <- file.path(output_path, paste0("ASV_sequences_", dir_suffix, ".fasta"))
    con <- file(fasta_file, open = "w")
    on.exit(close(con), add = TRUE)
    
    # tax_table_df ha rownames = sequenze ASV
    taxid_map <- setNames(as.character(tax_table_df$NCBI_taxid), rownames(tax_table_df))
    genus_map_hdr <- setNames(as.character(tax_table_df$Genus), rownames(tax_table_df))
    
    for (i in seq_along(seqs)) {
      s <- seqs[i]
      taxid_i <- taxid_map[[s]]
      genus_i <- genus_map_hdr[[s]]
      
      if (is.null(taxid_i) || is.na(taxid_i) || taxid_i == "") taxid_i <- "NA"
      if (is.null(genus_i) || is.na(genus_i) || genus_i == "") genus_i <- "NA"
      
      # sanitizza genus per header FASTA
      genus_i <- gsub("\\s+", "_", genus_i)
      
      # header con ID ASV + taxid + genus
      writeLines(paste0(">ASV_", i, " taxid=", taxid_i, " genus=", genus_i), con)
      writeLines(s, con)
    }
    
    message("FASTA written (with taxid/genus headers): ", fasta_file)
  }
}

# ---------- quick sanity previews ----------
cat("\nPreview OTU_table (ASV-level):\n")
print(utils::head(OTU_table[, seq_len(min(6, ncol(OTU_table))), drop = FALSE]))

cat("\nPreview Genus_Table (Genus-level aggregated):\n")
print(utils::head(Genus_Table[, seq_len(min(6, ncol(Genus_Table))), drop = FALSE]))

cat("\nPreview taxonomy_table (ASV-level):\n")
print(utils::head(tax_table_df))