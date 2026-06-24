#!/usr/bin/env Rscript
# =============================================================================
# Standalone test: measure the ASV reduction caused by removeBimeraDenovo
# =============================================================================
#
# Reuses the *already-filtered* fastq.gz files in DADA2 output/filtered/ so it
# skips the slow filterAndTrim step. Runs learnErrors + dada + sequence-table
# construction, then calls isBimeraDenovo (which only *flags* bimeras without
# removing them) so we can report both the pre- and post-chimera ASV counts in
# a single run.
#
# Reuses the same SEED as dada2_pipe.R so the dada results are byte-identical
# to your main pipeline run.
#
# Usage:
#   Rscript test_chimera_removal.R <path_to_DADA2_output_dir>
#
# Example:
#   Rscript test_chimera_removal.R "/Users/.../1_IBD_package/DADA2 output"
#
# Expected runtime: a few minutes (learnErrors + dada on 166 filtered samples).
# =============================================================================

suppressPackageStartupMessages({
  library(dada2)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) {
  stop("Usage: Rscript test_chimera_removal.R <path_to_DADA2_output_dir>")
}

dada_out_dir <- args[1]
filt_dir <- file.path(dada_out_dir, "filtered")

if (!dir.exists(filt_dir)) {
  stop(sprintf("Expected filtered/ subdir not found: %s\n  Run the main pipeline at least once first so DADA2 produces the filtered fastqs.",
               filt_dir))
}

# Match the SEED used in dada2_pipe.R
SEED <- 42L

# Find filtered fastq.gz files
filt_files <- list.files(filt_dir, pattern = "\\.fastq(\\.gz)?$",
                         full.names = TRUE, ignore.case = TRUE)
filt_files <- sort(filt_files)
if (length(filt_files) == 0) {
  stop("No fastq.gz files found in: ", filt_dir)
}

# Derive sample.names from filtered filenames: "ERR1464217_F_filt.fastq.gz" -> "ERR1464217"
sample.names <- sub("_(F|R|P)_filt\\.fastq(\\.gz)?$", "", basename(filt_files),
                    ignore.case = TRUE)
names(filt_files) <- sample.names

cat("================================================================\n")
cat("Chimera-removal effect test\n")
cat("================================================================\n")
cat(sprintf("Filtered files:  %d  (in %s)\n", length(filt_files), filt_dir))
cat(sprintf("Seed:            %d\n\n", SEED))

cat("[1/4] learnErrors ...\n")
set.seed(SEED)
err <- learnErrors(filt_files, multithread = TRUE)

cat("\n[2/4] dada (denoising) ...\n")
set.seed(SEED)
dada_res <- dada(filt_files, err = err, multithread = TRUE)

cat("\n[3/4] makeSequenceTable ...\n")
seqtab <- makeSequenceTable(dada_res)
n_pre <- ncol(seqtab)

cat("\n[4/4] isBimeraDenovoTable (consensus flag, no removal) ...\n")
# isBimeraDenovoTable applies isBimeraDenovo per-sample and returns a
# consensus flag per ASV — exactly what removeBimeraDenovo(method="consensus")
# uses internally, but without dropping the flagged sequences.
bimeras <- isBimeraDenovoTable(seqtab, multithread = TRUE, verbose = TRUE)
n_bimeras <- sum(bimeras)
n_post <- n_pre - n_bimeras

cat("\n================================================================\n")
cat("RESULTS\n")
cat("================================================================\n")
cat(sprintf("ASVs detected by DADA2:           %6d\n", n_pre))
cat(sprintf("ASVs flagged as bimeras (removed): %6d  (%.1f%%)\n",
            n_bimeras, 100 * n_bimeras / n_pre))
cat(sprintf("ASVs kept after chimera removal:  %6d  (%.1f%%)\n",
            n_post, 100 * n_post / n_pre))
cat("================================================================\n\n")

# Cross-check against pipeline output (if present)
pipeline_otu <- file.path(dada_out_dir, "OTU_table_F.csv")
if (file.exists(pipeline_otu)) {
  hdr <- readLines(pipeline_otu, n = 1)
  pipeline_cols <- length(strsplit(hdr, ",", fixed = TRUE)[[1]]) - 1L
  cat(sprintf("Pipeline OTU_table_F.csv columns: %d  (should equal post-chimera ASV count above)\n",
              pipeline_cols))
  if (pipeline_cols == n_post) {
    cat("=> MATCHES post-chimera count. Chimera removal is the cause of the reduction.\n")
  } else {
    cat(sprintf("=> MISMATCH: pipeline has %d, post-chimera computed %d.\n",
                pipeline_cols, n_post))
    cat("   Something other than removeBimeraDenovo is also affecting the count.\n")
  }
}
