#!/usr/bin/env Rscript
# ---------------------------------------------------------------------------
# DecontX ambient-RNA removal, one sample.
#
# Why this exists alongside CellBender: CellBender is a deep generative model
# trained for `--epochs` passes over every droplet, which on a CPU queue costs
# roughly six hours per library and, on 4-lane libraries, has taken days. DecontX
# is a variational-EM model over the CELL matrix, so it does not depend on a GPU
# to finish in reasonable time. The 2026 seven-method benchmark puts CellBender,
# DecontX and SoupX in the same performing group and tells you to choose on
# compute environment -- which is exactly the choice this script implements.
#
# Output contract: a CELLRANGER-layout .h5, because the QC notebook reads it with
# Seurat::Read10X_h5 and that reader treats every root-level group as a genome.
# write10xCounts(version = "3") writes the single `matrix` group it expects. The
# file is re-read and its dimensions checked before this script exits -- an
# unreadable matrix must fail HERE, not sixty lines into someone's QC render.
#
# Usage: decontx_run.R <filtered.h5> <sample_id> [background.h5|NO_FILE] [seed] [threads]
# ---------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(celda)
  library(DropletUtils)
  library(SingleCellExperiment)
  library(Matrix)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2L)
  stop("usage: decontx_run.R <filtered.h5> <sample_id> [background.h5|NO_FILE] [seed] [threads]")

filtered   <- args[[1]]
sample_id  <- args[[2]]
background <- if (length(args) >= 3L) args[[3]] else "NO_FILE"
seed       <- if (length(args) >= 4L) as.integer(args[[4]]) else 12345L
threads    <- if (length(args) >= 5L) as.integer(args[[5]]) else 1L

msg <- function(...) cat(sprintf(...), "\n", sep = "")

# ---- input ----------------------------------------------------------------
# read10xCounts rather than Read10X_h5: it keeps the Ensembl ID *and* the symbol
# in rowData, and both have to survive into the output file. Reading with
# Read10X_h5 collapses them to whichever one `use.names` picks, and the gene IDs
# would then be lost for every stage after QC.
sce <- DropletUtils::read10xCounts(filtered, col.names = TRUE)
msg("%s: read %d genes x %d cells", sample_id, nrow(sce), ncol(sce))

# read10xCounts returns an HDF5-BACKED TENxMatrix for an .h5 input, and DecontX
# then returns its result as one too -- whereupon round()/drop0() die with
# `invalid class "DelayedMatrix" in 'R_Matrix_as_Csparse'`, AFTER the EM has
# run. Materialise the cell matrix now: it is the filtered matrix, so it is
# small, and the EM is faster in memory. The background matrix is deliberately
# left delayed -- it holds every droplet and can be millions of columns.
SummarizedExperiment::assay(sce, "counts") <-
  as(SummarizedExperiment::assay(sce, "counts"), "dgCMatrix")

# ---- optional ambient profile from empty droplets --------------------------
# DecontX does not require empty droplets -- that is its advantage over
# CellBender and SoupX -- but it uses them when given, and a cellranger `outs`
# directory always has the raw matrix. Any failure here degrades to the
# self-contained estimate rather than taking the sample down.
bg <- NULL
if (!grepl("NO_FILE", background, fixed = TRUE) && file.exists(background)) {
  bg <- tryCatch({
    bg_sce <- DropletUtils::read10xCounts(background, col.names = TRUE)
    if (!identical(rownames(bg_sce), rownames(sce))) {
      common <- intersect(rownames(bg_sce), rownames(sce))
      if (length(common) < 0.5 * nrow(sce))
        stop(sprintf("background shares only %d of %d features", length(common), nrow(sce)))
      bg_sce <- bg_sce[rownames(sce), ]
    }
    msg("%s: background profile from %d droplets", sample_id, ncol(bg_sce))
    bg_sce
  }, error = function(e) {
    warning(sprintf("%s: ignoring background matrix (%s)", sample_id, conditionMessage(e)))
    NULL
  })
}

if (threads > 1L) {
  try({
    BiocParallel::register(BiocParallel::MulticoreParam(workers = threads))
  }, silent = TRUE)
}

# ---- decontaminate ---------------------------------------------------------
# z is left NULL: DecontX derives its own cell populations. Supplying Seurat
# clusters would mean clustering before QC filtering, which is the wrong order
# and a second thing to go wrong per sample.
t0  <- Sys.time()
sce <- celda::decontX(sce, background = bg, seed = seed)
secs <- as.numeric(difftime(Sys.time(), t0, units = "secs"))

contamination <- as.numeric(SummarizedExperiment::colData(sce)$decontX_contamination)
msg("%s: decontX finished in %.1f min, mean contamination %.1f%% (median %.1f%%)",
    sample_id, secs / 60, 100 * mean(contamination, na.rm = TRUE),
    100 * median(contamination, na.rm = TRUE))

# ---- corrected counts ------------------------------------------------------
# decontXcounts is non-integer by construction. Round, then drop0: without it the
# zeros stay as EXPLICIT entries and the file carries millions of stored zeros,
# which inflates it and slows every reader downstream.
corrected <- celda::decontXcounts(sce)
# Belt and braces: if a future celda hands back a delayed matrix anyway, coerce
# rather than fail after the expensive part is already done.
if (!inherits(corrected, "dgCMatrix")) {
  corrected <- tryCatch(as(corrected, "dgCMatrix"),
                        error = function(e) as(as(corrected, "sparseMatrix"), "dgCMatrix"))
}
corrected <- Matrix::drop0(round(corrected))

# A cell whose counts round away entirely cannot be analysed and makes Seurat
# fail later with an unrelated-looking error. Drop it here and say so.
depth    <- Matrix::colSums(corrected)
n_empty  <- sum(depth == 0)
if (n_empty > 0) {
  msg("%s: dropping %d cell(s) left with zero counts after correction", sample_id, n_empty)
  keep          <- depth > 0
  corrected     <- corrected[, keep, drop = FALSE]
  sce           <- sce[, keep]
  contamination <- contamination[keep]
}

out_h5 <- sprintf("decontx_%s_filtered.h5", sample_id)
if (file.exists(out_h5)) unlink(out_h5)
DropletUtils::write10xCounts(
  path        = out_h5,
  x           = corrected,
  barcodes    = colnames(sce),
  gene.id     = as.character(SummarizedExperiment::rowData(sce)$ID),
  gene.symbol = as.character(SummarizedExperiment::rowData(sce)$Symbol),
  type        = "HDF5",
  version     = "3",
  overwrite   = TRUE
)

# ---- prove the output is readable by the consumer that matters -------------
chk <- Seurat::Read10X_h5(out_h5)
if (!identical(dim(chk), dim(corrected)))
  stop(sprintf("wrote %s but Read10X_h5 returns %d x %d, expected %d x %d",
               out_h5, nrow(chk), ncol(chk), nrow(corrected), ncol(corrected)))
msg("%s: wrote %s (%d genes x %d cells), verified with Read10X_h5",
    sample_id, out_h5, nrow(chk), ncol(chk))

# ---- per-cell contamination and a run summary ------------------------------
write.csv(data.frame(barcode = colnames(sce), contamination = contamination),
          sprintf("decontx_%s_contamination.csv", sample_id), row.names = FALSE)

cat(sprintf(paste0('{\n  "sample": "%s",\n  "method": "decontx",\n',
                   '  "genes": %d,\n  "cells_in": %d,\n  "cells_out": %d,\n',
                   '  "cells_dropped_empty": %d,\n',
                   '  "mean_contamination": %.4f,\n  "median_contamination": %.4f,\n',
                   '  "background_used": %s,\n  "seed": %d,\n  "seconds": %.1f\n}\n'),
            sample_id, nrow(corrected), ncol(corrected) + n_empty, ncol(corrected),
            n_empty, mean(contamination, na.rm = TRUE), median(contamination, na.rm = TRUE),
            if (is.null(bg)) "false" else "true", seed, secs),
    file = sprintf("decontx_%s_summary.json", sample_id))
