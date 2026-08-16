#!/usr/bin/env Rscript
# ---------------------------------------------------------------------------
# Split the annotated object into the per-sample inputs Numbat needs.
#
# Numbat is run per sample. Two things have to be exactly right or the results
# cannot be joined back:
#
#   1. The BAM carries RAW cellranger barcodes ("AACGT...-1"), while the merged
#      Seurat object carries whatever Seurat's merge produced, typically a
#      prefixed or suffixed variant. This script writes an explicit
#      barcode -> cell_id map per sample so the join back is exact rather than
#      inferred from a string pattern.
#   2. The count matrix must contain the cells that survived QC, not the raw
#      cellranger set, otherwise Numbat models droplets that were already
#      discarded.
#
# Usage: numbat_prep.R <annotated.rds> <sample_col> <outdir>
# Writes per sample: <s>_counts.rds, <s>_barcodes.tsv, <s>_cellmap.tsv
# ---------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
  library(data.table)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3L) stop("usage: numbat_prep.R <annotated.rds> <sample_col> <outdir>")

obj_path   <- args[[1]]
sample_col <- args[[2]]
outdir     <- args[[3]]
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

obj <- readRDS(obj_path)
meta <- obj@meta.data

candidates <- c(sample_col, "sample_id", "orig.ident", "sample", "patient_id")
col <- candidates[candidates %in% colnames(meta)][1]
if (is.na(col)) stop("No sample column found; tried: ", paste(candidates, collapse = ", "))
message("Sample column: ", col)

counts <- tryCatch(
  SeuratObject::LayerData(obj, assay = DefaultAssay(obj), layer = "counts"),
  error = function(e) GetAssayData(obj, assay = DefaultAssay(obj), slot = "counts")
)

# Recover the raw 16bp+suffix cellranger barcode from whatever the merge produced.
# Handles "PREFIX_AACGT-1", "AACGT-1_1", and a bare "AACGT-1".
raw_barcode <- function(x) {
  b <- sub("^.*[_:]([ACGTN]{16}(-\\d+)?)$", "\\1", x)          # prefixed
  b <- sub("^([ACGTN]{16})(-\\d+)?_\\d+$", "\\1\\2", b)         # suffixed
  ifelse(grepl("^[ACGTN]{16}(-\\d+)?$", b), b, NA_character_)
}

samples <- unique(as.character(meta[[col]]))
message("Samples: ", length(samples))

written <- character(0)
for (s in samples) {
  cells <- rownames(meta)[as.character(meta[[col]]) == s]
  if (length(cells) < 20) {
    message(sprintf("[skip] %s: only %d cells", s, length(cells)))
    next
  }

  bc <- raw_barcode(cells)
  ok <- !is.na(bc)
  if (sum(ok) < 20) {
    message(sprintf("[skip] %s: could not recover cellranger barcodes for %d/%d cells",
                    s, sum(!ok), length(cells)))
    next
  }
  if (any(!ok)) {
    message(sprintf("[warn] %s: dropping %d cells with unrecognised barcode format",
                    s, sum(!ok)))
  }

  cells <- cells[ok]; bc <- bc[ok]
  dup <- duplicated(bc)
  if (any(dup)) {
    message(sprintf("[warn] %s: %d duplicate raw barcodes dropped", s, sum(dup)))
    cells <- cells[!dup]; bc <- bc[!dup]
  }

  mat <- counts[, cells, drop = FALSE]
  colnames(mat) <- bc            # Numbat matches the BAM's barcode space

  saveRDS(mat, file.path(outdir, paste0(s, "_counts.rds")))
  writeLines(bc, file.path(outdir, paste0(s, "_barcodes.tsv")))
  fwrite(data.table(cell_id = cells, barcode = bc, sample = s),
         file.path(outdir, paste0(s, "_cellmap.tsv")), sep = "\t")

  written <- c(written, s)
  message(sprintf("[ok] %s: %d genes x %d cells", s, nrow(mat), ncol(mat)))
}

if (!length(written)) stop("No sample produced usable Numbat input.")
writeLines(written, file.path(outdir, "samples.txt"))
