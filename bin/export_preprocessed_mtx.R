#!/usr/bin/env Rscript
# ---------------------------------------------------------------------------
# Export a preprocessed per-sample matrix from .rds to MatrixMarket + labels.
#
# This is the language-neutral seam between the R metaprogram preprocessing
# (Metaprog_NMFprepprocessing.qmd, which stays authoritative for the science)
# and the Python NMF engine. It deliberately performs NO transformation: the
# matrix is written out exactly as the R step produced it — log-CPM/10,
# gene-centred, clipped at zero — so the two engines factorise identical input
# and any difference in output is attributable to the solver alone.
#
# Usage: export_preprocessed_mtx.R <in.rds> <out_prefix>
# Writes: <out_prefix>.mtx, <out_prefix>_genes.tsv, <out_prefix>_cells.tsv
# ---------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(Matrix)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2L) stop("usage: export_preprocessed_mtx.R <in.rds> <out_prefix>")

in_path    <- args[[1]]
out_prefix <- args[[2]]

mat <- readRDS(in_path)

# Accept base matrix or any Matrix subclass; MatrixMarket needs a sparse form.
if (!inherits(mat, "Matrix") && !is.matrix(mat)) {
  stop(sprintf("Expected matrix/Matrix, got: %s", paste(class(mat), collapse = ", ")))
}
if (!inherits(mat, "dgCMatrix")) {
  mat <- as(as(as(mat, "dMatrix"), "generalMatrix"), "CsparseMatrix")
}

if (is.null(rownames(mat))) rownames(mat) <- paste0("g", seq_len(nrow(mat)))
if (is.null(colnames(mat))) colnames(mat) <- paste0("c", seq_len(ncol(mat)))

if (any(mat@x < 0)) {
  stop("Preprocessed matrix contains negative values; NMF requires a nonnegative input.")
}

Matrix::writeMM(mat, file = paste0(out_prefix, ".mtx"))
writeLines(as.character(rownames(mat)), paste0(out_prefix, "_genes.tsv"))
writeLines(as.character(colnames(mat)), paste0(out_prefix, "_cells.tsv"))

cat(sprintf("Exported %s: %d genes x %d cells -> %s.mtx\n",
            basename(in_path), nrow(mat), ncol(mat), out_prefix))
