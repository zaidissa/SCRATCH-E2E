#!/usr/bin/env Rscript
# ---------------------------------------------------------------------------
# Split the annotated object into the per-GROUP inputs Numbat needs.
#
# A "group" is the individual, normally patient_id. It is NOT the library.
#
# Numbat is designed to run per individual across that individual's libraries --
# its own pileup_and_phase.R says so:
#
#     --label    Individual label. One per run.
#     --samples  Sample name(s); comma delimited if multiple.
#                All samples must belong to the same individual.
#
# Running it per library instead breaks it in a way that is invisible in the
# output and wrong in a specific, reproducible direction. numbat:::mark_tumor_lineage
# assigns the tumour compartment as
#
#     tumor_root = mut_nodes[which.max(mut_burdens)]
#     compartment = ifelse(seq > 0, "tumor", "normal")
#
# -- it takes whatever segments survived the LLR filter, puts them on a tree,
# picks the heaviest node and calls everything below it tumour. There is no
# significance test and no "this sample has no tumour" branch; which.max always
# returns something. So a library with no malignant cells does not come back
# empty, it comes back with its own noise split into tumour and normal. That is
# what produced neoplastic calls in PBMC libraries, which should be ~0% tumour.
#
# Pooled with the rest of the individual's libraries, the real aneuploid clones
# carry overwhelming mutation burden, tumor_root lands on the genuine tumour
# clade, and the normal cells sit in the normal compartment where they belong.
# Joint genotyping across the libraries also gives more coverage at heterozygous
# sites, which is the whole advantage Numbat has over the expression-only callers.
#
# Three things have to be exactly right or the results cannot be joined back:
#
#   1. BARCODES COLLIDE ACROSS LIBRARIES. Two libraries of one patient both
#      contain AACGT...-1. Pooling them naively silently merges two different
#      cells. Every cell is therefore namespaced as "<library>_<raw barcode>",
#      which is the same shape Numbat's own bundled example uses
#      (df_allele_example$cell is "ATC2_TCCTAATCACAACCGC-1").
#   2. The BAM carries RAW cellranger barcodes, the merged Seurat object carries
#      whatever Seurat's merge produced. The cellmap written here is the explicit
#      join, rather than a string pattern inferred downstream.
#   3. The count matrix must hold the cells that survived QC, not the raw
#      cellranger set, or Numbat models droplets that were already discarded.
#
# Also SELECTS the cells for a patient-internal expression reference: that
# individual's own immune and stromal cells. numbat::ref_hca is 12 generic
# profiles from a different study, and deviation between it and the real cohort
# is read as expression CNV signal. Per library there were often no normal cells
# to build from -- per individual, the PBMC/normal library is exactly that
# reference.
#
# The selection happens here because this is where the annotation labels and the
# shared label matcher are. The AGGREGATION happens in numbat_run.R, because
# numbat::aggregate_counts() needs the numbat package and this script runs in
# image_cnv, which has Seurat but no numbat (the numbat image is the mirror
# image: numbat 1.4.2, no Seurat at all). Verified in both images rather than
# assumed.
#
# Usage:
#   numbat_prep.R <annotated.rds> <sample_col> <group_col> <celltype_col> \
#                 <reference_labels> <min_ref_cells> <outdir> [helpers.R]
#
# Writes:
#   <library>_barcodes.tsv   raw barcodes, one file per library (pileup input)
#   <group>_counts.rds       pooled counts, colnames "<library>_<barcode>"
#   <group>_cellmap.tsv      cell_id, barcode, cell, library, group
#   <group>_refcells.tsv     cells for the internal expression reference
#   groups.tsv               library -> group map (drives the Nextflow grouping)
#   grouping_report.tsv      per-group summary, for the report
# ---------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
  library(data.table)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 7L)
  stop("usage: numbat_prep.R <annotated.rds> <sample_col> <group_col> <celltype_col> ",
       "<reference_labels> <min_ref_cells> <outdir> [helpers.R]")

obj_path      <- args[[1]]
sample_col_in <- args[[2]]
group_col_in  <- args[[3]]
ct_col_in     <- args[[4]]
ref_labels_in <- args[[5]]
min_ref_cells <- suppressWarnings(as.integer(args[[6]]))
outdir        <- args[[7]]
helpers       <- if (length(args) >= 8) args[[8]] else ""
if (is.na(min_ref_cells)) min_ref_cells <- 50L

dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

# The label matcher is SHARED with stratify and cnv_concordance rather than
# copied: the reference lineages that define the null there must be the same
# cells that define the expression reference here, and two copies of one rule
# diverge silently. Falls back to exact, case-insensitive matching if the
# helper file was not staged.
if (nzchar(helpers) && file.exists(helpers)) {
  source(helpers, local = FALSE)
}
match_labels <- if (exists("scratch_match_labels")) {
  scratch_match_labels
} else {
  message("[warn] scratch_report.R not staged; falling back to exact label matching.")
  function(observed, requested) tolower(trimws(as.character(observed))) %in%
                                tolower(trimws(as.character(requested)))
}

obj  <- readRDS(obj_path)
meta <- obj@meta.data

pick_col <- function(candidates, what) {
  hit <- candidates[candidates %in% colnames(meta)]
  if (!length(hit)) return(NA_character_)
  message(sprintf("%s column: %s", what, hit[1]))
  hit[1]
}

sample_col <- pick_col(c(sample_col_in, "sample_id", "orig.ident", "sample"), "Library")
if (is.na(sample_col))
  stop("No library column found; tried: ",
       paste(c(sample_col_in, "sample_id", "orig.ident", "sample"), collapse = ", "))

# patient_id is a REQUIRED samplesheet column when starting at align, and the
# merge carries it through -- but a cohort entered further down the DAG with a
# hand-made metadata table may not have it. Falling back to the library column
# reproduces the OLD per-library behaviour, which is the behaviour this script
# exists to stop, so it is reported as a warning rather than taken silently.
group_col <- pick_col(c(group_col_in, "patient_id", "patient", "subject_id", "donor_id"),
                      "Group (individual)")
grouping_is_per_library <- FALSE
if (is.na(group_col)) {
  group_col <- sample_col
  grouping_is_per_library <- TRUE
  message(sprintf(paste(
    "[WARN] No group column found (tried: %s). Falling back to the library",
    "column '%s', which runs Numbat PER LIBRARY -- the configuration that",
    "over-calls tumour in normal libraries. Add patient_id to the metadata, or",
    "set --numbat_group_col, to group libraries by individual."),
    paste(c(group_col_in, "patient_id", "patient", "subject_id", "donor_id"),
          collapse = ", "), sample_col))
}

ct_col <- pick_col(c(ct_col_in, "scratch_celltype", "azimuth_labels", "sctype",
                     "cell_type", "celltype"), "Cell type")

ref_labels <- trimws(strsplit(ref_labels_in, "[;,]")[[1]])
ref_labels <- ref_labels[nzchar(ref_labels)]

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

meta$.cell_id <- rownames(meta)
meta$.library <- as.character(meta[[sample_col]])
meta$.group   <- as.character(meta[[group_col]])
meta$.ct      <- if (is.na(ct_col)) NA_character_ else as.character(meta[[ct_col]])
meta$.barcode <- raw_barcode(meta$.cell_id)

bad <- is.na(meta$.barcode) | is.na(meta$.library) | is.na(meta$.group)
if (any(bad)) {
  message(sprintf("[warn] dropping %d of %d cells with no recoverable barcode, library or group",
                  sum(bad), nrow(meta)))
  meta <- meta[!bad, , drop = FALSE]
}

dt <- as.data.table(meta[, c(".cell_id", ".library", ".group", ".ct", ".barcode")])
setnames(dt, sub("^\\.", "", names(dt)))

# The namespaced cell id. This is what Numbat sees, in both the count matrix
# and df_allele, and it is unique across the whole cohort.
dt[, cell := paste0(library, "_", barcode)]

# A duplicate raw barcode WITHIN one library is a genuine ambiguity (two merged
# cells claiming one droplet); across libraries it is expected and now harmless.
dup <- dt[, .N, by = .(library, barcode)][N > 1L]
if (nrow(dup)) {
  message(sprintf("[warn] %d barcode(s) duplicated within a library; dropping all copies",
                  nrow(dup)))
  dt <- dt[!dup, on = .(library, barcode)]
}

groups <- sort(unique(dt$group))
message(sprintf("Libraries: %d | Groups: %d", uniqueN(dt$library), length(groups)))

# ---------------------------------------------------------------------------
# One barcode file per LIBRARY: cellsnp-lite walks one BAM at a time and wants
# that library's raw barcodes, even though the phasing is joint per individual.
# ---------------------------------------------------------------------------
for (lb in sort(unique(dt$library))) {
  writeLines(dt[library == lb, barcode],
             file.path(outdir, paste0(lb, "_barcodes.tsv")))
}

report  <- list()
written <- character(0)

for (g in groups) {
  sub <- dt[group == g]
  libs <- sort(unique(sub$library))

  if (nrow(sub) < 20L) {
    message(sprintf("[skip] %s: only %d cells", g, nrow(sub)))
    next
  }

  mat <- counts[, sub$cell_id, drop = FALSE]
  colnames(mat) <- sub$cell

  # ----- cells for the individual's own expression reference --------------
  # Reference lineages ONLY: including tumour cells would flatten the very
  # signal Numbat is looking for. Cell types with too few cells are dropped
  # here -- a profile built from a handful of cells is noise wearing a cell-type
  # name, and Numbat picks per pseudobulk among whatever columns it is given.
  ref_cells <- data.table(cell = character(), celltype = character())
  n_ref <- 0L
  if (!is.na(ct_col) && length(ref_labels)) {
    is_ref <- match_labels(sub$ct, ref_labels)
    is_ref[is.na(is_ref)] <- FALSE
    n_ref <- sum(is_ref)
    if (n_ref > 0L) {
      cand <- data.table(cell = sub$cell[is_ref], celltype = as.character(sub$ct[is_ref]))
      keep_ct <- cand[, .N, by = celltype][N >= 20L, celltype]
      ref_cells <- cand[celltype %in% keep_ct]
    }
  }
  ref_basis <- if (nrow(ref_cells) >= min_ref_cells) "internal" else "ref_hca"
  if (ref_basis == "internal") {
    message(sprintf("[ref] %s: %d reference-lineage cell(s) over %d type(s): %s",
                    g, nrow(ref_cells), uniqueN(ref_cells$celltype),
                    paste(sort(unique(ref_cells$celltype)), collapse = ", ")))
  } else {
    message(sprintf(paste("[ref] %s: only %d usable reference-lineage cell(s) (need %d);",
                          "numbat_run.R will fall back to numbat::ref_hca. Tumour calls",
                          "for this individual then rest on an external reference."),
                    g, nrow(ref_cells), min_ref_cells))
  }

  saveRDS(mat, file.path(outdir, paste0(g, "_counts.rds")))
  # Always written, possibly header-only: an optional per-group file would have
  # to be joined with `remainder`, and a group silently missing its reference
  # would fall back to ref_hca with nothing saying so.
  fwrite(ref_cells, file.path(outdir, paste0(g, "_refcells.tsv")), sep = "\t")
  fwrite(sub[, .(cell_id, barcode, cell, library, group)],
         file.path(outdir, paste0(g, "_cellmap.tsv")), sep = "\t")

  written <- c(written, g)
  report[[length(report) + 1L]] <- data.table(
    group = g, n_libraries = length(libs), libraries = paste(libs, collapse = ","),
    n_cells = nrow(sub), n_reference_cells = n_ref,
    n_reference_usable = nrow(ref_cells), reference_basis = ref_basis)
  message(sprintf("[ok] %s: %d genes x %d cells from %d librar%s (%s)",
                  g, nrow(mat), ncol(mat), length(libs),
                  if (length(libs) == 1L) "y" else "ies", paste(libs, collapse = ", ")))
}

if (!length(written)) stop("No group produced usable Numbat input.")

rep_dt <- rbindlist(report)
fwrite(rep_dt, file.path(outdir, "grouping_report.tsv"), sep = "\t")

# Drives the Nextflow grouping: which libraries get piled up together.
# Only libraries belonging to a group that was actually written.
fwrite(unique(dt[group %in% written, .(library, group)])[order(group, library)],
       file.path(outdir, "groups.tsv"), sep = "\t")
writeLines(written, file.path(outdir, "samples.txt"))

# ---------------------------------------------------------------------------
# Say plainly whether the grouping did anything. A cohort where every group has
# exactly one library is still running per library -- the exact configuration
# that over-calls tumour in normal libraries -- and the only symptom otherwise
# is a number in a report that nobody has a reason to doubt.
# ---------------------------------------------------------------------------
multi <- rep_dt[n_libraries > 1L]
cat("\n==== Numbat grouping ====\n")
print(rep_dt)
if (grouping_is_per_library || !nrow(multi)) {
  cat(sprintf(paste(
    "\n[WARN] Numbat is still running PER LIBRARY: %d group(s), none with more",
    "than one library%s.\n",
    "       Numbat cannot report \"no tumour here\" -- mark_tumor_lineage() always\n",
    "       labels the heaviest clade tumour -- so a library with no malignant\n",
    "       cells will still come back with neoplastic calls.\n",
    "       Group libraries by individual via patient_id or --numbat_group_col.\n"),
    nrow(rep_dt),
    if (grouping_is_per_library) sprintf(" (no group column; used '%s')", group_col) else ""))
} else {
  cat(sprintf("\n%d of %d group(s) pool more than one library. Grouped by '%s'.\n",
              nrow(multi), nrow(rep_dt), group_col))
}
no_internal <- rep_dt[reference_basis != "internal"]
if (nrow(no_internal)) {
  cat(sprintf(paste("[note] %d individual(s) had too few reference-lineage cells for an",
                    "internal expression reference and use numbat::ref_hca: %s\n"),
              nrow(no_internal), paste(no_internal$group, collapse = ", ")))
}
