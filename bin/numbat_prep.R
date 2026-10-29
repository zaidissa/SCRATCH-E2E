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
# REFERENCE COLUMNS ARE CLUSTERS, NOT CELL TYPE LABELS.
#
# `lambdas_ref` is a gene x pseudobulk matrix and Numbat picks the best-fitting
# COLUMN per pseudobulk; the column names are labels for us, not for it. The
# docs say as much ("columns are pseudobulk names"), and the author's own worked
# example (vignettes/mouse.Rmd section 4) builds the reference from raw
# unsupervised clusters with no annotation at all:
#
#     clusters  = p2$clusters$PCA$multilevel
#     ref_annot = data.frame(cell = names(clusters), group = unname(clusters))
#     ref       = aggregate_counts(count_mat_ref,
#                                  ref_annot %>% group_by(group) %>% filter(n() > 100))
#
# So what a reference column needs is transcriptional homogeneity and enough
# cells -- not a name. scType assigns ONE label per cluster (see
# notebook_sctype_major.qmd: cluster-summed scores, top_n(1), stamped on every
# cell in the cluster), so aggregating by label COLLAPSES several clusters into
# one averaged column and throws away resolution we already paid for. Grouping
# by <label>__c<cluster> instead keeps the same eligibility rule and hands
# Numbat finer columns to match against.
#
# It also makes the per-column cell floor mean something. At label granularity
# every column has hundreds to thousands of cells and the floor never fires; at
# cluster granularity it is the actual quality control, which is why it is 100
# here (the author's number) rather than 20.
#
# The eligibility rule is unchanged and deliberately conservative -- the twelve
# reference lineages, which exclude Malignant and the four labels that Neftel
# malignant states also carry (Astrocyte, OPC, Radial_Glia, Neuron). Numbat's
# maintainer is explicit that "no tumor/malignant cells" may enter the diploid
# reference (kharchenkolab/numbat#224), and in glioma those four labels cannot
# distinguish a normal astrocyte from an AC-like tumour cell. The cost is that
# the GBM reference has no glial or neural column at all, which is in tension
# with the same maintainer's "include the relevant normal cell types present in
# your sample". Cluster granularity is what makes that resolvable: a single
# astrocyte CLUSTER judged normal can be admitted by hand, via
# --numbat_extra_reference_clusters, where the LABEL never could. It is empty by
# default, because deciding it from the annotation alone is circular.
#
# Usage:
#   numbat_prep.R <annotated.rds> <sample_col> <group_col> <celltype_col> \
#                 <reference_labels> <min_ref_cells> <outdir> [helpers.R] \
#                 [cluster_col] [min_ref_group_cells] [extra_ref_clusters]
#
# The last three are appended rather than inserted so an older caller keeps
# working: with no cluster column resolved the reference falls back to one
# column per cell type, which is what this script did before.
#
# Writes:
#   <library>_barcodes.tsv   raw barcodes, one file per library (pileup input)
#   <group>_counts.rds       pooled counts, colnames "<library>_<barcode>"
#   <group>_cellmap.tsv      cell_id, barcode, cell, library, group,
#                            in_reference_lineage
#   <group>_refcells.tsv     cells for the internal expression reference, with
#                            the ref_group column that becomes a lambdas_ref
#                            COLUMN (<label>__c<cluster>)
#   <group>_refgroups.tsv    one row per reference column: what the baseline is
#                            actually built from, and what the floor dropped
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
clust_col_in  <- if (length(args) >= 9) args[[9]] else ""
min_ref_group <- if (length(args) >= 10) suppressWarnings(as.integer(args[[10]])) else NA_integer_
extra_clust_in <- if (length(args) >= 11) args[[11]] else ""
if (is.na(min_ref_cells)) min_ref_cells <- 100L
# Per-COLUMN floor. 100 is the number used in the package author's own reference
# construction (vignettes/mouse.Rmd section 4: `filter(n() > 100)`).
if (is.na(min_ref_group)) min_ref_group <- 100L
# A column cannot be larger than the total, so a total floor below the per-column
# floor can never be satisfied by a single column and would silently mean "any
# one valid column is enough". Say so rather than letting the two disagree.
if (min_ref_cells < min_ref_group) {
  message(sprintf(paste("[warn] min_ref_cells (%d) is below the per-column floor (%d);",
                        "raising it to %d so the two are consistent."),
                  min_ref_cells, min_ref_group, min_ref_group))
  min_ref_cells <- min_ref_group
}
extra_clusters <- trimws(strsplit(extra_clust_in, "[;,]")[[1]])
extra_clusters <- extra_clusters[nzchar(extra_clusters)]

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

# The cluster the annotation was actually assigned on. Optional: without it the
# reference degrades to one column per cell type, which is the previous
# behaviour, so a cohort entered with a hand-made metadata table still runs.
clust_col <- pick_col(c(clust_col_in, "seurat_clusters", "clusters", "cluster",
                        "SCT_snn_res.0.8", "RNA_snn_res.0.8"), "Cluster")
if (is.na(clust_col))
  message(paste("[warn] No cluster column found. The expression reference will use ONE",
                "COLUMN PER CELL TYPE, which averages several clusters together and is",
                "coarser than Numbat's own worked example. Set --numbat_cluster_col."))

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
meta$.clust   <- if (is.na(clust_col)) NA_character_ else as.character(meta[[clust_col]])
meta$.barcode <- raw_barcode(meta$.cell_id)

bad <- is.na(meta$.barcode) | is.na(meta$.library) | is.na(meta$.group)
if (any(bad)) {
  message(sprintf("[warn] dropping %d of %d cells with no recoverable barcode, library or group",
                  sum(bad), nrow(meta)))
  meta <- meta[!bad, , drop = FALSE]
}

dt <- as.data.table(meta[, c(".cell_id", ".library", ".group", ".ct", ".clust", ".barcode")])
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
  # signal Numbat is looking for, and the maintainer is explicit about it
  # (kharchenkolab/numbat#224: "Make sure no tumor/malignant cells are included
  # in your diploid reference").
  #
  # Each surviving CLUSTER becomes its own lambdas_ref column -- see the header
  # for why that is finer than one column per label and not merely different.
  ref_cols <- c("cell", "celltype", "cluster", "ref_group")
  ref_cells <- setNames(data.table(character(), character(), character(), character()),
                        ref_cols)
  refgroups <- data.table(ref_group = character(), celltype = character(),
                          cluster = character(), n_cells = integer(),
                          admitted = character(), basis = character())
  n_ref <- 0L
  # Computed unconditionally and kept on `sub`: this is the SAME predicate
  # stratify_compartments.qmd gates on, so the gated malignant fraction reported
  # downstream means exactly what the final compartments mean. It is the label
  # test ALONE -- a reference-lineage cell whose cluster fell below the column
  # floor is still a reference-lineage cell, it just did not help build the
  # baseline, and conflating the two would understate the gate.
  is_ref <- if (!is.na(ct_col) && length(ref_labels)) {
    v <- match_labels(sub$ct, ref_labels); v[is.na(v)] <- FALSE; v
  } else rep(FALSE, nrow(sub))
  sub[, in_reference_lineage := is_ref]

  if (length(ref_labels) || length(extra_clusters)) {

    # Explicit per-cluster opt-in. This is the ONLY route by which a label
    # outside the reference list can reach the baseline, and it exists because
    # in glioma the normal glial lineages share their labels with malignant
    # states -- so the choice has to be made per cluster, by a human, with
    # independent evidence. Never set by default.
    is_extra <- !is.na(sub$clust) & as.character(sub$clust) %in% extra_clusters
    if (any(is_extra)) {
      ex <- sort(unique(as.character(sub$clust[is_extra])))
      message(sprintf(paste("[ref] %s: admitting %d cell(s) from opted-in cluster(s) %s",
                            "regardless of label (--numbat_extra_reference_clusters)."),
                      g, sum(is_extra), paste(ex, collapse = ", ")))
    }
    eligible <- is_ref | is_extra
    n_ref <- sum(eligible)

    if (n_ref > 0L) {
      cand <- data.table(
        cell     = sub$cell[eligible],
        celltype = as.character(sub$ct[eligible]),
        cluster  = if (is.na(clust_col)) NA_character_ else as.character(sub$clust[eligible]),
        # is_ref wins when both apply, so "opt_in" marks only the columns that
        # the label rule would NOT have admitted -- which is the whole point of
        # auditing them. Note the opt-in does NOT make a cell count as
        # in_reference_lineage: it changes what builds the baseline, never what
        # the lineage gate overrides.
        basis    = ifelse(is_ref[eligible], "reference_lineage", "opt_in"))
      cand[is.na(celltype) | !nzchar(celltype), celltype := "unlabelled"]
      # No cluster column resolved -> one column per cell type, the old shape.
      cand[, ref_group := if (is.na(clust_col)) celltype
                          else paste0(celltype, "__c", ifelse(is.na(cluster), "NA", cluster))]

      counts_by <- cand[, .(n_cells = .N,
                            celltype = celltype[1], cluster = cluster[1],
                            basis = paste(sort(unique(basis)), collapse = "+")),
                        by = ref_group]
      counts_by[, admitted := ifelse(n_cells >= min_ref_group, "yes", "no_below_floor")]
      refgroups <- counts_by[order(-n_cells),
                             .(ref_group, celltype, cluster, n_cells, admitted, basis)]
      ref_cells <- cand[ref_group %in% counts_by[admitted == "yes", ref_group],
                        ..ref_cols]

      dropped <- counts_by[admitted != "yes"]
      if (nrow(dropped))
        message(sprintf("[ref] %s: dropped %d column(s) below the %d-cell floor (%d cell(s)).",
                        g, nrow(dropped), min_ref_group, sum(dropped$n_cells)))
    }
  }
  n_ref_cols <- uniqueN(ref_cells$ref_group)
  ref_basis <- if (nrow(ref_cells) >= min_ref_cells) "internal" else "ref_hca"
  if (ref_basis == "internal") {
    message(sprintf("[ref] %s: %d cell(s) over %d reference column(s) (%d cell type(s)): %s",
                    g, nrow(ref_cells), n_ref_cols, uniqueN(ref_cells$celltype),
                    paste(sort(unique(ref_cells$ref_group)), collapse = ", ")))
  } else {
    message(sprintf(paste("[ref] %s: only %d usable reference cell(s) over %d column(s)",
                          "(need %d); numbat_run.R will fall back to numbat::ref_hca.",
                          "Tumour calls for this individual then rest on an external",
                          "reference from a different study."),
                    g, nrow(ref_cells), n_ref_cols, min_ref_cells))
  }

  saveRDS(mat, file.path(outdir, paste0(g, "_counts.rds")))
  # Always written, possibly header-only: an optional per-group file would have
  # to be joined with `remainder`, and a group silently missing its reference
  # would fall back to ref_hca with nothing saying so.
  fwrite(ref_cells, file.path(outdir, paste0(g, "_refcells.tsv")), sep = "\t")
  # What the baseline is actually made of, INCLUDING the columns the floor
  # dropped. Without this the only trace of a discarded column is a log line.
  fwrite(refgroups, file.path(outdir, paste0(g, "_refgroups.tsv")), sep = "\t")
  # in_reference_lineage travels with the cell map so numbat_run.R can report a
  # lineage-gated malignant fraction beside the raw one WITHOUT needing the
  # annotation, the label matcher or Seurat -- none of which exist in the numbat
  # image. The raw fraction is the one a reviewer sees first and the one that is
  # structurally guaranteed to be non-zero.
  fwrite(sub[, .(cell_id, barcode, cell, library, group, in_reference_lineage)],
         file.path(outdir, paste0(g, "_cellmap.tsv")), sep = "\t")

  written <- c(written, g)
  report[[length(report) + 1L]] <- data.table(
    group = g, n_libraries = length(libs), libraries = paste(libs, collapse = ","),
    n_cells = nrow(sub), n_reference_cells = n_ref,
    n_reference_usable = nrow(ref_cells), n_reference_columns = n_ref_cols,
    reference_basis = ref_basis)
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

# How many COLUMNS each individual's baseline has, which is the thing this
# script changed. Numbat matches per pseudobulk against these columns, so one
# averaged column per cell type is strictly less to match against than one per
# cluster. A group reporting a single column is a baseline with nothing to
# choose from.
cat("\n==== expression reference ====\n")
print(rep_dt[, .(group, n_reference_cells, n_reference_usable,
                 n_reference_columns, reference_basis)])
thin <- rep_dt[reference_basis == "internal" & n_reference_columns < 2L]
if (nrow(thin)) {
  cat(sprintf(paste("[note] %d individual(s) have an internal reference with fewer than two",
                    "columns (%s). Numbat picks the best-fitting column per pseudobulk, so a",
                    "single-column reference gives it no choice; raising",
                    "--numbat_ref_min_group_cells above the cluster sizes is the usual",
                    "cause.\n"),
              nrow(thin), paste(thin$group, collapse = ", ")))
}
