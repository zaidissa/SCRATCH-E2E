#!/usr/bin/env Rscript
# ---------------------------------------------------------------------------
# Run Numbat on one GROUP (one individual, all of their libraries) and emit
# per-cell malignancy calls keyed by the pipeline's own cell ids.
#
# Numbat combines expression-derived CNV signal with allele-specific signal from
# phased heterozygous SNPs. The allelic channel is what expression-only callers
# (inferCNV, SCEVAN, CopyKAT) structurally cannot see: a copy-neutral LOH event
# changes allele ratios without changing dosage, so it produces no expression
# signal at all.
#
# WHY PER INDIVIDUAL AND NOT PER LIBRARY. numbat:::mark_tumor_lineage() assigns
# the tumour compartment with
#     tumor_root = mut_nodes[which.max(mut_burdens)]
#     compartment = ifelse(seq > 0, "tumor", "normal")
# There is no significance test and no null: whatever segments survive min_LLR
# go on a tree, the heaviest clade becomes "tumor". Numbat cannot answer "there
# is no tumour in this sample". Run per library, a PBMC library therefore comes
# back with its own noise partitioned into tumour and normal rather than with
# nothing. Pooled with the individual's tumour libraries, the real aneuploid
# clones dominate the mutation burden, tumor_root lands on the real tumour, and
# the normal cells stay normal.
#
# The pileup is also joint: pileup_and_phase.R was given every library for this
# individual at once (--label <group> --samples lib1,lib2,...), so the
# heterozygous sites are genotyped and phased across all of them. It writes ONE
# allele table per library, with that library's raw barcodes; this script
# namespaces each to "<library>_<barcode>" so they can be concatenated and still
# line up with the pooled count matrix.
#
# Output: <group>_numbat_calls.tsv with cell_id, numbat_call, p_cnv, clone.
# STRATIFY reads cell_id and numbat_call; the rest is for the report.
#
# Usage: numbat_run.R <counts.rds> <cellmap.tsv> <refcells.tsv> <group> <outdir> \
#                     <ncores> <min_cells> <t> <max_iter> <min_LLR> <min_ref_cells> \
#                     <max_entropy> <alleles...>
# ---------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(numbat)
  library(Matrix)
  library(data.table)
})

a <- commandArgs(trailingOnly = TRUE)
if (length(a) < 13L)
  stop("usage: numbat_run.R <counts.rds> <cellmap.tsv> <refcells.tsv> <group> <outdir> ",
       "<ncores> <min_cells> <t> <max_iter> <min_LLR> <min_ref_cells> <max_entropy> ",
       "<allele_counts.tsv.gz...>")

counts_path <- a[[1]]; map_path <- a[[2]]; refcells_path <- a[[3]]
group_id    <- a[[4]]; outdir   <- a[[5]]
ncores      <- as.integer(a[[6]])
min_cells   <- as.integer(a[[7]])
t_param     <- as.numeric(a[[8]])
max_iter    <- as.integer(a[[9]])
# Numbat's own failure message says "Consider reducing min_LLR", but the value
# was not reachable from the pipeline. Exposed so that advice is actionable.
min_llr       <- as.numeric(a[[10]])
min_ref_cells <- as.integer(a[[11]])
# Single-cell filtering threshold. Numbat's maintainer names this and min_LLR as
# the two levers for "a tumor proportion much higher than orthogonal
# measurements suggest" -- i.e. CNV false positives in normal cells
# (kharchenkolab/numbat#224). Exposed at Numbat's own default of 0.5, so passing
# it changes nothing until it is deliberately lowered.
max_entropy   <- as.numeric(a[[12]])
if (is.na(max_entropy)) max_entropy <- 0.5
allele_paths  <- a[13:length(a)]

dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

count_mat <- readRDS(counts_path)
cellmap   <- fread(map_path)

# A group that yields no calls is a RESULT, not a pipeline failure. Emitting an
# empty-but-valid table keeps the output contract, lets STRATIFY see that this
# group contributed nothing, and -- the reason this exists -- stops one group
# from destroying the whole run. A 7h42m run died here with the OTHER sample's
# calls already complete, and threw them away with it.
emit_empty <- function(outdir, group_id, why) {
  message(sprintf("[%s] NO CALLS: %s", group_id, why))
  dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
  empty <- data.table::data.table(
    cell = character(), barcode = character(), cell_id = character(),
    sample = character(), group = character(),
    numbat_call = character(), p_cnv = numeric(),
    clone_opt = integer(), compartment_opt = character())
  data.table::fwrite(empty,
    file.path(outdir, paste0(group_id, "_numbat_calls.tsv")), sep = "\t")
  # Machine-readable reason, so a report can say WHY, not merely "0 rows".
  writeLines(why, file.path(outdir, paste0(group_id, "_numbat_nocalls.txt")))
  quit(save = "no", status = 0)
}

# ---------------------------------------------------------------------------
# Concatenate the per-library allele tables, namespacing each cell.
#
# pileup_and_phase.R writes {outdir}/{sample}_allele_counts.tsv.gz per library,
# and the `cell` column holds that library's RAW cellranger barcodes. Two
# libraries of one patient both contain AACGT...-1, so concatenating without a
# prefix silently merges distinct cells. The prefix is the library id recovered
# from the filename, which is the same string numbat_prep.R used to build the
# count matrix column names.
# ---------------------------------------------------------------------------
libs_seen <- character(0)
allele_list <- lapply(allele_paths, function(p) {
  lb <- sub("_allele_counts\\.tsv\\.gz$", "", basename(p))
  d  <- fread(p)
  if (!nrow(d)) {
    message(sprintf("[%s] %s: allele table is empty", group_id, lb))
    return(NULL)
  }
  if (!"cell" %in% names(d))
    stop(sprintf("[%s] %s: no 'cell' column in %s", group_id, lb, basename(p)))
  d[, cell := paste0(lb, "_", as.character(cell))]
  libs_seen <<- c(libs_seen, lb)
  d
})
allele_list <- Filter(Negate(is.null), allele_list)
if (!length(allele_list))
  emit_empty(outdir, group_id, sprintf(
    "All %d allele table(s) for this group were empty; the pileup produced no phased sites.",
    length(allele_paths)))

df_allele <- rbindlist(allele_list, use.names = TRUE, fill = TRUE)
rm(allele_list); invisible(gc())

message(sprintf("[%s] %d genes x %d cells | %d allele records from %d librar%s (%s)",
                group_id, nrow(count_mat), ncol(count_mat), nrow(df_allele),
                length(libs_seen), if (length(libs_seen) == 1L) "y" else "ies",
                paste(libs_seen, collapse = ", ")))

# The pooled count matrix and the concatenated allele table must address the
# same cells. A total mismatch means the prefixes disagree, which is silent and
# fatal to the result rather than to the run: Numbat would model expression with
# no allelic evidence at all and fall back to behaving like inferCNV.
shared <- length(intersect(colnames(count_mat), unique(df_allele$cell)))
message(sprintf("[%s] %d cell(s) carry both expression and allele data", group_id, shared))
if (shared == 0L)
  emit_empty(outdir, group_id, paste(
    "No cell appears in both the count matrix and the allele tables, so the",
    "allelic channel is empty. The two are namespaced '<library>_<barcode>' from",
    "the library id; check that the BAM directory names match the library column."))
if (shared < 0.5 * ncol(count_mat))
  message(sprintf(paste("[%s] WARNING: only %.1f%% of cells have allele data.",
                        "Numbat will model the rest from expression alone."),
                  group_id, 100 * shared / ncol(count_mat)))

if (ncol(count_mat) < min_cells)
  emit_empty(outdir, group_id,
             sprintf("%d cells is below --numbat_min_cells=%d; group skipped.",
                     ncol(count_mat), min_cells))

# ---------------------------------------------------------------------------
# Expression reference. Prefer the one numbat_prep.R built from THIS
# individual's own immune and stromal cells: numbat::ref_hca is 12 generic
# profiles from a different study, and the gap between it and the real cohort is
# read as expression CNV signal. The internal reference is why running per
# individual helps twice -- it is also what makes a normal library useful rather
# than merely harmless.
# ---------------------------------------------------------------------------
# numbat_prep.R SELECTS the cells (it has the annotation and the shared label
# matcher); the aggregation happens here because aggregate_counts() needs the
# numbat package, and prep runs in image_cnv, which has Seurat and no numbat.
# The two images are mirrors of each other -- verified, not assumed.
#
# The `group` column handed to aggregate_counts is the lambdas_ref COLUMN name.
# prep writes it as <celltype>__c<cluster>, one column per cluster, because
# Numbat picks the best-fitting column per pseudobulk and the author's own
# reference construction uses raw clusters (vignettes/mouse.Rmd section 4).
# `ref_group` is read in preference to `celltype`; falling back to `celltype`
# keeps an older refcells.tsv working.
#
# count_mat here is RAW counts and aggregate_counts does the normalisation
# (normalized = TRUE). That ordering is a requirement, not a preference: the
# maintainer's instruction is "created using numbat::aggregate_counts() from raw
# counts ... Don't use Seurat normalization or other methods"
# (kharchenkolab/numbat#224). numbat_prep.R reads layer = "counts" for the same
# reason.
ref <- NULL
ref_basis <- "ref_hca"
refcells <- tryCatch(fread(refcells_path), error = function(e) NULL)
grp_col <- if (!is.null(refcells) && "ref_group" %in% names(refcells)) "ref_group" else "celltype"
if (!is.null(refcells) && nrow(refcells) >= min_ref_cells &&
    all(c("cell", grp_col) %in% names(refcells))) {
  refcells <- refcells[cell %in% colnames(count_mat)]
  if (nrow(refcells) >= min_ref_cells) {
    annot <- data.frame(cell = refcells$cell,
                        group = as.character(refcells[[grp_col]]),
                        stringsAsFactors = FALSE)
    ref <- tryCatch(
      numbat::aggregate_counts(count_mat[, annot$cell, drop = FALSE], annot,
                               normalized = TRUE, verbose = FALSE),
      error = function(e) {
        message(sprintf("[%s] aggregate_counts failed (%s); falling back to ref_hca",
                        group_id, conditionMessage(e))); NULL })
    if (!is.null(ref)) {
      ref_basis <- "internal"
      message(sprintf(paste("[%s] internal reference from %d of this individual's own cells",
                            "over %d column(s), grouped by '%s'"),
                      group_id, nrow(annot), length(unique(annot$group)), grp_col))
      if (grp_col == "celltype")
        message(sprintf(paste("[%s] [warn] refcells.tsv has no ref_group column, so the",
                              "reference is ONE COLUMN PER CELL TYPE. That averages several",
                              "clusters into each column; set --numbat_cluster_col so prep",
                              "can write cluster-granular columns."), group_id))
    }
  }
}
if (is.null(ref)) {
  ref <- tryCatch(numbat::ref_hca, error = function(e) NULL)
  message(sprintf(paste("[%s] no internal reference (need %d reference-lineage cells);",
                        "using numbat::ref_hca"), group_id, min_ref_cells))
}
if (is.null(ref))
  stop(sprintf("[%s] no expression reference available (internal missing and numbat::ref_hca unreadable)",
               group_id))
message(sprintf("[%s] expression reference: %s (%d genes x %d profile(s): %s)",
                group_id, ref_basis, nrow(ref), ncol(ref),
                paste(colnames(ref), collapse = ", ")))

out_sub <- file.path(outdir, group_id)
dir.create(out_sub, recursive = TRUE, showWarnings = FALSE)

res <- try(
  numbat::run_numbat(
    count_mat = count_mat,
    lambdas_ref = ref,
    df_allele = df_allele,
    genome = Sys.getenv("NUMBAT_GENOME", "hg38"),
    t = t_param,
    ncores = ncores,
    max_iter = max_iter,
    min_LLR = min_llr,
    max_entropy = max_entropy,
    plot = FALSE,
    out_dir = out_sub
  ),
  silent = TRUE
)

if (inherits(res, "try-error")) {
  msg <- as.character(res)
  # Numbat drops cells with no usable SNP coverage BEFORE it checks the two
  # inputs against each other, so a pileup that produced nothing usable surfaces
  # as "No matching cell names between count_mat and df_allele" -- which reads
  # like the namespacing is broken when it is not. The explicit overlap check
  # above already proved the names line up, so reaching here means the coverage
  # is the problem. That is a RESULT for this individual, not a broken run, and
  # killing the task would take the whole cohort's other individuals with it.
  if (grepl("No matching cell names|0 coverage", msg)) {
    emit_empty(outdir, group_id, sprintf(paste(
      "%d cell(s) carried both expression and allele data, but Numbat found no",
      "usable SNP coverage in any of them and dropped them all. The pileup ran",
      "and produced %d allele record(s) across %d librar%s; none fell on a",
      "heterozygous site this individual's reads cover deeply enough. Check the",
      "phasing panel matches --numbat_genome (%s) and that the BAMs are the ones",
      "these barcodes came from."),
      shared, nrow(df_allele), length(libs_seen),
      if (length(libs_seen) == 1L) "y" else "ies",
      Sys.getenv("NUMBAT_GENOME", "hg38")))
  }
  stop(sprintf("[%s] run_numbat failed: %s", group_id, msg))
}

# run_numbat writes its per-cell posterior to clone_post_2.tsv. Later Numbat
# versions may emit a different iteration index, so take the highest available.
posts <- list.files(out_sub, pattern = "^clone_post_\\d+\\.tsv$", full.names = TRUE)
if (!length(posts)) {
  emit_empty(outdir, group_id, paste(
    "Numbat completed but no CNV survived filtering by log-likelihood ratio,",
    "so no clone posterior was produced. Either this individual carries no",
    "detectable clonal copy-number change, or the threshold is too strict for",
    "their coverage -- lower --numbat_min_llr to test the second possibility."))
}
post_path <- posts[order(as.integer(sub(".*clone_post_(\\d+)\\.tsv$", "\\1", posts)))][length(posts)]
message("[", group_id, "] reading ", basename(post_path))

post <- fread(post_path)

cell_col <- intersect(c("cell", "cell_id", "barcode"), names(post))[1]
if (is.na(cell_col)) stop(sprintf("[%s] no cell column in %s", group_id, basename(post_path)))
setnames(post, cell_col, "cell")
post[, cell := as.character(cell)]

# `compartment_opt` is Numbat's own tumour/normal assignment; p_cnv is the
# posterior probability of carrying any CNV. Prefer the explicit assignment and
# fall back to thresholding the posterior.
if ("compartment_opt" %in% names(post)) {
  post[, numbat_call := fifelse(grepl("tumou?r", tolower(compartment_opt)), "malignant",
                        fifelse(grepl("normal", tolower(compartment_opt)), "non_malignant",
                                NA_character_))]
} else if ("p_cnv" %in% names(post)) {
  post[, numbat_call := fifelse(is.na(p_cnv), NA_character_,
                        fifelse(p_cnv > 0.5, "malignant", "non_malignant"))]
} else {
  stop(sprintf("[%s] neither compartment_opt nor p_cnv present", group_id))
}

keep <- intersect(c("cell", "numbat_call", "p_cnv", "clone_opt", "compartment_opt"), names(post))
post <- post[, ..keep]

# Join back onto the pipeline's cell ids via the explicit map. `cell` is the
# namespaced "<library>_<barcode>" that both sides were built from.
# in_reference_lineage comes from prep and is the SAME predicate
# stratify_compartments.qmd gates on. Kept optional so an older cellmap still
# joins; the gated column is then simply absent rather than wrong.
map_cols <- c("cell", "cell_id", "barcode", "library", "group")
has_gate <- "in_reference_lineage" %in% names(cellmap)
if (has_gate) map_cols <- c(map_cols, "in_reference_lineage")
res_dt <- merge(cellmap[, ..map_cols], post, by = "cell", all.x = FALSE)
setnames(res_dt, "library", "sample")
if (!nrow(res_dt)) stop(sprintf("[%s] no Numbat cell matched the cell map", group_id))

fwrite(res_dt, file.path(outdir, paste0(group_id, "_numbat_calls.tsv")), sep = "\t")

# ---------------------------------------------------------------------------
# Per-library breakdown. This is the number the per-library configuration got
# wrong and the number a reviewer will look at first: a library that should be
# ~0% neoplastic (PBMC, matched normal) saying otherwise is visible here,
# in .command.log and in a published table, rather than only after someone
# plots the final object.
# ---------------------------------------------------------------------------
by_lib <- res_dt[, .(
  n_cells    = .N,
  malignant  = sum(numbat_call == "malignant", na.rm = TRUE),
  pct_malignant = round(100 * mean(numbat_call == "malignant", na.rm = TRUE), 1)
), by = .(sample)][order(sample)]

# The RAW percentage above is Numbat's compartment_opt, ungated, and for a
# library with no malignant cells it is structurally guaranteed to be non-zero:
# mark_tumor_lineage() always labels the heaviest clade tumour. The percentage
# that survives the lineage gate is what the final compartments actually use.
# Reporting only the first is how a PBMC library comes to look neoplastic in a
# table while the stratified object is fine -- two different numbers, one of
# which nobody had reason to doubt.
if (has_gate) {
  gated <- res_dt[, .(
    n_reference_lineage = sum(in_reference_lineage %in% c(TRUE, "TRUE", "True")),
    malignant_gated = sum(numbat_call == "malignant" &
                          !(in_reference_lineage %in% c(TRUE, "TRUE", "True")), na.rm = TRUE)
  ), by = .(sample)]
  by_lib <- merge(by_lib, gated, by = "sample", all.x = TRUE)
  by_lib[, pct_malignant_gated := round(100 * malignant_gated / pmax(n_cells, 1L), 1)]
  by_lib[, gate_overrode := malignant - malignant_gated]
  setnames(by_lib, "pct_malignant", "pct_malignant_raw")
  setcolorder(by_lib, c("sample", "n_cells", "n_reference_lineage", "malignant",
                        "pct_malignant_raw", "malignant_gated",
                        "pct_malignant_gated", "gate_overrode"))
}
by_lib[, reference_basis := ref_basis]
by_lib[, n_reference_columns := if (ref_basis == "internal") ncol(ref) else NA_integer_]
fwrite(by_lib, file.path(outdir, paste0(group_id, "_numbat_by_library.tsv")), sep = "\t")

message(sprintf("[%s] %d cells called | per library:", group_id, nrow(res_dt)))
pct_raw <- if (has_gate) by_lib$pct_malignant_raw else by_lib$pct_malignant
for (i in seq_len(nrow(by_lib))) {
  if (has_gate) {
    message(sprintf("  %-40s %6d cells  %6d malignant raw (%.1f%%)  %6d after lineage gate (%.1f%%)",
                    by_lib$sample[i], by_lib$n_cells[i],
                    by_lib$malignant[i], pct_raw[i],
                    by_lib$malignant_gated[i], by_lib$pct_malignant_gated[i]))
  } else {
    message(sprintf("  %-40s %6d cells  %6d malignant (%.1f%%)",
                    by_lib$sample[i], by_lib$n_cells[i],
                    by_lib$malignant[i], pct_raw[i]))
  }
}
if (has_gate)
  message(sprintf(paste("[%s] 'raw' is Numbat's ungated compartment_opt and is non-zero even",
                        "for a library with no malignant cells; 'after lineage gate' is what",
                        "the stratified compartments use."), group_id))
