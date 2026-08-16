#!/usr/bin/env Rscript
# ---------------------------------------------------------------------------
# Run Numbat on one sample and emit per-cell malignancy calls keyed by the
# pipeline's own cell ids.
#
# Numbat combines expression-derived CNV signal with allele-specific signal from
# phased heterozygous SNPs. The allelic channel is what expression-only callers
# (inferCNV, SCEVAN, CopyKAT) structurally cannot see: a copy-neutral LOH event
# changes allele ratios without changing dosage, so it produces no expression
# signal at all.
#
# Output: <sample>_numbat_calls.tsv with cell_id, numbat_call, p_cnv, clone.
# The STRATIFY step reads that file; nothing downstream needs the Numbat object.
#
# Usage: numbat_run.R <counts.rds> <allele_counts.tsv.gz> <cellmap.tsv> \
#                     <sample_id> <outdir> <ncores> <min_cells> <t> <max_iter>
# ---------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(numbat)
  library(Matrix)
  library(data.table)
})

a <- commandArgs(trailingOnly = TRUE)
if (length(a) < 6L) stop("usage: numbat_run.R <counts.rds> <alleles.tsv.gz> <cellmap.tsv> <sample> <outdir> <ncores> [min_cells] [t] [max_iter] [min_LLR]")

counts_path <- a[[1]]; allele_path <- a[[2]]; map_path <- a[[3]]
sample_id   <- a[[4]]; outdir      <- a[[5]]
ncores      <- as.integer(a[[6]])
min_cells   <- if (length(a) >= 7) as.integer(a[[7]]) else 50L
t_param     <- if (length(a) >= 8) as.numeric(a[[8]]) else 1e-5
max_iter    <- if (length(a) >= 9) as.integer(a[[9]]) else 2L
# Numbat's own failure message says "Consider reducing min_LLR", but the value
# was not reachable from the pipeline. Exposed so that advice is actionable.
min_llr     <- if (length(a) >= 10) as.numeric(a[[10]]) else 5

dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

count_mat <- readRDS(counts_path)
df_allele <- fread(allele_path)
cellmap   <- fread(map_path)

message(sprintf("[%s] %d genes x %d cells | %d allele records",
                sample_id, nrow(count_mat), ncol(count_mat), nrow(df_allele)))

# A sample that yields no calls is a RESULT, not a pipeline failure. Emitting an
# empty-but-valid table keeps the output contract, lets STRATIFY see that this
# sample contributed nothing, and -- the reason this exists -- stops one sample
# from destroying the whole run. A 7h42m run died here with the OTHER sample's
# calls already complete, and threw them away with it.
emit_empty <- function(outdir, sample_id, why) {
  message(sprintf("[%s] NO CALLS: %s", sample_id, why))
  dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
  empty <- data.table::data.table(
    barcode = character(), cell_id = character(), sample = character(),
    numbat_call = character(), p_cnv = numeric(),
    clone_opt = integer(), compartment_opt = character())
  data.table::fwrite(empty,
    file.path(outdir, paste0(sample_id, "_numbat_calls.tsv")), sep = "\t")
  # Machine-readable reason, so a report can say WHY, not merely "0 rows".
  writeLines(why, file.path(outdir, paste0(sample_id, "_numbat_nocalls.txt")))
  quit(save = "no", status = 0)
}

if (ncol(count_mat) < min_cells) {
  emit_empty(outdir, sample_id,
             sprintf("%d cells is below --numbat_min_cells=%d; sample skipped.",
                     ncol(count_mat), min_cells))
}

# Numbat needs an expression reference. `ref_hca` ships with the package and is
# the documented default; a cohort-internal reference built from this study's own
# immune cells would be preferable, but it is not available per-sample here.
ref <- tryCatch(numbat::ref_hca, error = function(e) NULL)
if (is.null(ref)) stop("numbat::ref_hca unavailable in this image; cannot supply an expression reference.")

out_sub <- file.path(outdir, sample_id)
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
    plot = FALSE,
    out_dir = out_sub
  ),
  silent = TRUE
)

if (inherits(res, "try-error")) {
  stop(sprintf("[%s] run_numbat failed: %s", sample_id, as.character(res)))
}

# run_numbat writes its per-cell posterior to clone_post_2.tsv. Later Numbat
# versions may emit a different iteration index, so take the highest available.
posts <- list.files(out_sub, pattern = "^clone_post_\\d+\\.tsv$", full.names = TRUE)
if (!length(posts)) {
  emit_empty(outdir, sample_id, paste(
    "Numbat completed but no CNV survived filtering by log-likelihood ratio,",
    "so no clone posterior was produced. Either this sample carries no detectable",
    "clonal copy-number change, or the threshold is too strict for its coverage --",
    "lower --numbat_min_llr to test the second possibility."))
}
post_path <- posts[order(as.integer(sub(".*clone_post_(\\d+)\\.tsv$", "\\1", posts)))][length(posts)]
message("[", sample_id, "] reading ", basename(post_path))

post <- fread(post_path)

cell_col <- intersect(c("cell", "cell_id", "barcode"), names(post))[1]
if (is.na(cell_col)) stop(sprintf("[%s] no cell column in %s", sample_id, basename(post_path)))
setnames(post, cell_col, "barcode")
post[, barcode := as.character(barcode)]

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
  stop(sprintf("[%s] neither compartment_opt nor p_cnv present", sample_id))
}

keep <- intersect(c("barcode", "numbat_call", "p_cnv", "clone_opt", "compartment_opt"), names(post))
post <- post[, ..keep]

# Join back onto the pipeline's cell ids via the explicit map.
res_dt <- merge(cellmap[, .(cell_id, barcode, sample)], post, by = "barcode", all.x = FALSE)
if (!nrow(res_dt)) stop(sprintf("[%s] no Numbat cell matched the barcode map", sample_id))

fwrite(res_dt, file.path(outdir, paste0(sample_id, "_numbat_calls.tsv")), sep = "\t")

message(sprintf("[%s] %d cells called: %s", sample_id, nrow(res_dt),
                paste(sprintf("%s=%d", names(table(res_dt$numbat_call)),
                              as.integer(table(res_dt$numbat_call))), collapse = " ")))
