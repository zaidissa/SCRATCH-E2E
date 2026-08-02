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
if (length(a) < 6L) stop("usage: numbat_run.R <counts.rds> <alleles.tsv.gz> <cellmap.tsv> <sample> <outdir> <ncores> [min_cells] [t] [max_iter]")

counts_path <- a[[1]]; allele_path <- a[[2]]; map_path <- a[[3]]
sample_id   <- a[[4]]; outdir      <- a[[5]]
ncores      <- as.integer(a[[6]])
min_cells   <- if (length(a) >= 7) as.integer(a[[7]]) else 50L
t_param     <- if (length(a) >= 8) as.numeric(a[[8]]) else 1e-5
max_iter    <- if (length(a) >= 9) as.integer(a[[9]]) else 2L

dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

count_mat <- readRDS(counts_path)
df_allele <- fread(allele_path)
cellmap   <- fread(map_path)

message(sprintf("[%s] %d genes x %d cells | %d allele records",
                sample_id, nrow(count_mat), ncol(count_mat), nrow(df_allele)))

if (ncol(count_mat) < min_cells) {
  stop(sprintf("[%s] %d cells < min_cells=%d", sample_id, ncol(count_mat), min_cells))
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
if (!length(posts)) stop(sprintf("[%s] no clone_post_*.tsv produced", sample_id))
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
