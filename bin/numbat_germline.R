#!/usr/bin/env Rscript
# ---------------------------------------------------------------------------
# One individual's germline genotypes, from a WGS VCF, in the two shapes the
# Numbat stage needs. Runs in the numbat image (data.table, no bcftools there).
#
# WHY. Without this, numbat:::genotype() INFERS the heterozygous sites from the
# scRNA itself: it pools the individual's reads at 1000G common SNPs and calls a
# site het when 0.1 <= AR <= 0.9 (hom-alt when AR == 1 and DP >= 10). That has
# two blind spots, and the second is the one that matters:
#   - only expressed, covered sites can be called at all;
#   - in a region of tumour LOH the tumour reads carry ONE allele, so a truly
#     heterozygous site looks homozygous and is dropped -- the sites that
#     carry Numbat's strongest LOH evidence are removed before it looks.
# A germline call from normal-tissue WGS (here PBMC) has neither problem.
#
# OUTPUTS (in <outdir>/):
#   <label>_germline_targets.vcf    sites for cellsnp-lite -R. 8 columns, no
#                                   "chr" -- the same shape as the 1000G SNP VCF
#                                   the pileup otherwise uses, which cellsnp-lite
#                                   already matches against cellranger's chr BAMs.
#   phasing/<label>_chr<N>.vcf.gz   Eagle2 targets, chr1-22, bgzipped + tabix.
#                                   The shape numbat:::make_vcf_chr() writes:
#                                   "chr" prefix, FORMAT GT only, one sample
#                                   column named <label>, 0/1 and 1/1 rows.
#   <label>_germline_sites.tsv.gz   CHROM POS REF ALT GT, for the identity check.
#   <label>_germline_report.tsv     what each filter kept, and the build check.
#
# 1/1 sites are kept as well as hets, as genotype() keeps them: Eagle2 gets
# them as target genotypes, and the identity check needs them (at a true 1/1
# site the scRNA reads should be ALT, whoever's tumour it is).
#
# Usage: numbat_germline.R <vcf> <sample|''> <label> <outdir> <min_dp> <min_gq> <snpvcf>
# ---------------------------------------------------------------------------
suppressPackageStartupMessages(library(data.table))

a <- commandArgs(trailingOnly = TRUE)
if (length(a) < 7) stop("usage: numbat_germline.R <vcf> <sample|''> <label> <outdir> <min_dp> <min_gq> <snpvcf>")
vcf_path <- a[1]; want_sample <- a[2]; label <- a[3]; outdir <- a[4]
min_dp <- as.integer(a[5]); min_gq <- as.integer(a[6]); snpvcf <- a[7]
min_build_overlap <- 0.30

dir.create(file.path(outdir, "phasing"), recursive = TRUE, showWarnings = FALSE)
report <- list()
note <- function(step, n, detail = "") {
  report[[length(report) + 1]] <<- data.table(step = step, n = n, detail = detail)
  message(sprintf("[germline %s] %-28s %10s  %s", label, step, format(n, big.mark = ","), detail))
}

# ---- read --------------------------------------------------------------------
# grep -v '^##' leaves the #CHROM header for fread. zcat -f reads both plain
# and gzipped/bgzipped VCFs.
v <- fread(cmd = sprintf("zcat -f '%s' | grep -v '^##'", vcf_path),
           sep = "\t", colClasses = "character", showProgress = FALSE)
setnames(v, 1, "CHROM")
fixed <- c("CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT")
if (!all(fixed %in% names(v))) stop("not a VCF with genotypes: missing ",
                                    paste(setdiff(fixed, names(v)), collapse = ", "))
samples <- setdiff(names(v), fixed)
if (nzchar(want_sample)) {
  if (!want_sample %in% samples)
    stop(sprintf("sample '%s' is not in %s. Samples there: %s",
                 want_sample, basename(vcf_path), paste(samples, collapse = ", ")))
  smp <- want_sample
} else {
  # A joint-called VCF holds several people. Guessing which column is this
  # individual would be guessing whose genotypes phase the tumour.
  if (length(samples) != 1)
    stop(sprintf("%s holds %d samples (%s); name the one for %s in the 'sample' column.",
                 basename(vcf_path), length(samples), paste(samples, collapse = ", "), label))
  smp <- samples
}
v <- v[, c(fixed, smp), with = FALSE]; setnames(v, smp, "SMP")
note("records", nrow(v), sprintf("sample %s", smp))

# ---- filters -----------------------------------------------------------------
v[, CHROM := sub("^chr", "", CHROM)]
v <- v[CHROM %in% as.character(1:22)]
note("chr1-22", nrow(v))

# An unfiltered HaplotypeCaller VCF has FILTER "." on every row. Requiring PASS
# there would keep nothing, so "." is accepted only when NO row carries a filter
# -- i.e. when no filtering was ever applied -- and that is said out loud.
if (all(v$FILTER %in% c(".", ""))) {
  note("FILTER PASS", nrow(v), "VCF is unfiltered (FILTER '.'): no filter applied; DP/GQ below still are")
} else {
  v <- v[FILTER == "PASS"]
  note("FILTER PASS", nrow(v))
}

v <- v[nchar(REF) == 1L & REF %in% c("A", "C", "G", "T") &
       nchar(ALT) == 1L & ALT %in% c("A", "C", "G", "T")]
note("biallelic SNPs", nrow(v))

# FORMAT differs row to row (GT:AD:DP:GQ:PL vs GT:AD:DP:GQ:PGT:PID:PL:PS), so
# each key is found by name, per distinct FORMAT string.
get_key <- function(fmt, val, key) {
  out <- rep(NA_character_, length(fmt))
  for (f in unique(fmt)) {
    i <- which(fmt == f); k <- match(key, strsplit(f, ":", fixed = TRUE)[[1]])
    if (is.na(k)) next
    parts <- tstrsplit(val[i], ":", fixed = TRUE, fill = NA_character_)
    if (k <= length(parts)) out[i] <- parts[[k]]
  }
  out
}
v[, GT := get_key(FORMAT, SMP, "GT")]
v[, DPg := suppressWarnings(as.integer(get_key(FORMAT, SMP, "DP")))]
v[, GQ := suppressWarnings(as.integer(get_key(FORMAT, SMP, "GQ")))]
v[, GT := chartr("|", "/", GT)]
v[GT == "1/0", GT := "0/1"]
v <- v[GT %in% c("0/1", "1/1")]
note("GT 0/1 or 1/1", nrow(v), sprintf("%d het, %d hom-alt", sum(v$GT == "0/1"), sum(v$GT == "1/1")))
v <- v[!is.na(DPg) & DPg >= min_dp & !is.na(GQ) & GQ >= min_gq]
note(sprintf("DP>=%d, GQ>=%d", min_dp, min_gq), nrow(v),
     sprintf("%d het, %d hom-alt", sum(v$GT == "0/1"), sum(v$GT == "1/1")))

v[, POS := as.integer(POS)]
v <- unique(v, by = c("CHROM", "POS"))
v[, CHROM := factor(CHROM, levels = as.character(1:22))]
setorder(v, CHROM, POS)
v[, CHROM := as.character(CHROM)]
if (!nrow(v)) stop("no germline sites survived the filters; see the report above.")

# ---- build check -------------------------------------------------------------
# An hg19 VCF against hg38 BAMs does not error: every position is simply wrong.
# Most of a person's het SNPs are common, so on the right build most of them
# sit on a 1000G common-SNP position with the same REF. On the wrong build
# almost none do.
snp <- fread(cmd = sprintf("zcat -f '%s' | grep -v '^#' | cut -f1,2,4", snpvcf),
             sep = "\t", header = FALSE, col.names = c("CHROM", "POS", "REF"),
             colClasses = c("character", "integer", "character"), showProgress = FALSE)
snp[, CHROM := sub("^chr", "", CHROM)]
snp <- unique(snp, by = c("CHROM", "POS", "REF"))  # multi-ALT sites appear once per ALT
ov <- nrow(merge(v[, .(CHROM, POS, REF)], snp, by = c("CHROM", "POS", "REF")))
frac <- ov / nrow(v); rm(snp); invisible(gc())
note("on 1000G common SNPs", ov, sprintf("%.1f%% of kept sites", 100 * frac))
if (frac < min_build_overlap)
  stop(sprintf(paste("only %.1f%% of germline sites match a 1000G hg38 common SNP (position + REF).",
                     "On the right build this is typically most of them. Check that the VCF is",
                     "hg38/GRCh38, the same build as the BAMs and --numbat_genome."), 100 * frac))

# ---- write -------------------------------------------------------------------
fwrite(rbindlist(report), file.path(outdir, sprintf("%s_germline_report.tsv", label)), sep = "\t")
fwrite(v[, .(CHROM, POS, REF, ALT, GT)],
       file.path(outdir, sprintf("%s_germline_sites.tsv.gz", label)), sep = "\t")

tg <- file.path(outdir, sprintf("%s_germline_targets.vcf", label))
writeLines(c("##fileformat=VCFv4.2", "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO"), tg)
fwrite(v[, .(CHROM, POS, ID = ".", REF, ALT, QUAL = ".", FILTER = "PASS", INFO = ".")],
       tg, sep = "\t", append = TRUE, col.names = FALSE)

for (chr in 1:22) {
  d <- v[CHROM == as.character(chr)]
  if (!nrow(d)) { message(sprintf("[germline %s] no sites on chr%d", label, chr)); next }
  f <- file.path(outdir, "phasing", sprintf("%s_chr%d.vcf", label, chr))
  writeLines(c("##fileformat=VCFv4.2",
               sprintf("##contig=<ID=chr%d>", chr),
               '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">',
               paste(c("#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT", label),
                     collapse = "\t")), f)
  fwrite(d[, .(paste0("chr", CHROM), POS, ".", REF, ALT, ".", "PASS", ".", "GT", GT)],
         f, sep = "\t", append = TRUE, col.names = FALSE)
  if (system2("bgzip", c("-f", shQuote(f))) != 0 ||
      system2("tabix", c("-f", "-p", "vcf", shQuote(paste0(f, ".gz")))) != 0)
    stop("bgzip/tabix failed for ", f)
}
message(sprintf("[germline %s] done: %s sites", label, format(nrow(v), big.mark = ",")))
