#!/usr/bin/env Rscript
# ---------------------------------------------------------------------------
# Is the germline VCF the same person as the scRNA? Checked BEFORE phasing.
#
# A WGS VCF from the wrong individual does not fail anywhere downstream. Eagle2
# phases it, Numbat reads allele imbalance off a genotype that is not this
# person's, and the result is a plausible-looking CNV profile built on noise.
# So this compares the two, using the pileup the scRNA has just produced at the
# germline sites (cellSNP.base.vcf INFO AD/DP, pooled over the libraries).
#
#   hom-alt concordance  at germline 1/1 sites with pooled DP >= 10, the share
#                        whose reads are >= 90% ALT. The primary test: no LOH
#                        or allele-specific expression can make a true 1/1 site
#                        show REF reads, while for another person many of those
#                        sites are 0/0 or 0/1.
#   het biallelic        at germline 0/1 sites with DP >= 20, the share showing
#                        >= 2 reads of each allele. Reported, not tested: LOH,
#                        imprinting and allele-specific expression all lower it
#                        legitimately.
#
# Exits 3 when hom-alt concordance is below <min_concordance> over at least
# <min_sites> covered sites. 0 disables the test (the numbers are still
# written). The threshold is a starting point, not a calibrated cut-off.
#
# Usage: numbat_germline_check.R <label> <sites.tsv.gz> <out.tsv> <min_concordance>
#                                <pileup_dir>...
# ---------------------------------------------------------------------------
suppressPackageStartupMessages(library(data.table))
a <- commandArgs(trailingOnly = TRUE)
label <- a[1]; sites_path <- a[2]; out <- a[3]; min_conc <- as.numeric(a[4])
dirs <- a[5:length(a)]
min_sites <- 200L

sites <- fread(sites_path, colClasses = c(CHROM = "character"))
pu <- rbindlist(lapply(dirs, function(d) {
  x <- fread(cmd = sprintf("grep -v '^##' '%s'", file.path(d, "cellSNP.base.vcf")),
             select = c(1, 2, 4, 5, 8), colClasses = "character", showProgress = FALSE)
  setnames(x, c("CHROM", "POS", "REF", "ALT", "INFO"))
  x[, AD := as.integer(sub(".*(^|;)AD=([0-9]+).*", "\\2", INFO))]
  x[, DP := as.integer(sub(".*(^|;)DP=([0-9]+).*", "\\2", INFO))]
  x[, .(CHROM = sub("^chr", "", CHROM), POS = as.integer(POS), REF, ALT, AD, DP)]
}))
pu <- pu[, .(AD = sum(AD), DP = sum(DP)), by = .(CHROM, POS, REF, ALT)]
m <- merge(sites, pu, by = c("CHROM", "POS", "REF", "ALT"))

ha <- m[GT == "1/1" & DP >= 10]
he <- m[GT == "0/1" & DP >= 20]
conc <- if (nrow(ha)) mean(ha$AD / ha$DP >= 0.9) else NA_real_
bial <- if (nrow(he)) mean(he$AD >= 2 & (he$DP - he$AD) >= 2) else NA_real_

res <- data.table(
  individual = label,
  germline_sites = nrow(sites), sites_with_reads = nrow(m),
  homalt_sites_dp10 = nrow(ha), homalt_concordance = round(conc, 4),
  het_sites_dp20 = nrow(he), het_biallelic = round(bial, 4),
  min_concordance = min_conc,
  verdict = if (min_conc <= 0) "not tested (threshold 0)"
            else if (nrow(ha) < min_sites) sprintf("not tested (< %d covered hom-alt sites)", min_sites)
            else if (conc >= min_conc) "pass" else "FAIL")
fwrite(res, out, sep = "\t")
message(sprintf(paste("[germline %s] identity: %d of %d germline sites have reads;",
                      "hom-alt concordance %.3f over %d sites; het biallelic %.3f over %d sites -> %s"),
                label, nrow(m), nrow(sites), conc, nrow(ha), bial, nrow(he), res$verdict))
# Exit 3, not stop()'s 1: a mismatch is deterministic, so conf/base.config
# gives up on the individual at once instead of retrying it.
if (res$verdict == "FAIL") {
  message(sprintf(paste("germline VCF for %s does not look like the same person as its scRNA:",
                     "only %.1f%% of covered hom-alt germline sites read as ALT (threshold %.0f%%).",
                     "Check the patient -> VCF / sample mapping. --numbat_germline_min_concordance 0",
                     "skips this test."), label, 100 * conc, 100 * min_conc))
  quit(status = 3)
}
