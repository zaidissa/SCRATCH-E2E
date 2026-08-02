set.seed(7)
samples <- c("GBM_DFCI1_CSF","GBM_DFCI2_CSF","GBM_DFCI3","GBM_DFCI4","GBM_DFCI5",
             "GBM_MSK1","GBM_MSK2","OV_TCGA1","OV_TCGA2","PDAC_01","PDAC_02","PDAC_03")
lin <- c("T_Cells","Myeloid","NK_Cells","B_Plasma_Cells","Endothelial_Cells","Fibroblast","Malignant")
qc <- do.call(rbind, lapply(seq_along(samples), function(i) {
  s <- samples[i]
  csf <- grepl("CSF", s)
  n  <- if (csf) round(runif(1, 380, 1400)) else round(runif(1, 4200, 14000))
  data.frame(sample=s, cohort=sub("_.*","",s), csf=csf, n_cells=n,
    median_genes = round(rnorm(1, if (csf) 780 else 2100, 250)),
    median_umi   = round(rnorm(1, if (csf) 1900 else 6400, 700)),
    pct_mito     = round(pmax(1, rnorm(1, if (csf) 14 else 6.5, 4)), 1),
    doublet_rate = round(pmax(0.5, rnorm(1, 5.2, 2)), 1))
}))
cells <- do.call(rbind, lapply(seq_len(nrow(qc)), function(i) {
  r <- qc[i,]; k <- min(r$n_cells, 900)
  p <- if (r$csf) c(.42,.30,.08,.05,.02,.03,.10) else c(.22,.20,.05,.06,.07,.10,.30)
  data.frame(sample=r$sample, cohort=r$cohort,
             lineage=sample(lin, k, TRUE, prob=p),
             genes=round(pmax(200, rnorm(k, r$median_genes, r$median_genes*.28))),
             mito =pmax(0, rnorm(k, r$pct_mito, 3)),
             umap1=rnorm(k), umap2=rnorm(k))
}))
saveRDS(list(qc=qc, cells=cells, samples=samples), "demo.rds")
