#!/usr/bin/env Rscript

# Build a smaller Seurat object for testing that still LOOKS like the cohort.
#
# The naive recipe -- "keep at most N cells per sample x lineage" -- silently
# rewrites the biology. Capping is uniform in absolute terms, so it truncates
# exactly the strata that are large and leaves small ones untouched: on the
# HGSOC test object it took Epithelial from 22.8% to 7.8%, which would make any
# downstream composition, differential-abundance or malignant-fraction result
# meaningless while looking perfectly healthy.
#
# Proportional allocation instead gives every stratum a share of the budget in
# proportion to its true size, so cohort composition survives the shrink. A
# floor keeps rare strata from vanishing entirely -- that is a deliberate,
# reported distortion (rare types end up slightly over-represented), which is
# the opposite of the silent one above.

suppressPackageStartupMessages({
  library(Seurat); library(SeuratObject); library(Matrix); library(optparse)
})

option_list <- list(
  make_option(c("-f", "--file"), type = "character",
              help = "Input .rds/.RDS Seurat object"),
  make_option(c("-o", "--out"), type = "character",
              help = "Output .rds path"),
  make_option(c("-n", "--n_cells"), type = "integer", default = 12000L,
              help = "Target cell count [default %default]"),
  make_option(c("--strata"), type = "character", default = "sample_id,celltype",
              help = "Comma-separated metadata columns defining strata [default %default]"),
  make_option(c("--floor"), type = "integer", default = 10L,
              help = "Minimum cells kept per non-empty stratum [default %default]"),
  make_option(c("--seed"), type = "integer", default = 1L, help = "RNG seed"),
  make_option(c("--report"), type = "character", default = NULL,
              help = "Optional CSV path for the composition-drift table"),
  make_option(c("--normalize"), action = "store_true", default = FALSE,
              help = "NormalizeData + FindVariableFeatures on the subset")
)
opt <- parse_args(OptionParser(option_list = option_list))
if (is.null(opt$file) || is.null(opt$out))
  stop("--file and --out are required", call. = FALSE)

# ---- allocation ------------------------------------------------------------
# Proportional (largest-remainder) allocation of n_target across strata, with a
# per-stratum floor and a per-stratum cap at what is actually available.
# Saturated strata hand their surplus back to the others, so the budget is spent
# rather than quietly lost.
allocate_proportional <- function(avail, n_target, floor_n) {
  if (sum(avail) <= n_target) return(avail)          # nothing to shrink

  base <- pmin(avail, floor_n)
  if (sum(base) >= n_target) {
    warning("floor (", floor_n, " x ", length(avail), " strata = ", sum(base),
            " cells) meets or exceeds the ", n_target, "-cell budget; returning ",
            "the floor allocation. Composition is NOT preserved -- lower --floor ",
            "or raise --n_cells.", immediate. = TRUE)
    return(base)
  }

  extra <- n_target - sum(base)
  cap   <- avail - base                 # room left above the floor
  w     <- avail / sum(avail)           # weights from TRUE sizes, not from `cap`
  alloc <- base

  repeat {
    open <- cap > 0
    if (!any(open) || extra <= 0L) break
    ww  <- w[open] / sum(w[open])
    raw <- extra * ww
    add <- pmin(floor(raw), cap[open])

    left <- extra - sum(add)            # largest-remainder for the rounding dust
    if (left > 0) {
      room <- cap[open] - add
      for (i in order(raw - floor(raw), decreasing = TRUE)) {
        if (left == 0) break
        if (room[i] > 0) { add[i] <- add[i] + 1; room[i] <- room[i] - 1; left <- left - 1 }
      }
    }

    alloc[open] <- alloc[open] + add
    cap[open]   <- cap[open]   - add
    gained <- sum(add)
    extra  <- extra - gained
    if (gained == 0) break              # every open stratum saturated; stop
  }
  alloc
}

# ---- load ------------------------------------------------------------------
cat("reading", opt$file, "\n")
obj <- readRDS(opt$file)
obj <- tryCatch(UpdateSeuratObject(obj), error = function(e) obj)
md  <- obj@meta.data

strata_cols <- trimws(strsplit(opt$strata, ",")[[1]])
missing <- setdiff(strata_cols, colnames(md))
if (length(missing))
  stop("stratum column(s) not in the object: ", paste(missing, collapse = ", "),
       "\n  available: ", paste(colnames(md), collapse = ", "), call. = FALSE)

key <- do.call(paste, c(lapply(strata_cols, function(c) as.character(md[[c]])), sep = "\r"))
idx_by <- split(seq_len(nrow(md)), key)
avail  <- lengths(idx_by)

cat(sprintf("%d cells in %d strata (%s)\n", nrow(md), length(avail),
            paste(strata_cols, collapse = " x ")))

# ---- allocate & sample -----------------------------------------------------
set.seed(opt$seed)
alloc <- allocate_proportional(avail, opt$n_cells, opt$floor)
keep  <- sort(unlist(Map(function(ii, k) if (k >= length(ii)) ii else sample(ii, k),
                         idx_by, alloc[names(idx_by)]), use.names = FALSE))
cat(sprintf("keeping %d cells (target %d)\n", length(keep), opt$n_cells))

# ---- composition drift -----------------------------------------------------
# The whole point of the exercise: show what moved, per stratum column.
drift_all <- list()
for (cl in strata_cols) {
  before <- prop.table(table(md[[cl]])) * 100
  after  <- prop.table(table(md[[cl]][keep])) * 100
  lev    <- union(names(before), names(after))
  d <- data.frame(
    column = cl, level = lev,
    n_before = as.integer(table(md[[cl]])[lev]),
    pct_before = round(as.numeric(before[lev]), 2),
    n_after = as.integer(table(md[[cl]][keep])[lev]),
    pct_after = round(as.numeric(after[lev]), 2), row.names = NULL
  )
  d[is.na(d)] <- 0
  d$delta_pp <- round(d$pct_after - d$pct_before, 2)
  drift_all[[cl]] <- d[order(-d$pct_before), ]
}
for (cl in strata_cols) {
  cat("\n-- composition drift:", cl, "--\n")
  print(drift_all[[cl]], row.names = FALSE)
  cat(sprintf("   max |drift| = %.2f pp\n", max(abs(drift_all[[cl]]$delta_pp))))
}
drift <- do.call(rbind, drift_all)
if (!is.null(opt$report)) { write.csv(drift, opt$report, row.names = FALSE)
                            cat("\nwrote", opt$report, "\n") }

# ---- rebuild ---------------------------------------------------------------
# subset() drops every cell on objects whose `data` layer is an empty
# placeholder (raw objects written by older Seurat), so rebuild from counts.
cts <- tryCatch(LayerData(obj, layer = "counts"),
                error = function(e) GetAssayData(obj, slot = "counts"))
sub <- CreateSeuratObject(counts = cts[, keep, drop = FALSE],
                          meta.data = md[keep, , drop = FALSE],
                          project = Project(obj))

# carry reductions over, restricted to the kept cells
for (rd in Reductions(obj)) {
  emb <- Embeddings(obj, rd)
  emb <- emb[colnames(sub)[colnames(sub) %in% rownames(emb)], , drop = FALSE]
  if (nrow(emb) == ncol(sub))
    sub[[rd]] <- CreateDimReducObject(embeddings = emb, key = Key(obj[[rd]]), assay = "RNA")
}

if (opt$normalize) {
  sub <- NormalizeData(sub, verbose = FALSE)
  sub <- FindVariableFeatures(sub, nfeatures = 2000, verbose = FALSE)
}

saveRDS(sub, opt$out)
cat("\nsaved", opt$out, ":", ncol(sub), "cells x", nrow(sub), "genes,",
    length(Reductions(sub)), "reductions\n")
