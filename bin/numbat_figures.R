#!/usr/bin/env Rscript
# ---------------------------------------------------------------------------
# Numbat's own diagnostic figures, drawn after the fact from a finished
# run_numbat output directory.
#
# Why a separate script rather than a few lines inside numbat_run.R:
#   run_numbat() is called with plot = FALSE, so the pipeline has never drawn
#   any of Numbat's figures -- the CNV report carried exactly one Numbat plot
#   (segment size by chromosome), which the concordance notebook builds itself
#   from segs_consensus. Everything in the Numbat vignette was missing.
#   The obvious fix is plot = TRUE, and it is the wrong one twice over: it
#   produces only 5 of the figures (no phylogeny heatmap, no single-cell tree,
#   no mutation history, no consensus plot), and Nextflow hashes a process's
#   script text, so changing NUMBAT_RUN re-runs the model fit for every sample.
#   On GBM-DFCI3 that fit cost 957 CPU-hours over 15 libraries.
#
#   Every figure below is derivable from files run_numbat has ALREADY written,
#   so they are drawn here instead: a separate, minutes-long task reading the
#   cached output directory. Nothing upstream is re-computed.
#
# What is on disk to draw from. run_numbat writes these OUTSIDE its `if (plot)`
# branch, so plot = FALSE costs no data -- only the pictures:
#   joint_post_<i>.tsv   exp_post_<i>.tsv    allele_post_<i>.tsv
#   clone_post_<i>.tsv   segs_consensus_<i>.tsv   geno_<i>.tsv
#   treeML_<i>.rds       tree_final_<i>.rds  mut_graph_<i>.rds
#   bulk_clones_final.tsv.gz   gexp_roll_wide.tsv.gz   hc.rds
# (gexp_roll_wide and hc exist only on the init_k > 1 path, which is the
#  default; bulk_subtrees is absent when segs_consensus_fix was supplied.)
#
# Every figure is drawn independently and a failure is RECORDED, not fatal.
# A sample whose phylogeny is degenerate must still yield its segment plots,
# and a reader must be able to tell "not drawn" from "nothing to draw": the
# manifest names each figure, its status, and the reason when there is one.
#
# Usage: numbat_figures.R <out_dir> <sample_id> <figdir> <datadir>
#                         [genome] [n_cut_max] [min_clone_cells] [min_llr]
#                         [exp_roll_k] [max_vector_cells]
# ---------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(numbat)
  library(ggplot2)
  library(data.table)
  library(dplyr)
  library(patchwork)
})

a <- commandArgs(trailingOnly = TRUE)
if (length(a) < 4L) {
  stop("usage: numbat_figures.R <out_dir> <sample_id> <figdir> <datadir> [genome] [n_cut_max] [min_clone_cells] [min_llr] [exp_roll_k] [max_vector_cells]")
}

out_dir    <- a[[1]]
sample_id  <- a[[2]]
figdir     <- a[[3]]
datadir    <- a[[4]]
genome     <- if (length(a) >= 5)  a[[5]]            else "hg38"
n_cut_max  <- if (length(a) >= 6)  as.integer(a[[6]])  else 4L
min_cells  <- if (length(a) >= 7)  as.integer(a[[7]])  else 50L
min_llr    <- if (length(a) >= 8)  as.numeric(a[[8]])  else 10
exp_roll_k <- if (length(a) >= 9)  as.integer(a[[9]])  else 3L
# Above this many cells no PDF is written for the dense heatmaps: one tile per
# cell per segment, so 10,000 cells x 60 segments is 600,000 vector rectangles
# -- a PDF no viewer will open. The PNG is still drawn at full size.
max_vec    <- if (length(a) >= 10) as.integer(a[[10]]) else 3000L

dir.create(figdir,  recursive = TRUE, showWarnings = FALSE)
dir.create(datadir, recursive = TRUE, showWarnings = FALSE)

# The transcript annotation the plotting functions place segments against. An
# hg19 gtf against hg38 coordinates does not error, it draws the wrong picture.
gtf_use <- switch(genome, hg38 = numbat::gtf_hg38, hg19 = numbat::gtf_hg19,
                  stop(sprintf("unsupported genome '%s' (expected hg38 or hg19)", genome)))

# ---- what is actually in the directory ------------------------------------
# Numbat suffixes most files with the iteration it converged on. numbat_run.R
# already takes the HIGHEST clone_post index rather than assuming 2; do the
# same here so the two always describe the same iteration.
files <- list.files(out_dir)
iters <- as.integer(sub("^joint_post_(\\d+)\\.tsv$", "\\1",
                        grep("^joint_post_\\d+\\.tsv$", files, value = TRUE)))
if (!length(iters)) {
  message(sprintf("[%s] no joint_post_<i>.tsv in %s -- Numbat produced no posterior, so there is nothing to draw.",
                  sample_id, out_dir))
  fwrite(data.table(sample = sample_id, figure = NA_character_, status = "no_input",
                    reason = "no joint_post_<i>.tsv in the Numbat output directory"),
         file.path(datadir, paste0(sample_id, "_numbat_figure_manifest.tsv")), sep = "\t")
  quit(save = "no", status = 0)
}
i <- max(iters)
message(sprintf("[%s] %d file(s) in %s; drawing iteration %d", sample_id, length(files), out_dir, i))

p <- function(...) file.path(out_dir, paste0(...))
rd_tsv <- function(f) if (file.exists(f)) tryCatch(fread(f), error = function(e) NULL) else NULL
rd_rds <- function(f) if (file.exists(f)) tryCatch(readRDS(f), error = function(e) NULL) else NULL

# `Numbat$new()` is the documented entry point and does two things worth having:
# it relevels CHROM to 1:22 (ggplot would otherwise order chromosomes
# lexically -- 1, 10, 11, ... 2) and it lazily loads gexp_roll_wide. It also
# dies outright on a missing geno_<i>.tsv, because that read is piped straight
# into column_to_rownames(). So: try it, and fall back to reading the files
# directly. The fallback is not dead code -- a sample whose fit stopped before
# genotyping reaches it.
nb <- tryCatch(numbat::Numbat$new(out_dir = out_dir, i = i, gtf = gtf_use, verbose = FALSE),
               error = function(e) { message(sprintf("[%s] Numbat$new() failed (%s); reading files directly.",
                                                     sample_id, conditionMessage(e))); NULL })

relevel_chrom <- function(d) {
  if (!is.null(d) && "CHROM" %in% names(d)) d[, CHROM := factor(as.character(CHROM), levels = as.character(1:22))]
  d
}

if (!is.null(nb)) {
  joint_post     <- nb$joint_post
  segs_consensus <- nb$segs_consensus
  clone_post     <- nb$clone_post
  bulk_clones    <- nb$bulk_clones
  gtree          <- nb$gtree
  mut_graph      <- nb$mut_graph
  hc             <- nb$hc
} else {
  joint_post     <- relevel_chrom(rd_tsv(p("joint_post_", i, ".tsv")))
  segs_consensus <- relevel_chrom(rd_tsv(p("segs_consensus_", i, ".tsv")))
  clone_post     <- rd_tsv(p("clone_post_", i, ".tsv"))
  bulk_clones    <- relevel_chrom(rd_tsv(p("bulk_clones_final.tsv.gz")))
  gtree          <- rd_rds(p("tree_final_", i, ".rds"))
  mut_graph      <- rd_rds(p("mut_graph_", i, ".rds"))
  hc             <- rd_rds(p("hc.rds"))
}

n_cells <- if (!is.null(joint_post) && "cell" %in% names(joint_post)) uniqueN(joint_post$cell) else NA_integer_
n_segs  <- if (!is.null(segs_consensus)) nrow(segs_consensus) else 0L
message(sprintf("[%s] %s cell(s), %d consensus segment(s)",
                sample_id, if (is.na(n_cells)) "?" else format(n_cells, big.mark = ","), n_segs))

# ---- figure harness -------------------------------------------------------
manifest <- list()

have <- function(x) !is.null(x) && (!is.data.frame(x) || nrow(x) > 0)

# A sample above max_vec cells gets PNG only. The tempting alternative --
# numbat's own `raster = TRUE`, which rasterises just the heatmap layer -- is
# NOT available: it calls ggrastr::rasterize(), and ggrastr is a Suggests of
# numbat that pkharchenkolab/numbat-rbase does not install (verified in the
# image). Passing raster = TRUE would therefore fail on exactly the large
# samples it exists for, so the PNG is drawn unrasterised and only the vector
# copy is dropped.
dense <- !is.na(n_cells) && n_cells > max_vec

save_fig <- function(name, plt, width, height, vector = TRUE) {
  png_path <- file.path(figdir, sprintf("%s_%s.png", sample_id, name))
  dev_png  <- if (requireNamespace("ragg", quietly = TRUE)) ragg::agg_png else "png"
  ggsave(png_path, plt, width = width, height = height, dpi = 200,
         device = dev_png, limitsize = FALSE)
  if (vector) {
    # A PDF that fails to write must not discard the PNG that already did.
    tryCatch(ggsave(file.path(figdir, sprintf("%s_%s.pdf", sample_id, name)),
                    plt, width = width, height = height, device = "pdf", limitsize = FALSE),
             error = function(e) message(sprintf("[%s] %s: PDF skipped (%s)",
                                                 sample_id, name, conditionMessage(e))))
  }
  invisible(png_path)
}

fig <- function(name, need, width, height, expr, vector = TRUE) {
  missing <- names(need)[!vapply(need, have, logical(1))]
  if (length(missing)) {
    manifest[[name]] <<- data.table(sample = sample_id, figure = name, status = "skipped",
                                    reason = paste("missing:", paste(missing, collapse = ", ")))
    message(sprintf("[%s] %-24s SKIPPED (missing %s)", sample_id, name, paste(missing, collapse = ", ")))
    return(invisible(NULL))
  }
  res <- try(expr, silent = TRUE)
  if (inherits(res, "try-error")) {
    manifest[[name]] <<- data.table(sample = sample_id, figure = name, status = "failed",
                                    reason = trimws(gsub("\\s+", " ", as.character(res))))
    message(sprintf("[%s] %-24s FAILED: %s", sample_id, name,
                    trimws(gsub("\\s+", " ", as.character(res)))))
    return(invisible(NULL))
  }
  ok <- try(save_fig(name, res, width, height, vector = vector), silent = TRUE)
  if (inherits(ok, "try-error")) {
    manifest[[name]] <<- data.table(sample = sample_id, figure = name, status = "failed",
                                    reason = trimws(gsub("\\s+", " ", as.character(ok))))
    message(sprintf("[%s] %-24s FAILED on write: %s", sample_id, name,
                    trimws(gsub("\\s+", " ", as.character(ok)))))
  } else {
    manifest[[name]] <<- data.table(sample = sample_id, figure = name, status = "drawn", reason = NA_character_)
    message(sprintf("[%s] %-24s drawn", sample_id, name))
  }
  invisible(NULL)
}

# The clone colour bar is drawn as a STACKED bar of per-clone probabilities,
# which plot_stack_bar() reads from clone_post's p_1, p_2, ... columns. Those
# come from get_clone_post(); a clone_post written by a different numbat
# version, or one carrying only the hard assignment, has none, and
# plot_phylo_heatmap then dies in str_remove() with "object 'clone' not
# found" -- losing the whole figure over its narrowest strip. Measured against
# a clone_post holding clone_opt but no p_<k>. Fall back to the hard
# assignment instead of failing.
clone_stack <- !is.null(clone_post) &&
               any(grepl("^p_[0-9]+$", names(clone_post)))

# 1. The signature figure: single-cell phylogeny beside the per-cell CNV calls,
#    with a clone colour bar and the tumour/normal divide marked. This is the
#    one plot that shows the CLONAL STRUCTURE -- whether the malignant cells
#    are one population or several, and which CNV separates them.
fig("phylo_heatmap",
    list(gtree = gtree, joint_post = joint_post, segs_consensus = segs_consensus),
    width = 12, height = max(5, min(12, 4 + n_segs / 15)),
    vector = !dense,
    plot_phylo_heatmap(gtree = gtree, joint_post = joint_post,
                       segs_consensus = segs_consensus, clone_post = clone_post,
                       clone_bar = TRUE, clone_stack = clone_stack,
                       p_min = 0.9))

# 2. The same heatmap re-cut at 1..n clones. Numbat's clone count is a CHOICE,
#    not a measurement: this is the figure that shows what that choice costs.
#    Needs the class -- cutree() recomputes the tree from treeML and the
#    genotype matrix, which the file-reading fallback does not carry.
if (n_cut_max >= 2L) {
  fig("phylo_heatmap_by_ncut",
      list(numbat_object = nb, joint_post = joint_post, segs_consensus = segs_consensus),
      width = 7 * min(2, n_cut_max), height = 5 * ceiling(n_cut_max / 2),
      vector = FALSE,
      {
        plots <- lapply(seq_len(n_cut_max), function(k) {
          nb$cutree(n_cut = k)
          nb$plot_phylo_heatmap() + ggtitle(paste0("n_cut = ", k))
        })
        # cutree() mutates the object, so put it back the way it was found;
        # anything drawn after this would otherwise silently use n_cut_max.
        nb$cutree(n_cut = 0)
        wrap_plots(plots, ncol = min(2, n_cut_max))
      })
}

# 3. Consensus segments along the genome: the EVENTS Numbat resolved, as
#    opposed to the cells carrying them.
fig("consensus_segments", list(segs_consensus = segs_consensus),
    width = 11, height = 3,
    plot_consensus(segs = segs_consensus))

# 4. The same segments as a genome-wide heatmap, one row per group, which is
#    how a reader compares two clones' profiles side by side.
fig("segment_heatmap", list(segs_consensus = segs_consensus),
    width = 11, height = max(2.5, 1 + 0.4 * uniqueN(segs_consensus[["group"]])),
    cnv_heatmap(segs = segs_consensus, genome = genome))

# 5. Pseudobulk profiles per clone: expression log-ratio on top, phased allele
#    frequency underneath. This is where the ALLELIC evidence is visible, and
#    it is the only figure that shows why Numbat can call copy-neutral LOH
#    when the expression-only callers cannot.
fig("bulk_clone_profiles", list(bulk_clones = bulk_clones),
    width = 12,
    height = {
      k <- if ("sample" %in% names(bulk_clones)) uniqueN(bulk_clones$sample) else 1L
      max(3, 2.6 * max(1L, k))
    },
    vector = FALSE,
    {
      b <- bulk_clones
      if ("n_cells" %in% names(b)) b <- b %>% filter(n_cells > min_cells)
      if (!nrow(b)) stop(sprintf("no clone has more than %d cells", min_cells))
      plot_bulks(b, min_LLR = min_llr, legend = TRUE)
    })

# 6. The single-cell phylogeny on its own, with each CNV placed on the branch
#    where it was acquired -- the ORDER of events, which the heatmap cannot show.
fig("sc_phylogeny", list(gtree = gtree),
    width = 9, height = 7,
    plot_sc_tree(gtree = gtree, label_size = 3, branch_width = 0.5, tip_length = 0.5, tip = TRUE))

# 7. The same history with cells of one genotype collapsed to a single node:
#    the clone tree, sized by how many cells each clone holds.
fig("mutation_history", list(mut_graph = mut_graph),
    width = 9, height = 5,
    plot_mut_history(G = mut_graph, clone_post = clone_post))

# 8. Window-smoothed expression, clustered. This is the EXPRESSION channel on
#    its own, so a reader can see how much of the call it carried -- and,
#    where it disagrees with the allelic channel, which one Numbat followed.
fig("expression_rolling",
    list(gexp_roll = if (!is.null(nb)) TRUE else NULL, hc = hc),
    width = 11, height = 6,
    vector = FALSE,
    nb$plot_exp_roll(k = exp_roll_k, n_sample = min(300L, if (is.na(n_cells)) 300L else n_cells)))

# ---- manifest -------------------------------------------------------------
man <- rbindlist(manifest, fill = TRUE)
if (!nrow(man)) man <- data.table(sample = sample_id, figure = NA_character_,
                                  status = "none", reason = "no figure attempted")
man[, cells := n_cells][, segments := n_segs][, iteration := i]
fwrite(man, file.path(datadir, paste0(sample_id, "_numbat_figure_manifest.tsv")), sep = "\t")

message(sprintf("[%s] %d drawn, %d skipped, %d failed", sample_id,
                man[status == "drawn", .N], man[status == "skipped", .N], man[status == "failed", .N]))

# A sample that drew nothing is a RESULT, not a task failure: it means Numbat
# found no clonal structure in it. Exiting non-zero here would take the whole
# cohort's figures down with it.
quit(save = "no", status = 0)
