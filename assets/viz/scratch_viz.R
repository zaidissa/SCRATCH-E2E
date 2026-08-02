# =============================================================================
# scratch_viz.R — publication figure system for SCRATCH-E2E
# =============================================================================
#
# Sourced by every notebook in the pipeline. It upgrades figures WITHOUT
# requiring the notebooks' plotting code to change, by setting three things
# globally:
#
#   1. theme_set()                     typography, grid, spacing, legends
#   2. ggplot2.discrete.*/continuous.* default colour scales
#   3. a ggsave() shim                 raster + vector, at journal sizes
#
# Any existing `ggplot(...) + geom_*()` therefore inherits the house style, and
# any existing `ggsave("fig.png")` also emits `fig.pdf`. Notebooks that want
# more control can call theme_scratch()/scratch_save() explicitly.
#
# Before this existed the pipeline had 45 ad-hoc theme calls across four theme
# families, DPI hardcoded at 49 sites ranging 150-300, and no vector output at
# all — nothing that could go into a manuscript without being redrawn.
# =============================================================================

suppressPackageStartupMessages({
  library(ggplot2)
  library(grid)
})

# -----------------------------------------------------------------------------
# Palette
#
# Categorical slots are ordered, not cycled: a series keeps its colour when
# other series are filtered out. This order is validated for colour-vision
# deficiency separation and for a normal-vision distinctness floor on both a
# light and a dark ground. Slots beyond 8 are NOT generated — excess levels fold
# to "Other", because an invented 9th hue is not separable from the first eight.
# -----------------------------------------------------------------------------

SCRATCH_PAL <- list(
  categorical = c("#2a78d6", "#eb6834", "#1baf7a", "#eda100",
                  "#e87ba4", "#008300", "#4a3aa7", "#e34948"),

  # Single hue, light to dark. For continuous magnitude.
  sequential  = c("#cde2fb", "#9ec5f4", "#6da7ec", "#3987e5",
                  "#256abf", "#184f95", "#0d366b"),

  # Two poles with a NEUTRAL midpoint — never a hue at the middle, or the
  # zero-crossing reads as a category rather than as "nothing".
  diverging   = list(low = "#2a78d6", mid = "#f0efec", high = "#e34948"),

  # Reserved. Never reused as a series colour, and always paired with a glyph
  # or label so state is never carried by hue alone.
  status      = c(good = "#0ca30c", warning = "#fab219",
                  serious = "#ec835a", critical = "#d03b3b"),

  # Compartment colours, fixed pipeline-wide so the tumour compartment is the
  # same colour in the stratification report and in the cross-stage report.
  compartment = c(malignant = "#eb6834", non_malignant = "#2a78d6",
                  ambiguous = "#8a9095"),

  ink   = "#0f1417", ink2 = "#4e5a60", muted = "#898781",
  grid  = "#e4e6e3", axis = "#b9bfc2", surface = "#ffffff"
)

scratch_colours <- function(n = NULL) {
  p <- SCRATCH_PAL$categorical
  if (is.null(n)) return(p)
  if (n <= length(p)) return(p[seq_len(n)])
  # Never fabricate hues. Recycling is visibly wrong, which is the point:
  # it signals that the data needs folding or facetting instead.
  warning(sprintf("%d levels requested but only %d validated slots exist; fold excess levels to 'Other' or facet.",
                  n, length(p)))
  rep_len(p, n)
}

# -----------------------------------------------------------------------------
# Theme
#
# Journals reproduce figures small. The defaults here assume the figure will be
# printed at 89 mm (single column) and still need to be legible: hairline axes,
# recessive grid, no chartjunk, text in ink tokens rather than series colours.
# -----------------------------------------------------------------------------

theme_scratch <- function(base_size = 9, base_family = "", grid = "y",
                          legend = "right") {
  gx <- grepl("x", grid); gy <- grepl("y", grid)

  theme_minimal(base_size = base_size, base_family = base_family) +
    theme(
      text  = element_text(colour = SCRATCH_PAL$ink2),
      plot.background  = element_rect(fill = SCRATCH_PAL$surface, colour = NA),
      panel.background = element_rect(fill = SCRATCH_PAL$surface, colour = NA),

      panel.grid.major.x = if (gx) element_line(colour = SCRATCH_PAL$grid, linewidth = 0.25) else element_blank(),
      panel.grid.major.y = if (gy) element_line(colour = SCRATCH_PAL$grid, linewidth = 0.25) else element_blank(),
      panel.grid.minor   = element_blank(),

      axis.line  = element_line(colour = SCRATCH_PAL$axis, linewidth = 0.3),
      axis.ticks = element_line(colour = SCRATCH_PAL$axis, linewidth = 0.3),
      axis.ticks.length = unit(2, "pt"),
      axis.text  = element_text(colour = SCRATCH_PAL$muted, size = rel(0.92)),
      axis.title = element_text(colour = SCRATCH_PAL$ink2, size = rel(1.0)),

      plot.title    = element_text(colour = SCRATCH_PAL$ink, face = "bold",
                                   size = rel(1.15), hjust = 0,
                                   margin = margin(b = 3)),
      plot.subtitle = element_text(colour = SCRATCH_PAL$ink2, size = rel(0.95),
                                   hjust = 0, margin = margin(b = 7)),
      plot.caption  = element_text(colour = SCRATCH_PAL$muted, size = rel(0.82),
                                   hjust = 0, margin = margin(t = 7)),
      plot.title.position   = "plot",
      plot.caption.position = "plot",

      legend.position   = legend,
      legend.title      = element_text(colour = SCRATCH_PAL$ink2, size = rel(0.92)),
      legend.text       = element_text(colour = SCRATCH_PAL$ink2, size = rel(0.88)),
      legend.key.size   = unit(9, "pt"),
      legend.background = element_blank(),
      legend.margin     = margin(0, 0, 0, 4),

      strip.background = element_rect(fill = "#f2f4f2", colour = NA),
      strip.text       = element_text(colour = SCRATCH_PAL$ink2, face = "bold",
                                      size = rel(0.92), margin = margin(3, 3, 3, 3)),

      plot.margin = margin(6, 8, 6, 6)
    )
}

# Applied globally: every existing ggplot() in every notebook picks this up.
theme_set(theme_scratch())

# Default colour scales, likewise global. ggplot2 >= 3.3 consults these options
# whenever a plot does not set an explicit scale, so `aes(fill = celltype)`
# gets the validated palette with no code change.
#
# ggplot2 calls these with arguments that differ between versions — 4.0 passes
# `na.value` and `call`, 3.x does not. Supplying our own `na.value` on top of
# that is a duplicate-argument error, and passing `call` to a scale that has no
# such formal is another. So build the argument list, fill defaults only where
# the caller left a gap, and drop anything the target function cannot accept.
.scratch_scale_factory <- function(fn, default_args) {
  force(fn); force(default_args)
  function(...) {
    args <- list(...)
    for (nm in names(default_args)) {
      if (is.null(args[[nm]])) args[[nm]] <- default_args[[nm]]
    }
    ok <- names(formals(fn))
    if (!"..." %in% ok) args <- args[names(args) %in% ok]
    # quote = TRUE is essential: ggplot2 4.0 passes `call` as a quoted language
    # object for error reporting, and do.call() would otherwise EVALUATE it,
    # trying to execute the caller's own expression.
    do.call(fn, args, quote = TRUE)
  }
}

options(
  ggplot2.discrete.fill = .scratch_scale_factory(
    scale_fill_manual,   list(values = SCRATCH_PAL$categorical, na.value = "#c9cdc9")),
  ggplot2.discrete.colour = .scratch_scale_factory(
    scale_colour_manual, list(values = SCRATCH_PAL$categorical, na.value = "#c9cdc9")),
  ggplot2.continuous.fill = .scratch_scale_factory(
    scale_fill_gradientn,   list(colours = SCRATCH_PAL$sequential, na.value = "#eeefee")),
  ggplot2.continuous.colour = .scratch_scale_factory(
    scale_colour_gradientn, list(colours = SCRATCH_PAL$sequential, na.value = "#eeefee"))
)

# -----------------------------------------------------------------------------
# Ready-made scales, for notebooks that want to be explicit
# -----------------------------------------------------------------------------

scale_fill_scratch    <- function(...) scale_fill_manual(..., values = SCRATCH_PAL$categorical, na.value = "#c9cdc9")
scale_colour_scratch  <- function(...) scale_colour_manual(..., values = SCRATCH_PAL$categorical, na.value = "#c9cdc9")
scale_color_scratch   <- scale_colour_scratch

scale_fill_scratch_c  <- function(...) scale_fill_gradientn(..., colours = SCRATCH_PAL$sequential)
scale_colour_scratch_c<- function(...) scale_colour_gradientn(..., colours = SCRATCH_PAL$sequential)

scale_fill_scratch_div <- function(midpoint = 0, ...) {
  d <- SCRATCH_PAL$diverging
  scale_fill_gradient2(low = d$low, mid = d$mid, high = d$high, midpoint = midpoint, ...)
}
scale_colour_scratch_div <- function(midpoint = 0, ...) {
  d <- SCRATCH_PAL$diverging
  scale_colour_gradient2(low = d$low, mid = d$mid, high = d$high, midpoint = midpoint, ...)
}

scale_fill_compartment   <- function(...) scale_fill_manual(..., values = SCRATCH_PAL$compartment)
scale_colour_compartment <- function(...) scale_colour_manual(..., values = SCRATCH_PAL$compartment)

# -----------------------------------------------------------------------------
# Sizing
#
# Journal column widths in millimetres. Sizing to these up front means the
# figure is never rescaled at submission, which is what breaks font sizes.
# -----------------------------------------------------------------------------

SCRATCH_WIDTH <- c(single = 89, onehalf = 120, double = 183)   # mm

scratch_size <- function(width = "single", ratio = 0.75) {
  w <- if (is.character(width)) unname(SCRATCH_WIDTH[width]) else width
  if (is.na(w)) w <- SCRATCH_WIDTH[["single"]]
  list(width = w / 25.4, height = (w * ratio) / 25.4)          # -> inches
}

# -----------------------------------------------------------------------------
# Saving: raster AND vector, always
#
# Journals want vector for line art and >=300 dpi for raster. Every figure is
# written both ways from one call, so nothing has to be regenerated later.
# -----------------------------------------------------------------------------

SCRATCH_FIG_LOG <- new.env(parent = emptyenv())
SCRATCH_FIG_LOG$rows <- list()

scratch_save <- function(plot, name, width = "single", ratio = 0.75,
                         dir = "figures", dpi = 400, vector = TRUE,
                         caption = NULL) {
  dir.create(dir, recursive = TRUE, showWarnings = FALSE)
  sz  <- scratch_size(width, ratio)
  stem <- file.path(dir, sub("\\.(png|pdf|svg|tiff)$", "", name))

  ok_png <- tryCatch({
    ggplot2::ggsave(paste0(stem, ".png"), plot, width = sz$width, height = sz$height,
                    dpi = dpi, units = "in", bg = SCRATCH_PAL$surface)
    TRUE
  }, error = function(e) { message("figure PNG failed: ", conditionMessage(e)); FALSE })

  ok_pdf <- FALSE
  if (isTRUE(vector)) {
    # cairo_pdf embeds fonts, which is what stops journals rejecting the PDF.
    dev <- if (capabilities("cairo")) grDevices::cairo_pdf else NULL
    ok_pdf <- tryCatch({
      if (is.null(dev)) {
        ggplot2::ggsave(paste0(stem, ".pdf"), plot, width = sz$width, height = sz$height, units = "in")
      } else {
        ggplot2::ggsave(paste0(stem, ".pdf"), plot, width = sz$width, height = sz$height,
                        units = "in", device = dev)
      }
      TRUE
    }, error = function(e) { message("figure PDF failed: ", conditionMessage(e)); FALSE })
  }

  SCRATCH_FIG_LOG$rows[[length(SCRATCH_FIG_LOG$rows) + 1]] <- data.frame(
    figure = basename(stem), png = ok_png, pdf = ok_pdf,
    width_mm = round(sz$width * 25.4, 1), height_mm = round(sz$height * 25.4, 1),
    dpi = dpi, caption = caption %||% NA_character_, stringsAsFactors = FALSE
  )
  invisible(plot)
}

`%||%` <- function(a, b) if (is.null(a)) b else a

# Shim: existing `ggsave("fig.png", p)` calls throughout the vendored notebooks
# now also emit `fig.pdf`, with no edit to those notebooks. The raster call is
# delegated verbatim so any arguments they pass still apply.
ggsave <- function(filename, plot = ggplot2::last_plot(), ...) {
  out <- ggplot2::ggsave(filename, plot, ...)
  if (grepl("\\.png$", filename, ignore.case = TRUE)) {
    try({
      pdfname <- sub("\\.png$", ".pdf", filename, ignore.case = TRUE)
      dots <- list(...)
      dots$dpi <- NULL                       # meaningless for vector output
      dev <- if (capabilities("cairo")) grDevices::cairo_pdf else NULL
      if (!is.null(dev)) dots$device <- dev
      do.call(ggplot2::ggsave, c(list(filename = pdfname, plot = plot), dots))
    }, silent = TRUE)
  }
  invisible(out)
}

scratch_figure_manifest <- function(path = "figures/figure_manifest.csv") {
  if (!length(SCRATCH_FIG_LOG$rows)) return(invisible(NULL))
  df <- do.call(rbind, SCRATCH_FIG_LOG$rows)
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  utils::write.csv(df, path, row.names = FALSE)
  invisible(df)
}

# -----------------------------------------------------------------------------
# Helpers used across stages
# -----------------------------------------------------------------------------

# Cap categorical levels at the number of validated slots; the rest fold to
# "Other" rather than generating hues that are not separable.
scratch_fold <- function(x, max_levels = 8L, other = "Other") {
  x <- as.character(x)
  keep <- names(sort(table(x), decreasing = TRUE))
  keep <- keep[seq_len(min(max_levels, length(keep)))]
  factor(ifelse(x %in% keep, x, other), levels = c(keep, other))
}

# A single-series bar chart needs no legend — the title names the measure.
# Direct labels are used up to `max_n`, beyond which they collide.
scratch_direct_labels <- function(p, mapping, max_n = 25, size = 2.4, ...) {
  n <- tryCatch(nrow(p$data), error = function(e) Inf)
  if (is.finite(n) && n <= max_n) p + geom_text(mapping, size = size, colour = SCRATCH_PAL$ink2, ...) else p
}

# Seurat reductions: strip the default heavy styling and apply the house look.
scratch_dimplot <- function(obj, group.by, reduction = "umap", label = TRUE,
                            title = NULL, pt.size = NULL, ...) {
  if (!requireNamespace("Seurat", quietly = TRUE)) stop("Seurat not available")
  n_cells <- ncol(obj)
  if (is.null(pt.size)) pt.size <- if (n_cells > 1e5) 0.05 else if (n_cells > 2e4) 0.2 else 0.5

  p <- Seurat::DimPlot(obj, group.by = group.by, reduction = reduction,
                       label = label, repel = TRUE, raster = n_cells > 1e5,
                       pt.size = pt.size, ...) +
    scale_colour_scratch() +
    labs(title = title %||% group.by, x = NULL, y = NULL) +
    theme_scratch(grid = "none") +
    theme(axis.text = element_blank(), axis.ticks = element_blank(),
          axis.line = element_blank())
  p
}

scratch_featureplot <- function(obj, features, reduction = "umap", ...) {
  if (!requireNamespace("Seurat", quietly = TRUE)) stop("Seurat not available")
  Seurat::FeaturePlot(obj, features = features, reduction = reduction,
                      raster = ncol(obj) > 1e5, ...) &
    scale_colour_scratch_c() &
    theme_scratch(grid = "none") &
    theme(axis.text = element_blank(), axis.ticks = element_blank(),
          axis.line = element_blank())
}

# Heatmap colours for pheatmap/ComplexHeatmap, which do not go through ggplot
# and therefore do not pick up the global scale options.
scratch_heat_colours <- function(n = 256, type = c("sequential", "diverging")) {
  type <- match.arg(type)
  if (type == "sequential") {
    grDevices::colorRampPalette(SCRATCH_PAL$sequential)(n)
  } else {
    d <- SCRATCH_PAL$diverging
    grDevices::colorRampPalette(c(d$low, d$mid, d$high))(n)
  }
}

message("[scratch_viz] publication figure system active — theme, palette, vector export")

# -----------------------------------------------------------------------------
# Locked category colours
#
# A per-plot scale assigns colours from the levels PRESENT in that plot, so the
# same cell type gets different colours in different panels — the single most
# common way a multi-panel report becomes unreadable. Locking the map once, from
# the full set of levels in the dataset, fixes every downstream figure.
# -----------------------------------------------------------------------------

SCRATCH_LOCK <- new.env(parent = emptyenv())

#' Establish (or fetch) the colour map for a categorical variable. Call once per
#' report with every level that exists anywhere in the data.
scratch_lock_levels <- function(levels, key = "default", palette = NULL) {
  levels <- unique(as.character(levels))
  levels <- levels[!is.na(levels)]
  ordered <- sort(levels)                       # stable regardless of input order
  pal <- palette %||% SCRATCH_PAL$categorical
  if (length(ordered) > length(pal)) {
    warning(sprintf("%d levels but %d validated slots; the surplus share a neutral grey — fold them.",
                    length(ordered), length(pal)))
  }
  cols <- rep("#c9cdc9", length(ordered))
  n <- min(length(ordered), length(pal))
  cols[seq_len(n)] <- pal[seq_len(n)]
  m <- setNames(cols, ordered)
  assign(key, m, envir = SCRATCH_LOCK)
  invisible(m)
}

scratch_locked <- function(key = "default") {
  if (!exists(key, envir = SCRATCH_LOCK)) return(NULL)
  get(key, envir = SCRATCH_LOCK)
}

#' Scales that use the locked map. `drop = FALSE` keeps absent levels in the
#' legend, so the legend itself is identical across panels too.
scale_colour_locked <- function(key = "default", drop = FALSE, ...) {
  m <- scratch_locked(key); if (is.null(m)) return(scale_colour_scratch(...))
  scale_colour_manual(values = m, limits = names(m), drop = drop,
                      na.value = "#c9cdc9", ...)
}
scale_color_locked <- scale_colour_locked
scale_fill_locked <- function(key = "default", drop = FALSE, ...) {
  m <- scratch_locked(key); if (is.null(m)) return(scale_fill_scratch(...))
  scale_fill_manual(values = m, limits = names(m), drop = drop,
                    na.value = "#c9cdc9", ...)
}

# -----------------------------------------------------------------------------
# Module scoring, without depending on Seurat::AddModuleScore
#
# Seurat 5.3.0's AddModuleScore calls GetAssayData(slot=), which SeuratObject
# 5.3.0 made defunct — so it fails on ANY object, including a freshly created
# one. Since the pipeline runs inside containers whose exact Seurat/SeuratObject
# pairing cannot be checked ahead of time, the safer course is not to depend on
# it at all.
#
# This implements the same method (Tirosh et al.): bin all genes by mean
# expression, draw `ctrl` control genes per signature gene from the matching
# bin, and score as mean(signature) - mean(controls). Deterministic given a
# seed, and identical in behaviour across Seurat versions.
# -----------------------------------------------------------------------------

scratch_module_score <- function(obj, features, name = "Score", nbin = 24,
                                 ctrl = 100, seed = 1234, assay = NULL,
                                 layer = "data") {
  stopifnot(is.list(features))
  set.seed(seed)
  assay <- assay %||% SeuratObject::DefaultAssay(obj)

  mat <- tryCatch(SeuratObject::LayerData(obj, assay = assay, layer = layer),
                  error = function(e)
                    SeuratObject::GetAssayData(obj, assay = assay, layer = layer))

  present <- rownames(mat)
  features <- lapply(features, function(g) intersect(unique(g), present))
  features <- features[vapply(features, length, integer(1)) > 0]
  if (!length(features)) {
    warning("scratch_module_score: no signature gene is present in the object")
    return(obj)
  }

  # Bin genes by mean expression. The tiny jitter reproduces Seurat's tie-break
  # so genes with identical means do not all land in one bin.
  gene_mean <- Matrix::rowMeans(mat)
  jitter    <- stats::rnorm(length(gene_mean)) / 1e30
  breaks    <- stats::quantile(gene_mean + jitter, probs = seq(0, 1, length.out = nbin + 1),
                               na.rm = TRUE)
  bins <- cut(gene_mean + jitter, breaks = unique(breaks), labels = FALSE,
              include.lowest = TRUE)
  names(bins) <- names(gene_mean)

  scores <- list()
  for (i in seq_along(features)) {
    fs <- features[[i]]
    ctrl_use <- unique(unlist(lapply(fs, function(g) {
      pool <- names(bins)[bins == bins[[g]]]
      pool <- setdiff(pool, fs)
      if (!length(pool)) return(character(0))
      sample(pool, min(ctrl, length(pool)))
    })))

    sig_mean  <- Matrix::colMeans(mat[fs, , drop = FALSE])
    ctrl_mean <- if (length(ctrl_use)) Matrix::colMeans(mat[ctrl_use, , drop = FALSE]) else 0
    scores[[paste0(name, i)]] <- as.numeric(sig_mean - ctrl_mean)
  }

  df <- as.data.frame(scores)
  rownames(df) <- colnames(mat)
  obj <- SeuratObject::AddMetaData(obj, df)
  attr(obj, "scratch_score_map") <- setNames(names(scores), names(features))
  obj
}
