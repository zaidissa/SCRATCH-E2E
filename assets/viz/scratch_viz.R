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

# Tracks which level-counts have already produced the palette advisory, so it is
# emitted once per session instead of once per plot.
.scratch_pal_warned <- new.env(parent = emptyenv())

scratch_colours <- function(n = NULL) {
  p <- SCRATCH_PAL$categorical
  if (is.null(n)) return(p)
  if (n <= length(p)) return(p[seq_len(n)])
  # Past the validated slots there is no good answer, only less-bad ones.
  # Recycling produces two levels with the SAME colour, which silently misreads
  # as "these are the same group" -- worse than an unvalidated hue. Interpolating
  # keeps every level distinguishable and the warning makes the loss of the
  # colour-vision guarantee explicit.
  # Once per level-count per session, NOT once per plot. A 7-panel VlnPlot fired
  # it seven times and knitr printed all seven INTO the report, above the figure
  # -- an advisory about colour-vision safety became the largest block of text on
  # the page. The advice is worth giving; it is not worth giving seven times.
  .key <- paste0("n", n)
  if (is.null(.scratch_pal_warned[[.key]])) {
    assign(.key, TRUE, envir = .scratch_pal_warned)
    warning(sprintf(paste0("%d levels requested but only %d colour-vision-validated slots exist. ",
                           "Extra hues are interpolated and NOT CVD-validated -- fold excess ",
                           "levels to 'Other' or facet for a publication figure."),
                    n, length(p)), call. = FALSE)
  }
  grDevices::colorRampPalette(p)(n)
}

# -----------------------------------------------------------------------------
# Theme
#
# Journals reproduce figures small. The defaults here assume the figure will be
# printed at 89 mm (single column) and still need to be legible: hairline axes,
# recessive grid, no chartjunk, text in ink tokens rather than series colours.
# -----------------------------------------------------------------------------

# base_size 9 was too small to read at full report width. Everything else in the
# theme is expressed with rel(), so raising this one number lifts titles, axis
# text and legend text together -- except the two absolute `unit(..., "pt")`
# values below, which are now tied to base_size so legend keys and tick marks
# scale with the type instead of shrinking against it.
theme_scratch <- function(base_size = 12, base_family = "", grid = "y",
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
      axis.ticks.length = unit(base_size / 4.5, "pt"),
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
      legend.key.size   = unit(base_size, "pt"),
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
  # A palette FUNCTION, so the scale adapts to however many levels the data has.
  # Using scale_*_manual(values = <fixed 8>) here made every plot with 9+ levels
  # fail outright with "Insufficient values in manual scale".
  # ggplot2 supplies na.value (and `call`) through `...`, so these must be set
  # only when absent -- passing them positionally as well gives
  # "formal argument na.value matched by multiple actual arguments".
  # quote = TRUE for the same reason as the factory above: `call` arrives as a
  # language object and must not be evaluated.
  ggplot2.discrete.fill = function(...) {
    a <- list(...)
    if (is.null(a$na.value)) a$na.value <- "#c9cdc9"
    a$aesthetics <- "fill"; a$palette <- scratch_colours
    do.call(ggplot2::discrete_scale, a, quote = TRUE)
  },
  ggplot2.discrete.colour = function(...) {
    a <- list(...)
    if (is.null(a$na.value)) a$na.value <- "#c9cdc9"
    a$aesthetics <- "colour"; a$palette <- scratch_colours
    do.call(ggplot2::discrete_scale, a, quote = TRUE)
  },
  ggplot2.continuous.fill = .scratch_scale_factory(
    scale_fill_gradientn,   list(colours = SCRATCH_PAL$sequential, na.value = "#eeefee")),
  ggplot2.continuous.colour = .scratch_scale_factory(
    scale_colour_gradientn, list(colours = SCRATCH_PAL$sequential, na.value = "#eeefee"))
)

# -----------------------------------------------------------------------------
# Ready-made scales, for notebooks that want to be explicit
# -----------------------------------------------------------------------------

# `scratch_colours(n)` already interpolates past the validated slots and warns
# once per level-count. These scales did not use it: they handed
# scale_*_manual a FIXED vector, so any plot with more than
# length(SCRATCH_PAL$categorical) levels died with
#   "Insufficient values in manual scale. 11 needed but only 8 provided."
# and took the whole notebook down with it. A palette that errors on a category
# count nobody chose is a worse failure than an unvalidated hue -- which is
# exactly the trade-off scratch_colours() was written to make.
#
# Using the discrete-scale palette FUNCTION form means the level count is known
# at draw time, so the interpolation happens only when it is actually needed.
scale_fill_scratch <- function(...) {
  ggplot2::discrete_scale("fill", palette = function(n) scratch_colours(n),
                          na.value = "#c9cdc9", ...)
}
scale_colour_scratch <- function(...) {
  ggplot2::discrete_scale("colour", palette = function(n) scratch_colours(n),
                          na.value = "#c9cdc9", ...)
}
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
    try(scratch_interactive_twin(plot, paste0(stem, ".png")), silent = TRUE)
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

# Interactive twin for EVERY figure.
#
# The reporting layer's interactive components only ever reached the handful of
# notebooks that were rewritten to call them; the other ~38 vendored notebooks
# kept emitting flat PNGs. Rather than edit 40 notebooks, convert at the point
# every one of them already goes through: ggsave(). Any `ggsave("fig.png", p)`
# now also writes `fig.html`, a hoverable plotly version of the same plot.
#
# plotly.js is written ONCE into a shared `_figlibs/` directory rather than
# inlined per figure -- self-contained widgets are ~3 MB each, which across a
# cohort's worth of figures is hundreds of MB of duplicated JavaScript.
scratch_interactive_twin <- function(plot, filename) {
  if (!requireNamespace("plotly", quietly = TRUE) ||
      !requireNamespace("htmlwidgets", quietly = TRUE)) return(invisible(FALSE))
  if (!inherits(plot, "ggplot")) return(invisible(FALSE))
  if (!isTRUE(getOption("scratch.interactive", TRUE))) return(invisible(FALSE))

  html <- sub("\\.(png|jpe?g|tiff?)$", ".html", filename, ignore.case = TRUE)
  ok <- try({
    w <- plotly::ggplotly(plot)
    w <- plotly::config(w, displaylogo = FALSE, responsive = TRUE,
                        modeBarButtonsToRemove = c("lasso2d", "select2d"))
    # libdir is resolved RELATIVE to the widget file, so pass a bare name --
    # passing dirname(html)/_figlibs nests it as figures/figures/_figlibs.
    htmlwidgets::saveWidget(w, html, selfcontained = FALSE,
                            libdir = "_figlibs", title = basename(html))
    TRUE
  }, silent = TRUE)

  if (inherits(ok, "try-error")) {
    # Never fail a render over a figure twin -- some geoms have no plotly
    # equivalent and that is acceptable; the static figure still exists.
    return(invisible(FALSE))
  }
  SCRATCH_FIG_LOG$interactive <- c(SCRATCH_FIG_LOG$interactive, html)
  invisible(TRUE)
}

# Shim: existing `ggsave("fig.png", p)` calls throughout the vendored notebooks
# now also emit `fig.pdf` AND `fig.html`, with no edit to those notebooks. The
# raster call is delegated verbatim so any arguments they pass still apply.
ggsave <- function(filename, plot = ggplot2::last_plot(), ...) {
  plot <- tryCatch(scratch_autoscale_text(plot), error = function(e) plot)
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
    try(scratch_interactive_twin(plot, filename), silent = TRUE)
  }
  invisible(out)
}

# Automatic text scaling.
#
# ggplot text sizes are absolute points, so a theme tuned for a 6-category bar
# chart renders unreadably crowded once the same code is handed 39 samples --
# the labels do not shrink, they just collide. Every figure in this pipeline is
# produced by vendored code that cannot know how many categories it will get at
# runtime, so the adjustment has to happen here, after the plot is built and the
# real category count is knowable.
#
# Only text size and x-label rotation change; geometry, data and colour are
# untouched.
scratch_autoscale_text <- function(plot) {
  if (!inherits(plot, "ggplot")) return(plot)
  # In current ggplot2 `panel_params$x$get_labels` is a FUNCTION, so taking
  # length() of it yields 1 and the whole adjustment silently never fires.
  # It has to be CALLED. Older versions expose a plain `x.labels` vector.
  n <- try({
    b  <- ggplot2::ggplot_build(plot)
    px <- b$layout$panel_params[[1]]
    labs <- NULL
    if (!is.null(px$x) && is.function(px$x$get_labels)) labs <- px$x$get_labels()
    if (is.null(labs)) labs <- px$x.labels
    length(labs[!is.na(labs)])
  }, silent = TRUE)
  if (inherits(n, "try-error") || !length(n) || is.na(n)) n <- 0L

  # Text layers at a panel boundary get sliced in half by the default clipping --
  # visible in the cluster composition plot, where the 0% and 100% labels are cut
  # by the axis line. Turning clipping off fixes it, but ONLY when the plot has
  # not set its own coord: overwriting coord_flip/coord_polar would silently
  # change the chart.
  has_text <- tryCatch(
    any(vapply(plot$layers, function(l) inherits(l$geom, c("GeomText", "GeomLabel")),
               logical(1))), error = function(e) FALSE)
  coord_is_default <- tryCatch(
    identical(class(plot$coordinates)[1], "CoordCartesian") &&
      isTRUE(plot$coordinates$clip == "on"),
    error = function(e) FALSE)
  if (has_text && coord_is_default)
    plot <- plot + ggplot2::coord_cartesian(clip = "off")

  # Below ~12 categories the default is fine; past that shrink toward a floor so
  # dense axes stay legible instead of overlapping.
  if (n > 12) {
    shrink <- max(0.55, 1 - (n - 12) * 0.012)
    plot <- plot + ggplot2::theme(
      axis.text.x  = ggplot2::element_text(size = ggplot2::rel(shrink),
                                           angle = 45, hjust = 1),
      axis.text.y  = ggplot2::element_text(size = ggplot2::rel(max(0.75, shrink))),
      legend.text  = ggplot2::element_text(size = ggplot2::rel(max(0.75, shrink))),
      legend.title = ggplot2::element_text(size = ggplot2::rel(max(0.8, shrink)))
    )
  }
  plot
}

# Shim: DimPlot label geometry.
#
# Seurat draws cluster/sample labels at a fixed point size that does not scale
# with the figure, and does not repel them unless asked. Across the vendored
# notebooks that produced UMAPs with labels rendered enormous and stacked on top
# of each other -- one had label.size = 12 against a Seurat default of 4, with no
# repel, so sample names overlapped into an unreadable smear.
#
# Fixing this per notebook means editing eight-plus vendored files; fixing it
# here reaches all of them. Only label geometry is touched -- the plot, its data
# and its colours are untouched, and an explicit small label.size is respected.
DimPlot <- function(object, ...) {
  dots <- list(...)
  if (isTRUE(dots$label)) {
    if (is.null(dots$repel)) dots$repel <- TRUE
    dots$label.size <- if (is.null(dots$label.size)) 3.6
                       else min(as.numeric(dots$label.size), 4.5)
  }
  do.call(Seurat::DimPlot, c(list(object), dots))
}

# -----------------------------------------------------------------------------
# Interactive plots INSIDE the rendered report
#
# The ggsave shim writes a .html twin beside each saved figure, but the report
# itself still embeds the static PNG -- so a reader opening the report sees
# nothing interactive. This hook makes knitr render a ggplot AS a plotly widget
# in the document.
#
# COVERAGE LIMIT, stated plainly: this only fires for plots that are the VALUE of
# a chunk (`p` or a bare `ggplot(...)` call). An explicit `print(p)` still
# rasterises, because ggplot2 4.0 objects are S7 (class includes S7_object) and
# print() dispatches through S7 rather than S3 -- print.gg / print.ggplot
# overrides do not intercept it. Converting those call sites is a per-notebook
# job: a bare `p` inside a for loop prints NOTHING, so it cannot be done blind.
#
# Disable with options(scratch.interactive = FALSE).
# -----------------------------------------------------------------------------
# Auto-printed plots never touch ggsave(), so the knit_print hook below renders
# them into the HTML and leaves figures/ EMPTY -- which is how a successful run
# published zero figures. Export here so every figure reaches figures/ whether
# it was ggsave()d explicitly or just printed by a chunk.
.scratch_export_figure <- function(x) {
  if (!isTRUE(getOption("scratch.export", TRUE))) return(invisible(NULL))
  lbl <- tryCatch(knitr::opts_current$get("label"), error = function(e) NULL)
  if (is.null(lbl) || !nzchar(lbl)) return(invisible(NULL))
  # Unnamed chunks get knitr's generated labels; those make useless filenames.
  if (grepl("^unnamed-chunk", lbl)) return(invisible(NULL))
  dir.create("figures", showWarnings = FALSE, recursive = TRUE)
  w <- tryCatch(knitr::opts_current$get("fig.width"), error = function(e) NULL)
  h <- tryCatch(knitr::opts_current$get("fig.height"), error = function(e) NULL)
  if (is.null(w) || !is.numeric(w)) w <- 8
  if (is.null(h) || !is.numeric(h)) h <- 5
  base <- file.path("figures", gsub("[^A-Za-z0-9_.-]+", "_", lbl))
  # Never let a failed export abort the render: a missing figure file is a far
  # smaller problem than a notebook that does not finish.
  for (ext in c("png", "pdf")) {
    try(suppressMessages(ggplot2::ggsave(paste0(base, ".", ext), plot = x,
                                         width = w, height = h, dpi = 300,
                                         limitsize = FALSE)), silent = TRUE)
  }
  invisible(NULL)
}

.scratch_as_widget <- function(x) {
  if (!isTRUE(getOption("scratch.interactive", TRUE))) return(NULL)
  if (!requireNamespace("plotly", quietly = TRUE)) return(NULL)

  # SIZE GUARD -- this is not an optimisation, it is a correctness fix.
  # ggplotly serialises EVERY data point into the document's JSON. A 24,691-cell
  # UMAP x ~9 plots blew past the maximum string length and quarto died with
  # "failed to allocate string; buffer exceeds maximum length", taking the whole
  # Azimuth stage with it. Even when it fits, a 24k-point widget is unusable.
  # Big plots stay static rasters; small ones become interactive.
  #
  # The count that matters is RENDERED MARKS, not input rows. A geom_bin2d over
  # 18,000 cells draws ~1,500 tiles and serialises to a small widget, but an
  # input-row count rejected it as if it were 18,000 points -- which is why the
  # QC density plots came out static while the 3-row waterfall beside them was
  # interactive. Statistical layers that aggregate get measured on their output.
  .binning <- c("StatBin2d", "StatBinhex", "StatBin", "StatCount",
                "StatBoxplot", "StatSummary", "StatDensity", "StatYdensity")
  n <- tryCatch({
    all_binned <- length(x$layers) > 0 && all(vapply(x$layers, function(l)
      class(l$stat)[1] %in% .binning, logical(1)))
    if (all_binned) {
      # Build it and count what actually gets drawn.
      b <- ggplot2::ggplot_build(x)
      max(vapply(b$data, function(d) if (is.data.frame(d)) nrow(d) else 0L, integer(1)))
    } else {
      d <- x$data
      rows <- if (is.data.frame(d)) nrow(d) else 0L
      for (l in x$layers) {
        ld <- tryCatch(l$data, error = function(e) NULL)
        if (is.data.frame(ld)) rows <- max(rows, nrow(ld))
      }
      rows
    }
  }, error = function(e) .Machine$integer.max)
  if (!is.numeric(n) || is.na(n)) n <- .Machine$integer.max

  # A flat 5,000-mark cap made every full-cohort UMAP static, which is the single
  # thing most complained about. Measured: one 24,439-point UMAP serialises to
  # 2.7 MB, or 2.2 MB with coordinates rounded -- perfectly viable. The original
  # crash ("failed to allocate string; buffer exceeds maximum length") came from
  # ~9 of them in ONE document, i.e. ~24 MB. So the limit that matters is a
  # per-DOCUMENT budget, not a per-plot cap.
  #
  # Pure point layers get the larger allowance; anything else keeps the old cap,
  # because a 24k-row heatmap or boxplot is unreadable interactively anyway.
  # "Scatter-like" rather than "point-only". SCP::CellDimPlot -- which draws
  # nearly every UMAP in the annotation reports -- adds cluster labels and leader
  # lines on top of the points. A strict point-only test excluded all of them, so
  # only the small ones (under the 5k mark cap) ever became interactive and every
  # full-cohort UMAP stayed a static raster. The annotation layers carry a handful
  # of rows; the points are what costs, so they should not veto interactivity.
  # GeomCustomAnn is the one that actually mattered and the one I guessed wrong.
  # A real SCP::CellDimPlot is exactly `GeomCustomAnn, GeomPoint, GeomTextRepel`
  # -- the CustomAnn being the little corner axis arrows that `theme_use =
  # "theme_blank"` draws. Omitting it disqualified EVERY annotation-stage UMAP,
  # while my synthetic test (which used geom_segment for the same decoration)
  # passed happily. Verified against the real 24,439-cell object, not a mock.
  .anno_geoms <- c("GeomText", "GeomLabel", "GeomTextRepel", "GeomLabelRepel",
                   "GeomSegment", "GeomCurve", "GeomPath", "GeomLine",
                   "GeomBlank", "GeomCustomAnn", "GeomRug")
  .geoms <- tryCatch(vapply(x$layers, function(l) class(l$geom)[1], character(1)),
                     error = function(e) character(0))
  .point_only <- length(.geoms) > 0 &&
                 any(.geoms %in% c("GeomPoint", "GeomJitter")) &&
                 all(.geoms %in% c("GeomPoint", "GeomJitter", .anno_geoms))

  # Charge the budget for the POINT layers only -- a 20-row label layer is free.
  if (.point_only) {
    n <- tryCatch({
      pl <- x$layers[.geoms %in% c("GeomPoint", "GeomJitter")]
      rows <- vapply(pl, function(l) {
        ld <- tryCatch(l$data, error = function(e) NULL)
        if (is.data.frame(ld)) nrow(ld)
        else if (is.data.frame(x$data)) nrow(x$data) else 0L
      }, numeric(1))
      max(c(0, rows))
    }, error = function(e) n)
  }

  cap <- if (.point_only) as.numeric(getOption("scratch.interactive.scattermax", 60000))
         else              as.numeric(getOption("scratch.interactive.maxpoints", 5000))
  if (n > cap) return(NULL)

  if (.point_only && n > 5000) {
    spent  <- as.numeric(getOption("scratch.interactive.spent", 0))
    # ~100k points is about 4 full-cohort UMAPs, ~9 MB of widget JSON -- roomy
    # for any one notebook and far below the ~24 MB that broke the string buffer.
    budget <- as.numeric(getOption("scratch.interactive.budget", 100000))
    # Falling back to a static raster is always safe; blowing the string buffer
    # kills the entire stage. When the budget is gone, stay static.
    if (spent + n > budget) return(NULL)
    options(scratch.interactive.spent = spent + n)

    # Rounding the positional aesthetics is free visually (a UMAP coordinate does
    # not mean anything past 2dp) and strips ~20% off the serialised widget.
    for (nm in c("x", "y")) {
      col <- tryCatch(rlang::as_name(x$mapping[[nm]]), error = function(e) NULL)
      if (!is.null(col) && is.data.frame(x$data) && is.numeric(x$data[[col]]))
        x$data[[col]] <- round(x$data[[col]], 2)
    }

    # The NON-positional numeric aesthetics matter more, and were being missed.
    # ggplotly builds the hover string from the data, so a colour column carrying
    # `pseudotime: 14.306830402` costs ten wasted bytes per point per trace. On
    # the trajectory report the hover text came to 8.4 MB of a 15.9 MB widget --
    # more than the x and y arrays put together, in a 39 MB document.
    #
    # signif() rather than round(): a p-value or an adjusted p-value collapses to
    # zero under round(v, 3), while signif(v, 6) leaves 1e-300 intact and still
    # turns 14.306830402 into 14.3068. Six significant digits is far beyond what
    # a hover label or a rendered pixel can express.
    for (nm in setdiff(names(x$mapping), c("x", "y"))) {
      col <- tryCatch(rlang::as_name(x$mapping[[nm]]), error = function(e) NULL)
      if (!is.null(col) && is.data.frame(x$data) &&
          !is.null(x$data[[col]]) && is.numeric(x$data[[col]]))
        x$data[[col]] <- signif(x$data[[col]], 6)
    }
  }

  # Honour an explicit `text` aesthetic when the plot supplies one; otherwise
  # ggplotly builds the hover string out of the raw aes expressions, e.g.
  # "reorder(sample, cells): HRS371755  cells: 10365  cells: 10365".
  .has_text <- tryCatch("text" %in% names(x$mapping) ||
                        any(vapply(x$layers, function(l)
                              "text" %in% names(if (is.null(l$mapping)) list() else l$mapping),
                            logical(1))),
                        error = function(e) FALSE)
  # GeomCustomAnn wraps a raw grob (SCP's corner axis arrows). plotly has no
  # equivalent and emits a malformed trace of type "gl" that carries none of the
  # attributes it is given -- verified in the built widget. It is pure decoration
  # and the axis titles say the same thing, so drop it rather than ship a broken
  # trace into the document.
  if (isTRUE(.point_only) && any(.geoms == "GeomCustomAnn")) {
    x$layers <- x$layers[.geoms != "GeomCustomAnn"]
  }

  # suppressWarnings: ggplotly reports every ggplot attribute a plotly trace type
  # does not accept, one line per trace. Eight cell types meant eight identical
  # notices printed into the report above the figure.
  w <- try(suppressWarnings(
         if (.has_text) plotly::ggplotly(x, tooltip = "text") else plotly::ggplotly(x)),
         silent = TRUE)
  if (inherits(w, "try-error")) return(NULL)

  # SVG chokes above ~10k points; WebGL renders 24k smoothly. It does not shrink
  # the payload (measured: identical), it makes the widget usable once open.
  # WebGL ONLY when every layer is points. ggplotly turns a text/label layer into
  # a trace with no `type`, and toWebGL renames that to the invalid type "gl" --
  # confirmed in the built widget for SCP::CellDimPlot, which carries repelled
  # cluster labels. The speed-up is not worth emitting a trace plotly.js does not
  # recognise, so labelled UMAPs stay SVG (interactive, just less brisk).
  .pure_points <- length(.geoms) > 0 &&
                  all(.geoms[.geoms != "GeomCustomAnn"] %in% c("GeomPoint", "GeomJitter"))
  if (isTRUE(.point_only) && isTRUE(.pure_points) && n > 5000) {
    # suppressWarnings: WebGL traces legitimately lack a few SVG-only attributes
    # ("'scattergl' objects don't have these attributes: 'hoveron'"), and plotly
    # reports each one. Nothing is lost, but the notices land in the rendered
    # report above the figure.
    wg <- try(suppressWarnings(plotly::toWebGL(w)), silent = TRUE)
    if (!inherits(wg, "try-error")) w <- wg
  }
  # ggplotly DROPS the subtitle. Every subtitle in this pipeline carries the
  # interpretation -- which line is the applied cutoff, what the colour means --
  # so losing it silently gutted the interactive figures. Fold it into the
  # plotly title as a second, smaller line.
  .sub <- tryCatch(x$labels$subtitle, error = function(e) NULL)
  .ttl <- tryCatch(x$labels$title, error = function(e) NULL)
  if (!is.null(.sub) && nzchar(.sub)) {
    esc <- function(z) { z <- gsub("&", "&amp;", z); z <- gsub("<", "&lt;", z)
                         gsub(">", "&gt;", z) }
    w <- try(plotly::layout(w, title = list(
      text = sprintf("<b>%s</b><br><sup>%s</sup>",
                     esc(if (is.null(.ttl)) "" else .ttl), esc(.sub)),
      x = 0, xanchor = "left", font = list(size = 15)),
      margin = list(t = 70)), silent = TRUE)
    if (inherits(w, "try-error")) return(NULL)
  }
  # Legend on the RIGHT, one entry per cell type, and clickable: in plotly a
  # legend entry toggles its trace, so a single UMAP coloured by label replaces a
  # grid of one-panel-per-label figures -- click to isolate a population instead
  # of hunting for its panel. Double-click isolates it outright.
  w <- try(plotly::layout(w,
             legend = list(orientation = "v", x = 1.02, xanchor = "left",
                           y = 1, yanchor = "top", itemsizing = "constant",
                           itemclick = "toggle", itemdoubleclick = "toggleothers"),
             margin = list(r = 10)),
           silent = TRUE)
  if (inherits(w, "try-error")) return(NULL)

  try(plotly::config(w, displaylogo = FALSE, responsive = TRUE), silent = TRUE)
}

if (requireNamespace("knitr", quietly = TRUE)) {
  knit_print.ggplot <- function(x, ...) {
    .scratch_export_figure(x)
    w <- .scratch_as_widget(x)
    if (is.null(w) || inherits(w, "try-error")) return(knitr::normal_print(x))
    knitr::knit_print(w, ...)
  }
  # Register for every class ggplot2 4.0 (S7) puts on a plot object, so the
  # method is found wherever knitr looks.
  for (.cls in c("ggplot", "gg", "ggplot2::ggplot", "ggplot2::gg")) {
    try(registerS3method("knit_print", .cls, knit_print.ggplot,
                         envir = asNamespace("knitr")), silent = TRUE)
  }
  rm(.cls)
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


# ---------------------------------------------------------------------------
# scratch_feature_selector() -- one UMAP, a dropdown, one score at a time.
#
# Replaces the FeaturePlot small-multiples grid. That grid drew N near-identical
# UMAPs at ~2 inches each: the cluster labels collided into an unreadable smear,
# every panel spent most of its ink redrawing the same grey background, and the
# colour bar was too small to read. The information is one embedding and N
# scores, so the honest form is one embedding and a way to choose the score.
#
# Built with plotly directly rather than ggplotly: this needs `updatemenus`,
# which has no ggplot equivalent, and one scattergl trace per feature keeps the
# widget usable at cohort scale.
# ---------------------------------------------------------------------------
scratch_feature_selector <- function(object, features, reduction = "umap",
                                     max_cells = getOption("scratch.selector.maxcells", 30000),
                                     title = "Signature score") {
  if (!requireNamespace("plotly", quietly = TRUE)) return(NULL)
  emb <- try(as.data.frame(SeuratObject::Embeddings(object, reduction))[, 1:2], silent = TRUE)
  if (inherits(emb, "try-error") || !nrow(emb)) return(NULL)
  names(emb) <- c("d1", "d2")
  md    <- object@meta.data
  feats <- intersect(features, names(md))
  feats <- feats[vapply(feats, function(f) is.numeric(md[[f]]), logical(1))]
  if (!length(feats)) return(NULL)

  idx <- seq_len(nrow(emb))
  if (length(idx) > max_cells) idx <- sort(sample(idx, max_cells))  # keeps the widget openable
  emb <- emb[idx, , drop = FALSE]
  md  <- md[idx, , drop = FALSE]

  # ONE trace, not one per feature. The x/y coordinates are identical for every
  # signature, so a trace-per-feature serialised the same embedding N times --
  # 4 MB for five signatures, and this notebook can have sixteen. The dropdown
  # restyles the colour and hover arrays of a single trace instead, so the
  # coordinates are paid for once.
  cols <- lapply(feats, function(f) as.numeric(md[[f]]))
  txts <- lapply(seq_along(feats), function(i)
    sprintf("%s\nscore: %.3f", feats[i], cols[[i]]))

  p <- plotly::plot_ly(
    x = round(emb$d1, 2), y = round(emb$d2, 2),
    type = "scattergl", mode = "markers",
    marker = list(size = 3, color = cols[[1]], colorscale = "Viridis",
                  showscale = TRUE,
                  colorbar = list(title = list(text = "score"), len = .75)),
    text = txts[[1]], hoverinfo = "text")

  ttl <- function(i) list(text = sprintf("<b>%s</b><br><sup>%s</sup>", title, feats[i]),
                          x = 0, xanchor = "left")
  btns <- lapply(seq_along(feats), function(i)
    list(method = "update", label = feats[i],
         args = list(list(`marker.color` = list(cols[[i]]), text = list(txts[[i]])),
                     list(title = ttl(i)))))

  plotly::config(
    plotly::layout(p,
      title = ttl(1),
      xaxis = list(title = paste0(reduction, "_1"), zeroline = FALSE),
      yaxis = list(title = paste0(reduction, "_2"), zeroline = FALSE),
      showlegend = FALSE,
      margin = list(t = 90, r = 10),
      updatemenus = list(list(
        type = "dropdown", direction = "down", showactive = TRUE,
        x = 1.02, xanchor = "left", y = 1, yanchor = "top",
        pad = list(r = 4, t = 4), buttons = btns))),
    displaylogo = FALSE, responsive = TRUE)
}

# =============================================================================
# Publication figure grammar
# =============================================================================
#
# Patterns taken from what journal figures actually do, and which the reports
# were missing. Three things separate a professional panel from a default
# ggplot, and none of them is styling:
#
#   1. SHOW THE DATA. A violin without its points hides n, and hides whether a
#      distribution is bimodal or zero-inflated. Reference figures overlay the
#      observations.
#   2. LABEL DIRECTLY. A volcano with a legend makes the reader look things up;
#      a volcano with gene names on the points does not.
#   3. PUT THE STATISTIC IN THE PANEL. "r = 0.46, P = 0.03" inside the plotting
#      area, not in prose somewhere below it.

#' Violin + jittered points + median crossbar.
#'
#' The default violin collapses to a hairline when the data are zero-inflated,
#' which is exactly what happened to the CNV burden panel: 11 groups rendered as
#' 11 dashes because most cells sit at zero. Points make that visible rather than
#' invisible, and the count annotation states n instead of implying it.
scratch_violin_jitter <- function(df, x, y, fill = NULL, max_points = 3000,
                                  log_y = FALSE, show_n = TRUE, ...) {
  stopifnot(is.data.frame(df))
  df <- as.data.frame(df)
  fill <- fill %||% x
  df <- df[is.finite(df[[y]]), , drop = FALSE]
  if (!nrow(df)) return(NULL)

  # Subsample only the points; the violin still uses every observation, so the
  # shape never changes with the plotting budget.
  pts <- if (nrow(df) > max_points) df[sample.int(nrow(df), max_points), , drop = FALSE] else df

  p <- ggplot2::ggplot(df, ggplot2::aes(x = .data[[x]], y = .data[[y]])) +
    ggplot2::geom_violin(ggplot2::aes(fill = .data[[fill]]),
                         scale = "width", colour = NA, alpha = .55, width = .9) +
    ggplot2::geom_jitter(data = pts, width = .18, height = 0, size = .35,
                         alpha = .35, colour = "grey20") +
    ggplot2::stat_summary(fun = stats::median, geom = "crossbar",
                          width = .55, linewidth = .45, colour = "grey12") +
    scale_fill_scratch(guide = "none") +
    theme_scratch()

  if (log_y) p <- p + ggplot2::scale_y_log10(labels = scales::label_log())

  if (show_n) {
    n_df <- as.data.frame(table(df[[x]]))
    names(n_df) <- c(x, "n")
    n_df <- n_df[n_df$n > 0, , drop = FALSE]
    ymax <- max(df[[y]], na.rm = TRUE)
    p <- p + ggplot2::geom_text(
      data = n_df, inherit.aes = FALSE,
      ggplot2::aes(x = .data[[x]], y = ymax, label = paste0("n=", format(n, big.mark = ","))),
      vjust = -0.4, size = 2.7, colour = "grey35") +
      ggplot2::coord_cartesian(clip = "off") +
      ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(.03, .14)))
  }
  p
}

#' Volcano plot with direction colouring and direct gene labels.
#'
#' Five notebooks in this pipeline compute differential expression and not one
#' of them draws a volcano. This is the canonical view: effect on x, evidence on
#' y, colour by direction, threshold lines drawn, and the extreme genes named on
#' the points so the reader is not decoding a legend.
#' @param highlight optional character vector of genes to call out as their own
#'   category. Both reference figures do this -- the Cancer Discovery volcano
#'   separates interferon-response genes from other upregulated genes, which
#'   turns a generic volcano into a statement about a specific programme.
#' @param direction_labels two-element vector annotated as arrows at the top,
#'   e.g. c("low ITMB", "high ITMB"). Reference figures label the CONTRAST, not
#'   just the axis, so the reader never has to work out which side is which.
scratch_volcano <- function(df, lfc = "avg_log2FC", p = "p_val_adj", label = NULL,
                            lfc_cut = 0.25, p_cut = 0.05, n_label = 20,
                            up_name = "up", down_name = "down",
                            highlight = NULL, highlight_name = "highlighted",
                            direction_labels = NULL) {
  df <- as.data.frame(df)
  df <- df[is.finite(df[[lfc]]) & !is.na(df[[p]]), , drop = FALSE]
  if (!nrow(df)) return(NULL)

  # A p-value of exactly 0 is -log10 = Inf and drops the point. Floor it at the
  # smallest representable double instead of silently losing the best hits.
  # A volcano of one-sided data is a misleading picture: an empty half implies
  # nothing moved that way, when in fact it was filtered out upstream. Seurat's
  # FindAllMarkers(only.pos = TRUE) and any `lfc > cut` filter both produce this.
  # Refuse to draw it silently -- return the diagnosis so the caller can report
  # the truth instead of rendering half a plot.
  n_neg <- sum(df[[lfc]] < 0, na.rm = TRUE)
  n_pos <- sum(df[[lfc]] > 0, na.rm = TRUE)
  if (n_neg == 0 || n_pos == 0) {
    # Render the DIAGNOSIS, not a half-empty volcano and not nothing. Returning
    # NULL would leave a silent gap in the report; returning a misleading plot is
    # worse. A panel that explains itself is the only option that cannot be
    # misread.
    msg <- sprintf(paste0("Volcano not drawn: this table is one-sided\n",
                          "(%s up, %s down).\n\n",
                          "The source was filtered to a single direction, so the\n",
                          "empty half would imply nothing moved that way.\n",
                          "Use scratch_marker_panel() for ranked one-sided data."),
                   format(n_pos, big.mark = ","), format(n_neg, big.mark = ","))
    return(
      ggplot2::ggplot() +
        ggplot2::annotate("text", x = 0, y = 0, label = msg, size = 3.4,
                          colour = "grey25", lineheight = 1.25) +
        ggplot2::xlim(-1, 1) + ggplot2::ylim(-1, 1) +
        theme_scratch() +
        ggplot2::theme(axis.text = ggplot2::element_blank(),
                       axis.title = ggplot2::element_blank(),
                       panel.grid = ggplot2::element_blank())
    )
  }

  df$.y <- -log10(pmax(df[[p]], .Machine$double.xmin))
  df$.dir <- ifelse(df[[p]] < p_cut & df[[lfc]] >=  lfc_cut, up_name,
             ifelse(df[[p]] < p_cut & df[[lfc]] <= -lfc_cut, down_name, "ns"))

  # A highlighted gene keeps its direction but gets its own colour, so the panel
  # answers "did THIS programme move" as well as "what moved".
  hl <- character(0)
  if (!is.null(highlight) && !is.null(label) && label %in% names(df)) {
    hl <- intersect(toupper(highlight), toupper(df[[label]]))
    is_hl <- toupper(df[[label]]) %in% hl & df$.dir != "ns"
    df$.dir[is_hl] <- highlight_name
  }
  lev <- c(down_name, "ns", up_name)
  if (length(hl)) lev <- c(lev, highlight_name)
  df$.dir <- factor(df$.dir, levels = lev)

  cols <- stats::setNames(c("#2c6fb0", "grey78", "#e79a9a"), c(down_name, "ns", up_name))
  if (length(hl)) cols[highlight_name] <- "#b32020"
  n_up <- sum(df$.dir %in% c(up_name, highlight_name)); n_dn <- sum(df$.dir == down_name)

  p_obj <- ggplot2::ggplot(df, ggplot2::aes(x = .data[[lfc]], y = .data$.y, colour = .data$.dir)) +
    ggplot2::geom_point(size = .8, alpha = .75) +
    ggplot2::geom_vline(xintercept = c(-lfc_cut, lfc_cut), linetype = "dashed",
                        colour = "grey45", linewidth = .3) +
    ggplot2::geom_hline(yintercept = -log10(p_cut), linetype = "dashed",
                        colour = "grey45", linewidth = .3) +
    ggplot2::scale_colour_manual(values = cols, name = NULL) +
    ggplot2::labs(x = expression(log[2]~fold~change), y = expression(-log[10]~adjusted~P)) +
    theme_scratch()

  # Direct labels on the most extreme genes in each direction.
  # `.have()` lives in scratch_report.R, which may not be loaded; use the base
  # check directly so this file stands alone.
  if (!is.null(label) && label %in% names(df) &&
      requireNamespace("ggrepel", quietly = TRUE)) {
    top <- do.call(rbind, lapply(c(up_name, down_name), function(dd) {
      s <- df[df$.dir == dd, , drop = FALSE]
      if (!nrow(s)) return(NULL)
      s <- s[order(-s$.y * abs(s[[lfc]])), , drop = FALSE]
      utils::head(s, ceiling(n_label / 2))
    }))
    # Prefer the highlighted genes for labelling when a set was supplied.
    if (length(hl)) {
      hset <- df[df$.dir == highlight_name, , drop = FALSE]
      if (nrow(hset)) {
        hset <- hset[order(-hset$.y * abs(hset[[lfc]])), , drop = FALSE]
        top <- rbind(utils::head(hset, n_label), top)
        top <- top[!duplicated(top[[label]]), , drop = FALSE]
        top <- utils::head(top, n_label)
      }
    }
    if (!is.null(top) && nrow(top)) {
      # Boxed labels, as in both reference figures: a white fill keeps the name
      # readable where points are dense, which is exactly where the interesting
      # genes are.
      p_obj <- p_obj + ggrepel::geom_label_repel(
        data = top, ggplot2::aes(label = .data[[label]], colour = .data$.dir),
        size = 2.5, max.overlaps = 30, min.segment.length = 0,
        label.padding = grid::unit(0.12, "lines"), label.size = 0.25,
        fill = "white", fontface = "italic",
        segment.size = .2, segment.colour = "grey55", show.legend = FALSE)
    }
  }

  # Name the contrast with arrows, so which side means what needs no caption.
  if (!is.null(direction_labels) && length(direction_labels) == 2) {
    rng <- range(df[[lfc]], na.rm = TRUE); ytop <- max(df$.y, na.rm = TRUE)
    p_obj <- p_obj +
      ggplot2::annotate("segment", x = rng[1] * .85, xend = rng[1],
                        y = ytop * 1.12, yend = ytop * 1.12,
                        arrow = grid::arrow(length = grid::unit(.16, "cm")),
                        colour = "#2c6fb0", linewidth = .5) +
      ggplot2::annotate("text", x = rng[1] * .8, y = ytop * 1.12,
                        label = direction_labels[1], hjust = 0, size = 3,
                        colour = "#2c6fb0") +
      ggplot2::annotate("segment", x = rng[2] * .85, xend = rng[2],
                        y = ytop * 1.12, yend = ytop * 1.12,
                        arrow = grid::arrow(length = grid::unit(.16, "cm")),
                        colour = "#b32020", linewidth = .5) +
      ggplot2::annotate("text", x = rng[2] * .8, y = ytop * 1.12,
                        label = direction_labels[2], hjust = 1, size = 3,
                        colour = "#b32020")
  }

  # Counts stated in the panel, so "how many changed" needs no second look.
  p_obj +
    ggplot2::annotate("text", x = Inf, y = Inf, hjust = 1.05, vjust = 1.4,
                      label = sprintf("%s up", format(n_up, big.mark = ",")),
                      colour = cols[[up_name]], size = 3.1, fontface = "bold") +
    ggplot2::annotate("text", x = -Inf, y = Inf, hjust = -0.05, vjust = 1.4,
                      label = sprintf("%s down", format(n_dn, big.mark = ",")),
                      colour = cols[[down_name]], size = 3.1, fontface = "bold") +
    ggplot2::coord_cartesian(clip = "off")
}

#' Enrichment dot plot: size = gene count, colour = adjusted P.
#'
#' The standard GO/pathway display. Two visual channels carry two quantities, so
#' a reader sees both how big a term is and how strong its evidence is without
#' cross-referencing a table.
scratch_enrich_dot <- function(df, term = "term", ratio = "gene_ratio",
                               count = "count", padj = "p_adjust", n_top = 15) {
  df <- as.data.frame(df)
  df <- df[is.finite(df[[ratio]]) & !is.na(df[[padj]]), , drop = FALSE]
  if (!nrow(df)) return(NULL)
  df <- df[order(df[[padj]]), , drop = FALSE]
  df <- utils::head(df, n_top)
  # Ordered by the quantity on the x axis, so the eye reads top-to-bottom.
  df[[term]] <- factor(df[[term]], levels = df[[term]][order(df[[ratio]])])

  ggplot2::ggplot(df, ggplot2::aes(x = .data[[ratio]], y = .data[[term]])) +
    ggplot2::geom_point(ggplot2::aes(size = .data[[count]], colour = .data[[padj]])) +
    ggplot2::scale_colour_gradient(low = "#b32020", high = "#4a7fb5",
                                   name = "adj P", trans = "log10") +
    ggplot2::scale_size_continuous(name = "count", range = c(2, 7)) +
    ggplot2::labs(x = "gene ratio", y = NULL) +
    theme_scratch() +
    ggplot2::theme(panel.grid.major.y = ggplot2::element_line(colour = "grey92", linewidth = .3))
}

#' Put a statistic inside the plotting area, journal style.
scratch_stat_label <- function(p, text, corner = c("topright", "topleft",
                                                   "bottomright", "bottomleft")) {
  corner <- match.arg(corner)
  pos <- switch(corner,
    topright    = list(x =  Inf, y =  Inf, h =  1.1, v =  1.5),
    topleft     = list(x = -Inf, y =  Inf, h = -0.1, v =  1.5),
    bottomright = list(x =  Inf, y = -Inf, h =  1.1, v = -0.7),
    bottomleft  = list(x = -Inf, y = -Inf, h = -0.1, v = -0.7))
  p + ggplot2::annotate("text", x = pos$x, y = pos$y, hjust = pos$h, vjust = pos$v,
                        label = text, size = 3.1, colour = "grey20") +
    ggplot2::coord_cartesian(clip = "off")
}

#' Heatmap with categorical annotation bars, as in a journal figure.
#'
#' The reference figures never show a bare heatmap: there is always a coloured
#' sidebar naming what each row or column IS (method, patient, GO group), so the
#' block structure can be read without cross-referencing an axis of 50 rotated
#' labels. Our heatmaps had the matrix and nothing else.
#'
#' Built with patchwork rather than ComplexHeatmap deliberately: it stays a
#' ggplot, so it keeps the house theme, the interactive shim and the PDF export
#' that every other figure here gets.
#'
#' @param df long-form: one row per cell of the matrix.
#' @param col_anno,row_anno optional data frames keyed by the x / y column,
#'   each remaining column drawn as one annotation stripe.
scratch_heatmap_annotated <- function(df, x, y, fill,
                                      col_anno = NULL, row_anno = NULL,
                                      fill_name = fill, diverging = FALSE,
                                      title = NULL, subtitle = NULL) {
  df <- as.data.frame(df)
  if (!nrow(df)) return(NULL)

  hm <- ggplot2::ggplot(df, ggplot2::aes(x = .data[[x]], y = .data[[y]], fill = .data[[fill]])) +
    ggplot2::geom_tile(colour = "white", linewidth = .3) +
    theme_scratch() +
    ggplot2::theme(panel.grid = ggplot2::element_blank(),
                   axis.text.x = ggplot2::element_text(angle = 45, hjust = 1)) +
    ggplot2::labs(x = NULL, y = NULL, title = title, subtitle = subtitle)

  hm <- hm + if (diverging) {
    ggplot2::scale_fill_gradient2(low = "#2c6fb0", mid = "#f7f7f7", high = "#b32020",
                                  midpoint = 0, name = fill_name)
  } else {
    ggplot2::scale_fill_gradient(low = "#f2f6fa", high = "#1f4e79", name = fill_name)
  }

  if (!requireNamespace("patchwork", quietly = TRUE) ||
      (is.null(col_anno) && is.null(row_anno))) return(hm)

  # One thin stripe per annotation column, sharing the heatmap's axis so the
  # categories line up with the matrix.
  strip <- function(anno, key, horizontal) {
    anno <- as.data.frame(anno)
    vars <- setdiff(names(anno), key)
    if (!length(vars)) return(NULL)
    long <- do.call(rbind, lapply(vars, function(v)
      data.frame(k = anno[[key]], var = v, val = as.character(anno[[v]]),
                 stringsAsFactors = FALSE)))
    g <- ggplot2::ggplot(long, ggplot2::aes(fill = val)) +
      scale_fill_scratch(name = NULL) +
      theme_scratch(base_size = 9) +
      ggplot2::theme(panel.grid = ggplot2::element_blank(),
                     axis.title = ggplot2::element_blank(),
                     axis.ticks = ggplot2::element_blank(),
                     legend.position = "bottom", legend.box = "vertical")
    if (horizontal) {
      g + ggplot2::geom_tile(ggplot2::aes(x = k, y = var), colour = "white", linewidth = .3) +
        ggplot2::theme(axis.text.x = ggplot2::element_blank())
    } else {
      g + ggplot2::geom_tile(ggplot2::aes(x = var, y = k), colour = "white", linewidth = .3) +
        ggplot2::theme(axis.text.y = ggplot2::element_blank())
    }
  }

  top  <- if (!is.null(col_anno)) strip(col_anno, x, TRUE)  else NULL
  side <- if (!is.null(row_anno)) strip(row_anno, y, FALSE) else NULL

  out <- hm
  if (!is.null(side)) out <- side + out + patchwork::plot_layout(widths = c(1, 14))
  if (!is.null(top))  out <- top / out + patchwork::plot_layout(heights = c(1, 10))
  out
}

#' Bubble plot: colour is one quantity, size is another.
#'
#' The single most common panel in the reference figures -- both use it, for
#' ligand-receptor pairs (colour = mean expression, size = significance) and for
#' pathway enrichment. Generalised here because a dot plot restricted to
#' enrichment could not serve the cell-communication reports.
scratch_bubble <- function(df, x, y, size, colour,
                           size_name = size, colour_name = colour,
                           diverging = FALSE, size_range = c(1.5, 8)) {
  df <- as.data.frame(df)
  df <- df[is.finite(df[[size]]) & is.finite(df[[colour]]), , drop = FALSE]
  if (!nrow(df)) return(NULL)

  p <- ggplot2::ggplot(df, ggplot2::aes(x = .data[[x]], y = .data[[y]])) +
    ggplot2::geom_point(ggplot2::aes(size = .data[[size]], colour = .data[[colour]])) +
    ggplot2::scale_size_continuous(name = size_name, range = size_range) +
    ggplot2::labs(x = NULL, y = NULL) +
    theme_scratch() +
    ggplot2::theme(panel.grid.major = ggplot2::element_line(colour = "grey93", linewidth = .3),
                   axis.text.x = ggplot2::element_text(angle = 40, hjust = 1))

  p + if (diverging) {
    ggplot2::scale_colour_gradient2(low = "#2c6fb0", mid = "#f7f7f7", high = "#b32020",
                                    midpoint = 0, name = colour_name)
  } else {
    ggplot2::scale_colour_gradient(low = "#4a7fb5", high = "#b32020", name = colour_name)
  }
}


#' Ranked marker panel: the honest figure for one-sided data.
#'
#' Seurat's FindAllMarkers(only.pos = TRUE) returns markers, not a two-sided
#' contrast -- there is no "down" to show, so a volcano is the wrong instrument.
#' This is the standard alternative: top genes per group, effect size on x,
#' evidence as colour, detection rate as size.
scratch_marker_panel <- function(df, group = "cluster", gene = "gene",
                                 lfc = "avg_log2FC", p = "p_val_adj",
                                 pct = NULL, n_top = 8, n_groups = 12) {
  df <- as.data.frame(df)
  df <- df[is.finite(df[[lfc]]) & !is.na(df[[p]]), , drop = FALSE]
  if (!nrow(df)) return(NULL)

  keep <- names(sort(table(df[[group]]), decreasing = TRUE))[seq_len(min(n_groups,
              length(unique(df[[group]]))))]
  df <- df[df[[group]] %in% keep, , drop = FALSE]

  top <- do.call(rbind, lapply(split(df, df[[group]]), function(s) {
    s <- s[order(-s[[lfc]]), , drop = FALSE]
    utils::head(s, n_top)
  }))
  if (is.null(top) || !nrow(top)) return(NULL)

  # Unique y label per group so the same gene can head two clusters without the
  # rows collapsing onto each other.
  top$.lab <- paste0(top[[gene]], "  (", top[[group]], ")")
  top <- top[order(top[[group]], top[[lfc]]), , drop = FALSE]
  top$.lab <- factor(top$.lab, levels = unique(top$.lab))
  top$.nlp <- -log10(pmax(top[[p]], .Machine$double.xmin))

  aes_args <- list(x = ggplot2::sym(lfc), y = ggplot2::sym(".lab"),
                   colour = ggplot2::sym(".nlp"))
  if (!is.null(pct) && pct %in% names(top)) aes_args$size <- ggplot2::sym(pct)

  p_obj <- ggplot2::ggplot(top, do.call(ggplot2::aes, aes_args)) +
    ggplot2::geom_segment(ggplot2::aes(x = 0, xend = .data[[lfc]], yend = .data$.lab),
                          colour = "grey85", linewidth = .5) +
    ggplot2::geom_point() +
    ggplot2::scale_colour_gradient(low = "#9ecae1", high = "#b32020",
                                   name = expression(-log[10]~adj~P)) +
    ggplot2::facet_wrap(~ .data[[group]], scales = "free_y", ncol = 3) +
    ggplot2::labs(x = expression(log[2]~fold~change), y = NULL) +
    theme_scratch() +
    ggplot2::theme(panel.grid.major.y = ggplot2::element_line(colour = "grey94", linewidth = .3),
                   strip.text = ggplot2::element_text(face = "bold"))
  if (!is.null(pct) && pct %in% names(top))
    p_obj <- p_obj + ggplot2::scale_size_continuous(name = "pct expressing", range = c(1.5, 5))
  p_obj
}
