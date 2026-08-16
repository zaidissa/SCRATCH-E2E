# =============================================================================
# scratch_report.R — smart-report components for SCRATCH-E2E
# =============================================================================
#
# The pipeline's reports were linear: every figure stacked under the previous
# one, so a 12-sample cohort produced 12 near-identical plots the reader had to
# scroll past. Across 40 notebooks there were zero tabsets, one interactive
# plot and five interactive tables.
#
# This provides the components to fix that, all degrading gracefully when the
# optional packages are absent from a container:
#
#   scratch_tiles()     summary-first KPI row      — the answer before the detail
#   scratch_finding()   a callout COMPUTED from the data, not a static caption
#   scratch_tabset()    one tab per sample/group   — instead of stacking
#   scratch_interact()  ggplot -> hover/zoom       — falls back to static
#   scratch_table()     searchable, sortable, CSV-exportable
#   scratch_details()   collapsed long tail
#   scratch_grid()      small multiples
#
# Design rule: the reader should get the conclusion in the first screen, and
# have to opt in to the detail. Everything below the fold is drill-down.
# =============================================================================

suppressPackageStartupMessages({ library(ggplot2) })

.have <- function(pkg) requireNamespace(pkg, quietly = TRUE)

# -----------------------------------------------------------------------------
# Summary-first: a KPI row
# -----------------------------------------------------------------------------

#' @param items named list: label -> value, or list(value=, sub=, state=)
#'   state is one of good/warning/critical and adds a glyph, never colour alone.
scratch_tiles <- function(items, title = NULL) {
  pal <- list(good = "#0ca30c", warning = "#fab219", critical = "#d03b3b")
  glyph <- list(good = "&#10003;", warning = "&#9888;", critical = "&#10007;")

  css <- '
<style>
.sx-tiles{display:grid;grid-template-columns:repeat(auto-fit,minmax(150px,1fr));gap:10px;margin:1rem 0 1.4rem}
.sx-tile{border:1px solid rgba(128,128,128,.28);border-radius:8px;padding:12px 14px;background:rgba(128,128,128,.04)}
.sx-tile .lab{font-size:.7rem;letter-spacing:.06em;text-transform:uppercase;opacity:.65}
.sx-tile .val{font-size:1.6rem;font-weight:650;line-height:1.2;margin-top:2px;font-variant-numeric:tabular-nums}
.sx-tile .sub{font-size:.76rem;opacity:.7}
.sx-tile .st{font-size:.76rem;font-weight:600}
</style>'

  out <- c(if (!is.null(title)) sprintf("\n**%s**\n", title), css, '<div class="sx-tiles">')
  for (nm in names(items)) {
    it <- items[[nm]]
    if (!is.list(it)) it <- list(value = it)
    st <- it$state %||% NULL
    stline <- if (!is.null(st) && st %in% names(pal))
      sprintf('<div class="st" style="color:%s">%s %s</div>', pal[[st]], glyph[[st]], st) else ""
    out <- c(out, sprintf(
      '<div class="sx-tile"><div class="lab">%s</div><div class="val">%s</div><div class="sub">%s</div>%s</div>',
      nm, it$value, it$sub %||% "", stline))
  }
  cat(paste(c(out, "</div>"), collapse = "\n"), "\n")
  invisible(NULL)
}

`%||%` <- function(a, b) if (is.null(a)) b else a

# -----------------------------------------------------------------------------
# A finding computed from the data
#
# The difference between a report and a slide deck: this states what the numbers
# actually show, generated at render time, so it cannot go stale relative to the
# figure beside it.
# -----------------------------------------------------------------------------

scratch_finding <- function(text, type = c("note", "tip", "important", "warning", "caution", "info"),
                            title = NULL, collapse = FALSE) {
  type <- match.arg(type)
  # "info" reads more naturally at a call site than "note"; quarto has no
  # callout-info, so alias it rather than emitting an unstyled callout.
  if (type == "info") type <- "note"
  hdr <- sprintf('::: {.callout-%s%s%s}', type,
                 if (!is.null(title)) sprintf(' title="%s"', title) else "",
                 if (collapse) ' collapse="true"' else "")
  cat("\n", hdr, "\n", paste(text, collapse = "\n"), "\n:::\n", sep = "")
  invisible(NULL)
}

# Convenience: phrase a proportion as a sentence with the number in it.
scratch_state <- function(n, total, what, threshold = NULL, unit = "") {
  pct <- if (total > 0) 100 * n / total else NA_real_
  s <- sprintf("**%d of %d** %s (%.1f%%)%s", n, total, what, pct, unit)
  if (!is.null(threshold)) {
    s <- paste0(s, if (pct >= threshold) " — above" else " — below",
                sprintf(" the %.0f%% threshold.", threshold))
  }
  s
}

# -----------------------------------------------------------------------------
# Tabsets: one tab per sample/group instead of a vertical stack
#
# Uses knit_child so the contents can be ggplots, htmlwidgets, tables or plain
# text — a plain `cat()` loop cannot emit widgets correctly.
# -----------------------------------------------------------------------------

.sx_render_one <- function(x) {
  if (inherits(x, c("gg", "ggplot", "patchwork"))) { print(x) }
  else if (inherits(x, c("htmlwidget", "datatables", "plotly"))) { print(x) }
  else if (is.data.frame(x)) { print(knitr::kable(x)) }
  else if (is.character(x)) { cat(x, sep = "\n") }
  else if (is.function(x)) { x() }
  else print(x)
}

#' @param items named list; names become tab labels
#' @param level heading level for the tabs (3 = ### under a ## section)
scratch_tabset <- function(items, level = 3, fig_width = 7.2, fig_height = 4.4,
                           pills = FALSE) {
  items <- items[!vapply(items, is.null, logical(1))]
  if (!length(items)) { cat("\n_Nothing to show._\n"); return(invisible(NULL)) }

  # One item does not need a tabset; a lone tab is UI noise.
  if (length(items) == 1L) {
    cat("\n"); .sx_render_one(items[[1]]); cat("\n")
    return(invisible(NULL))
  }

  env <- new.env(parent = globalenv())
  chunks <- character(0)
  for (i in seq_along(items)) {
    assign(paste0(".sx_item_", i), items[[i]], envir = env)
    chunks <- c(chunks,
      sprintf("%s %s", strrep("#", level), names(items)[i]), "",
      sprintf("```{r sx-tab-%s-%d, echo=FALSE, message=FALSE, warning=FALSE, fig.width=%s, fig.height=%s}",
              substr(digest_or_rand(), 1, 6), i, fig_width, fig_height),
      sprintf(".sx_render_one(.sx_item_%d)", i),
      "```", "")
  }
  assign(".sx_render_one", .sx_render_one, envir = env)

  body <- knitr::knit_child(text = chunks, envir = env, quiet = TRUE)
  cat(sprintf("\n::: {.panel-tabset%s}\n", if (pills) " .nav-pills" else ""),
      body, "\n:::\n", sep = "\n")
  invisible(NULL)
}

digest_or_rand <- function() {
  if (.have("digest")) digest::digest(runif(1)) else paste0(sample(letters, 8, TRUE), collapse = "")
}

# -----------------------------------------------------------------------------
# Interactive plots
#
# ggplot in, hover/zoom out. Falls back to the static plot when plotly is not
# installed, so a lean container still renders a complete report.
# -----------------------------------------------------------------------------

scratch_interact <- function(p, tooltip = "all", height = NULL, static_ok = TRUE) {
  if (!.have("plotly")) {
    if (static_ok) return(p)
    return(p)
  }
  out <- try({
    w <- plotly::ggplotly(p, tooltip = tooltip, height = height)
    w <- plotly::config(w, displaylogo = FALSE, responsive = TRUE,
                        modeBarButtonsToRemove = list("select2d", "lasso2d",
                                                      "hoverClosestCartesian",
                                                      "hoverCompareCartesian",
                                                      "toggleSpikelines"),
                        toImageButtonOptions = list(format = "svg", scale = 2))
    plotly::layout(w, legend = list(orientation = "v", x = 1.02, y = 1),
                   margin = list(l = 50, r = 20, t = 46, b = 50))
  }, silent = TRUE)
  if (inherits(out, "try-error")) p else out
}

# -----------------------------------------------------------------------------
# Tables: searchable, sortable, exportable
# -----------------------------------------------------------------------------

scratch_table <- function(df, caption = NULL, page = 10, searchable = TRUE,
                          max_static = 25, digits = 3) {
  if (is.null(df) || !nrow(df)) { cat("\n_No rows._\n"); return(invisible(NULL)) }

  # `df[num]` with a logical is data.frame idiom for COLUMN selection. On a
  # data.table the same expression selects ROWS, and errors outright when the
  # logical length does not match nrow -- which is exactly what happened when a
  # 7-row propeller result was handed a 9-element column mask. Normalise to a
  # plain data.frame first so the rounding below means the same thing for both.
  df <- as.data.frame(df)
  num <- vapply(df, is.numeric, logical(1))
  if (any(num)) df[num] <- lapply(df[num], function(x) round(x, digits))

  if (.have("DT")) {
    w <- DT::datatable(
      df, rownames = FALSE, caption = caption, filter = "top",
      extensions = if (.have("DT")) "Buttons" else NULL,
      options = list(
        pageLength = page, scrollX = TRUE, dom = "Bfrtip",
        buttons = list(list(extend = "csv",   text = "Download CSV"),
                       list(extend = "excel", text = "Download Excel")),
        searching = searchable,
        columnDefs = list(list(className = "dt-right", targets = which(num) - 1L))
      ))
    return(w)
  }
  # Fallback: a plain table, truncated so a lean container does not emit
  # thousands of rows into the HTML.
  knitr::kable(utils::head(df, max_static), caption = caption)
}

# -----------------------------------------------------------------------------
# Collapsed detail and small multiples
# -----------------------------------------------------------------------------

scratch_details <- function(summary_text, expr) {
  cat(sprintf("\n<details>\n<summary>%s</summary>\n\n", summary_text))
  force(expr)
  cat("\n\n</details>\n")
  invisible(NULL)
}

#' Small multiples: one panel per group, on a shared scale so panels are
#' comparable. Preferred over a tabset when the reader needs to compare across
#' groups rather than inspect one at a time.
scratch_grid <- function(plots, ncol = 3, guides = "collect") {
  plots <- plots[!vapply(plots, is.null, logical(1))]
  if (!length(plots)) return(invisible(NULL))
  if (.have("patchwork")) {
    p <- patchwork::wrap_plots(plots, ncol = ncol, guides = guides)
    if (.have("patchwork")) p <- p + patchwork::plot_annotation(tag_levels = "A")
    return(p)
  }
  for (x in plots) print(x)
  invisible(NULL)
}

# -----------------------------------------------------------------------------
# Report scaffolding
# -----------------------------------------------------------------------------

#' A standard header: what this stage did, on what, with which settings.
scratch_stage_header <- function(stage, n_cells = NULL, n_features = NULL,
                                 n_samples = NULL, settings = NULL) {
  tiles <- list()
  if (!is.null(n_cells))    tiles[["Cells"]]    <- list(value = format(n_cells, big.mark = ","))
  if (!is.null(n_features)) tiles[["Features"]] <- list(value = format(n_features, big.mark = ","))
  if (!is.null(n_samples))  tiles[["Samples"]]  <- list(value = n_samples)
  if (length(tiles)) scratch_tiles(tiles)
  if (!is.null(settings) && length(settings)) {
    scratch_details("Settings used for this stage",
      cat(paste(sprintf("- `%s` = `%s`", names(settings), unlist(settings)), collapse = "\n"), "\n"))
  }
  invisible(NULL)
}

message("[scratch_report] smart-report components active — tabsets, interactive plots, tables, findings")

# =============================================================================
# Report architecture
# =============================================================================
#
# A stage report is not a figure dump. It answers, in this order:
#
#   VERDICT    is anything wrong, and where          <- first screen, always
#   HEADLINE   what happened, in sentences           <- computed, never static
#   EVIDENCE   per-sample drill-down                 <- opt-in
#   COMPARE    cross-sample, shared scales           <- opt-in
#   DATA       every figure's underlying table       <- downloadable
#   METHODS    parameters, versions, session         <- collapsed
#
# Everything above EVIDENCE fits on one screen. The reader should be able to
# close the tab after 10 seconds knowing whether to care.
# =============================================================================

.SX_FINDINGS <- new.env(parent = emptyenv()); .SX_FINDINGS$rows <- list()

#' Register a finding. Severity ranks it; the report surfaces the worst first,
#' and everything registered here is also written to findings.json for
#' downstream consumers (the interpretation agent, dashboards, CI gates).
scratch_register <- function(id, severity = c("info","good","warning","critical"),
                             headline, detail = NULL, stage = NULL,
                             value = NULL, scope = NULL) {
  severity <- match.arg(severity)
  .SX_FINDINGS$rows[[length(.SX_FINDINGS$rows) + 1]] <- list(
    id = id, severity = severity, headline = headline,
    detail = detail %||% "", stage = stage %||% "",
    value = value, scope = scope %||% "")
  invisible(NULL)
}

scratch_findings_json <- function(path = "data/findings.json", stage = NULL) {
  rows <- .SX_FINDINGS$rows
  if (!length(rows)) return(invisible(NULL))
  if (!is.null(stage)) rows <- lapply(rows, function(r) { if (!nzchar(r$stage)) r$stage <- stage; r })
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  if (.have("jsonlite")) {
    jsonlite::write_json(rows, path, auto_unbox = TRUE, pretty = TRUE, null = "null")
  }
  invisible(path)
}

#' The verdict bar: per-unit pass/warn/fail, readable at a glance, before any
#' figure. Status carries a glyph and a label, never colour alone.
scratch_verdict <- function(status, labels = NULL, title = "Sample status",
                            legend = TRUE) {
  pal <- c(good = "#0ca30c", warning = "#fab219", critical = "#d03b3b", info = "#8a9095")
  gly <- c(good = "&#10003;", warning = "&#9888;", critical = "&#10007;", info = "&#8211;")
  if (is.null(labels)) labels <- names(status)

  css <- '
<style>
.sx-verdict{display:flex;flex-wrap:wrap;gap:6px;margin:.6rem 0 1rem}
.sx-v{display:inline-flex;align-items:center;gap:6px;border:1px solid rgba(128,128,128,.3);
      border-radius:6px;padding:5px 10px;font-size:.78rem;font-variant-numeric:tabular-nums}
.sx-v .g{font-weight:700}
.sx-legend{font-size:.72rem;opacity:.7;margin:.2rem 0 1rem}
</style>'
  out <- c(css, sprintf('<div style="font-size:.7rem;letter-spacing:.08em;text-transform:uppercase;opacity:.65">%s</div>', title),
           '<div class="sx-verdict">')
  for (i in seq_along(status)) {
    s <- as.character(status[i]); if (!s %in% names(pal)) s <- "info"
    out <- c(out, sprintf('<span class="sx-v"><span class="g" style="color:%s">%s</span>%s</span>',
                          pal[[s]], gly[[s]], labels[i]))
  }
  out <- c(out, '</div>')
  if (legend) {
    n <- table(factor(status, levels = names(pal)))
    out <- c(out, sprintf('<div class="sx-legend">%s</div>',
      paste(sprintf("%s %s", gly[names(n)][n > 0], paste0(n[n > 0], " ", names(n)[n > 0])), collapse = " &nbsp;·&nbsp; ")))
  }
  cat(paste(out, collapse = "\n"), "\n")
  invisible(NULL)
}

#' Print all registered findings, worst first. This is the report's headline.
scratch_headline <- function(max_show = 6) {
  rows <- .SX_FINDINGS$rows
  if (!length(rows)) { cat("\n_No findings registered._\n"); return(invisible(NULL)) }
  ord <- c(critical = 1, warning = 2, good = 3, info = 4)
  rows <- rows[order(vapply(rows, function(r) ord[[r$severity]], numeric(1)))]
  typ <- c(critical = "important", warning = "warning", good = "tip", info = "note")
  for (r in utils::head(rows, max_show)) {
    scratch_finding(c(paste0("**", r$headline, "**"),
                      if (nzchar(r$detail)) paste0("\n", r$detail) else NULL),
                    type = typ[[r$severity]])
  }
  if (length(rows) > max_show)
    cat(sprintf("\n_%d further finding(s) in the data section._\n", length(rows) - max_show))
  invisible(NULL)
}

#' Scale-aware point layer. Above `raster_above` cells the layer is rasterised,
#' which is what keeps a 200k-cell UMAP from producing a 90 MB vector file that
#' no viewer will open.
scratch_points <- function(n_cells, raster_above = 20000, ...) {
  if (n_cells > raster_above && .have("ggrastr")) ggrastr::rasterise(geom_point(...), dpi = 300)
  else geom_point(...)
}

#' Interactivity has a cost: every plotly widget embeds ~3 MB. Use it where
#' hovering answers a question, and stay static elsewhere.
scratch_should_interact <- function(n_rows, budget = 5000) {
  isTRUE(n_rows <= budget) && .have("plotly")
}

# -----------------------------------------------------------------------------
# Offline interpretation agent
#
# Embeds the report's facts, tables and findings into the page, plus a query
# engine that answers from them with no network. On a machine that happens to
# have a local model serving (Ollama, vLLM, llama.cpp) the agent uses it to
# phrase the answer — but the numbers always come from the deterministic
# computation, never from the model.
# -----------------------------------------------------------------------------

scratch_agent <- function(tables = list(), stage = NULL, project = NULL,
                          viz_dir = ".", height = NULL) {
  js  <- file.path(viz_dir, "scratch_agent.js")
  css <- file.path(viz_dir, "scratch_agent.css")
  for (alt in c("assets/viz", ".")) {
    if (!file.exists(js))  js  <- file.path(alt, "scratch_agent.js")
    if (!file.exists(css)) css <- file.path(alt, "scratch_agent.css")
  }
  if (!file.exists(js)) { cat("\n_Agent assets not staged._\n"); return(invisible(NULL)) }

  payload <- list(
    findings = .SX_FINDINGS$rows,
    tables   = lapply(tables, function(d) { d <- as.data.frame(d); d[] <- lapply(d, as.character); d }),
    meta     = list(stage = stage %||% "", project = project %||% "",
                    generated = format(Sys.time(), "%Y-%m-%d %H:%M"))
  )
  json <- if (.have("jsonlite")) jsonlite::toJSON(payload, auto_unbox = TRUE, dataframe = "rows", null = "null")
          else "{}"

  cat("\n```{=html}\n")
  cat("<style>\n", paste(readLines(css, warn = FALSE), collapse = "\n"), "\n</style>\n")
  cat('<div id="scratch-agent"></div>\n')
  cat("<script>window.SCRATCH_DATA = ", json, ";</script>\n", sep = "")
  cat("<script>\n", paste(readLines(js, warn = FALSE), collapse = "\n"), "\n</script>\n")
  cat("```\n")
  invisible(NULL)
}

# -----------------------------------------------------------------------------
# Linked filtering: a dropdown that drives the figures
#
# Preferred over a tabset once a cohort has more than a handful of samples. A
# tab strip with 40 labels wraps over several lines and is unusable; a select
# control stays one line at any cohort size, allows multi-select, and — because
# every linked figure reads the same SharedData — the page carries ONE copy of
# the data instead of one per tab.
# -----------------------------------------------------------------------------

#' Wrap a data frame for linked filtering.
#' @param key column uniquely identifying a row (usually the cell/clone id)
#' @param group name tying several views together; views sharing a group filter
#'   as one.
scratch_shared <- function(df, key, group = "scratch") {
  if (!.have("crosstalk")) return(df)
  crosstalk::SharedData$new(as.data.frame(df), key = as.formula(paste0("~", key)), group = group)
}

#' A dropdown (or checkbox/slider) bound to a shared dataset.
scratch_filter <- function(shared, column, label = NULL,
                           type = c("select", "checkbox", "slider"), multiple = TRUE) {
  type <- match.arg(type)
  if (!.have("crosstalk") || !inherits(shared, "SharedData")) return(invisible(NULL))
  f <- as.formula(paste0("~", column))
  lab <- label %||% gsub("_", " ", column)
  switch(type,
    select   = crosstalk::filter_select(paste0("flt_", column), lab, shared, f, multiple = multiple),
    checkbox = crosstalk::filter_checkbox(paste0("flt_", column), lab, shared, f),
    slider   = crosstalk::filter_slider(paste0("flt_", column), lab, shared, f))
}

#' Lay out filters above their figures.
scratch_filter_layout <- function(filters, views, widths = NULL) {
  if (!.have("crosstalk")) { for (v in views) print(v); return(invisible(NULL)) }
  htmltools::tagList(
    htmltools::div(class = "sx-filters",
                   crosstalk::bscols(widths = widths %||% rep(12 %/% max(length(filters), 1), length(filters)),
                                     !!!filters)),
    htmltools::div(class = "sx-views", htmltools::tagList(views))
  )
}

# -----------------------------------------------------------------------------
# Hover that actually reports values
#
# ggplotly's default tooltip echoes the aesthetic mapping, which on a bar chart
# reads "x: GBM_DFCI1, y: 13860" — the axis labels the reader can already see.
# A hovertemplate states the measure and formats the number.
# -----------------------------------------------------------------------------

scratch_hover <- function(p, template = NULL, ...) {
  if (!.have("plotly")) return(p)
  w <- if (inherits(p, "plotly")) p else scratch_interact(p, ...)
  if (!is.null(template)) {
    w <- try(plotly::style(w, hovertemplate = template), silent = TRUE)
    if (inherits(w, "try-error")) return(p)
  }
  plotly::layout(w, hovermode = "closest",
                 hoverlabel = list(bgcolor = "#ffffff", bordercolor = "#c9cdc9",
                                   font = list(size = 11)))
}

#' Prose that introduces a section: what the figure shows, how to read it, and
#' what would count as a problem. A figure without this is a picture, not a
#' result.
scratch_section <- function(what, how = NULL, watch = NULL) {
  cat("\n::: {.sx-lede}\n")
  cat(what, "\n")
  if (!is.null(how))   cat("\n**How to read it.** ", how, "\n", sep = "")
  if (!is.null(watch)) cat("\n**What to watch for.** ", watch, "\n", sep = "")
  cat(":::\n\n")
  invisible(NULL)
}

#' Render several tables with headings, keeping htmlwidget dependencies.
#'
#' A DT/plotly widget printed inside knit_child loses its JS dependency
#' registration and renders as nothing. Building an htmltools tagList in a normal
#' chunk keeps the dependencies attached to the parent document, which is why
#' this exists rather than reusing scratch_tabset() for tables.
scratch_table_set <- function(tables, level = 3, page = 10, collapse_after = 3) {
  tables <- tables[!vapply(tables, is.null, logical(1))]
  if (!length(tables)) return(htmltools::tagList())
  if (!.have("htmltools")) { for (d in tables) print(knitr::kable(utils::head(d, 20))); return(invisible(NULL)) }

  htmltools::tagList(lapply(seq_along(tables), function(i) {
    nm <- names(tables)[i]; d <- as.data.frame(tables[[i]])
    hdr <- htmltools::tags[[paste0("h", level)]](
      nm, htmltools::tags$span(sprintf("  %s rows x %s cols", format(nrow(d), big.mark = ","), ncol(d)),
                               style = "font-weight:400;opacity:.55;font-size:.8em"))
    body <- scratch_table(d, page = page)
    if (i > collapse_after) {
      htmltools::tags$details(htmltools::tags$summary(nm), body)
    } else {
      htmltools::tagList(hdr, body)
    }
  }))
}
