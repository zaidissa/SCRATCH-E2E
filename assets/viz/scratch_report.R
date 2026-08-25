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

  # Register as well as render. `.SX_FINDINGS` was only ever written by
  # scratch_register(), which no notebook calls — every notebook uses THIS
  # function. So the in-page agent received `findings: []` on every report and
  # could not summarise the one thing a reader most wants summarised. Recording
  # it here means all ~40 existing call sites are covered without touching one.
  .SX_FINDINGS$rows[[length(.SX_FINDINGS$rows) + 1L]] <- list(
    type  = type,
    title = if (is.null(title)) "" else title,
    text  = paste(text, collapse = " ")
  )
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

# Every table a report renders is also REGISTERED here, so the in-page agent can
# answer from the same numbers the reader is looking at. Previously the agent had
# to be handed its tables by hand, which is why it only ever appeared in
# stage_report.qmd and never in the per-module notebooks anyone actually opens.
.SX_TABLES <- new.env(parent = emptyenv()); .SX_TABLES$rows <- list()

scratch_table <- function(df, caption = NULL, page = 10, searchable = TRUE,
                          max_static = 25, digits = 3) {
  if (is.null(df) || !nrow(df)) { cat("\n_No rows._\n"); return(invisible(NULL)) }

  # Register before rendering, capped so a 24k-row cell table cannot inflate the
  # page: the agent needs enough rows to compute over, not the whole matrix.
  local({
    nm <- caption
    if (is.null(nm) || !nzchar(nm)) nm <- paste0("table_", length(.SX_TABLES$rows) + 1L)
    nm <- substr(gsub("[^A-Za-z0-9 _-]", "", nm), 1, 60)
    cap <- getOption("scratch.agent.maxrows", 2000L)
    d <- as.data.frame(df)
    if (nrow(d) > cap) d <- d[seq_len(cap), , drop = FALSE]
    .SX_TABLES$rows[[nm]] <- d
  })

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
  # DT is absent from FIVE of the eight analysis images (cnv, celltrajectory,
  # cellcommunication, tumorheterogenity, batchcorr), so this fallback is not an
  # edge case -- it is what most of the pipeline actually renders. A static kable
  # meant those stages lost search, sort and download without anything saying so.
  #
  # This emits a self-contained interactive table instead: no R package, no CDN,
  # ~2 KB of inline JS per table. Search, click-to-sort (numeric-aware) and CSV
  # download, working identically in every container.
  .sx_interactive_table(df, caption = caption, page = page)
}

.sx_interactive_table <- function(df, caption = NULL, page = 10) {
  esc <- function(x) {
    x <- as.character(x); x[is.na(x)] <- ""
    x <- gsub("&", "&amp;", x, fixed = TRUE)
    x <- gsub("<", "&lt;",  x, fixed = TRUE)
    gsub(">", "&gt;", x, fixed = TRUE)
  }
  id   <- paste0("sxt", substr(gsub("[^a-z0-9]", "", tolower(paste0(caption, runif(1)))), 1, 12),
                 sample(1000:9999, 1))
  cols <- names(df)
  head_html <- paste0("<th data-c='", seq_along(cols) - 1L, "'>", esc(cols), "</th>", collapse = "")
  body_html <- paste0(
    apply(df, 1, function(r) paste0("<tr>", paste0("<td>", esc(r), "</td>", collapse = ""), "</tr>")),
    collapse = "")
  csv <- paste0(paste(cols, collapse = ","), "\n",
                paste(apply(df, 1, function(r)
                  paste(ifelse(grepl("[,\"]", r), paste0('"', gsub('"', '""', r), '"'), r),
                        collapse = ",")), collapse = "\n"))
  csv_uri <- paste0("data:text/csv;charset=utf-8,", utils::URLencode(csv, reserved = TRUE))

  knitr::asis_output(paste0(
    "\n```{=html}\n",
    "<div class='sxt-wrap' id='", id, "'>",
    if (!is.null(caption)) paste0("<div class='sxt-cap'>", esc(caption), "</div>") else "",
    "<div class='sxt-bar'>",
      "<input class='sxt-q' placeholder='Filter rows…' aria-label='Filter rows'>",
      "<span class='sxt-n'></span>",
      "<a class='sxt-dl' download='", gsub("[^A-Za-z0-9]+", "_", caption %||% "table"),
        ".csv' href='", csv_uri, "'>Download CSV</a>",
    "</div>",
    "<div class='sxt-scroll'><table class='sxt'><thead><tr>", head_html,
    "</tr></thead><tbody>", body_html, "</tbody></table></div></div>",
    "<style>",
    ".sxt-wrap{margin:1rem 0;font-size:.86rem}",
    ".sxt-cap{font-weight:600;margin-bottom:.4rem}",
    ".sxt-bar{display:flex;gap:.6rem;align-items:center;margin-bottom:.35rem}",
    ".sxt-q{flex:0 1 240px;padding:4px 8px;border:1px solid rgba(128,128,128,.4);border-radius:6px;font:inherit}",
    ".sxt-n{opacity:.6;font-size:.8em}",
    ".sxt-dl{margin-left:auto;font-size:.8em;text-decoration:none;border:1px solid rgba(128,128,128,.4);",
      "padding:3px 9px;border-radius:6px}",
    ".sxt-scroll{max-height:", 34 * max(page, 5), "px;overflow:auto;",
      "border:1px solid rgba(128,128,128,.25);border-radius:8px}",
    "table.sxt{border-collapse:collapse;width:100%}",
    "table.sxt th{position:sticky;top:0;background:var(--bs-body-bg,#fff);cursor:pointer;",
      "text-align:left;padding:6px 9px;border-bottom:2px solid rgba(128,128,128,.35);white-space:nowrap}",
    "table.sxt th:hover{background:rgba(128,128,128,.12)}",
    "table.sxt td{padding:5px 9px;border-bottom:1px solid rgba(128,128,128,.15);white-space:nowrap}",
    "table.sxt tr:hover td{background:rgba(42,120,214,.07)}",
    "</style>",
    "<script>(function(){",
    "var w=document.getElementById('", id, "');if(!w)return;",
    "var tb=w.querySelector('tbody'),rows=[].slice.call(tb.rows),n=w.querySelector('.sxt-n');",
    "function count(v){n.textContent=v+' of ", nrow(df), " rows';}count(rows.length);",
    "w.querySelector('.sxt-q').addEventListener('input',function(e){",
      "var q=e.target.value.toLowerCase(),k=0;",
      "rows.forEach(function(r){var hit=r.textContent.toLowerCase().indexOf(q)>-1;",
        "r.style.display=hit?'':'none';if(hit)k++;});count(k);});",
    # Numeric-aware sort: a column of numbers must not sort as text, or 10 lands
    # between 1 and 2 -- which is exactly how a 'top N' table misleads.
    "var dir={};w.querySelectorAll('th').forEach(function(th){th.addEventListener('click',function(){",
      "var c=+th.dataset.c;dir[c]=!dir[c];var s=dir[c]?1:-1;",
      "rows.sort(function(a,b){var x=a.cells[c].textContent,y=b.cells[c].textContent;",
        "var nx=parseFloat(x.replace(/[, ]/g,'')),ny=parseFloat(y.replace(/[, ]/g,''));",
        "if(!isNaN(nx)&&!isNaN(ny))return (nx-ny)*s;",
        "return x.localeCompare(y)*s;});",
      "rows.forEach(function(r){tb.appendChild(r);});});});",
    "})();</script>",
    "\n```\n"))
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

#' Mount the agent using everything the report registered as it rendered.
#' Called automatically by the document hook below, so no notebook needs to know
#' the agent exists.
scratch_agent_auto <- function(stage = NULL, project = NULL) {
  if (isFALSE(getOption("scratch.agent", TRUE))) return(invisible(NULL))
  scratch_agent(tables  = .SX_TABLES$rows,
                stage   = stage   %||% getOption("scratch.stage", NULL),
                project = project %||% getOption("scratch.project", NULL))
}

# Appending via knitr's `document` hook is what makes this universal: it runs
# once per rendered notebook, after every chunk, so the panel carries the
# report's COMPLETE set of findings and tables. Wiring it per-notebook would
# mean 40 edits and would still miss whichever one was added next.
.sx_install_agent_hook <- function() {
  if (!requireNamespace("knitr", quietly = TRUE)) return(invisible(NULL))
  if (isFALSE(getOption("scratch.agent", TRUE))) return(invisible(NULL))
  knitr::knit_hooks$set(document = function(x) {
    html <- try(utils::capture.output(scratch_agent_auto()), silent = TRUE)
    if (inherits(html, "try-error") || !length(html)) return(x)
    c(x, "", html)
  })
  invisible(NULL)
}

scratch_agent <- function(tables = list(), stage = NULL, project = NULL,
                          viz_dir = ".", height = NULL) {
  js  <- file.path(viz_dir, "scratch_agent.js")
  css <- file.path(viz_dir, "scratch_agent.css")
  for (alt in c("assets/viz", ".")) {
    if (!file.exists(js))  js  <- file.path(alt, "scratch_agent.js")
    if (!file.exists(css)) css <- file.path(alt, "scratch_agent.css")
  }
  if (!file.exists(js)) { cat("\n_Agent assets not staged._\n"); return(invisible(NULL)) }

  # An empty R list serialises to JSON `[]` (an array), but the agent indexes
  # tables by name and expects `{}`. Naming the empty list forces the object
  # form, so a report with no registered tables still presents a valid shape
  # instead of an array the JS then treats as a table called "0".
  tbl <- lapply(tables, function(d) { d <- as.data.frame(d); d[] <- lapply(d, as.character); d })
  if (!length(tbl)) tbl <- stats::setNames(list(), character(0))

  payload <- list(
    findings = .SX_FINDINGS$rows,
    tables   = tbl,
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

# Install the auto-mount hook at source() time. Every notebook already sources
# this file in its scratch-viz-setup chunk, so every report gets the panel.
# Disable with options(scratch.agent = FALSE).
try(.sx_install_agent_hook(), silent = TRUE)

# =============================================================================
# Cell-type vocabulary matching, and the CNV reference-null cutoff
# =============================================================================
#
# These live here because TWO notebooks need them to agree exactly:
# stratify_compartments.qmd decides the compartment split, and
# notebook_cnv_concordance.qmd draws the confusion matrix that claims to
# describe that same decision. Two copies of one rule silently diverge; the
# comparison then describes a decision the pipeline never made.
#
# The matching problem is two problems wearing one symptom:
#
#   SPELLING     "T_Cells" (scType) vs "T cell" (Azimuth) — punctuation and
#                plurality. Normalisation handles it.
#   GRANULARITY  "Myeloid" vs "Macrophage" — a lineage and a cell type within
#                it. No string transform connects them; that needs a table.
#
# Under exact `%in%` matching, an scType reference list against an Azimuth
# column intersects on exactly ONE label (Fibroblast). That built a CNV null
# from 544 fibroblasts all scoring zero, so mean + 3*SD collapsed to a cutoff of
# 0 and every cell with any CNV signal was called malignant.

scratch_norm_label <- function(x) {
  x <- tolower(trimws(as.character(x)))
  x <- gsub("[^a-z0-9]+", " ", x)
  x <- gsub("\\b(cells|cell)\\b", " ", x)   # drop the "cell(s)" token
  x <- gsub("s\\b", "", x)                  # naive singular
  trimws(gsub(" +", " ", x))
}

# One-directional on purpose: a coarse lineage expands to the finer labels that
# fall under it, never the reverse. Expanding "Macrophage" up to "Myeloid" would
# make a macrophage match a request for the whole myeloid compartment.
.SCRATCH_LINEAGE_ALIASES <- list(
  "myeloid"     = c("macrophage","monocyte","dendritic","dc","pdc","mast",
                    "neutrophil","granulocyte","mono","mdc"),
  "b plasma"    = c("b","plasma","plasmablast"),
  "t"           = c("t","cd4 t","cd8 t","treg","regulatory t"),
  "nk"          = c("nk","nkt","innate lymphoid"),
  "endothelial" = c("endothelial","ec","lymphatic endothelial"),
  "fibroblast"  = c("fibroblast","caf","stromal","smooth muscle","pericyte"),
  "epithelial"  = c("epithelial","tumor","tumour","malignant","carcinoma")
)

scratch_expand_labels <- function(labels) {
  n <- scratch_norm_label(labels)
  unique(c(n, unlist(.SCRATCH_LINEAGE_ALIASES[intersect(n, names(.SCRATCH_LINEAGE_ALIASES))],
                     use.names = FALSE)))
}

#' TRUE for each observed label that falls under any requested lineage.
scratch_match_labels <- function(observed, requested) {
  scratch_norm_label(observed) %in% scratch_expand_labels(requested)
}

#' The CNV cutoff for one sample, from the reference-lineage null.
#'
#' Returns list(cutoff, n_ref, basis). `basis` names the rule that was used, so
#' a report can say which one applied rather than presenting every threshold as
#' if it came from the same place.
#'
#' A null with ZERO SPREAD is refused: mean + k*sd collapses to mean, and when
#' the reference cells all score 0 the cutoff is 0 and `score > cutoff` calls
#' everything malignant. The old `n_ref >= 30` guard never caught it because
#' count was never the problem — spread was.
scratch_cnv_cutoff <- function(scores, is_ref, k = 3, fallback_q = 0.75) {
  nulls <- scores[is_ref & is.finite(scores)]
  sd_n  <- if (length(nulls) >= 2) stats::sd(nulls) else NA_real_
  if (length(nulls) >= 30 && is.finite(sd_n) && sd_n > 0) {
    list(cutoff = mean(nulls) + k * sd_n, n_ref = length(nulls), basis = "per_sample")
  } else {
    # The null has no spread (or too few cells to have one). A plain quantile of
    # the whole distribution is NOT a safe fallback: when most cells score zero,
    # the 75th percentile is also zero and the cutoff is degenerate again under a
    # different name. Take the quantile of the scores that actually EXCEED the
    # null instead, so the threshold sits inside the signal rather than on the
    # floor. If nothing exceeds the null, there is no separation to find and the
    # null's own maximum is returned — every reference cell then stays below the
    # cutoff by construction, which is the one property that must hold.
    above <- scores[is.finite(scores) & scores > max(c(nulls, -Inf), na.rm = TRUE)]
    cut <- if (length(above)) as.numeric(stats::quantile(above, fallback_q, na.rm = TRUE))
           else max(c(nulls, 0), na.rm = TRUE)
    list(cutoff = cut, n_ref = length(nulls),
         basis = if (length(nulls) >= 30) "quantile_degenerate_null"
                 else "quantile_too_few_reference")
  }
}
