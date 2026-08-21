"""
scratch_viz.py — publication figure system for SCRATCH-E2E (python side)

The R counterpart (scratch_viz.R) covers the 37 knitr notebooks. This covers the
python ones — CellTypist, MetaTIME, GeneVector — plus the NMF engine's
diagnostics, so a figure produced on either side of the pipeline looks the same.

Import it and the defaults apply; nothing else in the notebook has to change:

    import scratch_viz as sv        # rcParams + scanpy defaults now set
    sv.save(fig, "umap_celltypist")  # -> .png (400 dpi) AND .pdf

The palette is identical to the R side, in the same slot order, so a cell type
is the same colour whichever engine drew it.
"""

from __future__ import annotations

import os
import csv
import warnings

import matplotlib

# Backend selection decides whether ANY figure appears in a jupyter-engine
# notebook. Forcing "Agg" unconditionally -- which this did -- replaces the
# IPython inline backend, and the inline backend is the thing that captures a
# figure into the rendered document. With Agg, plt.show() is a silent no-op:
# CellTypist declared six figures and emitted ONE, with the rest showing up as
# the repr of the Axes object. Setting sc.settings.autoshow = True was necessary
# but useless on its own, because there was nothing downstream to catch the draw.
#
# So: inside an IPython kernel, ask for inline; everywhere else (plain Rscript /
# python in a container, no display) fall back to Agg as before.
def _activate_backend():
    try:
        from IPython import get_ipython
        ip = get_ipython()
        if ip is not None and type(ip).__name__ == "ZMQInteractiveShell":
            ip.run_line_magic("matplotlib", "inline")
            return
    except Exception:
        pass
    matplotlib.use("Agg")


_activate_backend()
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap, to_hex
from cycler import cycler

# -----------------------------------------------------------------------------
# Palette — mirrors SCRATCH_PAL in scratch_viz.R, same order
# -----------------------------------------------------------------------------

CATEGORICAL = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100",
               "#e87ba4", "#008300", "#4a3aa7", "#e34948"]

SEQUENTIAL = ["#cde2fb", "#9ec5f4", "#6da7ec", "#3987e5",
              "#256abf", "#184f95", "#0d366b"]

DIVERGING = ["#2a78d6", "#f0efec", "#e34948"]

STATUS = {"good": "#0ca30c", "warning": "#fab219",
          "serious": "#ec835a", "critical": "#d03b3b"}

COMPARTMENT = {"malignant": "#eb6834", "non_malignant": "#2a78d6",
               "ambiguous": "#8a9095"}

INK, INK2, MUTED = "#0f1417", "#4e5a60", "#898781"
GRID, AXIS, SURFACE = "#e4e6e3", "#b9bfc2", "#ffffff"

cmap_sequential = LinearSegmentedColormap.from_list("scratch_seq", SEQUENTIAL)
cmap_diverging  = LinearSegmentedColormap.from_list("scratch_div", DIVERGING)
try:
    matplotlib.colormaps.register(cmap_sequential, name="scratch_seq", force=True)
    matplotlib.colormaps.register(cmap_diverging,  name="scratch_div", force=True)
except Exception:
    pass

# Journal column widths, millimetres -> inches
WIDTH_MM = {"single": 89.0, "onehalf": 120.0, "double": 183.0}


def size(width="single", ratio=0.75):
    w = WIDTH_MM.get(width, width) if isinstance(width, str) else float(width)
    return (w / 25.4, (w * ratio) / 25.4)


# -----------------------------------------------------------------------------
# Global defaults
#
# Set on import. Journals reproduce figures small, so the baseline assumes 89 mm
# and still-legible type: hairline spines, recessive grid, no chartjunk.
# -----------------------------------------------------------------------------

def apply_defaults():
    plt.rcParams.update({
        "figure.figsize":   size("single"),
        "figure.dpi":       120,
        "savefig.dpi":      400,
        "savefig.bbox":     "tight",
        "savefig.pad_inches": 0.02,
        "savefig.facecolor": SURFACE,
        "figure.facecolor":  SURFACE,
        "axes.facecolor":    SURFACE,

        "font.family":     "sans-serif",
        "font.sans-serif": ["Helvetica", "Arial", "DejaVu Sans"],
        # Kept in step with theme_scratch(base_size = 12) on the R side, so a
        # scanpy figure and a ggplot figure in the same report read at the same
        # size. Ratios below are the originals scaled from a base of 9.
        "font.size":       12,
        "axes.titlesize":  13,
        "axes.titleweight": "bold",
        "axes.titlelocation": "left",
        "axes.labelsize":  12,
        "xtick.labelsize": 11,
        "ytick.labelsize": 11,
        "legend.fontsize": 11,

        "axes.edgecolor":  AXIS,
        "axes.linewidth":  0.6,
        "axes.labelcolor": INK2,
        "axes.titlecolor": INK,
        "text.color":      INK2,
        "xtick.color":     MUTED,
        "ytick.color":     MUTED,
        "xtick.major.width": 0.6,
        "ytick.major.width": 0.6,

        "axes.grid":       True,
        "axes.axisbelow":  True,
        "grid.color":      GRID,
        "grid.linewidth":  0.5,
        "axes.spines.top":   False,
        "axes.spines.right": False,

        "legend.frameon":  False,
        "image.cmap":      "scratch_seq",
        "axes.prop_cycle": cycler(color=CATEGORICAL),

        # Vector output that journals accept: real text, embedded fonts.
        "pdf.fonttype": 42,
        "ps.fonttype":  42,
        "svg.fonttype": "none",
    })

    # scanpy's set_figure_params REWRITES rcParams — including prop_cycle, which
    # is how its default vega palette was winning over ours. Call it first, then
    # re-assert the settings we care about on top.
    try:
        import scanpy as sc
        sc.set_figure_params(dpi=120, dpi_save=400, frameon=False,
                             vector_friendly=True, format="png")
        # MUST stay True. With autoshow=False every `sc.pl.*` call returns its
        # Axes instead of calling plt.show(), so quarto's jupyter engine captured
        # only the repr -- `<Axes: xlabel='% of total counts'>` as literal text --
        # and FIVE of the six CellTypist figures silently vanished from the
        # report. Nothing here needs the return value: `save()` below is never
        # called by a notebook, so suppressing the draw bought nothing at all.
        sc.settings.autoshow = True
        # scanpy reads this for categorical obs colouring.
        sc.settings.categorical_palette = CATEGORICAL
    except ImportError:
        pass
    except Exception:
        pass

    plt.rcParams.update({
        "axes.prop_cycle": cycler(color=CATEGORICAL),
        "image.cmap":      "scratch_seq",
        "savefig.dpi":     400,
        "pdf.fonttype":    42,
        "ps.fonttype":     42,
        "svg.fonttype":    "none",
        # This block runs AFTER scanpy's set_figure_params, which rewrites
        # rcParams -- and scanpy sets the type sizes RELATIVELY ('medium',
        # 'large'), which then resolve against whatever font.size ends up being.
        # Setting them in the block above is therefore useless: scanpy clobbers
        # it, and ticks came out the same size as axis labels. Every explicit
        # size has to be re-asserted here, after scanpy, to stick.
        "font.size":       12,
        "axes.titlesize":  13,
        "axes.labelsize":  12,
        "xtick.labelsize": 11,
        "ytick.labelsize": 11,
        "legend.fontsize": 11,
        "axes.titlelocation": "left",
        "axes.titleweight":   "bold",
        "figure.facecolor":  SURFACE,
        "axes.facecolor":    SURFACE,
        "savefig.facecolor": SURFACE,
    })


apply_defaults()


# -----------------------------------------------------------------------------
# Saving: raster AND vector from one call
# -----------------------------------------------------------------------------

_MANIFEST = []


def save(fig, name, width="single", ratio=0.75, outdir="figures",
         dpi=400, vector=True, caption=None, tight=True):
    """Write <outdir>/<name>.png and .pdf, and record the figure."""
    os.makedirs(outdir, exist_ok=True)
    stem = os.path.join(outdir, os.path.splitext(name)[0])

    if fig is None:
        fig = plt.gcf()
    w, h = size(width, ratio)
    try:
        fig.set_size_inches(w, h)
    except Exception:
        pass
    if tight:
        try:
            fig.tight_layout()
        except Exception:
            pass

    ok_png = ok_pdf = False
    try:
        fig.savefig(stem + ".png", dpi=dpi, bbox_inches="tight",
                    facecolor=SURFACE)
        ok_png = True
    except Exception as exc:                      # a figure must never kill a run
        warnings.warn(f"PNG save failed for {name}: {exc}")

    if vector:
        try:
            fig.savefig(stem + ".pdf", bbox_inches="tight", facecolor=SURFACE)
            ok_pdf = True
        except Exception as exc:
            warnings.warn(f"PDF save failed for {name}: {exc}")

    _MANIFEST.append({
        "figure": os.path.basename(stem), "png": ok_png, "pdf": ok_pdf,
        "width_mm": round(w * 25.4, 1), "height_mm": round(h * 25.4, 1),
        "dpi": dpi, "caption": caption or "",
    })
    return fig


def manifest(path="figures/figure_manifest.csv"):
    if not _MANIFEST:
        return None
    os.makedirs(os.path.dirname(path) or ".", exist_ok=True)
    with open(path, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(_MANIFEST[0].keys()))
        w.writeheader()
        w.writerows(_MANIFEST)
    return path


# -----------------------------------------------------------------------------
# Helpers
# -----------------------------------------------------------------------------

def colours(n=None):
    """Categorical slots in fixed order. Never fabricates a 9th hue."""
    if n is None:
        return list(CATEGORICAL)
    if n <= len(CATEGORICAL):
        return CATEGORICAL[:n]
    warnings.warn(
        f"{n} levels requested but only {len(CATEGORICAL)} validated slots exist; "
        "fold excess levels to 'Other' or facet."
    )
    return [CATEGORICAL[i % len(CATEGORICAL)] for i in range(n)]


def fold_other(values, max_levels=8, other="Other"):
    """Cap categorical levels; the remainder folds rather than inventing hues."""
    from collections import Counter
    counts = Counter(str(v) for v in values)
    keep = {k for k, _ in counts.most_common(max_levels)}
    return [v if str(v) in keep else other for v in values]


def palette_for(categories, max_levels=8):
    """Stable category -> colour map, so a cell type keeps its colour."""
    cats = list(dict.fromkeys(str(c) for c in categories))
    pal = colours(min(len(cats), max_levels))
    out = {c: pal[i % len(pal)] for i, c in enumerate(cats[:max_levels])}
    for c in cats[max_levels:]:
        out[c] = "#c9cdc9"
    return out


def style_axes(ax, title=None, xlabel=None, ylabel=None, grid="y"):
    """Apply the house look to an axis produced by a library that ignores rcParams."""
    if title:  ax.set_title(title, loc="left", fontweight="bold", color=INK)
    if xlabel is not None: ax.set_xlabel(xlabel, color=INK2)
    if ylabel is not None: ax.set_ylabel(ylabel, color=INK2)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    for s in ("left", "bottom"):
        ax.spines[s].set_color(AXIS)
        ax.spines[s].set_linewidth(0.6)
    ax.grid(axis=grid if grid in ("x", "y", "both") else "y",
            color=GRID, linewidth=0.5)
    ax.set_axisbelow(True)
    return ax


print("[scratch_viz] publication figure system active — rcParams, palette, vector export")


def umap(adata, color, figsize=(9, 9), legend_ncol=1, legend_fontsize=7, **kw):
    """sc.pl.umap with a legend that does not eat the plot.

    scanpy picks the legend column count itself -- 2 columns above 14 categories,
    3 above 30 -- and draws it at the ambient font size. With 25 CellTypist
    labels that produced a two-column block of large text occupying two thirds of
    the figure, leaving the embedding squeezed into the left third and unreadable.

    There is no scanpy argument for this, so the legend has to be rebuilt after
    the fact: draw with show=False, replace the legend with a single narrow
    column at a smaller size, then show. Continuous colourings get a colorbar
    rather than a legend and are passed through untouched.
    """
    import scanpy as sc

    kw.setdefault("frameon", False)
    kw.setdefault("sort_order", False)
    with plt.rc_context({"figure.figsize": figsize}):
        axes = sc.pl.umap(adata, color=color, show=False, **kw)
        for ax in (axes if isinstance(axes, (list, tuple)) else [axes]):
            if ax is None or ax.get_legend() is None:
                continue
            handles, labels = ax.get_legend_handles_labels()
            if not handles:
                continue
            ax.legend(handles, labels, loc="center left",
                      bbox_to_anchor=(1.01, 0.5), ncol=legend_ncol,
                      fontsize=legend_fontsize, frameon=False,
                      markerscale=0.7, handletextpad=0.3,
                      labelspacing=0.25, borderaxespad=0.0)
        plt.show()


def violin_selector(adata, keys, groupby, title="Signature score",
                    height=520, points=False):
    """One violin panel with a dropdown, instead of N stacked static panels.

    The CellTypist notebook looped `sc.pl.violin` over ~16 signatures, emitting a
    separate full-size figure for each. They share an x axis and a purpose, so
    the reader is really choosing between them -- which is a dropdown, not a
    stack. One plotly trace per signature; the menu toggles which is visible.

    Returns a plotly Figure (display it as the cell's value), or None if plotly
    or the requested columns are unavailable -- callers should fall back to
    sc.pl.violin in that case rather than fail.
    """
    try:
        import plotly.graph_objects as go
    except Exception:
        return None
    if groupby not in adata.obs.columns:
        return None
    keys = [k for k in keys if k in adata.obs.columns]
    if not keys:
        return None

    grp = adata.obs[groupby].astype(str).values
    order = sorted(set(grp), key=lambda v: (len(v), v))

    fig = go.Figure()
    for i, k in enumerate(keys):
        fig.add_trace(go.Violin(
            x=grp, y=adata.obs[k].astype(float).values,
            name=k, visible=(i == 0), points="outliers" if points else False,
            box_visible=True, meanline_visible=False,
            line=dict(width=1), marker=dict(size=2),
            hovertemplate=f"{groupby}: %{{x}}<br>{k}: %{{y:.3f}}<extra></extra>"))

    fig.update_layout(
        title=dict(text=f"<b>{title}</b><br><sup>{keys[0]}</sup>",
                   x=0, xanchor="left"),
        xaxis=dict(title=groupby, categoryorder="array", categoryarray=order),
        yaxis=dict(title="score"), showlegend=False, height=height,
        margin=dict(t=90, r=10),
        updatemenus=[dict(
            type="dropdown", direction="down", showactive=True,
            x=1.02, xanchor="left", y=1, yanchor="top",
            buttons=[dict(method="update", label=k,
                          args=[{"visible": [j == i for j in range(len(keys))]},
                                {"title": {"text": f"<b>{title}</b><br><sup>{k}</sup>",
                                           "x": 0, "xanchor": "left"}}])
                     for i, k in enumerate(keys)])])
    return fig


# =============================================================================
# In-page agent, python side
# =============================================================================
#
# The R notebooks get the agent automatically: scratch_report.R installs a knitr
# `document` hook that mounts it once per render. The Jupyter-engine notebooks
# (CellTypist is the only one today) never source that file, so they rendered
# without a panel at all — verified: zero occurrences of `scratch_agent` in
# notebook_celltypist.html.
#
# This is the same payload contract, built from python. Call it from a cell with
# `#| output: asis` at the end of the notebook.

_SX_FINDINGS = []
_SX_TABLES = {}


def finding(text, type="note", title=""):
    """Render a Quarto callout AND register it for the agent."""
    _SX_FINDINGS.append({"type": type, "title": title, "text": str(text)})
    hdr = f'::: {{.callout-{type}' + (f' title="{title}"' if title else "") + "}"
    print("\n" + hdr + "\n" + str(text) + "\n:::\n")


def table(df, caption=None, max_rows=2000):
    """Register a dataframe for the agent. Returns it so the cell still displays."""
    name = caption or f"table_{len(_SX_TABLES) + 1}"
    try:
        d = df.head(max_rows)
        _SX_TABLES[name] = [
            {k: ("" if v is None else str(v)) for k, v in row.items()}
            for row in d.to_dict(orient="records")
        ]
    except Exception:
        pass
    return df


def agent(stage=None, project=None, viz_dir="."):
    """Emit the agent panel. Mirrors scratch_agent() in scratch_report.R."""
    import json
    import os

    js = css = None
    for root in (viz_dir, "assets/viz", "."):
        j, c = os.path.join(root, "scratch_agent.js"), os.path.join(root, "scratch_agent.css")
        if os.path.exists(j) and os.path.exists(c):
            js, css = j, c
            break
    if not js:
        print("\n_Agent assets not staged._\n")
        return

    payload = {
        "findings": _SX_FINDINGS,
        "tables": _SX_TABLES,
        "meta": {"stage": stage or "", "project": project or ""},
    }
    print("\n```{=html}")
    print("<style>\n" + open(css).read() + "\n</style>")
    print('<div id="scratch-agent"></div>')
    print("<script>window.SCRATCH_DATA = " + json.dumps(payload) + ";</script>")
    print("<script>\n" + open(js).read() + "\n</script>")
    print("```\n")
