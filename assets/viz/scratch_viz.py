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
matplotlib.use("Agg")                       # headless: notebooks render in containers
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
        "font.size":       9,
        "axes.titlesize":  10,
        "axes.titleweight": "bold",
        "axes.titlelocation": "left",
        "axes.labelsize":  9,
        "xtick.labelsize": 8,
        "ytick.labelsize": 8,
        "legend.fontsize": 8,

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
        sc.settings.autoshow = False
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
        "font.size":       9,
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
