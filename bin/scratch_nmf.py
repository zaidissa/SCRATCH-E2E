#!/usr/bin/env python3
"""
Per-sample NMF rank sweep for the SCRATCH metaprogram stage.

WHY THIS EXISTS
---------------
The original stage fits NMF with the R `NMF` package (`method = "snmf/r"`,
`nrun = 2`, ranks 3-5). That package is the dominant cost of the whole
metaprogram stage: it is pure-R/C hybrid, re-derives the full factorisation per
restart, and holds several dense gene x cell copies at once.

`Metaprog_PostNMF.qmd` consumes exactly one thing from each fit object:
`NMF::basis(fit)` — the gene x program loading matrix, per rank. Everything the
R fit carries beyond that (the mixture coefficients, the consensus matrices, the
residual traces) is discarded before the metaprograms are built. So the solver
can be replaced without touching the metaprogram science, provided the basis
matrices come out in the same layout.

That is what this script does.

INPUT
-----
The MatrixMarket export of the *same* `<sample>_preprocessed.rds` the R engine
would have read (see bin/export_preprocessed_mtx.R). No re-preprocessing is
performed here: the log-CPM/10, gene-centred, zero-clipped matrix produced by
Metaprog_NMFprepprocessing.qmd remains authoritative.

OUTPUT
------
  <sample>_basis_matrix.csv    genes x programs, columns "<rank>_<component>"
  <sample>_top_genes.csv       top-N genes per program (long format)
  <sample>_nmf_metrics.csv     per-rank reconstruction error, iterations, timing
  <sample>_nmf_rank.png        reconstruction error vs rank

The basis CSV is the R-compatible artefact: PostNMF reads it and prefixes the
columns with "<sample>." to reproduce the `sample.rank.component` naming the
metaprogram clustering expects.

NOTES ON THE REFERENCE SNIPPET
------------------------------
This implementation differs from the commonly-circulated scanpy/sklearn NMF
snippet in four ways that matter:

  * `fit()` followed by `transform()` solves for W twice. `fit_transform()`
    returns the W already computed during the fit, which is a straight saving of
    one full NNLS pass per rank per restart.
  * sklearn's default `max_iter=200` frequently stops before convergence on
    single-cell matrices and emits a ConvergenceWarning that is easy to miss.
    The default here is 1000 and non-convergence is recorded per rank in the
    metrics CSV rather than being silently accepted.
  * the snippet mutates `W` inside the per-component loop (`W['rank'] = ...`
    then `W.sort_values`), so column positions shift under the loop. Ranking is
    done here on an untouched copy.
  * R's `nrun` performs multiple restarts and keeps the best fit. sklearn has no
    equivalent, so restarts are run explicitly and the lowest-reconstruction-
    error solution is kept, keeping the two engines comparable.
"""

from __future__ import annotations

import argparse
import json
import os
import sys
import time
import warnings

# BLAS threads must be pinned BEFORE numpy/sklearn import, otherwise each
# process grabs every core and the per-sample fan-out oversubscribes the node.
_THREADS = os.environ.get("NMF_NUM_THREADS") or os.environ.get("OMP_NUM_THREADS") or "1"
for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS",
           "NUMEXPR_NUM_THREADS", "VECLIB_MAXIMUM_THREADS"):
    os.environ.setdefault(_v, _THREADS)

import numpy as np
import pandas as pd
from scipy.io import mmread
from sklearn.decomposition import NMF
from sklearn.exceptions import ConvergenceWarning


# ---------------------------------------------------------------------------
# IO
# ---------------------------------------------------------------------------

def load_matrix(prefix: str) -> pd.DataFrame:
    """Load the exported matrix as a dense genes x cells frame.

    NMF needs dense arrays for the coordinate-descent solver anyway, so the
    conversion happens once here rather than repeatedly inside the sweep.
    """
    mtx = mmread(f"{prefix}.mtx")
    arr = np.asarray(mtx.todense() if hasattr(mtx, "todense") else mtx, dtype=np.float32)

    with open(f"{prefix}_genes.tsv") as fh:
        genes = [l.rstrip("\n") for l in fh if l.strip()]
    with open(f"{prefix}_cells.tsv") as fh:
        cells = [l.rstrip("\n") for l in fh if l.strip()]

    if arr.shape != (len(genes), len(cells)):
        raise ValueError(
            f"Label/matrix shape mismatch: matrix {arr.shape} vs "
            f"{len(genes)} genes x {len(cells)} cells"
        )

    # The R preprocessing already clips at zero; guard against float drift
    # reintroducing tiny negatives, which NMF rejects outright.
    neg = arr < 0
    if neg.any():
        smallest = float(arr[neg].min())
        if smallest < -1e-6:
            raise ValueError(f"Input contains negative values (min {smallest}); "
                             "NMF requires a nonnegative matrix.")
        arr[neg] = 0.0

    return pd.DataFrame(arr, index=genes, columns=cells)


def select_hvgs(df: pd.DataFrame, n_keep: int) -> pd.DataFrame:
    """Variance-based gene selection, matching Metaprog_NMF.qmd.

    The R notebook computes the row variance with the sum-of-squares identity
    and keeps the top `hvg_keep`. Reproduced here so both engines factorise the
    same gene set.
    """
    if n_keep <= 0 or n_keep >= df.shape[0]:
        return df
    n = df.shape[1]
    s1 = df.values.sum(axis=1)
    s2 = (df.values ** 2).sum(axis=1)
    var = (s2 - (s1 ** 2) / n) / max(n - 1, 1)
    keep = np.argsort(-var)[:n_keep]
    keep.sort()
    return df.iloc[keep]


# ---------------------------------------------------------------------------
# Factorisation
# ---------------------------------------------------------------------------

def fit_rank(x: np.ndarray, rank: int, args, base_seed: int):
    """Fit one rank with `nrun` restarts; return the best by reconstruction error.

    Returns (W, H, info) where W is genes x rank.
    """
    best = None
    for run in range(max(1, args.nrun)):
        seed = base_seed + run
        # nndsvd* initialisations are deterministic, so extra restarts only add
        # value when the initialisation itself is random. Switch after the first.
        init = args.init if run == 0 else ("nndsvdar" if args.init.startswith("nndsvd") else "random")

        model = NMF(
            n_components=rank,
            init=init,
            solver=args.solver,
            beta_loss=args.beta_loss,
            max_iter=args.max_iter,
            tol=args.tol,
            random_state=seed,
        )

        t0 = time.time()
        with warnings.catch_warnings(record=True) as caught:
            warnings.simplefilter("always", ConvergenceWarning)
            # fit_transform returns the W solved during the fit; calling fit()
            # then transform() would solve for W a second time.
            w = model.fit_transform(x)
            converged = not any(issubclass(c.category, ConvergenceWarning) for c in caught)
        elapsed = time.time() - t0

        info = {
            "rank": rank,
            "run": run,
            "seed": seed,
            "init": init,
            "reconstruction_err": float(model.reconstruction_err_),
            "n_iter": int(model.n_iter_),
            "converged": bool(converged),
            "seconds": round(elapsed, 2),
        }

        if best is None or info["reconstruction_err"] < best[2]["reconstruction_err"]:
            best = (w, model.components_, info)

    return best


def main() -> int:
    p = argparse.ArgumentParser(description="Per-sample NMF rank sweep (SCRATCH metaprograms)")
    p.add_argument("--prefix", required=True, help="input prefix (expects <prefix>.mtx etc.)")
    p.add_argument("--sample", required=True, help="sample id, used for output names")
    p.add_argument("--outdir", default=".", help="output directory")
    p.add_argument("--rank-lb", type=int, default=3)
    p.add_argument("--rank-ub", type=int, default=5)
    p.add_argument("--nrun", type=int, default=2, help="restarts per rank; best kept")
    p.add_argument("--hvg-keep", type=int, default=7000)
    p.add_argument("--min-cells", type=int, default=100)
    p.add_argument("--max-cells", type=int, default=0, help="0 disables downsampling")
    p.add_argument("--genes-per-program", type=int, default=50)
    p.add_argument("--init", default="nndsvda",
                   choices=["nndsvd", "nndsvda", "nndsvdar", "random"])
    p.add_argument("--solver", default="cd", choices=["cd", "mu"])
    p.add_argument("--beta-loss", default="frobenius",
                   choices=["frobenius", "kullback-leibler", "itakura-saito"])
    p.add_argument("--max-iter", type=int, default=1000)
    p.add_argument("--tol", type=float, default=1e-4)
    p.add_argument("--seed", type=int, default=1234)
    args = p.parse_args()

    os.makedirs(args.outdir, exist_ok=True)
    rng = np.random.default_rng(args.seed)

    df = load_matrix(args.prefix)
    n_genes, n_cells = df.shape
    print(f"[{args.sample}] loaded {n_genes} genes x {n_cells} cells", flush=True)

    if n_cells < args.min_cells:
        raise SystemExit(
            f"[{args.sample}] {n_cells} cells < --min-cells {args.min_cells}; "
            "the preprocessing gate should have excluded this sample."
        )

    if args.max_cells and n_cells > args.max_cells:
        pick = rng.choice(n_cells, size=args.max_cells, replace=False)
        pick.sort()
        df = df.iloc[:, pick]
        print(f"[{args.sample}] downsampled to {df.shape[1]} cells", flush=True)

    df = select_hvgs(df, args.hvg_keep)
    print(f"[{args.sample}] {df.shape[0]} genes retained after variance selection", flush=True)

    x = np.ascontiguousarray(df.values, dtype=np.float64)
    genes = df.index.to_list()

    ranks = list(range(args.rank_lb, args.rank_ub + 1))
    if not ranks:
        raise SystemExit(f"Empty rank range: {args.rank_lb}..{args.rank_ub}")

    basis_blocks, metrics, top_rows = [], [], []

    for rank in ranks:
        w, _h, info = fit_rank(x, rank, args, base_seed=args.seed)
        metrics.append(info)

        status = "converged" if info["converged"] else f"NOT converged in {info['n_iter']} iter"
        print(f"[{args.sample}] rank {rank}: err={info['reconstruction_err']:.4f} "
              f"{status} ({info['seconds']}s)", flush=True)

        # genes x rank, columns "<rank>_<component>" — mirrors the R basis layout
        block = pd.DataFrame(
            w, index=genes,
            columns=[f"{rank}_{i}" for i in range(rank)],
        )
        basis_blocks.append(block)

        # Top-N genes per program, ranked on an untouched copy of the loadings.
        n_top = min(args.genes_per_program, block.shape[0])
        for col in block.columns:
            ordered = block[col].sort_values(ascending=False)
            for pos, (gene, score) in enumerate(ordered.head(n_top).items(), start=1):
                top_rows.append({
                    "sample": args.sample, "program": col, "rank": rank,
                    "position": pos, "gene": gene, "loading": float(score),
                })

    basis = pd.concat(basis_blocks, axis=1)

    out = lambda name: os.path.join(args.outdir, f"{args.sample}_{name}")
    basis.to_csv(out("basis_matrix.csv"))
    pd.DataFrame(top_rows).to_csv(out("top_genes.csv"), index=False)
    pd.DataFrame(metrics).to_csv(out("nmf_metrics.csv"), index=False)

    with open(out("nmf_meta.json"), "w") as fh:
        json.dump({
            "sample": args.sample, "engine": "python-sklearn",
            "sklearn_solver": args.solver, "init": args.init,
            "beta_loss": args.beta_loss, "max_iter": args.max_iter, "tol": args.tol,
            "nrun": args.nrun, "seed": args.seed,
            "ranks": ranks, "n_genes": int(basis.shape[0]),
            "n_cells": int(df.shape[1]), "n_programs": int(basis.shape[1]),
            "threads": _THREADS,
            "total_seconds": round(sum(m["seconds"] for m in metrics), 2),
            "any_non_converged": any(not m["converged"] for m in metrics),
        }, fh, indent=2)

    # Reconstruction error vs rank — the python analogue of the R RSS plot.
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt

        m = pd.DataFrame(metrics)
        fig, ax = plt.subplots(figsize=(6, 4))
        ax.plot(m["rank"], m["reconstruction_err"], marker="o")
        ax.set_xlabel("rank")
        ax.set_ylabel("reconstruction error")
        ax.set_title(f"NMF reconstruction error vs rank — {args.sample}")
        ax.set_xticks(ranks)
        ax.grid(alpha=0.3)
        fig.tight_layout()
        fig.savefig(out("nmf_rank.png"), dpi=150)
        plt.close(fig)
    except Exception as exc:  # plotting must never fail the factorisation
        print(f"[{args.sample}] rank plot skipped: {exc}", file=sys.stderr)

    print(f"[{args.sample}] wrote {basis.shape[0]} genes x {basis.shape[1]} programs", flush=True)
    return 0


if __name__ == "__main__":
    sys.exit(main())
