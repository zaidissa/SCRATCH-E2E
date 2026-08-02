# NMF engines for the metaprogram stage

The tumour-heterogeneity stage factorises each sample's expression matrix over a
rank sweep, then assembles recurrent metaprograms from the resulting gene
loadings. Two interchangeable engines are available.

```
--metaprog_nmf_engine python   # default
--metaprog_nmf_engine r        # the original NMF::nmf implementation
--metaprog_nmf_engine both     # run both; PostNMF keeps the two sets side by side
```

## Why the engine can be swapped safely

`Metaprog_PostNMF.qmd` consumes exactly one thing from each fit:

```r
get_basis_by_rank(fit)   # -> NMF::basis(fit), i.e. the gene x program matrix
```

Everything else the R `NMFfit` object carries — mixture coefficients, consensus
matrices, residual traces — is discarded before metaprograms are built. So the
solver is replaceable provided the gene × program loadings come out in the same
layout, which is what the Python engine produces.

Both engines also read the **same input**: the `<sample>_preprocessed.rds`
written by `METAPROG_NMF_PREP` (log-CPM/10, gene-centred, clipped at zero). No
re-preprocessing happens on the Python side. `bin/export_preprocessed_mtx.R`
only changes the container format (`.rds` → MatrixMarket + label files) and
performs no transformation, so a difference in output is attributable to the
solver alone. That is what makes `--metaprog_nmf_engine both` a meaningful
concordance check.

## Pipeline shape

```
METAPROG_NMF_PREP  ──>  <sample>_preprocessed.rds
                            │
        ┌───────────────────┴────────────────────┐
        │ engine = r                             │ engine = python
        ▼                                        ▼
   METAPROG_NMF                           METAPROG_EXPORT_MTX
   (NMF::nmf, snmf/r)                            │
        │                                        ▼
        │                                  METAPROG_NMF_PY
        │                                  (scanpy / scikit-learn)
        ▼                                        ▼
   <sample>_nmf_fit.rds                   <sample>_basis_matrix.csv
        └───────────────────┬────────────────────┘
                            ▼
                      METAPROG_POST
              (auto-detects either input format)
```

## Output contract

The Python engine writes, per sample, into
`<outdir>/<project>/metaprog/data/per_sample_mat/nmf_fit/`:

| File | Contents |
|---|---|
| `<sample>_basis_matrix.csv` | genes × programs, columns `"<rank>_<component>"` |
| `<sample>_top_genes.csv` | top-N genes per program, long format |
| `<sample>_nmf_metrics.csv` | per-rank reconstruction error, iterations, convergence, seconds |
| `<sample>_nmf_meta.json` | solver settings, thread count, totals |
| `<sample>_nmf_rank.png` | reconstruction error vs rank |

`PostNMF` renames the columns to `<sample>.<rank>.<component>` — the naming the
metaprogram clustering relies on to recover a program's sample of origin — and
converts Python's 0-based component index to R's 1-based one.

> A sample id containing a `.` breaks that lookup (the clustering splits program
> names on `.`). `PostNMF` now fails fast with a clear message instead of
> silently mis-assigning programs. This constraint exists in the R path too; it
> was previously unchecked.

## Differences from the reference scanpy/sklearn snippet

The implementation in `bin/scratch_nmf.py` differs in four ways that matter:

1. **`fit_transform` instead of `fit` then `transform`.** Calling `fit(df)`
   followed by `transform(df)` solves for W twice; `fit_transform` returns the W
   already computed during the fit. That is one full NNLS pass saved per rank,
   per restart.

2. **Convergence is controlled and recorded.** sklearn's default
   `max_iter=200` frequently stops early on single-cell matrices and emits a
   `ConvergenceWarning` that is easy to miss. The default here is
   `--metaprog_py_max_iter 1000`, and per-rank convergence is written to
   `<sample>_nmf_metrics.csv` rather than silently accepted.

3. **Ranking is done on an unmutated copy.** The circulated snippet does

   ```python
   W['rank'] = W[i].rank(...)
   W = W.sort_values(by='rank')
   ```

   inside the per-component loop, which adds a column to `W` and reorders its
   rows while the loop is still indexing into it. Top-gene extraction here works
   from an untouched loadings frame.

4. **Restarts approximate R's `nrun`.** The R engine runs `nrun` restarts and
   keeps the best; sklearn has no equivalent. `--metaprog_nrun` runs explicit
   seeded restarts and keeps the lowest reconstruction error. Because `nndsvd*`
   initialisations are deterministic, restarts after the first switch to
   `nndsvdar`, which is randomised — otherwise the extra restarts would recompute
   an identical answer.

## Parameters

| Parameter | Default | Applies to |
|---|---|---|
| `--metaprog_nmf_engine` | `python` | both |
| `--metaprog_rank_lb` / `--metaprog_rank_ub` | `3` / `5` | both |
| `--metaprog_nrun` | `2` | both |
| `--metaprog_hvg_keep` | `7000` | both |
| `--metaprog_max_cells` | `0` (no downsampling) | both |
| `--metaprog_genes_per_program` | `50` | both |
| `--metaprog_nmf_method` | `snmf/r` | R only |
| `--metaprog_py_init` | `nndsvda` | python only |
| `--metaprog_py_solver` | `cd` | python only |
| `--metaprog_py_beta_loss` | `frobenius` | python only |
| `--metaprog_py_max_iter` | `1000` | python only |
| `--metaprog_py_tol` | `1e-4` | python only |

Gene selection (top-variance, `hvg_keep`) is reimplemented in the Python engine
to match `Metaprog_NMF.qmd`'s sum-of-squares variance calculation, so both
engines factorise the same gene set.

## Validation performed

`bin/scratch_nmf.py` was run standalone on a synthetic matrix with three planted
gene programs (800 genes × 400 cells). Reconstruction error collapses from
114.0 at rank 2 to 7.31 at rank 3 and then flattens (6.97 at rank 4, 6.44 at
rank 5) — the elbow lands on the true rank. Output was confirmed to be
genes × programs, non-negative, with `"<rank>_<component>"` columns. The
multi-restart selection changed the retained fit on three of four ranks.

**Not yet validated:** a real side-by-side run of the two engines on your data.
The equivalence argument above is structural (both consume the same matrix and
both feed `NMF::basis`-shaped loadings into unchanged downstream code), but the
numerical agreement between `snmf/r` and sklearn's coordinate descent on real
tumour matrices has not been measured here. Run

```bash
nextflow run . --from metaprog --input_seurat_object <object>.RDS \
               --metaprog_nmf_engine both
```

on one cohort and compare the metaprogram sets before switching the default for
production analyses. The two program sets appear side by side in `PostNMF`
(python programs are suffixed `_py`), and the Jaccard overlap heatmap it already
produces is the natural place to read the concordance off.

## Performance expectation

The speed-up comes from the solver, not from doing less work: sklearn's
coordinate-descent NMF is a compiled implementation that avoids the repeated
dense copies the R `NMF` package makes per restart. Because the stage fans out
one task per sample, `conf/base.config` gives `METAPROG_NMF_PY` a smaller slot
(8 cpu / 24 GB) than `METAPROG_NMF` (12 cpu / 48 GB), so more samples run
concurrently for the same cluster allocation.

BLAS threads are pinned to `task.cpus` before numpy is imported. Without that,
every concurrent per-sample job would try to use every core on the node and the
fan-out would thrash.
