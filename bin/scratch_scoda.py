#!/usr/bin/env python3
"""Compositional differential abundance with scCODA.

Cell-type proportions live on a simplex: they are constrained to sum to 1, so
one population expanding mechanically deflates every other one. Rank tests on
proportions (Wilcoxon, Kruskal-Wallis) treat cell types as independent and will
report that arithmetic deflation as a real depletion.

scCODA models the per-sample count vector as Dirichlet-multinomial and applies
spike-and-slab selection, so "credible" means the effect survives the
compositional constraint rather than merely moving a proportion.

The reference cell type matters and is not a formality: the model has one fewer
degree of freedom than it has cell types, so effects are only identifiable
RELATIVE to a reference. "automatic" picks a type that is present everywhere and
has low dispersion. Every effect is a change relative to that type -- there is no
such thing as an absolute scCODA effect.

Reads a counts matrix (samples x cell types) and per-sample covariates; writes a
tidy CSV of effects. Exit code 2 means "scCODA unavailable", which the caller is
expected to treat as a degradation, not a failure.
"""

import argparse
import sys

import pandas as pd


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--counts", required=True,
                    help="CSV: first column sample id, remaining columns cell-type COUNTS")
    ap.add_argument("--covariates", required=True,
                    help="CSV: first column sample id, remaining columns covariates")
    ap.add_argument("--formula", default="condition",
                    help="patsy formula over covariate columns [default: %(default)s]")
    ap.add_argument("--reference", default="automatic",
                    help="reference cell type, or 'automatic' [default: %(default)s]")
    ap.add_argument("--fdr", type=float, default=0.05,
                    help="target FDR for the spike-and-slab threshold [default: %(default)s]")
    ap.add_argument("--num-results", type=int, default=20000,
                    help="HMC samples [default: %(default)s]")
    ap.add_argument("--seed", type=int, default=1234)
    ap.add_argument("--out", required=True, help="output CSV of effects")
    ap.add_argument("--out-summary", default=None, help="optional text summary")
    a = ap.parse_args()

    try:
        import numpy as np
        from sccoda.util import cell_composition_data as scdat
        from sccoda.util import comp_ana as scmod
    except ImportError as e:                       # noqa: PERF203 - explicit path
        sys.stderr.write(
            f"scCODA unavailable ({e}). Install with `pip install sccoda` "
            "(pulls tensorflow + tensorflow-probability).\n")
        return 2

    counts = pd.read_csv(a.counts, index_col=0)
    covs = pd.read_csv(a.covariates, index_col=0)

    common = counts.index.intersection(covs.index)
    if len(common) < len(counts.index):
        sys.stderr.write(f"warning: {len(counts.index) - len(common)} sample(s) "
                         "in counts have no covariates and were dropped\n")
    counts, covs = counts.loc[common], covs.loc[common]

    # scCODA is a SAMPLE-level model: each row is a replicate. Too few replicates
    # per group and the posterior is prior-dominated, so refuse rather than emit
    # confident-looking noise.
    terms = [t for t in covs.columns if t in a.formula]
    for t in terms:
        vc = covs[t].value_counts()
        if (vc < 2).any():
            sys.stderr.write(
                f"error: covariate '{t}' has level(s) with <2 samples "
                f"({vc[vc < 2].to_dict()}). scCODA needs replicates per group.\n")
            return 3
    if len(counts) < 4:
        sys.stderr.write(f"error: only {len(counts)} samples; scCODA needs more.\n")
        return 3

    # Counts must be integers -- passing proportions here silently destroys the
    # Dirichlet-multinomial's notion of sampling depth.
    if not np.allclose(counts.values, counts.values.round()):
        sys.stderr.write("error: --counts contains non-integers; pass COUNTS, not proportions.\n")
        return 3

    df = pd.concat([counts, covs], axis=1)
    data = scdat.from_pandas(df, covariate_columns=list(covs.columns))

    model = scmod.CompositionalAnalysis(
        data, formula=a.formula, reference_cell_type=a.reference)

    # scCODA 0.1.9's sample_hmc() takes no seed argument, so seed the TF global
    # RNG instead -- without this the HMC chain is not reproducible between runs
    # and the credible set can change.
    try:
        import tensorflow as tf
        tf.random.set_seed(a.seed)
    except Exception:                              # noqa: BLE001 - seeding is best-effort
        pass

    res = model.sample_hmc(num_results=a.num_results)
    res.set_fdr(est_fdr=a.fdr)

    eff = res.effect_df.reset_index()
    eff.columns = [str(c) for c in eff.columns]
    ren = {"Covariate": "covariate", "Cell Type": "cell_type",
           "Final Parameter": "effect", "Inclusion probability": "inclusion_prob",
           "log2-fold change": "log2fc", "Expected Sample": "expected_sample",
           "HDI 3%": "hdi_3", "HDI 97%": "hdi_97", "SD": "sd"}
    eff = eff.rename(columns={k: v for k, v in ren.items() if k in eff.columns})

    # scCODA's verdict is binary by construction: a non-zero Final Parameter is a
    # credible effect at the chosen FDR. There is no p-value to threshold.
    eff["credible"] = eff["effect"].ne(0)

    # Resolve what "automatic" actually chose. Recording the literal string
    # "automatic" would be worse than useless: every effect is relative to the
    # reference, so a reader who cannot see which type it was cannot interpret a
    # single number in this table.
    ref_name = a.reference
    for obj in (model, res):
        r = getattr(obj, "reference_cell_type", None)
        if r is None:
            continue
        if isinstance(r, (int, np.integer)):
            types = getattr(obj, "cell_types", None)
            if types is None:
                types = list(counts.columns)
            ref_name = str(types[int(r)])
        else:
            ref_name = str(r)
        break
    if ref_name == "automatic":
        sys.stderr.write("warning: could not resolve the automatic reference cell type\n")
    eff["reference_cell_type"] = ref_name
    eff["est_fdr"] = a.fdr
    eff.sort_values(["covariate", "inclusion_prob"], ascending=[True, False]).to_csv(
        a.out, index=False)

    n = int(eff["credible"].sum())
    sys.stderr.write(f"scCODA: {n} credible effect(s) at FDR {a.fdr} "
                     f"across {len(counts)} samples x {counts.shape[1]} cell types\n")

    if a.out_summary:
        with open(a.out_summary, "w") as fh:
            fh.write(f"reference cell type: {eff['reference_cell_type'].iloc[0]}\n")
            fh.write(f"target FDR: {a.fdr}\nsamples: {len(counts)}\n\n")
            fh.write(eff.to_string(index=False))
    return 0


if __name__ == "__main__":
    sys.exit(main())
