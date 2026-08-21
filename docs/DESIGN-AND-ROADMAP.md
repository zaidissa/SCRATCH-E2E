# SCRATCH-E2E — design, state, and what to do next

*Last updated 2026-08-15. Written to be read start-to-finish by someone picking
the project up cold. Claims marked **verified** were confirmed by running
something; claims marked **unverified** are read from code but not executed.*

---

## 1. What this pipeline is

Eight standalone SCRATCH pipelines — QC, Annotation, CNV, TumorHeterogeneity,
Batchcorrection, CellTrajectory, CellCommunication and TCR — merged into one
Nextflow DSL2 DAG with a single entry point and a cross-stage reporting layer.

The shape is a **linear backbone** that splits into **two arms**:

```
align → qc → cluster → annotation → cnv → stratify ─┬─ malignant      → metaprog
                                                     └─ non-malignant  → batchcorrect
                                                                          → trajectory
                                                                          → cellcomm
                                                                          → tme
```

TCR runs independently of the split. Numbat is an opt-in fourth CNV caller.
Reporting spans everything: one report per stage, then an integration layer.

`docs/pipeline-map.html` is the visual companion to this document — the DAG, the
per-stage process chains, the inter-stage artefact contracts, and entry recipes.

**Every stage is independently enterable** via `--from` / `--to` / `--only`. This
is a design commitment, not a convenience: it is what lets a stage be debugged
without re-running the ones before it, and it is why an intermediate object is a
first-class input everywhere.

---

## 2. Design principles

These are not aspirations. Each was learned from a specific failure, and a change
that violates one is a regression even if it passes the tests.

**1. A report must answer, not display.**
Verdict first, then computed findings, then evidence, then methods. Findings are
computed at render time from the same data the figures were drawn from, so they
cannot drift from the plot beside them. A static caption stating a number is a bug
waiting to happen.

**2. Independent methods must be cross-checked, and disagreement reported.**
The pipeline deliberately runs several annotators, several CNV callers, several
trajectory methods. The failure it exists to prevent is *N methods run, one used,
disagreement never mentioned*. Three separate instances of this class have been
found: the CNV callers, the annotators, and the doublet callers.

**3. Silent degradation is the enemy.**
A zero-match join, a dropped annotator, an empty reference, a caller that produced
nothing — each must fail loudly or be reported as a finding. The worst bugs in
this project have all been runs that reported success while producing nothing.

**4. Never silently change what a method does.**
When a better approach exists, the notebook *computes* the alternative, compares
it to the configured one, and reports the difference as an advisory. The user
decides. Established with the QC MAD-threshold advisory; applied since to PC
selection via knee detection.

**5. A null result is a result.**
A sample with no detectable CNV, an annotator with no matching column, a cohort
too small for a test: emit the empty output with a stated reason and carry on.
One sample's null must never destroy another sample's work.

**6. Everything runs in a container, declared centrally** in
`conf/containers.config`. A process with no container selector works on a laptop
and fails on HPC. Two such processes were found and fixed.

**7. Dependency cost is surfaced before building.**
Compute and image size are real costs on HPC/Cirro/AWS. When a feature pulls in a
heavy dependency, say so and offer the lighter default — propeller by default,
scCODA behind `--run_scoda`.

**8. Scale is a first-class constraint.**
Prefer the on-disk/lazy path even when a two-sample test would not need it. The
BPCells rewrite exists because an in-memory merge hits a hard ceiling at cohort
scale, long before anything looks wrong locally.

---

## 3. Where things stand

| stage | state |
|---|---|
| ALIGN | never executed. Needs CellRanger references. |
| QC → CLUSTER | **rebuilt on BPCells**, stub-verified only — see §4 |
| ANNOTATION | **verified** end to end on real data (24,439 cells) |
| ANNOTATION_CONCORDANCE | **verified**; 5 findings, 5 tables, 3 figures |
| CNV / inferCNV | subsetting bug fixed; **unverified** — not run since the fix |
| CNV / SCEVAN | **REMOVED 2026-08-20** — expression-only like inferCNV so it added no independent evidence, slow enough to dominate a run, a ~5,000-cell ceiling, and a silent-failure history |
| CNV / Numbat | **verified** end to end; see §5 |
| STRATIFY | not run since the CNV work. Concordance computed but not reported. |
| METAPROG, BATCHCORRECT, TRAJECTORY, CELLCOMM, TME | ran in earlier sessions |
| TCR, microbiome | never executed |

---

## 4. The BPCells rewrite (most important recent change)

QC, merge and clustering were rebuilt on an **on-disk BPCells counts store**,
ported from the `SCRATCH-QC` branch `multi-mode`.

**Why.** BPCells uses 64-bit indexing, so the joined matrix is not bound by the
~2.1 billion non-zero ceiling of an in-memory `dgCMatrix` — the ceiling that
caused a `vector::reserve` failure in `JoinLayers` at cohort scale. The old path
(`merge()` + `JoinLayers()`) works on two samples and fails on a real cohort.

**How it works.** `SEURAT_MERGE` writes each QC-approved sample to its own on-disk
matrix, aligns to a common gene set, lazily column-binds, and consolidates into
one store at `data/bpcells_counts`. Every downstream process stages that store at
**the same relative path**, because the object's lazy layers reference it by path.

**The subtle part.** `SCDBLFINDER` detects doublets per sample, then *rewrites the
store to contain only singlets* and re-emits it. The QC subworkflow therefore
emits the singlet store when doublets were removed, and the merge store otherwise.
Getting this backwards feeds doublets to every later stage.

**Status: stub-verified only (56/56).** The stub never invokes a notebook, so
BPCells has not executed once in this pipeline. **Proving it on real data is the
single highest-priority task.**

**Cost.** One visualization was lost: the doublet QC UMAP. It required an
embedding, and the new scdblfinder deliberately computes none — materialising an
embedding over the merged matrix is exactly the memory cost BPCells avoids. Its
quantitative companions (library-size ratio + Wilcoxon) survive.

Rollback: `.preport-backups/` holds the post-figure/pre-BPCells notebooks, since
the port was committed together with other work and cannot be reverted alone.

---

## 5. CNV calling and the Numbat result

**Numbat runs end to end.** On the two-sample HGSOC cohort:

| sample | cells | result |
|---|---|---|
| HRS371755 | 9,164 | **5,560 malignant / 3,604 non-malignant** at `min_LLR=5` |
| HRS371760 | 15,275 | **no calls** — no CNV survived the likelihood filter |

**The threshold hypothesis is dead.** At `min_LLR=2` (confirmed in Numbat's own
log, not merely the CLI flag) HRS371760 still produced nothing, and HRS371755
moved by 56 cells — 0.6%. A 2.5× threshold change is not the lever.

Two explanations remain untested. HGSOC is defined by copy-number instability, so
"no CNV" should prompt suspicion of the method before the biology:

1. **The expression reference.** `numbat_run.R` uses `numbat::ref_hca`, the
   package default. Its own comment notes a cohort-internal reference built from
   this study's immune cells would be preferable.
2. **Genuinely low tumour content** in that sample.

**Caller strategy.** inferCNV + Numbat is the recommended pairing because they
fail *differently* — expression-only versus allele-aware. SCEVAN and CopyKAT are
both expression-only, so pairing either with inferCNV cross-checks less than it
appears to. Numbat stays opt-in: it needs a genome BAM and ~10 GB of references
that a default user will not have, and a default must run out of the box.

---

## 6. What to do next, in order

**0. Prove the BPCells port on real data.**
`--from qc --to cluster` against the CellRanger input. Confirm `data/bpcells_counts`
is produced, and that the store reaching CLUSTER is the *singlet* store. Nothing
else should be built on an unproven merge.

**1. Chase HRS371760's null.**
Threshold is ruled out. Test a cohort-internal expression reference against
`ref_hca`, and check the sample's actual tumour content.

**2. Fix the CNV defaults.**
- ~~`cnv_scevan_max_cells`~~ — moot: SCEVAN is gone.
- ~~`skip_scevan`~~ — superseded by removal. Slow, cell-capped, and a history of
  silent failure.
- Add a **caution finding in STRATIFY when only one caller ran.** With SCEVAN off
  and Numbat opt-in, the default becomes a single caller with no cross-check —
  which violates principle 2 unless the report says so.

**3. Run inferCNV** to prove its subsetting fix on real data.

**4. STRATIFY reporting — the biggest gap.**
The caller-concordance number is computed and then discarded into a `cat()`.
The stage has 0 findings, 0 tables and 2 figures. It decides which cells are
malignant, so it carries more weight than any other report. This is the original
instance of the class that ANNOTATION_CONCORDANCE was built to fix.

**5. DGE visualisation.**
Five notebooks compute differential expression; there are no volcano plots, no
p-value histograms and no pseudobulk PCA anywhere in 42 notebooks. Two caveats
that shape the design:
- `FindAllMarkers` pre-filters on `logfc.threshold` and `min.pct`, so a volcano
  from that output is truncated, and a **p-value histogram would be meaningless** —
  it would show the filter, not the test. An unfiltered run is required.
- Pseudobulk PCA is not honest at n=2 samples. Build it with a refusal guard, the
  way `quality_table_report` declines to compute a 3-MAD outlier test below three
  samples.

**6. The report agent.** `assets/viz/scratch_agent.js` already exists —
deterministic answers computed from data embedded in the page, with an optional
local model (Ollama or OpenAI-compatible, no internet) as a phrasing layer only.
It is wired into stage reports, which have never rendered, so it has never
appeared. Verify it renders, then consider moving it to the per-module notebooks
people actually read.

---

## 7. Open questions

**About the science**
- Is HRS371760 genuinely CNV-free, or is `ref_hca` the wrong reference?
- scType assigns no cell to `Epithelial`, leaving the epithelium in `Unknown`,
  while Azimuth calls 4,764 cells epithelial with high confidence. Which
  annotation should downstream analysis trust?
- In a full run, `annotation.nf` prefers Azimuth, but STRATIFY's reference labels
  are scType's vocabulary — matching 1 of 6 rather than 6 of 6. **Unverified**,
  but it likely explains a previously observed high-ambiguity result.

**About the project — needed to prioritise well**
- What does this feed: a manuscript, a lab tool, a shared resource?
- Expected cohort size — tens of samples, or hundreds?
- Who else runs it? That determines how defensive the defaults must be.
- Any deadline shaping the order of the work above?

---

## 8. Operating notes

- `-profile docker,local` forces `linux/amd64`, so everything runs under emulation
  on Apple silicon. Numbat's pileup over 157 GB of BAM took ~8 hours.
- Nextflow budgets memory from the **host**, but tasks run in the **Docker VM**,
  which is smaller. Cap with `--max_memory` when running heavy stages.
- Nextflow does **not** hash `bin/` script contents. Editing a script there will
  not re-run a cached task; change the process script to bust the hash.
- The stub suite (`-profile test,stub -stub`) proves DAG wiring and output
  contracts. It never invokes a notebook, so a green stub says nothing about
  whether a stage produces correct results.
