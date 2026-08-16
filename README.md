# SCRATCH-E2E

End-to-end single-cell tumour analysis. Merges eight previously standalone
SCRATCH pipelines into one Nextflow DAG with a unified configuration, an explicit malignant/non-malignant
split, a cross-stage integration layer, and a single navigable report.

```
ALIGN ──┬── GEX ──> QC ──> CLUSTER ──> ANNOTATION ──> CNV ──> STRATIFY
        │                                 │                      │
        │                                 │      malignant ──────┼──> METAPROG
        │                                 │                      │
        │                                 │  non-malignant ──────┴──> BATCHCORRECT ──┬──> TRAJECTORY
        │                                 │                                          ├──> CELLCOMM
        │                                 │                                          └──> TME
        ├── VDJ ──────────────────────────┴───────────────────────────────────────────────> TCR
                                                                                             │
                      all stages ────────────────────────────────────────────────────────────┴──> INTEGRATION
```


CNV calling runs inferCNV and SCEVAN by default, with CopyKAT and
[Numbat](docs/NUMBAT.md) available. Numbat is opt-in because it needs the BAM
plus a phasing panel — but it is the only caller that can see copy-neutral LOH,
which the expression-based three are structurally blind to.

**The compartment split is the pivot of the design.** CNV calling runs on the
annotated object; `STRATIFY` uses those calls — with the annotation label as a
tie-break — to separate malignant from non-malignant cells. Metaprograms are a
property of the tumour compartment. Batch correction, trajectory, cell
communication and TME characterisation belong to the microenvironment:
correcting batch across a mixed object would regress out genuine tumour
heterogeneity along with the batch effect.

## Quick start

```bash
# Dry-run the entire DAG (no containers, no computation) — good first check
nextflow run . -profile test,stub -stub

# Real run
nextflow run . -profile docker \
    --project_name MYSTUDY \
    --input_samplesheet  samplesheet.csv \
    --input_exp_table    metadata.csv \
    --outdir             results
```

The front door of the results is
`results/<project_name>/master_report/index.html`.

## Running part of the pipeline

Any stage can be entered directly by supplying an existing Seurat object.

```bash
nextflow run . --from align --to stratify         # through the compartment split
nextflow run . --from annotation --input_seurat_object cluster_object.RDS
nextflow run . --from stratify   --input_seurat_object annotated.RDS
nextflow run . --from metaprog   --input_tumor_object malignant.RDS
nextflow run . --from tme        --input_nonmalignant_object nonmalignant.RDS
nextflow run . --run_metaprog false                # skip one arm
```

`--from`/`--to` select an inclusive range over the backbone
(`align, qc, cluster, annotation, cnv, stratify`); reaching the end also enables
both compartment arms (`metaprog`; `batchcorrect, trajectory, cellcomm, tme`),
`tcr`, and `integration`. Naming a branch in `--from` runs *that* branch only.
Any `--run_<stage>` flag overrides the range.

## Execution environments

| Profile | Use |
|---|---|
| `docker`, `singularity`, `apptainer`, `podman` | container engine |
| `local` | single workstation, bounded by `--max_cpus` / `--max_memory` |
| `lsf`, `slurm` | HPC schedulers, with per-task queue selection by wall time |
| `awsbatch`, `cirro` | cloud |
| `seadragon` | MD Anderson site profile (LSF + Singularity, 28 cpu / 900 GB caps) |
| `test`, `test_full`, `stub`, `debug` | development |

Combine them: `-profile singularity,lsf` or `-profile docker,local`.

## Output layout

Every stage publishes into the same three-way structure, which is what allows
the integration layer to find artefacts by convention:

```
results/
├── pipeline_info/                 timeline, report, trace, DAG
└── <project_name>/
    ├── align/  qc/  cluster/  annotation/  cnv/  stratify/
    ├── metaprog/                       malignant arm
    ├── batchcorrection/  trajectory/  cellcomm/  tme/    non-malignant arm
    ├── tcr/
    │   ├── report/                rendered .html
    │   ├── figures/               .png / .pdf
    │   └── data/                  .RDS / .csv / .h5ad
    ├── integration/
    │   ├── data/<project>_master_object.RDS    every stage's annotation, joined
    │   ├── data/<project>_cell_table.csv       the same, as a tidy table
    │   └── manifest/run_manifest.json          provenance
    └── master_report/index.html   ← start here
```

## What the merge changed

Beyond wiring, the following defects in the source pipelines were fixed. Each is
also commented at the site of the fix.

| Area | Problem | Fix |
|---|---|---|
| `conf/modules.config` | `paramas.input_integration_dimension` typo made `SEURAT_CLUSTER` fail at config evaluation | corrected |
| Parameters | `organism` was `"Human"` in Annotation but `"human"` in CNV/CellComm; `thr_n_features`, `thr_resolution`, `seed`, `project_name`, `outdir` all collided across pipelines | one canonical `organism` with derived `organism_capitalised`/`organism_lower`; all stage params namespaced (`qc_`, `cluster_`, `annot_`, `bc_`, `cnv_`, `metaprog_`, `traj_`, `cc_`, `tcr_`) |
| Containers | absolute `.sif` paths (`/home/sazaidi/…`, `/path/to/…`), local-only tags (`nf-quarto:latest`), `SEURAT_MERGE` with its container commented out | all images centralised in `conf/containers.config`, individually overridable, with a `--container_registry` mirror hook |
| QC | `SCDBLFINDER` ran but its result was discarded (`ch_final_object = ch_merge_object`) — because the notebook writes **two** objects matching one output glob | `--qc_doublet_filter sample\|cluster` selects one, and it propagates |
| Annotation | subworkflow emitted `Channel.empty()`, so nothing chained downstream | emits the annotated object (Azimuth → scType → input fallback) |
| Batch correction | no chainable output; `cpus 4 / memory 16.GB` hardcoded well below need | explicit `seurat_rds` emit; resources from `conf/base.config` |
| Trajectory | `shell:` block used `${var}` (a *shell* variable), so the staged object never reached Quarto; self-referential `cpus`/`memory`; output glob on an absolute path | rewritten as `script:`; resources from config; outputs captured under `data/` |
| Cell communication | NicheNet's inputs were hardcoded output paths built inside the `params` block, ordered by an artificial `_done` tick | real channel dependencies, so `-resume` is correct |
| Metaprograms | an artificial `LEI_TICK` barrier serialised Leiden and NMF preprocessing, which are independent | run concurrently |
| Compartments | the malignant/non-malignant split was implicit and duplicated — `notebook_batch_correctionR.qmd` built its own `obj_non_tumor` by excluding the `Epithelial` label, and the metaprogram stage separately subset to `Epithelial`; two notebooks deciding what a tumour cell is, from a lineage proxy, with no record | explicit `STRATIFY` stage using CNV evidence with the label as tie-break; writes `malignant_status`, an assignment table and caller-concordance figures. `bc_exclude_labels` and `metaprog_subset_*` retuned so cells are not filtered twice |
| TCR | eight optional inputs shared one `assets/NO_FILE`, causing `MASTER_SUMMARY` "input file name collision" whenever more than one was absent | one placeholder per slot |
| Alignment | `GEX+TCR` overwrote the GEX handle with the TCR one, so VDJ contigs were unreachable | separate `gex_outs` / `vdj_outs` / `bam` emissions |
| Stubs | 21 processes had none; several existing ones wrote filenames the real code never produces, and `SCEVAN`'s "stub" ran `quarto render` | full stub coverage; `-stub` now dry-runs all 47 tasks |
| Repo | 1.4 GB of committed macOS/aarch64 R libraries under `modules/local/scdblfinder/renv/` | not vendored |
| Test data | samplesheet pointed at another developer's absolute paths and had a UTF-8 BOM on the header | self-contained fixture in `assets/test_data/` |

## Speed

Three sources of improvement over running the eight pipelines by hand:

1. **Branch parallelism.** The two compartment arms are independent of each
   other, and within the non-malignant arm trajectory, cell communication and
   TME characterisation all run concurrently — rather than as separate
   sequential pipeline invocations.
2. **Staging.** The downstream pipelines all used `stageInMode 'copy'`, which
   duplicates every multi-GB `.RDS` into each task directory. The default is now
   `symlink` (`--stage_mode copy` restores the old behaviour).
3. **NMF engine.** The metaprogram factorisation defaults to
   scanpy/scikit-learn instead of the R `NMF` package — see
   [docs/NMF_ENGINES.md](docs/NMF_ENGINES.md).

## Documentation

- [docs/pipeline-map.html](docs/pipeline-map.html) — the full DAG as a rendered page: stage graph, per-stage process chains, inter-stage artefact contracts, concurrency, gating rules and entry recipes
- [docs/USAGE.md](docs/USAGE.md) — parameter reference and worked examples
- [docs/NMF_ENGINES.md](docs/NMF_ENGINES.md) — the R and Python NMF engines
- [docs/NUMBAT.md](docs/NUMBAT.md) — the opt-in allele-aware CNV caller
- [docs/MIGRATION.md](docs/MIGRATION.md) — mapping from the eight standalone pipelines
- [docs/MICROBIOME.md](docs/MICROBIOME.md) — why the microbiome stage is disabled
