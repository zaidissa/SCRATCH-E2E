# Migrating from the eight standalone pipelines

## Pipeline → stage

| Standalone pipeline | Entry point | Stage in SCRATCH-E2E |
|---|---|---|
| SCRATCH-QC-multi-mode | `scratch_align_entry.nf` | `align` |
| SCRATCH-QC-multi-mode | `scratch_qc_entry.nf` | `qc` |
| SCRATCH-QC-multi-mode | `scratch_cluster_entry.nf` | `cluster` |
| SCRATCH-QC-multi-mode | `scratch_microbiome_entry.nf` | `microbiome` (not implemented — see MICROBIOME.md) |
| SCRATCH-Annotation | `main.nf` | `annotation` |
| SCRATCH-Batchcorrection_Evaluation | `main.nf` | `batchcorrect` |
| SCRATCH-CNV | `main.nf` | `cnv` (now on the backbone) |
| — | — | `stratify` (new) |
| SCRATCH-TumorHeterogeneity | `main.nf` | `metaprog` |
| SCRATCH-CellTrajectory | `main.nf` | `trajectory` |
| SCRATCH_CellCommunication | `main.nf` | `cellcomm` |
| SCRATCH-TCR | `main.nf` | `tcr` |
| — | — | `tme` (new) |
| — | — | `integration` (new) |

## Parameter renames

Stage-specific parameters are namespaced because the same name meant different
things in different pipelines. Notebook-facing `-P` names are unchanged;
`conf/modules.config` maps between them.

### Collided across pipelines — renamed

| Old | New | Was in |
|---|---|---|
| `organism` (`"Human"`) | `organism` = `human`; stages read `organism_capitalised` | Annotation |
| `organism` (`"human"`) | `organism` = `human`; stages read `organism_lower` | CNV, CellComm |
| `thr_n_features` | `cluster_thr_n_features` | QC-multi-mode (cluster) |
| `thr_n_features` | `annot_thr_n_features` | Annotation |
| `thr_n_dimensions` | `cluster_thr_n_dimensions` | QC-multi-mode |
| `thr_resolution` | `cluster_thr_resolution` | QC-multi-mode |
| `thr_resolution` | `annot_thr_resolution` | Annotation |
| `thr_npc` | `annot_thr_npc` | Annotation |
| `thr_proportion` | `cluster_thr_proportion` | QC-multi-mode |
| `input_integration_dimension` | `cluster_integration_dimension` | QC-multi-mode |
| `input_group_plot` | `cluster_group_plot` | QC-multi-mode |
| `thr_estimate_n_cells` etc. | `qc_thr_*` | QC-multi-mode |
| `annotation_db` | `annot_db` | Annotation |
| `input_integration_method` | `bc_integration_method` | BatchCorrection |
| `input_target_variables` | `bc_target_variables` | BatchCorrection |
| `input_batch_step` | `bc_batch_step` | BatchCorrection |
| `exclude_labels` | `bc_exclude_labels` | BatchCorrection |
| `celltype_column` | `bc_celltype_column` | BatchCorrection |
| `label_candidates` | `bc_label_candidates` | BatchCorrection |
| `n_hvgs` / `n_pcs` | `bc_n_hvgs` / `bc_n_pcs` | BatchCorrection |
| `input_analysis_mode` | `cnv_analysis_mode` | CNV |
| `ct_species` | `traj_species` | CellTrajectory |
| `subset_col`, `subset_values`, `hvg_keep`, `rank_lb`, … | `metaprog_*` | TumorHeterogeneity |
| `celltype_col`, `sample_col`, `assay`, `slot`, `min_cells_per_group`, `liana_*`, `receiver_mode`, … | `cc_*` | CellCommunication |
| `expected_cells`, `total_droplets`, `fpr`, `epochs` | `cellbender_*` | QC-multi-mode |
| `skip_BatchCorr` | `skip_batchcorrect` | BatchCorrection |
| `skip_celltraj` | `skip_trajectory` | CellTrajectory |
| `n_threads` / `n_memory` | removed — processes use `task.cpus` / `task.memory` | BatchCorrection, CellTrajectory |

### Kept as-is

TCR's tool parameters were already uniquely prefixed (`gliph_`, `tcri_`,
`tcrdist_`, `giana_`, `conga_`, `repertoire_`, `consensus_`, `vdj_`) and merge
without collision — 212 of them live in `conf/tcr.config` under their original
names, so the ten TCR notebooks run unmodified. TCR's own `run_*` module toggles
were renamed to `tcr_run_*` to avoid clashing with the stage selectors.

### New

| Parameter | Purpose |
|---|---|
| `--from` / `--to` | stage range selection |
| `--run_<stage>` | per-stage override |
| `--metaprog_nmf_engine` | `python`, `r`, or `both` |
| `--metaprog_py_*` | Python NMF solver settings |
| `--container_*`, `--container_registry` | per-image and mirror overrides |
| `--stage_mode` | task input staging (`symlink` default, was `copy`) |
| `--publish_mode` | publish mode for all stages |
| `--integration_*`, `--report_*` | integration layer and master report |

## Topology change

The largest difference from the standalone pipelines. Previously every
downstream analysis ran on one batch-corrected object containing both tumour and
microenvironment. Now:

```
ANNOTATION -> CNV -> STRATIFY -+- malignant     -> METAPROG
                               +- non-malignant -> BATCHCORRECT -> TRAJECTORY
                                                                   CELLCOMM
                                                                   TME
```

* **CNV moved onto the backbone.** It was a terminal branch; STRATIFY consumes
  its per-cell calls, so it must run first.
* **Batch correction moved off the backbone** onto the non-malignant arm.
  Correcting batch across a mixed object regresses out genuine tumour
  heterogeneity along with the batch effect.
* **The split became explicit.** It was previously happening twice, implicitly
  and inconsistently: `notebook_batch_correctionR.qmd` built its own
  `obj_non_tumor` by excluding the `Epithelial` label, and the metaprogram stage
  separately subset *to* `Epithelial`. Two notebooks deciding what a tumour cell
  is, from a lineage proxy rather than the CNV evidence, with no record of the
  decision. `bc_exclude_labels` now defaults to empty and
  `metaprog_subset_col`/`_values` point at `malignant_status`, so cells are not
  filtered twice.
* **TME characterisation is new.** No source pipeline had an equivalent.

## Behavioural changes to be aware of

1. **Doublet filtering now takes effect.** The QC stage previously discarded
   `SCDBLFINDER`'s output and passed the pre-filter merged object downstream.
   It now propagates the doublet-filtered object. Cell counts in every
   downstream stage will be lower than in the old runs. Set
   `--skip_scdblfinder true` to reproduce the old behaviour.

2. **`stageInMode` is `symlink`, not `copy`.** No process mutates its inputs in
   place, so this is safe and much faster, but if a site's shared filesystem
   misbehaves with symlinks use `--stage_mode copy`.

3. **The metaprogram NMF engine defaults to Python.** Use
   `--metaprog_nmf_engine r` to keep the original solver. See
   [NMF_ENGINES.md](NMF_ENGINES.md) — validate concordance on one cohort before
   switching production analyses.

4. **Cells can now be dropped as `ambiguous`.** When CNV callers disagree,
   STRATIFY assigns `ambiguous` and, by default, excludes those cells from both
   compartments so neither is contaminated. They remain in
   `*_malignancy_assignment.csv` and in the master object. Use
   `--stratify_ambiguous malignant|non_malignant` to route them instead.

5. **Output paths changed.** Each pipeline used its own layout
   (`SCRATCH-LeidenMetaprogram/`, `SCRATCH_NMF-prep/`, `data/<task.process>/`,
   `<project>/cellcommunication/report/`). Everything is now
   `<outdir>/<project>/<stage>/{report,figures,data}/`.

6. **Resource ceilings are unified.** Defaults are `--max_cpus 16`,
   `--max_memory 128.GB`. The standalone pipelines ranged from 10 to 32 cpus and
   60 GB to 650 GB. Use `-profile seadragon` or raise the ceilings for large
   cohorts; per-process requests in `conf/base.config` scale on retry.

## Equivalence check

To confirm a stage behaves as it did standalone, run it in isolation from the
same input object and diff the outputs:

```bash
nextflow run . --from cnv --to cnv \
    --input_seurat_object <the object you fed SCRATCH-CNV> \
    --project_name equivalence_check \
    --outdir /tmp/e2e_cnv
```
