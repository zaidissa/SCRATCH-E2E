# Usage

## Inputs

Only the inputs required by your chosen starting stage need to be set.

| Starting at | Required |
|---|---|
| `align` | `--input_samplesheet`, `--input_exp_table`, `--genome`, `--modality` |
| `qc` | `--input_gex_matrices_path`, `--input_exp_table` |
| `cluster`, `annotation`, `cnv`, `stratify` | `--input_seurat_object` |
| `metaprog` | `--input_tumor_object` (or `--input_seurat_object`) |
| `batchcorrect`, `trajectory`, `cellcomm`, `tme` | `--input_nonmalignant_object` (or `--input_seurat_object`) |
| `tcr` | `--input_seurat_object` (annotated), optionally `--input_tcr_sample_sheet` |

Entering `trajectory`, `cellcomm` or `tme` without `--input_nonmalignant_object`
automatically pulls in `batchcorrect`, since those three read its output.

### Samplesheets

Three alignment modes, selected by `--multi` / `--demux`:

**Standard** (both false) — pre-demultiplexed FASTQs, one row per sample:

```csv
sample,fastq_1,fastq_2,modality
S1,/data/S1_R1.fastq.gz,/data/S1_R2.fastq.gz,GEX
S2,/data/S2_R1.fastq.gz,/data/S2_R2.fastq.gz,GEX
T1,/data/T1_R1.fastq.gz,/data/T1_R2.fastq.gz,TCR
```

**Demux** (`--demux true`) — pooled FASTQs, barcodes auto-detected from I1:

```csv
run_id,r1,r2,i1,i2,modality
RUN1,/data/pool_R1.fastq.gz,/data/pool_R2.fastq.gz,/data/pool_I1.fastq.gz,,GEX
```

**Multi** (`--multi true`) — hashtag-multiplexed, CellRanger multi demultiplexes:

```csv
sample_id,multi_config
GBM_DFCI1,/data/cellranger_multi_config_DFCI1.csv
```

Do not write the samplesheet with a UTF-8 BOM — it corrupts the header name and
`splitCsv` will not find the first column.

### Metadata table

Joined onto the merged object during QC:

```csv
sample,patient_id,timepoint,condition,batch
S1,pat01,pre,Control,b1
S2,pat02,post,Case,b2
```

The column names here are what you point `--cc_sample_col`,
`--cc_condition_col` and `--bc_target_variables` at.

## Worked examples

```bash
# Whole pipeline, HPC
nextflow run . -profile singularity,lsf \
    --project_name GBM_COHORT \
    --input_samplesheet sheet.csv --input_exp_table meta.csv \
    --genome GRCh38 --organism human --modality GEX+TCR \
    --outdir /scratch/gbm/results

# Through the compartment split, then inspect before continuing
nextflow run . -profile docker --from align --to stratify \
    --input_samplesheet sheet.csv --input_exp_table meta.csv

# Re-split with a different rule, without redoing CNV
nextflow run . -profile docker --from stratify --to stratify \
    --input_seurat_object annotated.RDS --stratify_method cnv --stratify_min_callers 2

# Only the microenvironment arm, from an existing non-malignant object
nextflow run . -profile docker --from tme \
    --input_nonmalignant_object nonmalignant.RDS

# Continue from the QC object
nextflow run . -profile docker --from cluster \
    --input_seurat_object results/GBM_COHORT/qc/data/GBM_COHORT_qc_dbl_singlet_object.RDS

# Just the tumour-heterogeneity stage, comparing both NMF engines
nextflow run . -profile docker --from metaprog --to metaprog \
    --input_seurat_object corrected.RDS \
    --metaprog_nmf_engine both

# Everything except CNV and TCR
nextflow run . -profile docker --run_cnv false --run_tcr false \
    --input_samplesheet sheet.csv --input_exp_table meta.csv

# Dry-run the full DAG without containers
nextflow run . -profile test,stub -stub
```

## Stage selection semantics

`--from` / `--to` describe an inclusive range over the backbone:

```
align → qc → cluster → annotation → cnv → stratify
```

* `--to end` (the default) additionally enables both compartment arms
  (`metaprog`; `batchcorrect → trajectory, cellcomm, tme`), `tcr` and
  `integration`.
* `--to <backbone stage>` stops there; no arms, no integration.
* `--from <branch name>` skips the backbone and runs **only that branch**
  (plus its prerequisites and `integration`). `--from trajectory` does not
  also launch metaprograms and TCR.
* Any `--run_<stage> true|false` overrides the range for that stage.

Note that params set inside a profile take precedence over the command line in
Nextflow, so a custom profile should not pin `from`/`to` unless you intend them
to be unoverridable.

## Selected parameters

Full defaults are in `nextflow.config`; TCR's 212 tool parameters are in
`conf/tcr.config`. Print the resolved set with:

```bash
nextflow config -profile docker .
```

### Global

| Parameter | Default | Notes |
|---|---|---|
| `--project_name` | `project` | names the output directory and most artefacts |
| `--outdir` | `./results` | |
| `--organism` | `human` | canonical; stages derive `Human`/`human` as needed |
| `--genome` | `GRCh38` | `GRCh38` or `GRCm39` |
| `--seed` | `1234` | applied across stages |
| `--max_cpus` / `--max_memory` / `--max_time` | `16` / `128.GB` / `240.h` | hard ceiling via `resourceLimits` |
| `--stage_mode` | `symlink` | `copy` if your filesystem needs it |
| `--publish_mode` | `copy` | |

### QC

| Parameter | Default |
|---|---|
| `--skip_cellbender` | `true` |
| `--qc_thr_n_feature_rna_min` / `_max` | `300` / `7500` |
| `--qc_thr_percent_mito` | `25` |
| `--qc_thr_median_genes_per_cell` | `900` |
| `--skip_scdblfinder` | `false` |

### Clustering and annotation

| Parameter | Default |
|---|---|
| `--cluster_thr_n_features` | `2000` |
| `--cluster_thr_n_dimensions` | `100` |
| `--cluster_thr_resolution` | `0.5` |
| `--cluster_group_plot` | `patient_id;timepoint` |
| `--tumor_type` | `ovarian` (or `glioblastoma`) — presets the four settings below that must agree |
| `--annot_db` | per tumour type: `assets/cell_markers_database.csv` / `assets/gbm_cell_markers_database.csv` |
| `--annot_state_exclude` | per tumour type: major classes with no state-level pass |
| `--skip_celltypist` / `--skip_sctype` / `--skip_azimuth` | `false` |

Azimuth only runs when `--input_reference_object` is set to a real file;
scType's marker database in this build is human-only.

### CNV and the compartment split

| Parameter | Default | Notes |
|---|---|---|
| `--skip_infercnv` | `false` | |
| `--run_numbat` | `true` (auto-skips with a warning when no BAM is available) | |
| `--cnv_caller_priority` | `numbat` (or `infercnv`) — which caller wins where they disagree | |
| `--skip_cnv_concordance` | `false` — the inferCNV vs Numbat comparison report | |
| `--skip_copykat` | `true` | |
| `--run_numbat` | `false` | allele-aware; needs BAM + phasing panel, see [NUMBAT.md](NUMBAT.md) |
| `--stratify_gate_by_lineage` | `true` | immune/stromal lineages are non-malignant regardless of CNV score |
| `--stratify_method` | `consensus` | `consensus` / `cnv` / `annotation` |
| `--stratify_min_callers` | `1` | callers that must agree before CNV decides |
| `--stratify_tumor_labels` | per tumour type: `Epithelial` / `Malignant` | labels treated as the tumour lineage |
| `--stratify_reference_labels` | per tumour type: immune + stromal / TAMs, T, NK, oligodendrocyte, endothelial | the normal baseline for inferCNV, stratification and CNV concordance |
| `--stratify_infercnv_sd` | `3` | SDs above the reference-lineage mean |
| `--stratify_ambiguous` | `exclude` | `exclude` / `malignant` / `non_malignant` |
| `--stratify_min_cells` | `50` | below this, that arm is skipped |

`STRATIFY` writes `malignant_status` (`malignant` / `non_malignant` /
`ambiguous`) onto the object, plus an assignment table and caller-concordance
figures. Ambiguous cells are excluded from **both** compartments by default so
neither is contaminated; they remain in the assignment table and the master
object.

### Non-malignant arm

| Parameter | Default |
|---|---|
| `--bc_integration_method` | `all` (`cca`/`rpca`/`harmony`/`mnn`) |
| `--bc_target_variables` | `batch` |
| `--bc_exclude_labels` | `''` — empty, because STRATIFY already removed malignant cells |
| `--traj_species` | `human` |
| `--skip_tme` | `false` |
| `--tme_celltype_col` | `azimuth_labels` |
| `--tme_top_n_celltypes` | `12` |
| `--tme_min_cells_per_group` | `30` |

### Cell communication

| Parameter | Default |
|---|---|
| `--cc_celltype_col` | `celltype` |
| `--cc_sample_col` | `sample_id` |
| `--cc_condition_col` | `condition` |
| `--cc_condition_ref` / `--cc_condition_test` | `Control` / `Case` |
| `--cc_min_cells_per_group` | `30` |

### Containers

| Parameter | Purpose |
|---|---|
| `--container_registry` | prefix every SCRATCH image with a private mirror |
| `--container_qc`, `--container_cnv`, … | override one image |

Third-party images (CellBender, CellRanger) are never re-prefixed by
`--container_registry`.

## Troubleshooting

**"No samples matched QC expectations."** Every sample failed the `qc_thr_*`
gates. Open `<outdir>/<project>/qc/report/` for the per-sample summary before
loosening thresholds.

**"Starting from 'X' requires --input_seurat_object."** Stages after QC need an
object to start from.

**"input file name collision"** in a TCR process — an optional input is being
supplied twice under the same basename. The eight VDJ-QC optional slots each
have their own placeholder in `assets/NO_FILE_*`.

**Metaprogram stage fails on sample ids containing `.`** — program names encode
`sample.rank.component`, so a dot in a sample id is ambiguous. Rename the
samples.

**Out-of-memory kills.** `conf/base.config` retries with proportionally more
memory (up to `maxRetries = 2`) for exit codes 104/130/134/137/139/140/143/247.
For inferCNV or CellChat on very large objects, raise `--max_memory` or use a
site profile.

**A stage's report is missing but the run succeeded.** `RUN_MANIFEST` and
`MASTER_REPORT` use `errorStrategy 'ignore'` deliberately: a report formatting
failure must not discard a completed multi-hour analysis. Check
`<outdir>/pipeline_info/` for the task trace.
