# Metadata tables for `--input_exp_table`

`--input_exp_table` is the per-sample **metadata** table QC's merge step reads.
It is not the alignment samplesheet. The contract, from
`modules/local/seurat/merge/notebook_seurat_merge.qmd`:

- A **`sample`** column is required; the notebook stops without it.
- The join is `merge(meta, by.x = "orig.ident", by.y = "sample")` — **an inner
  join, with no `all.x`**. Every library present in the object must appear here,
  or `AddMetaData` fails with *"Cannot add more or less meta data without cell
  names"*.
- `orig.ident` is the **cellranger output directory name**. cellranger runs with
  `--id="${sample}"`, so these match end to end — but they are compared as
  strings, and strings that differ join nothing.
- `sampleIndex`, `fastq_1..4` and `modality` are dropped and `distinct()` is
  applied, so an alignment samplesheet also works; the extra columns just go.

## Why these files exist at all

Starting at `--from align`, `preprocess.py` builds the samplesheet and
`build_samplesheet()` asserts `sample`, `modality` and `patient_id`.

Starting at **`--from qc`** or later, `--input_exp_table` points at a *previous
run's* `pipeline_info/samplesheet.valid.csv` and **nothing checks it**. If that
run was launched before `patient_id` was carried through the pivot, the column is
simply absent, the merged object never gains it, and every stage keyed on it
falls back to the library column in silence.

That is not cosmetic when Numbat is on: `numbat_prep.R` groups by `patient_id`,
so without it Numbat runs **per library**, and
`numbat:::mark_tumor_lineage()` has no null — it always labels its heaviest clade
`tumor`, so a library with no malignant cells still returns neoplastic calls.
`describe_exp_table()` now prints a loud block for exactly this combination.

## `gbm1_dfci3_exp_table.csv`

17 GEX libraries, **one patient**, six timepoints.

| | |
|---|---|
| `patient_id` | `DFCI3` for all 17 |
| `timepoint` | `S1`–`S6` |
| libraries per timepoint | S1:3, S2:1, S3:3, S4:4, S5:3, S6:3 |

**Derivation.** The `sample` strings are the 17 distinct values logged from
`datasets/f99458a6-.../pipeline_info/samplesheet.valid.csv` by the 2026-10-22 and
2026-10-23 runs. `patient_id` is the `DFCI<n>` token and `timepoint` is the first
`S<n>` token — the naming scheme confirmed by the user, and consistent with the
`grouping` column of the DFCI6/DFCI7 cohort
(`GBM1.DFCI6.d123.R.TS.L3_C2.v1`), where `DFCI<n>` tracks the individual. The
generator asserts that the repeated S-number in the `BTC-*-S<n>-72-S<n>-*` family
agrees with itself rather than taking the first match on trust.

**Verified** against the merge notebook's own join on a 17-identity object:
inner join keeps 680 of 680 cells, `patient_id` and `timepoint` are non-NA on
every cell, and `numbat_prep.R` then reports `Libraries: 17 | Groups: 1`. Without
the table it reports `Groups: 17` and the per-library warning.

### Two things this table does NOT say

- **No `sample_type`.** The original Cirro FASTQ samplesheet has one — the April
  logs show `['sample','modality','sample_type','patient_id','readgroup','grouping']`
  — and in the DFCI6/DFCI7 cohort every value was `TS` (tumour sample). For
  DFCI3 that column never reached QC and is not recoverable from any log here, so
  **which libraries are PBMC is unknown**. It is left out rather than guessed.
  Add it from the FASTQ dataset if you want the tumour-vs-PBMC split labelled.
- **`S4-L1` and `S4-L1Deeper` are treated as two libraries.** If "Deeper" is a
  deeper resequencing of the same GEM well, the same physical cell carries the
  same barcode in both, and namespacing (`<library>_<barcode>`) then makes it two
  cells. Nothing here merges them. Decide whether that is what you want before
  reading cell counts off the pooled run.

## Resource consequence

All 17 libraries in one group is the worst case for `NUMBAT_PILEUP`, because
Numbat piles up samples **sequentially** — `pileup_and_phase.R` writes one
`cellsnp-lite` command per sample into `run_pileup.sh` and runs the script — so
wall time is the sum, not the max. `conf/base.config` raises PILEUP to 120 h and
`NUMBAT_RUN` to 192 GB / 96 h. Both numbers are **estimates from the 2026-10
failure pattern, not measurements**; replace them with the peak RSS and runtimes
from that run's `report.html`.
