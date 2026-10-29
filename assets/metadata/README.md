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
| `sample_type` | `TP` (tumour) ×12, `NB` (normal blood / PBMC) ×5 |
| libraries per timepoint | S1:3, S2:1, S3:3, S4:4, S5:3, S6:3 |

`sample_type` is the split the design actually has: only S1, S3, S5 and S6 have
both compartments. **S2 is blood only** (its one library is the PBMC draw) and
**S4 is tumour only** (all four libraries, including both `Deeper` ones). The
`grouping` strings carry the same fact — `.S.PB.` for every `NB` library and
`.S.TS.` for every `TP` one — and the day offsets in them (`d0, d16, d30, d58,
d91, d121`) map S1–S6 in order, matching `nominal_timepoint` 1–6 on this
subject's bulkDNA biospecimens.

**Derivation.** The `sample` strings are the 17 distinct values logged from
`datasets/f99458a6-.../pipeline_info/samplesheet.valid.csv` by the 2026-10-22 and
2026-10-23 runs. `patient_id` is the `DFCI<n>` token and `timepoint` is the first
`S<n>` token — the naming scheme confirmed by the user, and consistent with the
`grouping` column of the DFCI6/DFCI7 cohort
(`GBM1.DFCI6.d123.R.TS.L3_C2.v1`), where `DFCI<n>` tracks the individual. The
generator asserts that the repeated S-number in the `BTC-*-S<n>-72-S<n>-*` family
agrees with itself rather than taking the first match on trust.

`sample_type` is copied verbatim from the **FASTQ-level samplesheet**, which is
the only place it exists — the QC-stage samplesheet has the column blank for
every library, and the BTC file registry's `biospecimen.sample_type` is the
constant `CC` on all 3,038 rows, so neither can supply it. The source vocabulary
(`TP`/`NB`) is kept unchanged rather than renamed to something friendlier, so a
run that is handed the real Cirro samplesheet instead of this file behaves
identically.

**Four independent FASTQ samplesheets agree on all 17 assignments**, with no
disagreement and nothing unmatched; the `-LIB` suffix that the QC stage uses is
the only difference between their `sample` values and these. The 12 `TP`
libraries are corroborated a second time, from a different source, by a prior
per-cell export of this cohort (`gex_sample` → `sample_type`, 12/12 agree, 0
disagree). The 5 `NB` libraries are absent from that export, which is consistent
with it being a tumour-and-CSF object.

**Verified** against the merge notebook's own join on a 17-identity object,
10/10 assertions in `scratch-cnv`: the inner join keeps 680 of 680 cells;
`patient_id`, `timepoint` and `sample_type` are non-NA on every cell;
`sample_type` is constant within every library and splits 12 `TP` / 5 `NB`; and
the notebook's own design-factor detection picks up `timepoint` and
`sample_type`. `numbat_prep.R` then reports `Libraries: 17 | Groups: 1`. Without
the table it reports `Groups: 17` and the per-library warning.

### What `sample_type` does and does not reach

Nothing in the pipeline reads a column called `sample_type`; the merge is
generic, so it lands in `@meta.data` and is detected as a design factor, and
that is all. In particular **the Numbat per-library table does not carry it** —
`numbat_run.R` keys on `library` and `group`, so the tumour-vs-PBMC malignant
fraction still has to be obtained by joining this file to
`cnv/data/<individual>_numbat_by_library.tsv` by hand. Wiring it through `numbat_prep.R`'s
cellmap would make that automatic.

Also note that the 5 `NB` libraries stay in the same Numbat group as the 12
tumour ones, which is deliberate: internal normals from the same experiment are
the reference the maintainer recommends over an external dataset (see
`docs/NUMBAT.md`). They are reference cells, not held-out cells, so they are
still scored — which is exactly why the lineage-gated percentage, not the raw
one, is the number to read for them.

### One thing this table does NOT say

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
