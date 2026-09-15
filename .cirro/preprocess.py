#!/usr/bin/env python3
"""
Cirro preprocessing for SCRATCH-E2E.

Three jobs, in order of how badly their absence bites:

1. Build the cellranger samplesheet when the run starts at `align`. Cirro
   annotates a FASTQ dataset one row per file; cellranger wants one row per
   sample with R1 and R2 in columns. Without this, `--from align` fails on a
   samplesheet that was never written.

2. Translate `alignment_mode` into the two booleans main.nf actually reads, and
   drop entry-point paths left blank in the form. The form offers three entry
   objects because the right one depends on where you start, so at most one is
   ever filled; passing the others as empty strings defeats the `?:` chain in
   main.nf's entry-point resolution and produces a file-not-found for a path
   nobody typed.

3. Log enough that a failed launch can be diagnosed without re-running it.
"""

import os

import pandas as pd
from cirro.helpers.preprocess_dataset import PreprocessDataset

# Blank in a form field means "not supplied", not "supplied as empty".
BLANK_MEANS_UNSET = [
    "input_seurat_object",
    "input_tumor_object",
    "input_nonmalignant_object",
    "input_reference_object",
    "input_bam_path",
    "input_vdj_contigs",
    # Blank means "use the --tumor_type preset". nextflow.config resolves these
    # four from tumor_type, and any value that reaches Nextflow — including "" —
    # overrides the preset: "" for annot_db reaches Channel.fromPath(...,
    # checkIfExists: true) and fails on a path nobody typed, and "" for
    # annot_state_exclude would silently run a state pass for every class.
    "annot_db",
    "annot_state_exclude",
]


def build_samplesheet(ds: PreprocessDataset) -> pd.DataFrame:
    """One row per sample, R1/R2 in columns — cellranger's layout, not Cirro's.

    Mirrors SCRATCH-QC-multi-mode/.cirro/align/preprocess.py, which is the
    version already known to work against a real Cirro FASTQ dataset.
    """
    REQUIRED = ["sample", "modality", "patient_id"]
    OPTIONAL = ["timepoint", "batch"]

    have_optional = [c for c in OPTIONAL if c in ds.samplesheet.columns]
    missing = [c for c in REQUIRED if c not in ds.samplesheet.columns]
    if missing:
        raise ValueError(
            f"Dataset samplesheet is missing required column(s): {', '.join(missing)}. "
            f"Present: {', '.join(ds.samplesheet.columns)}"
        )

    table = ds.pivot_samplesheet(
        index=["sampleIndex", "sample", "lane"],
        pivot_columns="read",
        metadata_columns=REQUIRED + have_optional,
        column_prefix="fastq_",
    ).sort_values(by="sample")

    ds.logger.info("Pivoted samplesheet:")
    ds.logger.info(table.to_csv(index=None))
    return table


def setup_parameters(ds: PreprocessDataset):
    start = str(ds.params.get("from", "align")).lower()

    # ---- alignment ------------------------------------------------------
    if start == "align":
        table = build_samplesheet(ds)
        table.to_csv("samplesheet.csv", index=None)
        # ${launchDir} rather than an absolute path: the head job and the tasks
        # do not share a filesystem view on Batch.
        ds.add_param("input_samplesheet", "${launchDir}/samplesheet.csv")
    else:
        ds.logger.info(f"Starting at '{start}', so no cellranger samplesheet is built.")

    # main.nf reads two booleans, not a mode string.
    mode = ds.params.get("alignment_mode", "standard")
    ds.add_param("multi", mode == "multi")
    ds.add_param("demux", mode == "demux")
    ds.remove_param("alignment_mode")

    # ---- blanks ---------------------------------------------------------
    dropped = [
        k for k in BLANK_MEANS_UNSET
        if k in ds.params and not str(ds.params[k]).strip()
    ]
    for k in dropped:
        ds.remove_param(k)
    if dropped:
        ds.logger.info(f"Dropped blank optional paths: {', '.join(dropped)}")


if __name__ == "__main__":

    ds = PreprocessDataset.from_running()

    ds.logger.info("Files annotated in the dataset:")
    ds.logger.info(ds.files)

    try:
        ds.logger.info("Samplesheet columns:")
        ds.logger.info(ds.samplesheet.columns)
    except Exception as e:                      # a dataset need not carry one
        ds.logger.info(f"  no samplesheet available: {e}")

    ds.logger.info(f"Launch directory: {os.getcwd()}")

    setup_parameters(ds)

    # "Which stages did this run actually do?" is the first question asked of any
    # output, and --from/--to is the only place the answer lives.
    ds.logger.info(
        f"Stage range: --from {ds.params.get('from', 'align')} "
        f"--to {ds.params.get('to', 'end')}"
    )

    ds.logger.info("Final parameters:")
    ds.logger.info(ds.params)
