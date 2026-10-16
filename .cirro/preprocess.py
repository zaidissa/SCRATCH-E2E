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

3. Find Numbat's BAM in the input dataset instead of making the user retype a
   path Cirro already knows. Selecting a dataset maps only the matrices and the
   metadata table, so a cellranger dataset full of possorted_genome_bam.bam
   files still produced "Numbat skipped: no BAM" unless the glob was typed into
   the form by hand. Twice, on real launches.

4. Log enough that a failed launch can be diagnosed without re-running it.
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


BAM_NAME = "possorted_genome_bam.bam"


def _truthy(value) -> bool:
    return str(value).strip().lower() in ("true", "yes", "y", "on", "1")


def _bam_glob_from_matrices(gex_path: str):
    """Turn the matrices glob into the BAM glob beside it.

    process-input.json builds input_gex_matrices_path as
        <dataset>/data/data/SCRATCH_ALIGN:CELLRANGER_COUNT/**/outs/*
    and cellranger writes possorted_genome_bam.bam into that same `outs`, so the
    BAM glob is the matrices glob with its last segment replaced. Returns None
    for any path that does not look like a cellranger layout, so an unusual
    input simply falls back to the old behaviour.
    """
    gex_path = (gex_path or "").strip()
    if not gex_path.startswith("s3://") or "/outs/" not in gex_path:
        return None
    return gex_path.rsplit("/outs/", 1)[0] + "/outs/" + BAM_NAME


def _files_table_lists_bam(ds: PreprocessDataset) -> bool:
    """True when the dataset's own file table already names a BAM."""
    try:
        files = ds.files
        if files is None or len(files) == 0:
            return False
        return bool(
            files.astype(str)
            .apply(lambda col: col.str.endswith(BAM_NAME))
            .to_numpy()
            .any()
        )
    except Exception:
        return False


def _s3_lists_bam(bam_glob: str, ds: PreprocessDataset, max_pages: int = 5) -> bool:
    """Confirm at least one BAM exists under the glob's fixed prefix.

    This check is the point of the whole function. main.nf reads the path with
    `Channel.fromPath(..., checkIfExists: true)`, so handing it a glob that
    matches nothing converts today's graceful "Numbat skipped" warning into a
    hard failure at DAG build. Deriving a path is only safe once something has
    been seen at the other end of it.
    """
    try:
        import boto3

        bucket, _, key = bam_glob[len("s3://"):].partition("/")
        prefix = key.split("**", 1)[0].split("*", 1)[0]
        pages = 0
        for page in boto3.client("s3").get_paginator("list_objects_v2").paginate(
            Bucket=bucket, Prefix=prefix
        ):
            pages += 1
            for obj in page.get("Contents", []):
                if obj["Key"].endswith(BAM_NAME):
                    return True
            if pages >= max_pages:
                break
        return False
    except Exception as exc:
        ds.logger.info(f"Numbat: could not list S3 for BAMs ({exc})")
        return False


def resolve_numbat_bam(ds: PreprocessDataset, start: str) -> None:
    """Fill --input_bam_path from the input dataset when the form left it blank.

    Starting at `align` the pipeline produces its own BAMs, and the input is a
    FASTQ dataset that has none, so nothing is derived there.
    """
    if start == "align":
        return
    if not _truthy(ds.params.get("run_numbat", False)):
        return
    if str(ds.params.get("input_bam_path", "")).strip():
        return  # an explicit path always wins

    bam_glob = _bam_glob_from_matrices(ds.params.get("input_gex_matrices_path", ""))
    if not bam_glob:
        ds.logger.info(
            "Numbat: input_gex_matrices_path is not a cellranger layout, so no BAM "
            "glob could be derived. Numbat will self-disable."
        )
        return

    if not (_files_table_lists_bam(ds) or _s3_lists_bam(bam_glob, ds)):
        ds.logger.info(
            f"Numbat: no {BAM_NAME} found under the input dataset, so "
            "--input_bam_path stays unset and Numbat will self-disable with a warning."
        )
        return

    ds.add_param("input_bam_path", bam_glob, overwrite=True)
    ds.logger.info(f"Numbat: derived --input_bam_path from the input dataset: {bam_glob}")


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
        # ...and the per-sample METADATA table that QC's merge step reads.
        #
        # process-input.json maps input_exp_table to
        #   <input dataset>/data/pipeline_info/samplesheet.valid.csv
        # which exists only inside a PREVIOUS PIPELINE'S OUTPUT dataset. Starting
        # at `align` the input is a FASTQ dataset, which has no such file, so the
        # run died at DAG build — before a single task — with
        #   No such file or directory: .../data/pipeline_info/samplesheet.valid.csv
        #
        # The samplesheet just written is a valid metadata table: the merge
        # notebook requires a `sample` column, drops fastq_1/fastq_2/modality and
        # calls distinct(), so the one-row-per-lane layout collapses to one row
        # per sample carrying patient_id, timepoint and batch. The sample names
        # match end to end — cellranger is run with --id="${sample}", QC sets
        # orig.ident from it, and merge joins orig.ident to metadata$sample.
        #
        # overwrite=True is required: add_param asserts when the key already
        # exists, and process-input.json always supplies this one.
        ds.add_param("input_exp_table", "${launchDir}/samplesheet.csv", overwrite=True)
    else:
        ds.logger.info(f"Starting at '{start}', so no cellranger samplesheet is built.")

    # main.nf reads two booleans, not a mode string.
    mode = ds.params.get("alignment_mode", "standard")
    ds.add_param("multi", mode == "multi")
    ds.add_param("demux", mode == "demux")
    ds.remove_param("alignment_mode")

    # TCR analysis is a TRI-STATE, not a checkbox, and it has to stay one.
    # main.nf reads run_tcr as "null = decide from --from/--to" and applies the
    # override only when it is non-null, so a boolean field defaulting to false
    # would send false on EVERY launch and silently switch the TCR stage off in
    # full runs that should fan out to it.
    #
    # 'yes' is what makes a TCR-only run possible: with Stop after = align the
    # range enables no stage that needs gene expression, and run_tcr = true adds
    # the TCR stage back on top of alignment.
    tcr = str(ds.params.get("tcr_analysis", "auto")).strip().lower()
    if tcr == "yes":
        ds.add_param("run_tcr", True, overwrite=True)
    elif tcr == "no":
        ds.add_param("run_tcr", False, overwrite=True)
    # 'auto' leaves run_tcr unset, which is the null the stage logic expects.
    if "tcr_analysis" in ds.params:
        ds.remove_param("tcr_analysis")

    # ---- Numbat's BAM ----------------------------------------------------
    # Must run BEFORE the blank-dropping below: a derived path is non-blank and
    # therefore survives it, while a dataset with no BAMs leaves the field blank
    # and is dropped exactly as before.
    resolve_numbat_bam(ds, start)

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
