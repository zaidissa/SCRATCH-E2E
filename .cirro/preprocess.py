#!/usr/bin/env python3
"""
Cirro preprocessing for SCRATCH-E2E.

Runs before the workflow is launched. Two jobs: log enough about the dataset
that a failed launch can be diagnosed without re-running, and drop the
entry-point object parameters that were left blank in the form.

That second one matters. The form offers three entry-point objects because the
right one depends on where you start, so at most one is ever filled. Passing the
other two through as empty strings makes Nextflow receive
`--input_tumor_object ''`, which is truthy enough to defeat the `?:` chain in
main.nf's entry-point resolution and produces a confusing "file not found" for a
path nobody typed.
"""

import os

from cirro.helpers.preprocess_dataset import PreprocessDataset

# Blank in the form means "not supplied", not "supplied as empty".
OPTIONAL_PATHS = [
    "input_seurat_object",
    "input_tumor_object",
    "input_nonmalignant_object",
    "input_reference_object",
]

if __name__ == "__main__":

    ds = PreprocessDataset.from_running()

    ds.logger.info("Files annotated in the dataset:")
    ds.logger.info(ds.files)

    ds.logger.info("Samplesheet columns:")
    try:
        ds.logger.info(ds.samplesheet.columns)
    except Exception as e:                      # a dataset need not carry one
        ds.logger.info(f"  no samplesheet available: {e}")

    ds.logger.info(f"Launch directory: {os.getcwd()}")

    dropped = [k for k in OPTIONAL_PATHS
               if k in ds.params and not str(ds.params[k]).strip()]
    for k in dropped:
        ds.remove_param(k)
    if dropped:
        ds.logger.info(f"Dropped blank optional paths: {', '.join(dropped)}")

    # Record the stage range explicitly. "Which stages did this run actually
    # do?" is the first question asked of any output, and --from/--to is the
    # only place the answer lives.
    ds.logger.info(
        f"Stage range: --from {ds.params.get('from', 'align')} "
        f"--to {ds.params.get('to', 'end')}"
    )

    ds.logger.info("Final parameters:")
    ds.logger.info(ds.params)
