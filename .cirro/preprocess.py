#!/usr/bin/env python3
"""
Cirro preprocessing for SCRATCH-E2E.

Three jobs, in order of how badly their absence bites:

1. Build the cellranger samplesheet when the run starts at `align`. Cirro
   annotates a FASTQ dataset one row per file; cellranger wants one row per
   sample with R1 and R2 in columns. Without this, `--from align` fails on a
   samplesheet that was never written. R1/R2 are taken from the file names
   (also SRR_1/SRR_2 from a fetchngs "Public Dataset (FASTQ)"), never from
   table order -- see normalise_files().

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


REPO_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))


def resolve_exp_table_override(value) -> str:
    """Turn the form's "Sample metadata table" into an --input_exp_table value.

    Pure, so it can be tested off-platform. Returns "" for blank (keep the
    default). Accepts two shapes and refuses everything else, because a path
    Nextflow cannot open fails only at DAG build -- after the launch looks
    accepted:

      s3://bucket/key.csv           passed through unchanged
      assets/metadata/<name>.csv    a file shipped in this repository,
                                    returned as ${projectDir}/assets/...
                                    (Nextflow interpolates projectDir in the
                                    params file, as it does launchDir for the
                                    align samplesheet).

    Why the field exists: from `qc` onward process-input.json points
    input_exp_table at the input dataset's pipeline_info/samplesheet.valid.csv,
    and for GBM1-DFCI3 that file has no patient_id. Nothing on the form could
    replace it, so per-patient Numbat silently fell back to per library.
    """
    v = str(value if value is not None else "").strip()
    if not v:
        return ""
    if v.startswith("s3://"):
        if not v.lower().endswith(".csv"):
            raise ValueError(f"Sample metadata table must be a .csv file; got {v!r}")
        return v
    rel = v[2:] if v.startswith("./") else v
    if rel.startswith("assets/") and ".." not in rel.split("/"):
        if not rel.lower().endswith(".csv"):
            raise ValueError(f"Sample metadata table must be a .csv file; got {v!r}")
        return "${projectDir}/" + rel
    raise ValueError(
        f"Sample metadata table {v!r} is neither an s3:// path nor a path inside "
        "this repository starting with assets/ (e.g. "
        "assets/metadata/gbm1_dfci3_exp_table.csv)."
    )


def apply_exp_table_override(ds: PreprocessDataset) -> bool:
    """Replace input_exp_table with the form's table, if one was given.

    Runs for every start stage, align included: a table typed into the form is
    an explicit choice and wins over both defaults. Returns whether it applied.
    """
    raw = ds.params.get("exp_table_override", "")
    if "exp_table_override" in ds.params:
        ds.remove_param("exp_table_override")   # not a pipeline parameter
    resolved = resolve_exp_table_override(raw)
    if not resolved:
        return False
    previous = ds.params.get("input_exp_table", "(unset)")
    ds.add_param("input_exp_table", resolved, overwrite=True)
    ds.logger.info(f"--input_exp_table from the form: {resolved}")
    ds.logger.info(f"  (replaces {previous})")
    return True


def describe_exp_table(body: str, run_numbat: bool = False) -> list:
    """Lines describing a metadata CSV. Pure, so it can be tested off-platform.

    `run_numbat` escalates the missing-`patient_id` case from a note to a
    blocking-looking warning, because that one combination is not cosmetic:
    Numbat groups libraries by patient_id, and without it numbat_prep.R falls
    back to the library column and runs PER LIBRARY. That is the configuration
    numbat:::mark_tumor_lineage() cannot survive -- it has no null, so a library
    with no malignant cells still returns neoplastic calls.

    This happened. The GBM1-DFCI3 runs of 2026-10-22 and 2026-10-23 both logged
    "no `patient_id` ... column" as one INFO line among forty, enabled Numbat,
    and spent 1,166.7 CPU hours producing per-library calls -- while the
    pipeline's own per-individual support sat there unused. The information was
    present and not loud enough to act on.
    """
    lines = [ln for ln in body.splitlines() if ln.strip()]
    if not lines:
        return ["the table is empty."]

    columns = [c.strip().strip('"') for c in lines[0].split(",")]
    out = [f"columns: {', '.join(columns)}"]

    if "sample" in columns:
        idx = columns.index("sample")
        values = []
        for ln in lines[1:]:
            cells = [c.strip().strip('"') for c in ln.split(",")]
            if len(cells) > idx and cells[idx] and cells[idx] not in values:
                values.append(cells[idx])
        out.append(
            f"sample values ({len(values)} distinct, possibly truncated): "
            f"{', '.join(values)}"
        )
        out.append(
            "QC's merge joins these to orig.ident, which is the cellranger "
            "output DIRECTORY name. Strings that differ join nothing."
        )
    else:
        out.append(
            "WARNING: no `sample` column. The merge notebook stops on this, but "
            "only after every per-sample QC task has already run."
        )

    absent = [c for c in ("patient_id", "timepoint", "condition", "batch")
              if c not in columns]
    if absent:
        out.append(
            f"no {', '.join('`' + c + '`' for c in absent)} column, so the merged "
            "object will not carry it. Stages keyed on it fall back to another "
            "sample column, and condition-level tests have no design factor."
        )

    if run_numbat and "patient_id" not in columns:
        out.append("")
        out.append("*" * 72)
        out.append("WARNING: --run_numbat is ON and this table has no `patient_id`.")
        out.append("")
        out.append("  Numbat groups libraries by patient_id. Without it, numbat_prep.R")
        out.append("  falls back to the library column and runs PER LIBRARY --")
        out.append("  numbat:::mark_tumor_lineage() has no null, so it always labels its")
        out.append("  heaviest clade 'tumor' and a library with no malignant cells comes")
        out.append("  back with neoplastic calls anyway. Pooling an individual's")
        out.append("  libraries is what gives the phylogeny a real normal anchor.")
        out.append("")
        out.append("  Fix it WITHOUT re-running alignment: fill the form's 'Sample metadata")
        out.append("  table' field (Entry point object) -- it replaces --input_exp_table -- with a")
        out.append("  CSV of `sample,patient_id[,timepoint]` whose `sample` values are")
        out.append("  exactly the strings listed above. QC's merge joins them to")
        out.append("  orig.ident with an INNER join, so every library must appear.")
        out.append("  assets/metadata/gbm1_dfci3_exp_table.csv is the one for this cohort.")
        out.append("")
        out.append("  Launching anyway is legitimate -- a genuinely one-library-per-patient")
        out.append("  cohort has nothing to pool -- but then expect per-library calls and")
        out.append("  read pct_malignant_gated rather than pct_malignant_raw.")
        out.append("*" * 72)

    return out


def log_exp_table_columns(ds: PreprocessDataset, start: str,
                          overridden: bool = False) -> None:
    """Print the columns of the metadata table QC's merge step will read.

    Starting at `align` the table is the samplesheet written above, and
    build_samplesheet() has already asserted its required columns. Starting
    anywhere later, input_exp_table points at a PREVIOUS run's
    pipeline_info/samplesheet.valid.csv and nothing checks it at all.

    That gap cost a full GBM run. SAMPLESHEET_CHECK copies the align
    samplesheet verbatim, the merge notebook drops the fastq columns and joins
    `orig.ident` to `sample`, and the join is an INNER join. Measured in
    scratch-cluster:latest: a TOTAL key mismatch makes AddMetaData fail loudly
    ("Cannot add more or less meta data without cell names"), but a column that
    is simply ABSENT from the table passes through in silence -- the object is
    built without it, and the first thing to notice is a notebook three stages
    later asking for `patient_id`.

    Logging the header is enough to answer it at launch: the column names say
    whether the design factors survived, and the distinct `sample` values are
    the exact strings the merge join needs `orig.ident` to equal. Diagnostic
    only -- it never changes a parameter and never fails a launch, because a
    cohort with no patient_id is still a legitimate run.
    """
    if start == "align" and not overridden:
        return

    path = str(ds.params.get("input_exp_table", "")).strip()
    if path.startswith("${projectDir}/"):
        local = os.path.join(REPO_ROOT, path[len("${projectDir}/"):])
        if not os.path.isfile(os.path.join(REPO_ROOT, "main.nf")):
            # preprocess.py is not running inside a checkout, so the repo
            # cannot be looked at from here. Nextflow will still open it.
            ds.logger.info(f"--input_exp_table: {path} (no repo checkout beside "
                           "preprocess.py; not inspected)")
            return
        if not os.path.isfile(local):
            # A repo path that is not in the repo fails at DAG build with a bare
            # "No such file"; say which file and where it was looked for, now.
            raise FileNotFoundError(
                f"Sample metadata table {path} not found in this checkout "
                f"(looked for {local}). Is it committed and pushed?"
            )
        with open(local, encoding="utf-8", errors="replace") as fh:
            body = fh.read(65536)
        _log_exp_table(ds, path, body)
        return
    if not path.startswith("s3://"):
        ds.logger.info(f"--input_exp_table is not an S3 path ({path!r}); not inspected.")
        return

    try:
        import boto3

        bucket, _, key = path[len("s3://"):].partition("/")
        # The header and a few rows are all that is needed, and a metadata table
        # for a 17-library cohort is a few kB in total.
        body = (
            boto3.client("s3")
            .get_object(Bucket=bucket, Key=key, Range="bytes=0-65535")["Body"]
            .read()
            .decode("utf-8", "replace")
        )
    except Exception as exc:
        ds.logger.info(f"--input_exp_table could not be read ({exc}); not inspected.")
        return

    _log_exp_table(ds, path, body)


def _log_exp_table(ds: PreprocessDataset, path: str, body: str) -> None:
    # run_numbat is a form checkbox, so it arrives as a bool or as a string.
    numbat_on = _truthy(ds.params.get("run_numbat", False))
    ds.logger.info(f"--input_exp_table: {path}")
    for line in describe_exp_table(body, run_numbat=numbat_on):
        ds.logger.info(f"  {line}")


import re as _re

_R_TOKEN = _re.compile(r'_R([12])(?=[_.])')
_END_DIGIT = _re.compile(r'[_.]([12])\.f(?:ast)?q\.gz$')
_LANE = _re.compile(r'_L(\d{3})_')


def read_number_from_name(path: str):
    """1 or 2 from a FASTQ name, or None. Same rule as bin/cellranger_renaming.py:
    the LAST _R1/_R2 token, else a _1/_2 ending (fetchngs: SRR123_1.fastq.gz)."""
    name = os.path.basename(str(path))
    tokens = _R_TOKEN.findall(name)
    if tokens:
        return int(tokens[-1])
    m = _END_DIGIT.search(name)
    return int(m.group(1)) if m else None


def normalise_files(files: pd.DataFrame) -> pd.DataFrame:
    """Give every file a trustworthy `read`, `lane` and `sampleIndex`. Pure.

    Cirro's pivot pairs files on these. A FASTQ dataset annotated from
    Illumina names has them; a fetchngs "Public Dataset (FASTQ)" may not, and
    then pivot_files() numbers each sample's files in TABLE ORDER -- which can
    put R2 in fastq_1 -- and pivots on a blank lane. So:

      read   kept when it is 1 or 2 AND agrees with the file name; filled from
             the name when missing; a disagreement or an unreadable name stops
             the launch, because guessing is how reads get swapped.
      lane   kept when present, else parsed from _L00n_, else 1.
      sampleIndex  kept when present, else the sample name.
    """
    f = files.copy()
    from_name = f["file"].map(read_number_from_name)

    if "read" in f.columns:
        given = pd.to_numeric(f["read"], errors="coerce")
        clash = given.notna() & from_name.notna() & (given != from_name)
        if clash.any():
            bad = f.loc[clash, "file"].map(os.path.basename).tolist()
            raise ValueError(f"Read number in the dataset disagrees with the file name for: {bad}")
        f["read"] = given.where(given.isin([1, 2]), from_name)
    else:
        f["read"] = from_name

    unreadable = f["read"].isna()
    if unreadable.any():
        bad = f.loc[unreadable, "file"].map(os.path.basename).tolist()
        raise ValueError(
            "Cannot tell R1 from R2 for: " + ", ".join(bad) + ". Expected names like "
            "..._R1_001.fastq.gz / ..._R2_001.fastq.gz or ..._1.fastq.gz / ..._2.fastq.gz."
        )
    f["read"] = f["read"].astype(int)

    lane_from_name = f["file"].map(
        lambda p: int(m.group(1)) if (m := _LANE.search(os.path.basename(str(p)))) else 1)
    if "lane" in f.columns:
        given = pd.to_numeric(f["lane"], errors="coerce")
        f["lane"] = given.where(given.notna(), lane_from_name).astype(int)
    else:
        f["lane"] = lane_from_name

    if "sampleIndex" not in f.columns:
        f["sampleIndex"] = f["sample"]
    else:
        f["sampleIndex"] = f["sampleIndex"].where(f["sampleIndex"].notna(), f["sample"])
    return f


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
            f"Present: {', '.join(ds.samplesheet.columns)}. For a public (fetchngs) "
            f"dataset, add them to the samples' metadata in Cirro."
        )

    # A fetchngs "Public Dataset (FASTQ)" names files SRR_1/SRR_2 with no lane,
    # and may carry no `read` column -- then Cirro pairs files in table order.
    ds.logger.info(f"File table columns: {', '.join(map(str, ds.files.columns))}")
    ds.files = normalise_files(ds.files)

    table = ds.pivot_samplesheet(
        index=["sampleIndex", "sample", "lane"],
        pivot_columns="read",
        metadata_columns=REQUIRED + have_optional,
        column_prefix="fastq_",
    ).sort_values(by="sample")

    # Every row must be a complete pair. A missing mate would otherwise reach
    # cellranger as a path of "nan".
    for col in ("fastq_1", "fastq_2"):
        if col not in table.columns:
            raise ValueError(f"No {col} column after pairing files; columns: {list(table.columns)}")
    incomplete = table[table["fastq_1"].isna() | table["fastq_2"].isna()]
    if len(incomplete):
        raise ValueError("Samples without both R1 and R2: "
                         + ", ".join(incomplete["sample"].astype(str).unique()))

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

    # ---- both CNV callers, always ----------------------------------------
    # inferCNV and Numbat are not optional on Cirro: the form no longer offers
    # "Run Numbat", "Skip inferCNV" or "Download Numbat references". Set here
    # rather than left to nextflow.config so a launch made with "Copy
    # parameters" from an older run cannot carry an old off-switch through.
    # Numbat still self-disables, with a warning, when there is no BAM. The
    # references are read from the Numbat image (/data); the upstream download
    # host now answers with an HTML page instead of the phasing panel.
    ds.add_param("run_numbat", True, overwrite=True)
    ds.add_param("skip_infercnv", False, overwrite=True)
    ds.add_param("numbat_fetch_refs", False, overwrite=True)

    # ---- Numbat's BAM ----------------------------------------------------
    # Must run BEFORE the blank-dropping below: a derived path is non-blank and
    # therefore survives it, while a dataset with no BAMs leaves the field blank
    # and is dropped exactly as before.
    resolve_numbat_bam(ds, start)

    # ---- the metadata table ---------------------------------------------
    # After the align block, so a table typed into the form also replaces the
    # samplesheet align just wrote.
    overridden = apply_exp_table_override(ds)
    log_exp_table_columns(ds, start, overridden=overridden)

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
