/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    CELLBENDER — ambient RNA removal

    Reads the RAW matrix, not the filtered one. CellBender infers the ambient
    profile from empty droplets, so the filtered matrix — which is exactly the
    barcodes cellranger already called as cells — contains none of the evidence
    the model needs.

    `--expected-cells` and `--total-droplets-included` are derived PER SAMPLE.
    They were previously two global constants (5,000 and 15,000) applied to every
    sample, which cannot be right for more than one library. On the HGSOC pair
    they were wrong in both directions: cellranger called 11,551 and 18,128 cells
    against an expectation of 5,000, and for HRS371760 the 15,000 droplet cap sat
    BELOW its own cell count — CellBender would have been told that 3,128 real
    cells were empty droplets and used them to build the background profile.
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

process CELLBENDER {

    tag "Removing ambient RNA from ${sample_id}"
    label 'process_high'

    input:
        // `raw` is raw_feature_bc_matrix.h5; qc.nf guarantees it is non-null.
        tuple val(sample_id), path(csv_metrics), path(raw)

    output:
        // The FILTERED output is what continues downstream: CellBender re-calls
        // cells after denoising, and the unfiltered file holds every droplet it
        // was asked to model (tens of thousands of empties), which is not a cell
        // matrix. Emitted under the tuple shape SEURAT_QUALITY already expects.
        tuple val(sample_id), path(csv_metrics), path("cellbender_${sample_id}_matrix_filtered.h5"), emit: corrected
        path("cellbender_${sample_id}_matrix.h5")                     , emit: full
        path("cellbender_${sample_id}_matrix_cell_barcodes.csv")      , emit: cell_barcodes, optional: true
        path("cellbender_${sample_id}_matrix.pdf")                    , emit: report, optional: true
        path("cellbender_${sample_id}_params.json")                   , emit: params_used

    when:
        task.ext.when == null || task.ext.when

    script:
        def args = task.ext.args ?: ''
        def out  = "cellbender_${sample_id}_matrix.h5"
        // null/0 means "derive from the data"; an explicit value overrides.
        def exp_override = (params.cellbender_expected_cells ?: 0) as int
        def tot_override = (params.cellbender_total_droplets ?: 0) as int
        """
        # ------------------------------------------------------------------
        # Derive the two sizing parameters from THIS sample.
        #
        # expected-cells    cellranger's own call, read from metrics_summary.csv
        # total-droplets    number of barcodes with >= 100 UMI, i.e. the point at
        #                   which the rank plot has clearly entered the ambient
        #                   plateau, clamped to [2x, 6x] expected cells so a
        #                   pathological curve cannot produce a value below the
        #                   cell count (too few empties to learn from) or an
        #                   enormous one (runtime grows with this number).
        # ------------------------------------------------------------------
        python3 <<'PYEOF'
import csv, json, re
import h5py, numpy as np

exp_override, tot_override = ${exp_override}, ${tot_override}

expected = None
try:
    row = next(iter(csv.DictReader(open("${csv_metrics}"))))
    for k, v in row.items():
        if k.strip().lower() == "estimated number of cells":
            expected = int(re.sub(r"[^0-9]", "", v))
            break
except Exception as e:
    print("could not read metrics_summary.csv:", e)

with h5py.File("${raw}", "r") as f:
    g = f["matrix"]
    indptr = g["indptr"][:]
    data   = g["data"][:]
    # Cumulative-sum differencing gives per-barcode UMI totals without
    # materialising the matrix; np.add.reduceat overruns on the final segment.
    cs  = np.concatenate(([0], np.cumsum(data, dtype=np.int64)))
    umi = cs[indptr[1:]] - cs[indptr[:-1]]

n_ge100 = int((umi >= 100).sum())
if expected is None:                       # no metrics file: fall back to the curve
    expected = int((umi >= 500).sum()) or 5000

total = tot_override or int(min(max(n_ge100, 2 * expected), 6 * expected))
expected = exp_override or expected

# Both must hold for the model to be identifiable at all.
if total <= expected:
    total = 3 * expected

json.dump({"sample": "${sample_id}", "expected_cells": int(expected),
           "total_droplets_included": int(total), "barcodes_ge_100_umi": n_ge100,
           "expected_source": "override" if exp_override else "metrics_summary.csv",
           "total_source": "override" if tot_override else "umi_curve"},
          open("cellbender_${sample_id}_params.json", "w"), indent=2)

with open("cb_args.txt", "w") as fh:
    fh.write(f"--expected-cells {expected} --total-droplets-included {total}")
print(f"${sample_id}: expected-cells={expected} total-droplets-included={total} "
      f"(barcodes >=100 UMI: {n_ge100})")
PYEOF

        cellbender remove-background \\
            --input ${raw} \\
            --output ${out} \\
            --cpu-threads ${task.cpus} \\
            \$(cat cb_args.txt) \\
            ${args}
        """

    stub:
        def out = "cellbender_${sample_id}_matrix.h5"
        """
        # Must match the real output names exactly, or -stub proves nothing about
        # what qc.nf consumes downstream.
        touch ${out}
        touch cellbender_${sample_id}_matrix_filtered.h5
        touch cellbender_${sample_id}_matrix_cell_barcodes.csv
        touch cellbender_${sample_id}_matrix.pdf
        echo '{"sample":"${sample_id}","expected_cells":0,"total_droplets_included":0}' \\
            > cellbender_${sample_id}_params.json
        """
}
