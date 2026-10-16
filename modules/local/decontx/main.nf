/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    DECONTX — ambient RNA removal without a GPU

    The alternative to CELLBENDER, selected with `--ambient_method decontx`.
    Same slot in QC, same output contract: a cellranger-layout .h5 that
    SEURAT_QUALITY reads with Read10X_h5, carried in the tuple shape that
    process already expects.

    Why it is here. CellBender trains for `--cellbender_epochs` passes over every
    droplet; on the Cirro CPU queue that is ~6 h per library, and the 4-lane GBM
    libraries retried for four days before the run was terminated. `--cuda`
    alone does not help: nothing in this pipeline requests an accelerator, so the
    task lands on a CPU instance regardless. DecontX is variational EM over the
    CELL matrix and has no such dependency.

    This is a substitution, not a downgrade to be hidden: the 2026 benchmark of
    seven decontamination methods places CellBender, DecontX and SoupX in the
    same performing band and recommends choosing between them on the compute
    environment available. Which method ran is recorded per sample and stated in
    that sample's QC report, so a run that mixes them is never silent.

    Unlike CELLBENDER this process reads the FILTERED matrix. DecontX does not
    need empty droplets; it uses them only to sharpen the ambient profile, so the
    raw matrix is optional here and its absence is not an error.
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

process DECONTX {

    tag "Removing ambient RNA from ${sample_id} (DecontX)"
    label 'process_medium'

    input:
        // `background` is raw_feature_bc_matrix.h5 when the sample has one, and
        // assets/NO_FILE when it does not — a path input cannot be null.
        tuple val(sample_id), path(csv_metrics), path(filtered), path(background)

    output:
        tuple val(sample_id), path(csv_metrics), path("decontx_${sample_id}_filtered.h5"), emit: corrected
        path("decontx_${sample_id}_contamination.csv"), emit: contamination
        path("decontx_${sample_id}_summary.json")     , emit: summary

    when:
        task.ext.when == null || task.ext.when

    script:
        def seed = (params.decontx_seed ?: 12345) as int
        """
        decontx_run.R ${filtered} ${sample_id} ${background} ${seed} ${task.cpus}
        """

    stub:
        // Output names must match the real ones exactly, or a -stub run proves
        // nothing about what qc.nf consumes downstream.
        """
        touch decontx_${sample_id}_filtered.h5
        echo "barcode,contamination" > decontx_${sample_id}_contamination.csv
        echo '{"sample":"${sample_id}","method":"decontx","mean_contamination":0}' \\
            > decontx_${sample_id}_summary.json
        """
}
