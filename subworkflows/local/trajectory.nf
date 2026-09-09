/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    TRAJECTORY — pseudotime / lineage inference
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { asBool } from '../../lib/booleans.nf'
include { CELLTRAJECTORY } from '../../modules/local/cell_trajectory/main.nf'

workflow TRAJECTORY {

    take:
        ch_seurat_object
        ch_page_config

    main:
        ch_rds = Channel.empty()
        ch_obj = Channel.empty()

        if (!asBool(params.skip_trajectory)) {
            CELLTRAJECTORY(ch_seurat_object,
                           Channel.fromPath(params.notebook_celltraj, checkIfExists: true),
                           ch_page_config)
            ch_rds = CELLTRAJECTORY.out.rds
            ch_obj = CELLTRAJECTORY.out.seurat_rds
        }

    emit:
        rds = ch_rds
        // The non-malignant object with pseudotime written into it. Optional at
        // the module level, so callers must fall back to their input when
        // trajectory did not produce one.
        seurat_rds = ch_obj
}
