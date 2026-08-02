/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    TRAJECTORY — pseudotime / lineage inference
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { CELLTRAJECTORY } from '../../modules/local/cell_trajectory/main.nf'

workflow TRAJECTORY {

    take:
        ch_seurat_object
        ch_page_config

    main:
        ch_rds = Channel.empty()

        if (!params.skip_trajectory) {
            CELLTRAJECTORY(ch_seurat_object,
                           Channel.fromPath(params.notebook_celltraj, checkIfExists: true),
                           ch_page_config)
            ch_rds = CELLTRAJECTORY.out.rds
        }

    emit:
        rds = ch_rds
}
