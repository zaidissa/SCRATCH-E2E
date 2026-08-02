/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    CELL COMMUNICATION — LIANA, CellChat, NicheNet
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    LIANA and CellChat are independent and run concurrently. NicheNet consumes
    both of their outputs.

    Changed on merge: the standalone pipeline declared NicheNet's inputs as
    hardcoded absolute paths built inside the `params` block —

        liana_csv    = "${params.outdir}/${params.project_name}/cellcommunication/liana/data/liana_results.csv"
        cellchat_rds = "${params.outdir}/.../cellchat_object.rds"

    which pointed at files that only exist after the upstream tasks publish, and
    additionally self-referenced `params.outdir` during params evaluation. It
    worked only because a separate `_done` tick channel imposed ordering. Those
    two params are gone; NicheNet now takes the real output channels, so the
    dependency is expressed in the DAG and `-resume` behaves correctly.
----------------------------------------------------------------------------------------
*/

include { CELLCOMM_LIANA    } from '../../modules/local/Cell_Communication/main.nf'
include { CELLCOMM_CELLCHAT } from '../../modules/local/Cell_Communication/main.nf'
include { CELLCOMM_NICHENET } from '../../modules/local/Cell_Communication/main.nf'

workflow CELL_COMMUNICATION {

    take:
        ch_seurat_object   // value channel

    main:
        ch_liana_csv    = Channel.empty()
        ch_cellchat_rds = Channel.empty()
        ch_nichenet     = Channel.empty()

        if (!params.skip_cellcomm) {

            nb_liana    = Channel.fromPath(params.liana_qmd,    checkIfExists: true)
            nb_cellchat = Channel.fromPath(params.cellchat_qmd, checkIfExists: true)
            nb_nichenet = Channel.fromPath(params.nichenet_qmd, checkIfExists: true)

            CELLCOMM_LIANA(ch_seurat_object.combine(nb_liana))
            CELLCOMM_CELLCHAT(ch_seurat_object.combine(nb_cellchat))

            ch_liana_csv    = CELLCOMM_LIANA.out.liana_csv
            ch_cellchat_rds = CELLCOMM_CELLCHAT.out.cellchat_rds

            // Real data dependency: NicheNet waits on the actual artefacts.
            CELLCOMM_NICHENET(
                ch_seurat_object
                    .combine(nb_nichenet)
                    .combine(ch_liana_csv)
                    .combine(ch_cellchat_rds)
                    .map { seurat, nb, liana_csv, cellchat_rds -> tuple(seurat, nb, liana_csv, cellchat_rds) },
                file(params.nichenet_assets_dir, type: 'dir', checkIfExists: true)
            )

            // NicheNet's `summary_csv` emit is commented out in the module, so
            // the generic `data` channel is the available handle.
            ch_nichenet = CELLCOMM_NICHENET.out.data
        }

    emit:
        liana_csv    = ch_liana_csv
        cellchat_rds = ch_cellchat_rds
        nichenet     = ch_nichenet
}
