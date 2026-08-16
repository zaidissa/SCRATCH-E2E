/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    CLUSTER — normalisation, dimensionality reduction, clustering
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { SEURAT_NORMALIZE } from '../../modules/local/seurat/normalization/main.nf'
include { SEURAT_CLUSTER   } from '../../modules/local/seurat/cluster/main.nf'

workflow CLUSTER {

    take:
        ch_seurat_object
        // BPCells store from QC; staged at data/bpcells_counts in both processes.
        ch_bpcells_store
        ch_page_config

    main:

        ch_notebook_normalize  = Channel.fromPath(params.notebook_normalize,  checkIfExists: true)
        ch_notebook_clustering = Channel.fromPath(params.notebook_clustering, checkIfExists: true)

        SEURAT_NORMALIZE(ch_seurat_object, ch_bpcells_store, ch_notebook_normalize, ch_page_config)
        SEURAT_CLUSTER(SEURAT_NORMALIZE.out.seurat_rds, ch_bpcells_store, ch_notebook_clustering, ch_page_config)

    emit:
        seurat_rds     = SEURAT_CLUSTER.out.seurat_rds
        normalized_rds = SEURAT_NORMALIZE.out.seurat_rds
        figures        = SEURAT_CLUSTER.out.figures
}
