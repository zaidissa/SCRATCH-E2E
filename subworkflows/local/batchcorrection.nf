/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    BATCH CORRECTION — CCA / RPCA / Harmony / MNN with side-by-side evaluation
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    When skipped, the input object is passed through unchanged so downstream
    branches still receive an object rather than an empty channel.
----------------------------------------------------------------------------------------
*/

include { BATCHCORRECTION } from '../../modules/local/batch_correction/main.nf'

workflow BATCH_CORRECTION {

    take:
        ch_seurat_object
        ch_page_config

    main:
        ch_notebook = Channel.fromPath(params.notebook_batchcorr, checkIfExists: true)

        if (!params.skip_batchcorrect) {
            BATCHCORRECTION(ch_seurat_object, ch_notebook, ch_page_config)
            ch_out     = BATCHCORRECTION.out.seurat_rds.ifEmpty { null }
            ch_report  = BATCHCORRECTION.out.report
            ch_figures = BATCHCORRECTION.out.figures
            ch_data    = BATCHCORRECTION.out.data
        } else {
            ch_out     = ch_seurat_object
            ch_report  = Channel.empty()
            ch_figures = Channel.empty()
            ch_data    = Channel.empty()
        }

        // Fall back to the uncorrected object if the notebook's autosave block
        // did not run (auto_save = false), so the DAG never stalls here.
        ch_final = ch_out.filter { it != null }.concat(ch_seurat_object).first()

    emit:
        seurat_rds = ch_final
        report     = ch_report
        figures    = ch_figures
        data       = ch_data
}
