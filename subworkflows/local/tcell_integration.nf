#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

include { TCELL_INTEGRATION } from '../../modules/local/TCELL_INTEGRATION/main.nf'

workflow TCELL_INTEGRATION_SW {

    take:
        contigs_after_qc
        annotated_object
        project_name

    main:
        ch_notebook_tcell = Channel.fromPath(
            params.notebook_tcell_integration,
            checkIfExists: true
        )

        ch_tcell = TCELL_INTEGRATION(
            contigs_after_qc,
            annotated_object,
            ch_notebook_tcell,
            project_name
        )

    emit:
        report_html            = ch_tcell.report_html
        seurat_tcells_with_tcr = ch_tcell.seurat_tcells_with_tcr
        export_cells           = ch_tcell.export_cells
}