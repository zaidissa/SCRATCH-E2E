#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

include { VDJ_QC } from '../../modules/local/VDJ_QC/main.nf'

workflow VDJ_QC_SW {

    take:
        ch_sample_sheet
        ch_project_name
        ch_input_annotated_object

    main:
        ch_notebook_vdj_qc = Channel.fromPath(params.notebook_vdj_qc, checkIfExists: true)

        ch_vdj_qc = VDJ_QC(
            ch_notebook_vdj_qc,
            ch_sample_sheet,
            ch_input_annotated_object,
            ch_project_name
        )

    emit:
        report_html      = ch_vdj_qc.report_html
        contigs_after_qc = ch_vdj_qc.contigs_after_qc
        qc_tables        = ch_vdj_qc.tables
        
        qc_figures       = ch_vdj_qc.figures
}
// #!/usr/bin/env nextflow
// nextflow.enable.dsl = 2

// include { VDJ_QC } from '../../modules/local/VDJ_QC/main.nf'

// workflow VDJ_QC_SW {

//     take:
//         ch_sample_sheet
//         ch_project_name

//     main:
//         ch_notebook_vdj_qc = Channel.fromPath(params.notebook_vdj_qc, checkIfExists: true)

//         ch_vdj_qc = VDJ_QC(
//             ch_notebook_vdj_qc,
//             ch_sample_sheet,
//             ch_project_name
//         )

//     emit:
//         report_html      = ch_vdj_qc.report_html
//         contigs_after_qc = ch_vdj_qc.contigs_after_qc
//         qc_tables        = ch_vdj_qc.tables
// }