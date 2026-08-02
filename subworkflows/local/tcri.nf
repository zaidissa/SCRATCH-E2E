#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

include { TCRI } from '../../modules/local/TCRI/main.nf'

workflow TCRI_SW {

    take:
        ch_seurat_tcells_with_tcr
        ch_export_cells
        ch_project_name

    main:
        ch_notebook_tcri = Channel.fromPath(params.notebook_tcri, checkIfExists: true)

        ch_tcri = TCRI(
            ch_notebook_tcri,
            ch_seurat_tcells_with_tcr,
            ch_export_cells,
            ch_project_name
        )

    emit:
        report_html      = ch_tcri.report_html
        seurat_with_tcri = ch_tcri.seurat_with_tcri
        export_cells     = ch_tcri.export_cells
}



// #!/usr/bin/env nextflow
// nextflow.enable.dsl = 2

// include { TCRI } from '../../modules/local/TCRI/main.nf'

// workflow TCRI_SW {

//     take:
//         ch_seurat_tcells_with_tcr
//         ch_export_cells
//         ch_project_name

//     main:
//         ch_notebook_tcri = Channel.fromPath(params.notebook_tcri, checkIfExists: true)

//         ch_tcri = TCRI(
//             ch_notebook_tcri,
//             ch_seurat_tcells_with_tcr,
//             ch_export_cells,
//             ch_project_name
//         )

//     emit:
//         report_html      = ch_tcri.report_html
//         seurat_with_tcri = ch_tcri.seurat_with_tcri
//         export_cells     = ch_tcri.export_cells
// }