#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

include { CONGA } from '../../modules/local/CONGA/main.nf'

workflow CONGA_SW {

    take:
        seurat_rds
        export_cells
        project_name

    main:
        qmd = Channel.fromPath("${projectDir}/modules/local/CONGA/CoNGA_Report.qmd", checkIfExists: true)

        CONGA(
            seurat_rds,
            export_cells,
            qmd,
            project_name
        )

    emit:
        report_html       = CONGA.out.report_html
        seurat_with_conga = CONGA.out.seurat_with_conga
        export_cells      = CONGA.out.export_cells
        data              = CONGA.out.data
        tables            = CONGA.out.tables
        figures           = CONGA.out.figures
}



// include { CONGA } from '../../modules/local/CONGA/main'

// workflow CONGA_SW {
//     take:
//     seurat_rds
//     export_cells
//     project_name

//     main:
//     qmd = Channel.fromPath("${projectDir}/modules/local/CONGA/CoNGA_Report.qmd", checkIfExists: true)
//     CONGA(seurat_rds, export_cells, qmd, project_name)

//     emit:
//     report_html       = CONGA.out.report_html
//     seurat_with_conga = CONGA.out.seurat_with_conga
//     export_cells      = CONGA.out.export_cells
// }