#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

include { GIANA } from '../../modules/local/GIANA/main.nf'

workflow GIANA_SW {

    take:
        seurat_rds
        export_cells
        project_name

    main:
        qmd = Channel.fromPath("${projectDir}/modules/local/GIANA/GIANA_Report.qmd", checkIfExists: true)

        GIANA(
            seurat_rds,
            export_cells,
            qmd,
            project_name
        )

    emit:
        report_html       = GIANA.out.report_html
        seurat_with_giana = GIANA.out.seurat_with_giana
        export_cells      = GIANA.out.export_cells
        tables            = GIANA.out.tables
        figures           = GIANA.out.figures
}

// include { GIANA } from '../../modules/local/GIANA/main'

// workflow GIANA_SW {
//     take:
//     seurat_rds
//     export_cells
//     project_name

//     main:
//     qmd = Channel.fromPath("${projectDir}/modules/local/GIANA/GIANA_Report.qmd", checkIfExists: true)
//     GIANA(seurat_rds, export_cells, qmd, project_name)

//     emit:
//     report_html       = GIANA.out.report_html
//     seurat_with_giana = GIANA.out.seurat_with_giana
//     export_cells      = GIANA.out.export_cells
// }