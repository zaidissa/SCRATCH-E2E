#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

include { GLIPH2 } from '../../modules/local/GLIPH2/main.nf'

workflow GLIPH2_SW {

    take:
        seurat_rds
        export_cells
        project_name

    main:
        qmd = Channel.fromPath("${projectDir}/modules/local/GLIPH2/GLIPH2_Report.qmd", checkIfExists: true)
        gliph_reference_bundle = Channel.fromPath(params.gliph_reference_bundle, checkIfExists: true)

        GLIPH2(
            seurat_rds,
            export_cells,
            qmd,
            gliph_reference_bundle,
            project_name
        )

    emit:
        report_html        = GLIPH2.out.report_html
        seurat_with_gliph2 = GLIPH2.out.seurat_with_gliph2
        export_cells       = GLIPH2.out.export_cells
        tables             = GLIPH2.out.tables
        figures            = GLIPH2.out.figures
}


// include { GLIPH2 } from '../../modules/local/GLIPH2/main.nf'

// workflow GLIPH2_SW {
//     take:
//     seurat_rds
//     export_cells
//     project_name

//     main:
//     qmd = Channel.fromPath("${projectDir}/modules/local/GLIPH2/GLIPH2_Report.qmd", checkIfExists: true)
//     gliph_reference_bundle = Channel.fromPath(params.gliph_reference_bundle, checkIfExists: true)
//     // GLIPH2(seurat_rds, export_cells, qmd, project_name)

//     GLIPH2(seurat_rds, export_cells, qmd, gliph_reference_bundle, project_name)

//     emit:
//     report_html        = GLIPH2.out.report_html
//     seurat_with_gliph2 = GLIPH2.out.seurat_with_gliph2
//     export_cells       = GLIPH2.out.export_cells
// }