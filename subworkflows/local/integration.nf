/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    INTEGRATION — cross-stage synthesis and the master report
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    Runs last. Every input is `.collect().ifEmpty([])` at the call site, so a
    stage that was skipped contributes an empty list rather than stalling the
    DAG — the report simply records that stage as not run.
----------------------------------------------------------------------------------------
*/

include { INTEGRATE_OBJECTS } from '../../modules/local/integration/main.nf'
include { CROSS_ANALYSIS    } from '../../modules/local/integration/main.nf'
include { RUN_MANIFEST      } from '../../modules/local/integration/main.nf'
include { MASTER_REPORT     } from '../../modules/local/integration/main.nf'

workflow INTEGRATION {

    take:
        ch_seurat_object
        ch_annotated_object
        ch_cnv_infercnv
        ch_cnv_scevan
        ch_metaprog
        ch_trajectory
        ch_liana
        ch_cellchat
        ch_tcr
        ch_compartment   // stratification assignment + TME tables
        ch_page_config

    main:

        ch_nb_integrate = Channel.fromPath("${projectDir}/modules/local/integration/integrate_objects.qmd", checkIfExists: true)
        ch_nb_cross     = Channel.fromPath("${projectDir}/modules/local/integration/cross_analysis.qmd",    checkIfExists: true)
        ch_nb_master    = Channel.fromPath("${projectDir}/modules/local/integration/master_report.qmd",     checkIfExists: true)

        // ---- 1. Join every stage's per-cell annotation onto one object ----
        INTEGRATE_OBJECTS(
            ch_seurat_object.first(),
            ch_nb_integrate,
            ch_page_config,
            ch_cnv_infercnv,
            ch_cnv_scevan,
            ch_metaprog,
            ch_trajectory,
            ch_tcr,
            ch_compartment
        )

        // ---- 2. Analyses that need two or more stages at once -------------
        // Cell-communication outputs are read directly by the notebook, so they
        // are staged in alongside the other stage artefacts.
        ch_cross_inputs = ch_liana
            .mix(ch_cellchat)
            .mix(ch_metaprog)
            .mix(ch_compartment)
            .collect()
            .ifEmpty([])

        CROSS_ANALYSIS(
            INTEGRATE_OBJECTS.out.cell_table,
            ch_nb_cross,
            ch_page_config,
            ch_cross_inputs
        )

        // ---- 3. Provenance -------------------------------------------------
        // Serialise the resolved parameters so the report shows what actually
        // ran rather than the defaults.
        ch_params_json = Channel.of(groovy.json.JsonOutput.prettyPrint(
                                        groovy.json.JsonOutput.toJson(params.findAll { k, v ->
                                            !(k.toString().startsWith('genomes'))
                                        })))
            .collectFile(name: 'resolved_params.json')

        // `ready` is an ordering barrier only: the manifest is written after the
        // cross analysis so the trace file it copies is as complete as possible.
        RUN_MANIFEST(
            CROSS_ANALYSIS.out.report.map { 'ready' }.ifEmpty('ready').first(),
            ch_params_json
        )

        // ---- 4. One navigable site over everything -------------------------
        ch_all_reports = INTEGRATE_OBJECTS.out.report
            .mix(CROSS_ANALYSIS.out.report)
            .collect()
            .ifEmpty([])

        ch_all_figures = INTEGRATE_OBJECTS.out.figures
            .mix(CROSS_ANALYSIS.out.figures)
            .flatten()
            .collect()
            .ifEmpty([])

        MASTER_REPORT(
            ch_nb_master,
            ch_page_config,
            RUN_MANIFEST.out.manifest.mix(RUN_MANIFEST.out.tables).collect().ifEmpty([]),
            ch_all_reports,
            ch_all_figures,
            CROSS_ANALYSIS.out.data.flatten().collect().ifEmpty([])
        )

    emit:
        master_rds = INTEGRATE_OBJECTS.out.master_rds
        cell_table = INTEGRATE_OBJECTS.out.cell_table
        report     = MASTER_REPORT.out.report
}
