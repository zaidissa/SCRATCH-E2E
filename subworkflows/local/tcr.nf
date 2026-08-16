/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    TCR — VDJ QC, T-cell integration, clonotype clustering, consensus, summary
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    Wraps the ten TCR subworkflows that previously lived in their own pipeline's
    `main.nf`, so the whole thing becomes one stage of the merged DAG that can be
    fed directly from ALIGN's VDJ output and ANNOTATION's labelled object.

    The optional-input convention here (a NO_FILE placeholder carried forward
    when a module is disabled) is the pattern the other stages were rebuilt
    around, because it lets the master summary run against whatever subset of
    tools actually executed.
----------------------------------------------------------------------------------------
*/

include { VDJ_QC_SW            } from './vdj_qc.nf'
include { TCELL_INTEGRATION_SW } from './tcell_integration.nf'
include { TCRI_SW              } from './tcri.nf'
include { CONGA_SW             } from './conga.nf'
include { GLIPH2_SW            } from './gliph2.nf'
include { TCRDIST3_SW          } from './tcrdist3.nf'
include { GIANA_SW             } from './giana.nf'
include { CONSENSUS_SW         } from './consensus_clustering.nf'
include { REPERTOIRE_SW        } from './repertoire.nf'
include { MASTER_SUMMARY_SW    } from './master_summary.nf'

// Each optional VDJ artefact gets its OWN NO_FILE placeholder. The standalone
// pipeline shared one for all eight, which makes MASTER_SUMMARY fail with an
// "input file name collision" whenever more than one is absent.
def pickQcTable(ch, name, slot) {
    ch.flatten().filter { it.getName() == name }
      .ifEmpty(file("${projectDir}/assets/NO_FILE_${slot}"))
}

def pickQcFigure(ch, name, slot) {
    ch.flatten().filter { it.getName() == name }
      .ifEmpty(file("${projectDir}/assets/NO_FILE_${slot}"))
}

workflow TCR {

    take:
        ch_annotated_object   // annotated Seurat object
        ch_sample_sheet       // TCR sample sheet (may be the NO_FILE placeholder)

    main:

        def nofile      = file("${projectDir}/assets/NO_FILE")
        ch_project_name = Channel.value(params.project_name)

        // ---- Step 1: VDJ QC ---------------------------------------------
        vdj_qc_out = VDJ_QC_SW(ch_sample_sheet, ch_project_name, ch_annotated_object)

        // Curated QC artefacts for the master summary.
        //
        // Each optional slot falls back to its OWN placeholder file. The
        // standalone pipeline used one shared `assets/NO_FILE` for all eight,
        // which makes MASTER_SUMMARY fail with "input file name collision"
        // whenever more than one of them is absent — i.e. on any run where VDJ
        // QC did not emit the full set.
        // Closure variables inside a workflow body are rejected by Nextflow
        // 26.x ("`pick_table` is not defined"); top-level functions resolve.

        ch_compact        = pickQcTable(vdj_qc_out.qc_tables, 'vdj_qc_per_sample_compact.tsv',      'compact')
        ch_before_after   = pickQcTable(vdj_qc_out.qc_tables, 'qc_contigs_before_after_summary.tsv','before_after')
        ch_sheet_resolved = pickQcTable(vdj_qc_out.qc_tables, 'sample_sheet_resolved.tsv',          'sheet_resolved')
        ch_rank_abundance = pickQcTable(vdj_qc_out.qc_tables, 'clone_rank_abundance.tsv',           'rank_abundance')

        ch_fig_retention  = pickQcFigure(vdj_qc_out.qc_figures, 'qc_before_after_retention.png',     'fig_retention')
        ch_fig_pairing    = pickQcFigure(vdj_qc_out.qc_figures, 'pairing_bar_by_sample.png',         'fig_pairing')
        ch_fig_rank       = pickQcFigure(vdj_qc_out.qc_figures, 'clone_rank_abundance.png',          'fig_rank')
        ch_fig_chains     = pickQcFigure(vdj_qc_out.qc_figures, 'multiple_chains_by_sample.png',     'fig_chains')

        // ---- Step 2: T-cell integration (baseline for everything below) --
        tcell_out = TCELL_INTEGRATION_SW(
            vdj_qc_out.contigs_after_qc,
            ch_annotated_object,
            ch_project_name
        )

        ch_tcells  = tcell_out.seurat_tcells_with_tcr
        ch_export  = tcell_out.export_cells

        // ---- Step 3: independent clustering / analysis modules -----------
        ch_tcri_done       = Channel.empty()
        ch_conga_done      = Channel.empty()
        ch_gliph_done      = Channel.empty()
        ch_tcrdist_done    = Channel.empty()
        ch_giana_done      = Channel.empty()
        ch_repertoire_done = Channel.empty()

        ch_gliph_export   = Channel.value(nofile)
        ch_tcrdist_export = Channel.value(nofile)
        ch_giana_export   = Channel.value(nofile)

        if (params.tcr_run_tcri) {
            TCRI_SW(ch_tcells, ch_export, ch_project_name)
            ch_tcri_done = TCRI_SW.out.report_html
        }
        if (params.tcr_run_conga) {
            CONGA_SW(ch_tcells, ch_export, ch_project_name)
            ch_conga_done = CONGA_SW.out.report_html
        }
        if (params.tcr_run_gliph2) {
            GLIPH2_SW(ch_tcells, ch_export, ch_project_name)
            ch_gliph_done   = GLIPH2_SW.out.report_html
            ch_gliph_export = GLIPH2_SW.out.export_cells
        }
        if (params.tcr_run_tcrdist3) {
            TCRDIST3_SW(ch_tcells, ch_export, ch_project_name)
            ch_tcrdist_done   = TCRDIST3_SW.out.report_html
            ch_tcrdist_export = TCRDIST3_SW.out.export_cells
        }
        if (params.tcr_run_giana) {
            GIANA_SW(ch_tcells, ch_export, ch_project_name)
            ch_giana_done   = GIANA_SW.out.report_html
            ch_giana_export = GIANA_SW.out.export_cells
        }
        if (params.tcr_run_repertoire) {
            REPERTOIRE_SW(ch_tcells, ch_export, ch_project_name)
            ch_repertoire_done = REPERTOIRE_SW.out.report_html
        }

        // ---- Step 4: consensus across clustering methods ------------------
        if (params.tcr_run_consensus) {
            CONSENSUS_SW(
                ch_tcells, ch_export,
                ch_gliph_export, ch_tcrdist_export, ch_giana_export,
                ch_project_name
            )
        }

        // ---- Step 5: master summary (waits for every enabled module) ------
        ch_downstream_done = ch_tcri_done
            .mix(ch_conga_done)
            .mix(ch_gliph_done)
            .mix(ch_tcrdist_done)
            .mix(ch_giana_done)
            .mix(ch_repertoire_done)

        ch_summary = Channel.empty()
        if (params.tcr_run_master_summary) {
            MASTER_SUMMARY_SW(
                ch_tcells, ch_export,
                ch_compact, ch_before_after, ch_sheet_resolved, ch_rank_abundance,
                ch_fig_retention, ch_fig_pairing, ch_fig_rank, ch_fig_chains,
                ch_downstream_done.collect().ifEmpty([]),
                ch_project_name
            )
            // The subworkflow emits three named channels; the summary tables are
            // what the integration layer joins on.
            ch_summary = MASTER_SUMMARY_SW.out.tables
        }

    emit:
        tcells      = ch_tcells
        export      = ch_export
        summary     = ch_summary
        contigs_qc  = vdj_qc_out.contigs_after_qc
}
