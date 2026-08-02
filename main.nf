#!/usr/bin/env nextflow
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    SCRATCH-E2E
    End-to-end single-cell tumour analysis
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    ALIGN --+-- GEX --> QC --> CLUSTER --> ANNOTATION --> CNV --> STRATIFY
            |                                  |                     |
            |                                  |     malignant ------+-----> METAPROG
            |                                  |                     |
            |                                  | non-malignant ------+-----> BATCHCORRECT
            |                                  |                                  |
            |                                  |                                  +--> TRAJECTORY
            |                                  |                                  +--> CELLCOMM
            |                                  |                                  +--> TME
            |                                  |
            +-- VDJ ---------------------------+--------------------------------------> TCR
            +-- BAM -------------------------------------------------------------------> MICROBIOME
                                                                                          |
                                     all stages ------------------------------------------+--> INTEGRATION

    The split is the pivot of the whole design. CNV calling runs on the
    annotated object; STRATIFY then uses those calls (with the annotation label
    as a tie-break) to separate malignant from non-malignant cells. Metaprograms
    are a property of the tumour compartment; batch correction, trajectory, cell
    communication and TME characterisation belong to the microenvironment.
    Correcting batch across a mixed object would regress out real tumour biology,
    which is why batch correction sits on the non-malignant arm.

    Run everything:
        nextflow run . -profile docker \
            --input_samplesheet sheet.csv --input_exp_table meta.csv \
            --project_name MYSTUDY

    Run part of it:
        nextflow run . --to cluster                     stop after clustering
        nextflow run . --from annotation --input_seurat_object cluster_object.RDS
        nextflow run . --only cnv        --input_seurat_object annotated.RDS
        nextflow run . --only cnv,stratify --input_seurat_object annotated.RDS

    See docs/USAGE.md for the full parameter reference.
----------------------------------------------------------------------------------------
*/

nextflow.enable.dsl = 2

include { ALIGN              } from './subworkflows/local/align.nf'
include { QC                 } from './subworkflows/local/qc.nf'
include { CLUSTER            } from './subworkflows/local/cluster.nf'
include { ANNOTATION         } from './subworkflows/local/annotation.nf'
include { CNV                } from './subworkflows/local/cnv.nf'
include { NUMBAT             } from './subworkflows/local/numbat.nf'
include { STRATIFY           } from './subworkflows/local/stratify.nf'
include { METAPROG           } from './subworkflows/local/metaprog.nf'
include { BATCH_CORRECTION   } from './subworkflows/local/batchcorrection.nf'
include { TRAJECTORY         } from './subworkflows/local/trajectory.nf'
include { CELL_COMMUNICATION } from './subworkflows/local/cellcomm.nf'
include { TME                } from './subworkflows/local/tme.nf'
include { TCR                } from './subworkflows/local/tcr.nf'
include { STAGE_REPORTS      } from './subworkflows/local/stage_reports.nf'
include { INTEGRATION        } from './subworkflows/local/integration.nf'


/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Stage selection

    Two mechanisms, in precedence order:
      1. --from / --to      inclusive range over the linear backbone
      2. --run_<stage>      explicit override; always wins when set
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// Self-contained: a Nextflow script's top-level `def` values are not visible
// from inside a `def` function, so the stage lists live in here.
def resolveStages() {

    // Linear backbone, in execution order, up to the compartment split.
    // CNV is on the backbone because STRATIFY consumes its calls.
    def backbone = ['align', 'qc', 'cluster', 'annotation', 'cnv', 'stratify']

    // The malignant arm.
    def tumour_arm = ['metaprog']

    // The non-malignant arm. Batch correction gates the three that follow it.
    def stromal_arm = ['batchcorrect', 'trajectory', 'cellcomm', 'tme']

    // Independent of the split.
    def independent = ['tcr', 'microbiome']

    def branches = tumour_arm + stromal_arm + independent
    def all = backbone + branches + ['integration']

    // --only is the ergonomic form of "just this step": it takes one stage or a
    // comma-separated set, runs exactly those, and nothing else — no range
    // arithmetic, no prerequisite pull-in, no integration tacked on the end.
    if (params.only) {
        def wanted = params.only.toString().toLowerCase()
                        .split('[,;]').collect { it.trim() }.findAll { it }
        def unknown = wanted.findAll { !(it in all) }
        if (unknown) error "--only ${unknown.join(', ')} is not a stage. Valid: ${all.join(', ')}"

        def picked = [:]
        all.each { picked[it] = (it in wanted) }
        // An explicit --run_<stage> still wins, so --only cnv --run_integration true works.
        all.each { name ->
            def ov = params."run_${name}"
            if (ov != null) picked[name] = (ov.toString().toLowerCase() in ['true','yes','1'])
        }
        return picked
    }

    def from = params.from?.toString()?.toLowerCase() ?: 'align'
    def to   = params.to?.toString()?.toLowerCase()   ?: 'end'

    if (!(from in all))
        error "--from '${from}' is not a stage. Valid: ${all.join(', ')}"
    if (to != 'end' && !(to in all))
        error "--to '${to}' is not a stage. Valid: end, ${all.join(', ')}"

    def enabled = [:]

    def fromIdx = backbone.indexOf(from)
    def toIdx   = (to == 'end') ? backbone.size() - 1 : backbone.indexOf(to)

    // Backbone: a branch name in --from (fromIdx < 0) skips the backbone entirely.
    backbone.eachWithIndex { name, i ->
        enabled[name] = (fromIdx >= 0) && (i >= fromIdx) && (toIdx < 0 || i <= toIdx)
    }

    // Branches.
    //
    // Naming a branch in --from means "run that branch", not "run every branch":
    // `--from trajectory` should not also launch metaprograms and TCR. Entering
    // from a BACKBONE stage with --to end is what fans out to all of them.
    def entering_at_branch = (fromIdx < 0) && (from in branches)

    branches.each { name ->
        enabled[name] = entering_at_branch ? (from == name)
                                           : ((to == 'end') || (to == name))
    }

    // The three stages on the non-malignant arm read the corrected object, so
    // entering at one of them still needs batch correction — unless that object
    // is supplied directly.
    if (entering_at_branch && (from in ['trajectory', 'cellcomm', 'tme'])
        && !params.input_nonmalignant_object) {
        enabled['batchcorrect'] = true
    }

    // Integration synthesises across stages, so it only makes sense at the end
    // of a full range (an explicit --run_integration true still forces it).
    enabled['integration'] = (to == 'end')

    // Explicit overrides always win.
    all.each { name ->
        def override = params."run_${name}"
        if (override != null)
            enabled[name] = (override.toString().toLowerCase() in ['true', 'yes', '1'])
    }

    return enabled
}


/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Main workflow
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow {

    def run  = resolveStages()
    def line = '=' * 74

    log.info """
    ${line}
     SCRATCH-E2E  v${workflow.manifest.version}
    ${line}
     Project        : ${params.project_name}
     Organism       : ${params.organism}  (genome ${params.genome})
     Output         : ${params.outdir}
     Range          : ${params.only ? "--only ${params.only}" : "--from ${params.from} --to ${params.to}"}
     Stages enabled : ${run.findAll { k, v -> v }.keySet().join(', ') ?: 'none'}
     NMF engine     : ${params.metaprog_nmf_engine}
     Profile        : ${workflow.profile}
    ${line}
    """.stripIndent()

    if (!run.any { k, v -> v })
        error "No stages selected. Check --from/--to and the run_* flags."

    // ---- Shared Quarto assets -------------------------------------------
    // Every notebook-rendering process takes the same template bundle. Built
    // once here rather than re-derived inside each of the eleven subworkflows,
    // which is what the standalone pipelines did.
    ch_template = Channel.fromPath(params.template, checkIfExists: true).collect()

    // The publication figure system travels with the Quarto assets, so every
    // notebook-rendering process receives it and its setup chunk can source it.
    // Without this the chunk is a guarded no-op and figures fall back to each
    // notebook's own ad-hoc styling.
    ch_viz = Channel.fromPath(params.viz_library, checkIfExists: true).collect()

    ch_page_config = ch_template
        .map { files -> files.find { it.toString().endsWith('.png') } }
        .combine(Channel.fromPath(params.page_config, checkIfExists: true).collect())
        .combine(ch_viz)
        .collect()

    // ---- Entry point resolution -----------------------------------------
    ch_seurat = Channel.empty()
    ch_gex    = Channel.empty()
    ch_vdj    = Channel.empty()
    ch_bam    = Channel.empty()

    // Any stage downstream of QC can be entered directly with an existing object.
    // A compartment object counts: entering the malignant arm with
    // --input_tumor_object should not also demand the pre-split object.
    if (!run['align'] && !run['qc']) {
        def entry_object = params.input_seurat_object
                        ?: params.input_nonmalignant_object
                        ?: params.input_tumor_object
        if (!entry_object)
            error "Starting from '${params.from}' requires an input object: " +
                  "--input_seurat_object, or --input_tumor_object / --input_nonmalignant_object " +
                  "to enter one arm directly."
        // Also becomes the base object the integration layer annotates.
        ch_seurat = Channel.value(file(entry_object, checkIfExists: true))
    }

    // ---- ALIGN -----------------------------------------------------------
    if (run['align']) {
        if (!params.input_samplesheet)
            error "Stage 'align' requires --input_samplesheet <PATH/TO/samplesheet.csv>."
        if (params.multi && params.demux)
            error "--multi and --demux are mutually exclusive. Pick one alignment mode."

        ALIGN(
            Channel.fromPath(params.input_samplesheet, checkIfExists: true),
            params.modality,
            params.genome,
            params.demux as boolean,
            params.multi as boolean
        )
        ch_gex = ALIGN.out.gex_outs.map { sample, files -> files }.flatten()
        ch_vdj = ALIGN.out.vdj_outs
        ch_bam = ALIGN.out.bam
    }
    else if (params.input_gex_matrices_path) {
        // Alignment already done elsewhere; start from the cellranger outputs.
        ch_gex = Channel.fromPath(params.input_gex_matrices_path, checkIfExists: true)
    }

    // ---- QC --------------------------------------------------------------
    if (run['qc']) {
        if (!params.input_exp_table)
            error "Stage 'qc' requires --input_exp_table <PATH/TO/metadata.csv>."
        QC(ch_gex, Channel.fromPath(params.input_exp_table, checkIfExists: true), ch_page_config)
        ch_seurat = QC.out.seurat_rds
        rep_qc = QC.out.qc_metrics.mix(QC.out.seurat_rds)
    }

    // ---- CLUSTER ---------------------------------------------------------
    if (run['cluster']) {
        CLUSTER(ch_seurat, ch_page_config)
        ch_seurat = CLUSTER.out.seurat_rds
        rep_cluster = CLUSTER.out.figures.flatten()
    }

    // ---- ANNOTATION ------------------------------------------------------
    ch_annotated = ch_seurat
    if (run['annotation']) {
        ANNOTATION(
            ch_seurat,
            Channel.fromPath(params.input_cell_mask),
            Channel.fromPath(params.annot_db, checkIfExists: true),
            Channel.fromPath(params.input_reference_object),
            ch_page_config
        )
        ch_annotated = ANNOTATION.out.seurat_rds
        ch_seurat    = ANNOTATION.out.seurat_rds
        rep_annot    = ANNOTATION.out.seurat_rds
    }

    // ---- CNV -------------------------------------------------------------
    // On the backbone, not a terminal branch: STRATIFY consumes its calls.
    ch_cnv_infercnv = Channel.empty()
    ch_cnv_scevan   = Channel.empty()
    ch_cnv_evidence = Channel.empty()

    // Artefacts feeding the per-stage smart reports.
    rep_qc = Channel.empty(); rep_cluster = Channel.empty(); rep_annot = Channel.empty()
    rep_cnv = Channel.empty(); rep_strat = Channel.empty(); rep_meta = Channel.empty()
    rep_bc = Channel.empty();  rep_traj = Channel.empty();  rep_cc = Channel.empty()
    rep_tme = Channel.empty(); rep_tcr = Channel.empty()

    if (run['cnv']) {
        CNV(ch_annotated, ch_page_config)
        ch_cnv_infercnv = CNV.out.infercnv
        ch_cnv_scevan   = CNV.out.scevan
        // Per-cell calls for the split; the raw caller directories stay behind.
        ch_cnv_evidence = CNV.out.infercnv_meta
            .mix(CNV.out.scevan_meta)
            .mix(CNV.out.copykat)
        rep_cnv = CNV.out.infercnv.mix(CNV.out.scevan).flatten()
    }

    // ---- NUMBAT (opt-in, allele-aware) -------------------------------------
    // Needs the BAM, so it can only run when ALIGN produced one (or when
    // --input_bam_path supplies them). Its calls join the same evidence pool.
    if (params.run_numbat) {
        ch_numbat_bam = ch_bam
        if (params.input_bam_path) {
            ch_numbat_bam = Channel.fromPath(params.input_bam_path, checkIfExists: true)
                .map { f -> tuple(f.parent.parent.name, f) }
        }

        NUMBAT(ch_annotated, ch_numbat_bam)
        ch_cnv_evidence = ch_cnv_evidence.mix(NUMBAT.out.calls)
    }

    // ---- STRATIFY ---------------------------------------------------------
    ch_tumor        = Channel.empty()
    ch_nonmalignant = Channel.empty()
    ch_assignment   = Channel.empty()
    ch_stratify_fig = Channel.empty()

    if (run['stratify']) {
        STRATIFY(ch_annotated, ch_cnv_evidence.collect().ifEmpty([]), ch_page_config)
        ch_tumor        = STRATIFY.out.tumor_rds
        ch_nonmalignant = STRATIFY.out.nonmalignant_rds
        ch_assignment   = STRATIFY.out.assignment
        ch_stratify_fig = STRATIFY.out.figures
        rep_strat = STRATIFY.out.assignment
            .mix(STRATIFY.out.data.flatten())
            .mix(STRATIFY.out.figures.flatten())
    }

    // Compartment objects can also be supplied directly, to enter either arm
    // without re-running the split.
    if (params.input_tumor_object) {
        ch_tumor = Channel.value(file(params.input_tumor_object, checkIfExists: true))
    }
    if (params.input_nonmalignant_object) {
        ch_nonmalignant = Channel.value(file(params.input_nonmalignant_object, checkIfExists: true))
    }
    // Falling back to the unstratified object keeps a partial run working, e.g.
    // `--from metaprog --input_seurat_object <obj>` with no split performed.
    if (!run['stratify'] && !params.input_tumor_object)        ch_tumor        = ch_seurat
    if (!run['stratify'] && !params.input_nonmalignant_object) ch_nonmalignant = ch_seurat

    // ---- Malignant arm ----------------------------------------------------
    ch_metaprog = Channel.empty()
    if (run['metaprog']) {
        METAPROG(ch_tumor.first())
        ch_metaprog = METAPROG.out.post_data
        rep_meta = METAPROG.out.post_data.flatten().mix(METAPROG.out.fits.flatten())
    }

    // ---- Non-malignant arm ------------------------------------------------
    // Batch correction is applied here rather than on the backbone: correcting
    // across a mixed object would regress out genuine tumour heterogeneity
    // along with the batch effect.
    ch_corrected  = ch_nonmalignant
    ch_trajectory = Channel.empty()
    ch_liana      = Channel.empty()
    ch_cellchat   = Channel.empty()
    ch_tme        = Channel.empty()

    if (run['batchcorrect']) {
        BATCH_CORRECTION(ch_nonmalignant, ch_page_config)
        ch_corrected = BATCH_CORRECTION.out.seurat_rds
        rep_bc = BATCH_CORRECTION.out.data.flatten().mix(BATCH_CORRECTION.out.figures.flatten())
    }

    if (run['trajectory']) {
        TRAJECTORY(ch_corrected, ch_page_config)
        ch_trajectory = TRAJECTORY.out.rds
        rep_traj = TRAJECTORY.out.rds.flatten()
    }

    if (run['cellcomm']) {
        CELL_COMMUNICATION(ch_corrected.first())
        ch_liana    = CELL_COMMUNICATION.out.liana_csv
        ch_cellchat = CELL_COMMUNICATION.out.cellchat_rds
        rep_cc = CELL_COMMUNICATION.out.liana_csv.mix(CELL_COMMUNICATION.out.nichenet.flatten())
    }

    if (run['tme']) {
        TME(ch_corrected, ch_page_config)
        ch_tme = TME.out.data
        rep_tme = TME.out.data.flatten().mix(TME.out.figures.flatten())
    }

    // ---- Independent of the split -----------------------------------------
    ch_tcr_summary = Channel.empty()
    if (run['tcr']) {
        ch_tcr_sheet = Channel.fromPath(params.input_tcr_sample_sheet ?: "${projectDir}/assets/NO_FILE")
        // T cells are non-malignant, so the corrected compartment is the better
        // input when it exists; --tcr_input annotated keeps the older behaviour.
        ch_tcr_object = (params.tcr_input == 'nonmalignant') ? ch_corrected : ch_annotated
        TCR(ch_tcr_object.first(), ch_tcr_sheet)
        ch_tcr_summary = TCR.out.summary
        rep_tcr = TCR.out.summary.flatten()
    }

    if (run['microbiome']) {
        // The eight microbiome modules its subworkflow includes were never
        // present in the source repository. See docs/MICROBIOME.md.
        error "Stage 'microbiome' is not implemented: modules/local/microbiome/* is " +
              "absent from the source pipelines. Re-run with --run_microbiome false."
    }

    // ---- Per-stage smart reports ------------------------------------------
    // Each stage's published tables and figures become a summary-first report:
    // computed verdict, findings ranked worst-first, searchable tables, linked
    // interactive views, figure browser and the offline agent. Built from the
    // published outputs, so the validated analysis notebooks are untouched.
    if (params.stage_reports) {
        // Stages emit the same file on more than one channel (STRATIFY publishes
        // the assignment table via both `assignment` and `data`), which would be
        // staged twice and collide. Deduplicate by basename before packing.
        def stagePack = { ch, id, title ->
            ch.collect().ifEmpty([]).map { f ->
                def files = (f instanceof List) ? f : [f]
                tuple(id, title, files.unique { it.name })
            }
        }

        ch_stage_outputs = stagePack(rep_qc,      'qc',             'Quality control')
            .mix( stagePack(rep_cluster, 'cluster',        'Clustering') )
            .mix( stagePack(rep_annot,   'annotation',     'Cell-type annotation') )
            .mix( stagePack(rep_cnv,     'cnv',            'Copy number') )
            .mix( stagePack(rep_strat,   'stratify',       'Compartment split') )
            .mix( stagePack(rep_meta,    'metaprog',       'Tumour heterogeneity') )
            .mix( stagePack(rep_bc,      'batchcorrection','Batch correction') )
            .mix( stagePack(rep_traj,    'trajectory',     'Trajectory') )
            .mix( stagePack(rep_cc,      'cellcomm',       'Cell communication') )
            .mix( stagePack(rep_tme,     'tme',            'Microenvironment') )
            .mix( stagePack(rep_tcr,     'tcr',            'TCR repertoire') )

        STAGE_REPORTS(ch_stage_outputs, ch_page_config)
    }

    // ---- Cross-stage integration and master report ------------------------
    // The master object is built on the annotated object, which is the last
    // point at which both compartments are present in one place. The two arms'
    // per-cell results are joined back onto it.
    if (run['integration']) {
        INTEGRATION(
            ch_annotated,
            ch_annotated,
            ch_cnv_infercnv.collect().ifEmpty([]),
            ch_cnv_scevan.collect().ifEmpty([]),
            ch_metaprog.collect().ifEmpty([]),
            ch_trajectory.collect().ifEmpty([]),
            ch_liana.collect().ifEmpty([]),
            ch_cellchat.collect().ifEmpty([]),
            ch_tcr_summary.collect().ifEmpty([]),
            ch_assignment.mix(ch_tme).collect().ifEmpty([]),
            ch_page_config
        )
    }
}


workflow.onComplete {
    def line = '=' * 74
    log.info(
        workflow.success
            ? """
            ${line}
             SCRATCH-E2E completed
            ${line}
             Duration : ${workflow.duration}
             Results  : ${params.outdir}/${params.project_name}
             Report   : ${params.outdir}/${params.project_name}/master_report/index.html
             Run info : ${params.outdir}/pipeline_info
            ${line}
            """.stripIndent()
            : """
            ${line}
             SCRATCH-E2E failed
            ${line}
             Exit status : ${workflow.exitStatus}
             Error       : ${workflow.errorMessage}
             Resume      : nextflow run . -resume
            ${line}
            """.stripIndent()
    )
}
