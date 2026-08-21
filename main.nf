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
            def ov = params["run_${name}"]
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

    // TCR is "independent of the split", but it is NOT independent of having any
    // VDJ data. Fanning out from a backbone stage used to switch it on
    // unconditionally, so a GEX-only run reached VDJ_QC with the NO_FILE
    // placeholder as its sample sheet and died on "missing required columns".
    // A stage with no possible input should not be scheduled at all.
    if (enabled['tcr']) {
        def has_contigs = params.input_vdj_contigs
        def sheet       = params.input_tcr_sample_sheet?.toString()
        def has_sheet   = sheet && !sheet.endsWith('NO_FILE')
        // ALIGN can produce VDJ itself from a samplesheet carrying TCR rows.
        def align_makes_vdj = enabled['align']
        if (!(has_contigs || has_sheet || align_makes_vdj)) {
            enabled['tcr'] = false
            log.info "Stage 'tcr' skipped: no VDJ input " +
                     "(--input_vdj_contigs / --input_tcr_sample_sheet unset and ALIGN not running). " +
                     "Pass --run_tcr true to force it."
        }
    }

    // ...and the mirror case. Entering AT batchcorrect with --to end means
    // "batch-correct, then everything that consumes the corrected object". The
    // rule above only handled the opposite direction, so `--from batchcorrect`
    // ran batch correction and then stopped, silently skipping the three stages
    // it exists to feed.
    if (entering_at_branch && from == 'batchcorrect' && to == 'end') {
        ['trajectory', 'cellcomm', 'tme'].each { enabled[it] = true }
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
        def override = params["run_${name}"]
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

// Declared as a top-level FUNCTION, not a closure variable inside the workflow
// body: Nextflow 26.x's strict parser reports "`stagePack` is not defined" for
// the latter. Functions are declarations and resolve normally.
def stagePack(ch, id, title) {
            ch.collect().ifEmpty([]).map { f ->
                def files = (f instanceof List) ? f : [f]
                files = files.findAll { !(it.toString().contains('_figlibs')) }
                tuple(id, title, files.unique { it.name })
            }
}

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
        // Entering at/after cluster: no BPCells store travels with a plain RDS.
        // The notebooks fall back to the object's own layers when the staged
        // store is the placeholder.
        ch_bpcells_store = Channel.value(file("${projectDir}/assets/NO_FILE_bpcells"))
    }

    // ---- Artefacts feeding the per-stage smart reports ---------------------
    // Declared HERE, before the first stage that fills one. They used to sit
    // between ANNOTATION and CNV, which meant `rep_qc`, `rep_cluster` and
    // `rep_annot` were assigned by their stages and then immediately reset to
    // Channel.empty() — so the QC, clustering and annotation stage reports were
    // built from nothing, every run, however much those stages published.
    rep_qc = Channel.empty(); rep_cluster = Channel.empty(); rep_annot = Channel.empty()
    rep_cnv = Channel.empty(); rep_strat = Channel.empty(); rep_meta = Channel.empty()
    rep_bc = Channel.empty();  rep_traj = Channel.empty();  rep_cc = Channel.empty()
    rep_tme = Channel.empty(); rep_tcr = Channel.empty()

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
        ch_bpcells_store = QC.out.bpcells_store
        rep_qc = QC.out.qc_metrics.mix(QC.out.seurat_rds)
    }

    // ---- CLUSTER ---------------------------------------------------------
    if (run['cluster']) {
        CLUSTER(ch_seurat, ch_bpcells_store, ch_page_config)
        ch_seurat = CLUSTER.out.seurat_rds
        rep_cluster = CLUSTER.out.figures.flatten()
    }

    // ---- ANNOTATION ------------------------------------------------------
    ch_annotated = ch_seurat
    // Per-cell annotation tables for the integration join. The master object is
    // built on the AZIMUTH output, so without this it carries azimuth_labels and
    // nothing from scType — this pipeline's own annotator.
    ch_annotation_tables = Channel.empty()
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
        ch_annotation_tables = ANNOTATION.out.concordance
    }

    // ---- CNV -------------------------------------------------------------
    // On the backbone, not a terminal branch: STRATIFY consumes its calls.
    // inferCNV and Numbat both live inside the CNV subworkflow and start from the
    // same annotated object, so they run CONCURRENTLY; CNV_CONCORDANCE compares
    // them at the end of the stage.
    ch_cnv_infercnv    = Channel.empty()
    ch_cnv_evidence    = Channel.empty()
    ch_cnv_concordance = Channel.empty()
    ch_cnv_object      = Channel.empty()

    // Numbat needs the BAM. Resolved here because --input_bam_path is a pipeline
    // input, not a CNV-stage concern; ALIGN's own emission is already a
    // (sample, bam, bai) tuple.
    //
    // Numbat is on by default but cannot run without a BAM. Switching it off with
    // a warning — rather than letting it reach the pileup and die on "no sample
    // matched between the supplied BAMs and the annotated object" — is the same
    // treatment TCR gets when there is no VDJ input.
    def numbat_has_bam = params.input_bam_path || run['align']
    if (params.run_numbat && !numbat_has_bam) {
        log.warn "Numbat skipped: no BAM (--input_bam_path unset and ALIGN not running). " +
                 "It phases reads at heterozygous SNPs, so it needs the aligned BAM; " +
                 "inferCNV does not. CNV will run inferCNV alone and the caller " +
                 "comparison will report that only one caller was available."
    }

    ch_numbat_bam = Channel.empty()
    if (params.run_numbat && numbat_has_bam) {
        ch_numbat_bam = ch_bam
        if (params.input_bam_path) {
            // Carry the .bai alongside the BAM. Without it the pileup re-indexes an
            // 80 GB file on every attempt, which costs more than the pileup.
            ch_numbat_bam = Channel.fromPath(params.input_bam_path, checkIfExists: true)
                .map { f ->
                    def idx = file("${f}.bai")
                    if (!idx.exists()) idx = file("${f.toString().replaceAll(/\.bam$/, '.bai')}")
                    tuple(f.parent.parent.name, f,
                          idx.exists() ? idx : file("${projectDir}/assets/NO_FILE_bai"))
                    // A DEDICATED placeholder: reusing NO_FILE_gmap here would put two
                    // inputs with the same basename in the task dir and recreate the
                    // "input file name collision" failure.
                }
        }
    }

    if (run['cnv']) {
        CNV(ch_annotated, ch_numbat_bam, (params.run_numbat && numbat_has_bam), ch_page_config)
        ch_cnv_infercnv = CNV.out.infercnv
        // Per-cell calls for the split; the raw caller directories stay behind.
        ch_cnv_evidence = CNV.out.infercnv_meta
            .mix(CNV.out.numbat)
            .mix(CNV.out.copykat)
        ch_cnv_concordance = CNV.out.concordance_csv

        // The annotated object with BOTH callers' per-cell calls written into it.
        // This is what keeps the backbone cumulative: STRATIFY reads this rather
        // than the pre-CNV object, so nothing upstream is dropped on the way.
        ch_cnv_object = CNV.out.seurat_rds

        rep_cnv = CNV.out.infercnv
            .mix(CNV.out.numbat)
            .mix(CNV.out.concordance)
            .mix(CNV.out.concordance_csv)
            .flatten()
    }

    // ---- STRATIFY ---------------------------------------------------------
    ch_tumor        = Channel.empty()
    ch_nonmalignant = Channel.empty()
    ch_assignment   = Channel.empty()
    ch_stratify_fig = Channel.empty()

    if (run['stratify']) {
        // Feed the CNV-enriched object when CNV ran, so the compartment objects
        // inherit both callers' per-cell calls instead of STRATIFY re-deriving
        // them and dropping the rest. Falls back to the annotated object when the
        // stage is entered directly.
        def ch_strat_in = run['cnv'] ? ch_cnv_object : ch_annotated
        STRATIFY(ch_strat_in, ch_cnv_evidence.collect().ifEmpty([]), ch_page_config)
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
    // Falls back through the arm: trajectory's object if it ran, else the
    // corrected object, else the raw compartment. Whichever it ends up being,
    // this is the object that represents the microenvironment compartment.
    ch_nonmalignant_final = ch_nonmalignant
    ch_trajectory = Channel.empty()
    ch_liana      = Channel.empty()
    ch_cellchat   = Channel.empty()
    ch_tme        = Channel.empty()

    if (run['batchcorrect']) {
        BATCH_CORRECTION(ch_nonmalignant, ch_page_config)
        ch_corrected = BATCH_CORRECTION.out.seurat_rds
        ch_nonmalignant_final = BATCH_CORRECTION.out.seurat_rds
        rep_bc = BATCH_CORRECTION.out.data.flatten().mix(BATCH_CORRECTION.out.figures.flatten())
    }

    if (run['trajectory']) {
        TRAJECTORY(ch_corrected, ch_page_config)
        ch_trajectory = TRAJECTORY.out.rds
        rep_traj = TRAJECTORY.out.rds.flatten()

        // Object 3 of the three the pipeline produces: the NON-MALIGNANT object,
        // now carrying its inherited history (QC, clustering, annotation, both CNV
        // callers, the compartment call) plus the batch-corrected reduction and
        // per-cell pseudotime.
        //
        // Cell communication and TME results are deliberately NOT written onto it:
        // LIANA/CellChat produce per-cell-type-PAIR interactions and TME produces
        // per-cell-type composition and signature scores (89 rows, not 19,657).
        // Neither is a per-cell quantity, so attaching them would invent a value
        // for each cell that the analysis never computed.
        ch_nonmalignant_final = TRAJECTORY.out.seurat_rds.ifEmpty(null).filter { it != null }
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
        // Also drop the interactive figures' JavaScript payload. Each figure
        // directory carries a `_figlibs/` holding plotly.js, crosstalk.js,
        // jquery.js and friends; staging several of those into one flat
        // directory collides on every filename, and `unique { it.name }` cannot
        // resolve it because they are genuinely different files that merely
        // share a name. The integration stage reads analysis outputs, not the
        // widget runtime, so these have no business being staged at all.
        //
        // Inlined deliberately: a separate `def` closure referenced from inside
        // this one does not resolve reliably here — the same trap that broke
        // publishStage() and the container helper.

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
        // INTEGRATE_OBJECTS stages each of these with `stageAs:
        // 'stage_inputs/<stage>/*'`, which FLATTENS the directory tree. The
        // interactive figures write a `_figlibs/` beside every figure directory
        // (plotly.js, crosstalk.js, jquery.js ...), and trajectory alone has
        // four of them under data/ — 64 files that all flatten onto the same
        // handful of names and abort the task with a filename collision.
        //
        // These are the widget runtime, not analysis outputs, so drop them
        // before staging. Applied to every stage channel, not just trajectory,
        // because any stage that saves a figure now carries the same payload.
        // Inlined rather than factored into a `def` closure: a helper closure
        // referenced from inside another closure does not resolve reliably here.
        // `.flatten()` BEFORE `.filter{}` is load-bearing, not tidying.
        //
        // An output declared `path "data/**"` emits ONE channel item holding a
        // LIST of every matched file. Calling `.filter{}` on that emission tests
        // the LIST's toString(), which is the concatenation of all the paths — so
        // a single `_figlibs` entry anywhere in it discarded the ENTIRE stage.
        // That is what happened to TRAJECTORY: it ran, published four figlib
        // directories among its data, and reached INTEGRATE_OBJECTS with nothing
        // at all — `stage_inputs/trajectory/` was never even created, while the
        // trajectory STAGE REPORT got its files because it flattens first.
        // Flattening makes the filter a per-file test, which is what it meant.
        INTEGRATION(
            ch_annotated,
            ch_annotated,
            ch_cnv_infercnv.flatten().filter { !it.toString().contains('_figlibs') }.collect().ifEmpty([]),
            ch_metaprog.flatten().filter     { !it.toString().contains('_figlibs') }.collect().ifEmpty([]),
            ch_trajectory.flatten().filter   { !it.toString().contains('_figlibs') }.collect().ifEmpty([]),
            ch_liana.flatten().filter        { !it.toString().contains('_figlibs') }.collect().ifEmpty([]),
            ch_cellchat.flatten().filter     { !it.toString().contains('_figlibs') }.collect().ifEmpty([]),
            ch_tcr_summary.flatten().filter  { !it.toString().contains('_figlibs') }.collect().ifEmpty([]),
            ch_assignment.mix(ch_tme).flatten().filter { !it.toString().contains('_figlibs') }.collect().ifEmpty([]),
            ch_annotation_tables.flatten().filter { !it.toString().contains('_figlibs') }.collect().ifEmpty([]),
            ch_page_config
        )
    }
}
