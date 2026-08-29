/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    CNV — inferCNV and Numbat (CopyKAT optional)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    The callers are independent and run in parallel — they share one input
    (the annotated object) and nothing else, so Nextflow schedules them
    concurrently. Each is individually skippable; per-tool `ext.when` in
    conf/modules.config does the gating so a disabled tool leaves an empty
    channel rather than an error.

    NUMBAT lives here rather than in main.nf so the two primary callers are
    siblings in one place: same input, same stage, compared at the end of it by
    CNV_CONCORDANCE.

    SCEVAN was REMOVED (2026-08-20). It was expression-only like inferCNV, so it
    added no independent evidence; it was slow enough to dominate a run (7 h
    under emulation); its `classifyTumorCells` breaks above ~5,000 cells/sample;
    and it had a silent-failure history — a doubled sample prefix that matched 0
    of 3000 cells while still being counted as a caller in the "consensus".
----------------------------------------------------------------------------------------
*/

include { INFERCNV        } from '../../modules/local/infercnv/main.nf'
include { COPYKAT_PROCESS } from '../../modules/local/copykat/main.nf'
include { NUMBAT          } from './numbat.nf'
include { CNV_CONCORDANCE } from '../../modules/local/cnv_concordance/main.nf'

workflow CNV {

    take:
        ch_seurat_object
        ch_bam            // tuple(sample_id, bam, bai); empty when there is none
        run_numbat        // val boolean: params.run_numbat AND a BAM actually exists
        ch_page_config

    main:
        ch_infercnv = Channel.empty()
        ch_copykat  = Channel.empty()
        ch_numbat   = Channel.empty()
        // Per-SEGMENT CNV events, as opposed to per-cell calls. Both go to the
        // concordance module: "how many cells did Numbat call" and "how many
        // CNVs did Numbat find" are different questions with answers two orders
        // of magnitude apart, and a report that can only answer the first
        // invites the second being answered wrongly.
        ch_numbat_segments = Channel.empty()

        // Per-cell calls, kept separate from the raw caller output directories.
        // These are what the stratification step consumes.
        ch_infercnv_meta = Channel.empty()

        if (!params.skip_infercnv) {
            INFERCNV(ch_seurat_object,
                     Channel.fromPath(params.notebook_infercnv, checkIfExists: true),
                     ch_page_config)
            ch_infercnv      = INFERCNV.out.results
            ch_infercnv_meta = INFERCNV.out.meta_object
        }

        // Allele-aware, and the only caller that can see copy-neutral LOH.
        // Started from the same input as inferCNV above, so the two run together.
        //
        // Gated on the resolved boolean, NOT on params.run_numbat: with the
        // parameter on by default and no BAM present, the subworkflow would still
        // be invoked, its BAM channel would be empty, and the join's `ifEmpty`
        // would abort the whole run with "no sample matched between the supplied
        // BAMs" — an error about mismatched sample names when the real situation
        // is that there are no BAMs at all.
        if (run_numbat) {
            NUMBAT(ch_seurat_object, ch_bam)
            ch_numbat = NUMBAT.out.calls

            // Stage the per-sample DIRECTORY, not the segment files themselves.
            //
            // Numbat names its consensus files segs_consensus_<iteration>.tsv and
            // distinguishes samples only by the enclosing directory, so collecting
            // the files flattens two different samples onto one basename and the
            // task dies with "input file name collision". That stayed hidden until
            // now because only HRS371755 ever produced segments — HRS371760 is
            // tumour-poor and Numbat correctly returned nothing for it. With
            // ambient RNA removed it produces 24, and the collision appeared.
            //
            // The notebook already expects this layout: it recovers the sample with
            // `sub("^.*/([^/]+)/segs_consensus_.*$", ...)`, i.e. from the parent
            // directory, which a flattened stage had also silently broken.
            // Directories stage as symlinks, so the 3.4 GB behind them costs nothing.
            ch_numbat_segments = NUMBAT.out.segments.flatten().map { it.parent }.unique()
        }

        if (!params.skip_copykat) {
            COPYKAT_PROCESS(ch_seurat_object,
                            Channel.fromPath(params.notebook_copykat, checkIfExists: true),
                            ch_page_config)
            ch_copykat = COPYKAT_PROCESS.out.results
        }

        // ---- Compare the callers before anything downstream consumes them ----
        // Both callers reduce to a per-cell malignant/non-malignant call, and
        // until now nothing ever put those two side by side: STRATIFY merged them
        // into a consensus and the disagreement disappeared into it. This runs
        // only when both produced output, and it is read-only — it changes no
        // call, it just reports where they differ.
        ch_concordance      = Channel.empty()
        ch_concordance_csv  = Channel.empty()
        // Object 1 of the three the pipeline produces: the PRE-STRATIFICATION
        // object, carrying QC, clustering, annotation and now both callers'
        // per-cell calls. STRATIFY subsets it into objects 2 (malignant) and 3
        // (non-malignant), so whatever is attached here is inherited by both.
        // Defaults to the stage's input so the chain never breaks when the
        // comparison is skipped — it is then simply the annotated object.
        ch_cnv_object       = ch_seurat_object
        if (!params.skip_cnv_concordance) {
            CNV_CONCORDANCE(
                ch_seurat_object,
                ch_infercnv_meta.collect().ifEmpty([]),
                ch_numbat.mix(ch_numbat_segments).collect().ifEmpty([]),
                Channel.fromPath(params.notebook_cnv_concordance, checkIfExists: true),
                ch_page_config
            )
            ch_concordance     = CNV_CONCORDANCE.out.figures
            ch_concordance_csv = CNV_CONCORDANCE.out.comparison
            ch_cnv_object      = CNV_CONCORDANCE.out.seurat_rds
        }

    emit:
        seurat_rds      = ch_cnv_object
        infercnv        = ch_infercnv
        copykat         = ch_copykat
        infercnv_meta   = ch_infercnv_meta
        numbat          = ch_numbat
        concordance     = ch_concordance
        concordance_csv = ch_concordance_csv
}
