/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ANNOTATION — CellTypist, scType (major -> state -> aggregate), Azimuth
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    Changed on merge: the standalone subworkflow ended with

        emit:
            ch_dumps = Channel.empty()

    so it produced no chainable object at all — every downstream pipeline had to
    be pointed by hand at a file path on disk. It now emits the annotated object,
    preferring Azimuth (which is what batch correction and the metaprogram stage
    were configured to consume, via `azimuth_labels`) and falling back through
    the scType aggregate to the unannotated input.
----------------------------------------------------------------------------------------
*/

include { HELPER_SEURAT_SUBSET                            } from '../../modules/local/helpers/subset/main.nf'
include { HELPER_SCEASY_CONVERTER as SCEASY_CONVERTER_ONE } from '../../modules/local/helpers/convert/main.nf'
include { CELLTYPIST_ANNOTATION                           } from '../../modules/local/celltypist/main.nf'
include { SCYTPE_MAJOR_ANNOTATION                         } from '../../modules/local/sctype/major/main.nf'
include { SCYTPE_STATE_ANNOTATION                         } from '../../modules/local/sctype/state/main.nf'
include { SCYTPE_AGGREGATE_ANNOTATION                     } from '../../modules/local/sctype/aggregate/main.nf'
include { AZIMUTH_ANNOTATION                              } from '../../modules/local/azimuth/main.nf'
include { ANNOTATION_CONCORDANCE                          } from '../../modules/local/annotation_concordance/main.nf'
include { HELPER_SCEASY_CONVERTER as SCEASY_CONVERTER_TWO } from '../../modules/local/helpers/convert/main.nf'

workflow ANNOTATION {

    take:
        ch_seurat_object
        ch_cell_malignancy
        ch_database
        ch_reference_object
        ch_page_config

    main:

        ch_notebook_celltypist = Channel.fromPath(params.notebook_celltypist,   checkIfExists: true)
        ch_notebook_sctype_mj  = Channel.fromPath(params.notebook_sctype_major, checkIfExists: true)
        ch_notebook_sctype_st  = Channel.fromPath(params.notebook_sctype_state, checkIfExists: true).collect()
        ch_notebook_sctype_ag  = Channel.fromPath(params.notebook_sctype_agg,   checkIfExists: true).collect()
        ch_notebook_azimuth    = Channel.fromPath(params.notebook_azimuth,      checkIfExists: true)

        // Optional malignancy-based subsetting.
        def has_mask = params.input_cell_mask && !params.input_cell_mask.toString().contains('NO_FILE')
        if (has_mask) {
            ch_working = HELPER_SEURAT_SUBSET(ch_seurat_object, ch_cell_malignancy)
        } else {
            ch_working = ch_seurat_object
        }

        SCEASY_CONVERTER_ONE(ch_working)

        // Track the best available annotated object as stages succeed.
        ch_annotated = ch_working
        ch_sctype_rds = Channel.empty()
        ch_azimuth_rds = Channel.empty()
        ch_celltypist_obs = Channel.empty()

        if (!params.skip_celltypist) {
            CELLTYPIST_ANNOTATION(ch_notebook_celltypist, SCEASY_CONVERTER_ONE.out.project_rds, ch_page_config)
            ch_celltypist_obs = CELLTYPIST_ANNOTATION.out.csv_file
        }

        // scType's marker database is human-only in this build.
        if (params.organism_capitalised == 'Human' && !params.skip_sctype) {

            SCYTPE_MAJOR_ANNOTATION(ch_notebook_sctype_mj, ch_working, ch_database, ch_page_config)
            ch_major_object = SCYTPE_MAJOR_ANNOTATION.out.seurat_rds

            // Populations that get a second, state-level annotation pass.
            ch_major_list = SCYTPE_MAJOR_ANNOTATION.out.major_list
                .splitText()
                .map    { it.split(':') }
                .filter { !(it[0] =~ /Unknown|Epithelial|Fibroblast|NK_Cells/) }
                .map    { it[0].trim() }

            SCYTPE_STATE_ANNOTATION(ch_notebook_sctype_st, ch_major_object, ch_database, ch_major_list, ch_page_config)

            SCYTPE_AGGREGATE_ANNOTATION(
                ch_notebook_sctype_ag,
                ch_major_object,
                SCYTPE_STATE_ANNOTATION.out.annotation.collect(),
                ch_page_config
            )

            ch_sctype_rds = SCYTPE_AGGREGATE_ANNOTATION.out.seurat_rds
            SCEASY_CONVERTER_TWO(ch_sctype_rds)
        }

        // Azimuth requires a reference object.
        def has_reference = params.input_reference_object && !params.input_reference_object.toString().contains('NO_FILE')
        if (has_reference && !params.skip_azimuth) {
            AZIMUTH_ANNOTATION(ch_notebook_azimuth, ch_working, ch_reference_object, ch_page_config)
            ch_azimuth_rds = AZIMUTH_ANNOTATION.out.seurat_rds
        }

        // The three annotators above run in PARALLEL on the same input and never
        // meet, so nothing in the pipeline ever checked whether they agree. This
        // joins their labels per cell and reports the disagreement rate. Optional
        // inputs are filled with the NO_FILE placeholder; the notebook renders a
        // "not enough annotators" finding when fewer than two are present.
        // Initialised outside the guard: an emit that only exists on one branch
        // is an "is not defined" at parse time when the stage is skipped.
        ch_concordance = Channel.empty()
        if (!params.skip_annotation_concordance) {
            // A DEDICATED placeholder per input, not one shared NO_FILE.
            // Nextflow stages inputs by basename, so two absent annotators both
            // resolving to `NO_FILE` put two files of that name in the task
            // directory and the run dies with "input file name collision". That
            // needs only two of the three to be missing, which is now the DEFAULT
            // combination: CellTypist is off unless asked for, and Azimuth is off
            // unless --input_reference_object is given. The same trap is recorded
            // against NO_FILE_gmap in numbat.nf.
            def nf_ct  = file("${projectDir}/assets/NO_FILE_celltypist")
            def nf_st  = file("${projectDir}/assets/NO_FILE_sctype")
            def nf_az  = file("${projectDir}/assets/NO_FILE_azimuth")
            ANNOTATION_CONCORDANCE(
                Channel.fromPath(params.notebook_annotation_concordance, checkIfExists: true),
                ch_celltypist_obs.ifEmpty(nf_ct),
                ch_sctype_rds.ifEmpty(nf_st),
                ch_azimuth_rds.ifEmpty(nf_az),
                ch_page_config
            )
            ch_concordance = ANNOTATION_CONCORDANCE.out.concordance
        }

        // Preference order: Azimuth -> scType aggregate -> unannotated input.
        // `first()` on a mixed channel would be non-deterministic, so the
        // fallback is resolved by concatenation and taking the head.
        ch_annotated = ch_azimuth_rds
            .concat(ch_sctype_rds)
            .concat(ch_working)
            .first()

    emit:
        seurat_rds  = ch_annotated
        azimuth_rds = ch_azimuth_rds
        sctype_rds  = ch_sctype_rds
        // Per-cell table of every annotator side by side. INTEGRATION joins this
        // so the master object is not Azimuth-only.
        concordance = ch_concordance
}
