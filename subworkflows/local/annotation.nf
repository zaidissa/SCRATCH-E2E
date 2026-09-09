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

include { asBool } from '../../lib/booleans.nf'
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

        // A typo here would silently fall through to scType and look like it
        // worked, so name the valid values instead.
        if (!(params.annot_primary?.toString()?.toLowerCase() in ['sctype', 'azimuth']))
            error "--annot_primary '${params.annot_primary}' is not valid. Use 'sctype' or 'azimuth'."


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

        if (!asBool(params.skip_celltypist)) {
            CELLTYPIST_ANNOTATION(ch_notebook_celltypist, SCEASY_CONVERTER_ONE.out.project_rds, ch_page_config)
            ch_celltypist_obs = CELLTYPIST_ANNOTATION.out.csv_file
        }

        // scType's marker database is human-only in this build.
        if (params.organism_capitalised == 'Human' && !asBool(params.skip_sctype)) {

            SCYTPE_MAJOR_ANNOTATION(ch_notebook_sctype_mj, ch_working, ch_database, ch_page_config)
            ch_major_object = SCYTPE_MAJOR_ANNOTATION.out.seurat_rds

            // Populations that get a second, state-level annotation pass.
            // The exclusion list is a property of the TUMOUR TYPE, not of the
            // pipeline: it was hardcoded to Unknown|Epithelial|Fibroblast|NK_Cells,
            // which is an ovarian-carcinoma assumption. A sarcoma wants fibroblast
            // states; a glioma has no epithelial compartment to exclude.
            def state_exclude = params.annot_state_exclude
                                    .toString()
                                    .split('[;,|]')
                                    .collect { it.trim() }
                                    .findAll { it }
                                    .join('|')
            ch_major_list = SCYTPE_MAJOR_ANNOTATION.out.major_list
                .splitText()
                .map    { it.split(':') }
                // `=~ state_exclude` with a plain String, NOT /${state_exclude}/:
                // Nextflow 26.x's parser rejects an interpolated slashy regex
                // with "Unexpected input: '/'". Groovy coerces the String itself.
                .filter { state_exclude ? !(it[0] =~ state_exclude) : true }
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
        if (has_reference && !asBool(params.skip_azimuth)) {
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
        if (!asBool(params.skip_annotation_concordance)) {
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

        // Which annotator's object goes downstream. Both still run and both are
        // reported; this only decides whose labels CNV, stratification and cell
        // communication read.
        //
        // `first()` on a mixed channel would be non-deterministic, so the
        // fallback is resolved by concatenation and taking the head — the
        // preferred annotator is simply concatenated first, and an annotator that
        // did not run contributes an empty channel rather than a null.
        //
        // Previously Azimuth always won. That made the choice depend on whether
        // --input_reference_object happened to be supplied, rather than on
        // anything anyone decided.
        ch_annotated = (params.annot_primary?.toString()?.toLowerCase() == 'azimuth'
                            ? ch_azimuth_rds.concat(ch_sctype_rds)
                            : ch_sctype_rds.concat(ch_azimuth_rds))
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
