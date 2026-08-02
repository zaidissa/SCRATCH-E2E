/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    STRATIFY — split into malignant and non-malignant compartments
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Consumes the annotated object plus whatever CNV callers produced, and emits
    the two compartment objects that drive the rest of the pipeline.
----------------------------------------------------------------------------------------
*/

include { STRATIFY_COMPARTMENTS } from '../../modules/local/stratify/main.nf'

workflow STRATIFY {

    take:
        ch_annotated_object
        ch_cnv_evidence     // collected CNV meta objects; may be empty
        ch_page_config

    main:
        ch_notebook = Channel.fromPath(
            "${projectDir}/modules/local/stratify/stratify_compartments.qmd", checkIfExists: true)

        STRATIFY_COMPARTMENTS(
            ch_annotated_object.first(),
            ch_notebook,
            ch_page_config,
            ch_cnv_evidence
        )

    emit:
        tumor_rds        = STRATIFY_COMPARTMENTS.out.tumor_rds
        nonmalignant_rds = STRATIFY_COMPARTMENTS.out.nonmalignant_rds
        assignment       = STRATIFY_COMPARTMENTS.out.assignment
        data             = STRATIFY_COMPARTMENTS.out.data
        figures          = STRATIFY_COMPARTMENTS.out.figures
}
