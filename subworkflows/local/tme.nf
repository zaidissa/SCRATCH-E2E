/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    TME — microenvironment characterisation of the non-malignant compartment
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { TME_CHARACTERISATION } from '../../modules/local/tme/main.nf'

workflow TME {

    take:
        ch_nonmalignant_object
        ch_page_config

    main:
        ch_data    = Channel.empty()
        ch_figures = Channel.empty()

        if (!params.skip_tme) {
            ch_notebook = Channel.fromPath(
                "${projectDir}/modules/local/tme/tme_characterisation.qmd", checkIfExists: true)

            TME_CHARACTERISATION(ch_nonmalignant_object.first(), ch_notebook, ch_page_config)
            ch_data    = TME_CHARACTERISATION.out.data
            ch_figures = TME_CHARACTERISATION.out.figures
        }

    emit:
        data    = ch_data
        figures = ch_figures
}
