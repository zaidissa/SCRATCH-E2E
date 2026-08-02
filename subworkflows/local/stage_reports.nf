/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    STAGE_REPORTS — build the smart report for every stage that produced output
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { STAGE_REPORT } from '../../modules/local/report/stage_report.nf'

workflow STAGE_REPORTS {

    take:
        ch_stage_outputs   // tuple(stage, stage_title, [files])
        ch_page_config

    main:
        ch_nb = Channel.fromPath("${projectDir}/modules/local/report/stage_report.qmd", checkIfExists: true)

        // Only stages that actually published something get a report.
        ch_in = ch_stage_outputs.filter { stage, title, files ->
            files && (files instanceof List ? files.size() > 0 : true)
        }

        STAGE_REPORT(ch_in, ch_nb.first(), ch_page_config)

    emit:
        reports  = STAGE_REPORT.out.report
        findings = STAGE_REPORT.out.findings
}
