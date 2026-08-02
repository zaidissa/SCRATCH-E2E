/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    NUMBAT — allele-aware CNV calling, opt-in
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Joins per-sample BAMs (from ALIGN) to per-sample expression (from the
    annotated object) and emits per-cell calls for STRATIFY's evidence pool.
----------------------------------------------------------------------------------------
*/

include { NUMBAT_PREP   } from '../../modules/local/numbat/main.nf'
include { NUMBAT_PILEUP } from '../../modules/local/numbat/main.nf'
include { NUMBAT_RUN    } from '../../modules/local/numbat/main.nf'

workflow NUMBAT {

    take:
        ch_annotated_object
        ch_bam               // tuple(sample_id, possorted_genome_bam.bam)

    main:

        ch_calls = Channel.empty()

        if (params.run_numbat) {

            // Fail early and legibly: the reference panel is a large external
            // download, and discovering it is missing after the pileup has
            // burned hours on a BAM is the worst time to find out.
            [ 'numbat_gmap'    : 'genetic map (Eagle2)',
              'numbat_snpvcf'  : 'SNP VCF',
              'numbat_paneldir': '1000G phasing panel directory' ].each { key, what ->
                if (!params[key]) {
                    error "--run_numbat requires --${key} (${what}). See docs/NUMBAT.md for how to obtain them."
                }
            }

            NUMBAT_PREP(ch_annotated_object.first())

            // Key the per-sample artefacts by sample id so they join the BAMs.
            ch_counts  = NUMBAT_PREP.out.counts.flatten()
                .map { f -> tuple(f.name.replaceFirst(/_counts\.rds$/, ''), f) }
            ch_barcodes = NUMBAT_PREP.out.barcodes.flatten()
                .map { f -> tuple(f.name.replaceFirst(/_barcodes\.tsv$/, ''), f) }
            ch_cellmap  = NUMBAT_PREP.out.cellmap.flatten()
                .map { f -> tuple(f.name.replaceFirst(/_cellmap\.tsv$/, ''), f) }

            // An inner join: only samples with BOTH a BAM and post-QC cells are
            // run. A sample dropped at QC contributes no cells; a sample whose
            // BAM was never produced cannot be phased.
            ch_pileup_in = ch_bam.join(ch_barcodes, failOnMismatch: false, remainder: false)

            NUMBAT_PILEUP(
                ch_pileup_in,
                file(params.numbat_gmap,     checkIfExists: true),
                file(params.numbat_snpvcf,   checkIfExists: true),
                file(params.numbat_paneldir, type: 'dir', checkIfExists: true)
            )

            ch_run_in = NUMBAT_PILEUP.out.alleles
                .join(ch_counts,  remainder: false)
                .join(ch_cellmap, remainder: false)

            NUMBAT_RUN(ch_run_in)
            ch_calls = NUMBAT_RUN.out.calls
        }

    emit:
        calls = ch_calls
}
