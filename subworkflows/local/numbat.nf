/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    NUMBAT — allele-aware CNV calling, opt-in
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Joins per-sample BAMs (from ALIGN) to per-sample expression (from the
    annotated object) and emits per-cell calls for STRATIFY's evidence pool.
----------------------------------------------------------------------------------------
*/

include { NUMBAT_FETCH_REFS } from '../../modules/local/numbat/main.nf'
include { NUMBAT_PREP   } from '../../modules/local/numbat/main.nf'
include { NUMBAT_PILEUP } from '../../modules/local/numbat/main.nf'
include { NUMBAT_RUN    } from '../../modules/local/numbat/main.nf'

workflow NUMBAT {

    take:
        ch_annotated_object
        ch_bam               // tuple(sample_id, possorted_genome_bam.bam, .bam.bai)

    main:

        ch_calls = Channel.empty()

        if (params.run_numbat) {

            // Fail early and legibly: the reference panel is a large external
            // download, and discovering it is missing after the pileup has
            // burned hours on a BAM is the worst time to find out.
            // The numbat image ships all three references, so "not supplied" means
            // "use the image copy" -- NOT "go and download 10 GB". The download
            // path stays available via --numbat_fetch_refs for anyone running a
            // different image, but it is no longer the default, because the
            // upstream panel URL is behind a WAF that serves a 212-byte HTML
            // challenge with HTTP 200. curl saves that AS the zip and unzip fails.
            // THREE distinct placeholders, not one file passed three times:
            // Nextflow stages each `path` input by its own basename, so reusing a
            // single NO_FILE makes all three land as "NO_FILE" in the task dir and
            // the run dies with "input file name collision". They still all match
            // the *NO_FILE* test the pileup script uses to detect "not supplied".
            def nf_gmap  = file("${projectDir}/assets/NO_FILE_gmap")
            def nf_vcf   = file("${projectDir}/assets/NO_FILE_snpvcf")
            def nf_panel = file("${projectDir}/assets/NO_FILE_paneldir")
            def ref_keys = [ 'numbat_gmap', 'numbat_snpvcf', 'numbat_paneldir' ]
            def supplied = ref_keys.findAll { params[it] }

            if (supplied && supplied.size() < ref_keys.size()) {
                error "Numbat references are partially specified (given: " +
                      "${supplied.join(', ')}). Supply all three for a matching " +
                      "${params.numbat_genome} set, or none and use the copies " +
                      "shipped in the Numbat image."
            }

            if (supplied) {
                ch_gmap     = Channel.value(file(params.numbat_gmap,     checkIfExists: true))
                ch_snpvcf   = Channel.value(file(params.numbat_snpvcf,   checkIfExists: true))
                ch_paneldir = Channel.value(file(params.numbat_paneldir, type: 'dir', checkIfExists: true))
            } else if (params.numbat_fetch_refs) {
                NUMBAT_FETCH_REFS()
                ch_gmap     = NUMBAT_FETCH_REFS.out.gmap
                ch_snpvcf   = NUMBAT_FETCH_REFS.out.snpvcf
                ch_paneldir = NUMBAT_FETCH_REFS.out.paneldir
            } else {
                // Placeholders: NUMBAT_PILEUP resolves these to the in-image paths.
                ch_gmap     = Channel.value(nf_gmap)
                ch_snpvcf   = Channel.value(nf_vcf)
                ch_paneldir = Channel.value(nf_panel)
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
            // `failOnMismatch: false` + `remainder: false` means a key that does not
            // match is DROPPED IN SILENCE. If no key matches at all the channel is
            // simply empty, PILEUP and RUN never run, and the pipeline reports
            // success having produced zero Numbat calls -- verified in a stub run,
            // which reported "Succeeded: 58, Failed: 0" with no pileup at all.
            //
            // That matters more for Numbat than for any other caller, because
            // cnv_caller_priority = numbat means STRATIFY is meant to trust it
            // over the others. Silently losing it changes the malignant split with
            // nothing in the report to say so. This is the same failure that let
            // the SCEVAN doubled-prefix bug run a "two-caller consensus" on one
            // caller, so it fails loudly here instead.
            //
            // BAM key   = the grandparent directory name (…/<sample>/outs/x.bam)
            // Sample key = --numbat_sample_col, falling back to orig.ident
            ch_pileup_in = ch_bam.join(ch_barcodes, failOnMismatch: false, remainder: false)
                .ifEmpty {
                    error "Numbat: no sample matched between the supplied BAMs and the " +
                          "annotated object, so there is nothing to pile up.\n" +
                          "  BAM key    = the directory two levels above the .bam " +
                          "(--input_bam_path '<...>/<SAMPLE>/outs/*.bam')\n" +
                          "  Sample key = --numbat_sample_col ('${params.numbat_sample_col}'), " +
                          "falling back to orig.ident/sample/patient_id\n" +
                          "These must be the same strings. Rename the BAM directories or " +
                          "set --numbat_sample_col to the column that matches them."
                }

            NUMBAT_PILEUP(ch_pileup_in, ch_gmap, ch_snpvcf, ch_paneldir)

            ch_run_in = NUMBAT_PILEUP.out.alleles
                .join(ch_counts,  remainder: false)
                .join(ch_cellmap, remainder: false)

            NUMBAT_RUN(ch_run_in)
            ch_calls = NUMBAT_RUN.out.calls
        }

    emit:
        calls = ch_calls
}
