/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    NUMBAT — allele-aware CNV calling, opt-in
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Joins per-sample BAMs (from ALIGN) to per-sample expression (from the
    annotated object) and emits per-cell calls for STRATIFY's evidence pool.
----------------------------------------------------------------------------------------
*/

include { asBool } from '../../lib/booleans.nf'
include { NUMBAT_FETCH_REFS } from '../../modules/local/numbat/main.nf'
include { NUMBAT_PREP   } from '../../modules/local/numbat/main.nf'
include { NUMBAT_PILEUP } from '../../modules/local/numbat/main.nf'
include { NUMBAT_RUN    } from '../../modules/local/numbat/main.nf'
include { NUMBAT_REPORT } from '../../modules/local/numbat/main.nf'

workflow NUMBAT {

    take:
        ch_annotated_object
        ch_bam               // tuple(sample_id, possorted_genome_bam.bam, .bam.bai)

    main:

        ch_calls = Channel.empty()
        ch_segments = Channel.empty()
        ch_figures = Channel.empty()

        if (asBool(params.run_numbat)) {

            // The pileup image must be PULLABLE BY THE EXECUTOR, not merely
            // present on the machine that built it. `scratch-numbat:procps` is a
            // bare local tag, and Docker reads a reference with no '/' as a Docker
            // Hub official image, so AWS Batch asked the Hub for
            // library/scratch-numbat and got
            //     CannotPullImageManifestError ... denied: requested access to the
            //     resource is denied / unauthorized: authentication required
            // That is a task-START failure, which `errorStrategy retry` cannot
            // help: all 15 pileups were submitted, each retried three times, and
            // the run died after 1h50m and 33 CPU-hours having produced no calls
            // (Cirro 2026-10-20, dataset 9c9d99d4). The error names neither the
            // image nor the parameter that chose it.
            //
            // conf/containers.config exists to prevent exactly this -- its header
            // lists "local-only tags never pushed to a registry" as one of the
            // three portability bugs it was written to remove -- and Numbat
            // reintroduced one. Check it here, before NUMBAT_PREP spends an hour
            // building per-sample inputs for a stage that cannot start.
            //
            // The test is just "does the reference contain a '/'", which is the
            // same rule Docker applies when deciding whether a name is
            // official-library shorthand. `syedsazaidi/scratch-cnv:latest` and
            // `us.gcr.io/broad-dsde-methods/cellbender:0.3.0` pass;
            // `nf-quarto:latest` and `scratch-numbat:procps` do not.
            //
            // Exempt only the local executor, where a locally built tag is the
            // whole point. `workflow.profile` is NOT the discriminator: the
            // documented laptop invocation is `-profile docker`, which leaves the
            // executor at its default, and Cirro reports profile 'standard' while
            // setting the executor from its own compute config. Read the executor
            // itself -- null means Nextflow's default, which is local. Measured on
            // both 25.04.8 (Cirro's version) and 26.04.6: `navigate` returns null
            // when unset and 'awsbatch' when .cirro/process-compute.config is in
            // play. A bare `session` is undefined under the strict parser, hence
            // the fully qualified field.
            def numbat_image = params.image_numbat?.toString()
            def nxf_executor = nextflow.Global.session.config
                                   .navigate('process.executor')?.toString() ?: 'local'
            if (nxf_executor != 'local' && numbat_image && !numbat_image.contains('/')) {
                error "Numbat's image '${numbat_image}' is a local-only Docker tag, so the " +
                      "'${nxf_executor}' executor cannot pull it.\n" +
                      "  A reference with no '/' is read as a Docker Hub official image " +
                      "(library/${numbat_image.tokenize(':').first()}), which does not " +
                      "exist, so every pileup task fails at image pull with \"denied: " +
                      "requested access to the resource is denied\" -- after the whole " +
                      "fan-out has already been submitted.\n" +
                      "  Push it somewhere the executor can reach, then point this there:\n" +
                      "      docker tag  ${numbat_image} <registry>/${numbat_image}\n" +
                      "      docker push <registry>/${numbat_image}\n" +
                      "      nextflow run ... --container_numbat '<registry>/${numbat_image}'\n" +
                      "  The image is containers/numbat/Dockerfile = " +
                      "pkharchenkolab/numbat-rbase:latest + procps. Almost all of its 43.8 GB " +
                      "is the bundled 1000G panel, which the public base image already " +
                      "carries, so a Docker Hub push cross-mounts those blobs and uploads " +
                      "only the 4.8 MB procps layer.\n" +
                      "  To reach STRATIFY without allele-aware CNV instead, set " +
                      "--run_numbat false; the malignant split then rests on inferCNV and " +
                      "the annotation label alone."
            }

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
            // A PARTIAL mismatch is still silent. `ifEmpty` below only fires when
            // NOTHING matches; with 17 libraries, three BAM directories named
            // differently from their sample column means three samples get no
            // Numbat call and the run says nothing at all. That is the caller
            // STRATIFY is told to trust over the others
            // (cnv_caller_priority = numbat), so those samples quietly fall back
            // to inferCNV or to the annotation label.
            //
            // Report the join itself, both sides, by name. `println` and not
            // `log.warn`: inside a `subscribe` closure the logger produces no
            // output at all (same reason `workflow.onComplete` has to use
            // `println`), so a log call here would look correct and print
            // nothing — verified against 26.04.6 on a three-key harness.
            //
            // The keys are collected through a TAGGED mix rather than paired with
            // `combine` or `merge`: both of those FLATTEN two list-emitting
            // channels into one list, which loses the boundary between the sides.
            // `combine` failed with "Unknown method invocation `intersect` on
            // String type" and `merge` handed the closure a single
            // [A, B, C, A, B] — measured, not assumed.
            ch_bam.map { tuple('bam', it[0]) }
                .mix(ch_barcodes.map { tuple('obj', it[0]) })
                .collect(flat: false)
                .subscribe { rows ->
                    def bam_keys = rows.findAll { it[0] == 'bam' }.collect { it[1] }.unique().sort()
                    def obj_keys = rows.findAll { it[0] == 'obj' }.collect { it[1] }.unique().sort()
                    def matched  = bam_keys.intersect(obj_keys)
                    def bam_only = bam_keys - obj_keys
                    def obj_only = obj_keys - bam_keys
                    if (bam_only || obj_only) {
                        println "WARNING: Numbat will phase ${matched.size()} of " +
                            "${bam_keys.size()} BAM sample(s). " +
                            (bam_only ? "No post-QC cells for: ${bam_only.join(', ')}. " : '') +
                            (obj_only ? "No BAM for: ${obj_only.join(', ')}. " : '') +
                            "Those samples get NO Numbat call, and STRATIFY falls back to " +
                            "the other callers or the annotation label for their cells. " +
                            "The BAM key is the directory two levels above the .bam; the " +
                            "sample key is --numbat_sample_col ('${params.numbat_sample_col}')."
                    } else {
                        println "Numbat: all ${matched.size()} sample(s) matched a BAM."
                    }
                }

            ch_pileup_in = ch_bam.join(ch_barcodes, failOnMismatch: false, remainder: false)
                .ifEmpty {
                    // This fires whenever the join is empty, which has TWO causes, and
                    // the common one is not a key mismatch. `ifEmpty` is a runtime
                    // check: if any earlier stage failed, NUMBAT_PREP never ran, the
                    // right-hand side of the join is empty, and this message is printed
                    // as the run's `Error :` even though nothing is wrong with the keys
                    // (2026-10-16 on Cirro: SEURAT_CLUSTER died and this is what the
                    // summary blamed). Say so first, so the message cannot send the
                    // next reader looking for a naming problem that does not exist.
                    error "Numbat: nothing to pile up -- the BAM channel and the " +
                          "annotated object had no sample in common.\n" +
                          "  If an earlier error appears above, THAT is the cause: " +
                          "NUMBAT_PREP never produced per-sample barcodes, so this join " +
                          "had nothing to match against. Fix that failure first; this " +
                          "line is a consequence of it, not an independent problem.\n" +
                          "  If the rest of the run succeeded, the two keys genuinely " +
                          "disagree:\n" +
                          "    BAM key    = the directory two levels above the .bam " +
                          "(--input_bam_path '<...>/<SAMPLE>/outs/*.bam')\n" +
                          "    Sample key = --numbat_sample_col " +
                          "('${params.numbat_sample_col}'), falling back to " +
                          "orig.ident/sample/patient_id\n" +
                          "  These must be the same strings. Rename the BAM directories " +
                          "or set --numbat_sample_col to the column that matches them."
                }

            NUMBAT_PILEUP(ch_pileup_in, ch_gmap, ch_snpvcf, ch_paneldir)

            ch_run_in = NUMBAT_PILEUP.out.alleles
                .join(ch_counts,  remainder: false)
                .join(ch_cellmap, remainder: false)

            NUMBAT_RUN(ch_run_in)
            ch_calls = NUMBAT_RUN.out.calls
            ch_segments = NUMBAT_RUN.out.segments

            // Numbat's own figures. Drawn in a SEPARATE task from the finished
            // output directory rather than inside NUMBAT_RUN, because the task
            // hash covers the process script: one plotting line added to
            // NUMBAT_RUN re-runs the model fit for all 15 libraries (957 CPU
            // hours, measured on GBM-DFCI3) to redraw pictures from numbers
            // that are already on disk.
            //
            // `detail` is declared as numbat/<sample>/**, so it emits the
            // FILES; the directory is recovered from them and the sample id
            // from the directory's name -- the same shape cnv.nf already uses
            // to stage the segment consensus. The grandparent test is what
            // keeps a nested file from contributing its own subdirectory as a
            // second, wrong "sample": only a file sitting directly in
            // numbat/<sample>/ has 'numbat' two levels up.
            if (!asBool(params.skip_numbat_figures)) {
                ch_numbat_dirs = NUMBAT_RUN.out.detail
                    .flatten()
                    .map { f -> f.parent }
                    .filter { d -> d.parent.name == 'numbat' }
                    .unique()
                    .map { d -> tuple(d.name, d) }

                NUMBAT_REPORT(ch_numbat_dirs)
                ch_figures = NUMBAT_REPORT.out.figures
            }
        }

    emit:
        calls = ch_calls
        // Per-SEGMENT consensus calls, keyed by sample. Distinct from `calls`,
        // which is per cell.
        segments = ch_segments
        // Numbat's own diagnostic figures, one set per sample.
        figures = ch_figures
}
