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
include { NUMBAT_GERMLINE } from '../../modules/local/numbat/main.nf'
include { NUMBAT_CELLSNP } from '../../modules/local/numbat/main.nf'
include { NUMBAT_PHASE   } from '../../modules/local/numbat/main.nf'
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
                // Placeholders: NUMBAT_CELLSNP / NUMBAT_PHASE resolve these to the in-image paths.
                ch_gmap     = Channel.value(nf_gmap)
                ch_snpvcf   = Channel.value(nf_vcf)
                ch_paneldir = Channel.value(nf_panel)
            }

            // The label matcher is shared with STRATIFY and CNV_CONCORDANCE:
            // the reference lineages that define the CNV null there have to be
            // the same cells that define Numbat's expression reference here.
            ch_helpers = Channel.value(
                file("${projectDir}/assets/viz/scratch_report.R", checkIfExists: true))

            NUMBAT_PREP(ch_annotated_object.first(), ch_helpers)

            // ---------------------------------------------------------------
            // GROUPING. Numbat runs per INDIVIDUAL, over all of that
            // individual's libraries -- see the header of bin/numbat_prep.R.
            // numbat_prep.R decides the grouping (it is the thing that can read
            // the object's metadata) and writes it out as groups.tsv; this
            // reads that back to drive the channel topology, so there is one
            // decision about what a group is and it lives in one place.
            // ---------------------------------------------------------------
            ch_lib2group = NUMBAT_PREP.out.groups
                .splitCsv(header: true, sep: '\t')
                .map { r -> tuple(r.library.toString(), r.group.toString()) }

            // Per-library artefacts, keyed by LIBRARY so they still join the BAMs.
            ch_barcodes = NUMBAT_PREP.out.barcodes.flatten()
                .map { f -> tuple(f.name.replaceFirst(/_barcodes\.tsv$/, ''), f) }

            // Per-group artefacts, keyed by GROUP.
            ch_counts  = NUMBAT_PREP.out.counts.flatten()
                .map { f -> tuple(f.name.replaceFirst(/_counts\.rds$/, ''), f) }
            ch_cellmap = NUMBAT_PREP.out.cellmap.flatten()
                .map { f -> tuple(f.name.replaceFirst(/_cellmap\.tsv$/, ''), f) }
            // The cells for the individual's own expression reference. PREP
            // selects them (it has the annotation and the shared matcher), RUN
            // aggregates them -- aggregate_counts() lives in the numbat package
            // and PREP runs in image_cnv, which has Seurat and no numbat.
            // Always written, possibly header-only: an optional per-group file
            // would have to be joined with `remainder`, and a group silently
            // missing its reference would fall back to ref_hca with nothing
            // saying so. numbat_run.R reports which basis it used either way.
            ch_refcells = NUMBAT_PREP.out.refcells.flatten()
                .map { f -> tuple(f.name.replaceFirst(/_refcells\.tsv$/, ''), f) }

            // An inner join: only libraries with BOTH a BAM and post-QC cells are
            // run. A library dropped at QC contributes no cells; a library whose
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
            // BAM key    = the grandparent directory name (.../<library>/outs/x.bam)
            // Library key = --numbat_sample_col, falling back to orig.ident
            // A PARTIAL mismatch is still silent. `ifEmpty` below only fires when
            // NOTHING matches; with 17 libraries, three BAM directories named
            // differently from their sample column means three libraries get no
            // Numbat call and the run says nothing at all.
            //
            // Report the join itself, both sides, by name. `println` and not
            // `log.warn`: inside a `subscribe` closure the logger produces no
            // output at all (same reason `workflow.onComplete` has to use
            // `println`), so a log call here would look correct and print
            // nothing -- verified against 26.04.6 on a three-key harness.
            //
            // The keys are collected through a TAGGED mix rather than paired with
            // `combine` or `merge`: both of those FLATTEN two list-emitting
            // channels into one list, which loses the boundary between the sides.
            // `combine` failed with "Unknown method invocation `intersect` on
            // String type" and `merge` handed the closure a single
            // [A, B, C, A, B] -- measured, not assumed.
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
                            "${bam_keys.size()} BAM librar(ies). " +
                            (bam_only ? "No post-QC cells for: ${bam_only.join(', ')}. " : '') +
                            (obj_only ? "No BAM for: ${obj_only.join(', ')}. " : '') +
                            "Those libraries get NO Numbat call, and STRATIFY falls back to " +
                            "the other callers or the annotation label for their cells. " +
                            "The BAM key is the directory two levels above the .bam; the " +
                            "library key is --numbat_sample_col ('${params.numbat_sample_col}')."
                    } else {
                        println "Numbat: all ${matched.size()} librar(ies) matched a BAM."
                    }
                }

            // Say what the grouping actually produced. A cohort where every
            // group holds one library is STILL RUNNING PER LIBRARY -- the
            // configuration that over-calls tumour in normal libraries -- and
            // without this line the only symptom is a plausible-looking number
            // in a report. numbat_prep.R warns in its own log; this surfaces it
            // in the run's console output where it will actually be read.
            NUMBAT_PREP.out.grouping_report
                .splitCsv(header: true, sep: '\t')
                .collect(flat: false)
                .subscribe { rows ->
                    def n      = rows.size()
                    def multi  = rows.count { (it.n_libraries as Integer) > 1 }
                    def ext    = rows.count { it.reference_basis != 'internal' }
                    println "Numbat: ${n} individual(s), ${multi} pooling more than one library."
                    if (multi == 0) {
                        println "WARNING: Numbat is running PER LIBRARY. numbat:::mark_tumor_lineage() " +
                            "has no null -- it always labels its heaviest clade \"tumor\" -- so a " +
                            "library with no malignant cells will still return neoplastic calls. " +
                            "Group libraries by individual: add patient_id to the metadata, or set " +
                            "--numbat_group_col to the column that identifies the individual."
                    }
                    if (ext > 0) {
                        println "Numbat: ${ext} individual(s) had too few reference-lineage cells " +
                            "for an internal expression reference and fall back to numbat::ref_hca " +
                            "(--numbat_min_ref_cells = ${params.numbat_min_ref_cells})."
                    }
                    // Whether the compartment split will exist at all. Saying
                    // this BEFORE the pileup matters: the run is hours long, and
                    // discovering afterwards that the one table a reviewer asked
                    // for was never written is a second run, not a second look.
                    def st = rows.collect { it.sample_types }
                                 .findAll { it && it != 'NA' && it != '' }
                    if (st) {
                        println "Numbat: sample types present: " +
                            "${st.join(' / ')}. The malignant fraction will also be " +
                            "published per compartment (*_numbat_by_sample_type.tsv)."
                    } else {
                        println "WARNING: no sample-type column, so Numbat cannot report " +
                            "the malignant fraction per compartment. \"There should be almost " +
                            "zero neoplastic cells in PBMC\" is a claim about a compartment, " +
                            "and without this column it can only be checked by joining the " +
                            "per-library table to the metadata by hand. Add sample_type to " +
                            "--input_exp_table, or set --numbat_sample_type_col."
                    }

                    // How many COLUMNS each internal reference has. Numbat picks
                    // the best-fitting column per pseudobulk, so a single-column
                    // baseline gives it nothing to choose from -- and that is
                    // invisible in the calls themselves.
                    def cols = rows.findAll { it.reference_basis == 'internal' && it.n_reference_columns }
                    if (cols) {
                        def ncol = cols.collect { it.n_reference_columns as Integer }
                        println "Numbat: internal expression reference has ${ncol.min()}-${ncol.max()} " +
                            "column(s) per individual (grouped by " +
                            "${params.numbat_cluster_col}, floor " +
                            "${params.numbat_ref_min_group_cells} cells/column)."
                        def thin = ncol.count { it < 2 }
                        if (thin > 0) {
                            println "WARNING: ${thin} individual(s) have a single-column expression " +
                                "reference. Numbat matches each pseudobulk against the columns of " +
                                "lambdas_ref, so one column means no choice at all. Lower " +
                                "--numbat_ref_min_group_cells, or check that " +
                                "--numbat_cluster_col names a real cluster column."
                        }
                    }
                }

            // library -> (group, library, bam, bai, barcodes): one pileup task
            // per LIBRARY. The regrouping onto the individual happens after the
            // pileup, in front of NUMBAT_PHASE.
            ch_cellsnp_in = ch_bam
                .join(ch_barcodes,  failOnMismatch: false, remainder: false)
                .join(ch_lib2group, failOnMismatch: false, remainder: false)
                .map { lib, bam, bai, bc, grp -> tuple(grp, lib, bam, bai, bc) }
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
                          "annotated object had no library in common.\n" +
                          "  If an earlier error appears above, THAT is the cause: " +
                          "NUMBAT_PREP never produced per-library barcodes, so this join " +
                          "had nothing to match against. Fix that failure first; this " +
                          "line is a consequence of it, not an independent problem.\n" +
                          "  If the rest of the run succeeded, the two keys genuinely " +
                          "disagree:\n" +
                          "    BAM key     = the directory two levels above the .bam " +
                          "(--input_bam_path '<...>/<LIBRARY>/outs/*.bam')\n" +
                          "    Library key = --numbat_sample_col " +
                          "('${params.numbat_sample_col}'), falling back to " +
                          "orig.ident/sample\n" +
                          "  These must be the same strings. Rename the BAM directories " +
                          "or set --numbat_sample_col to the column that matches them."
                }

            // ---- Germline genotypes (optional, per individual) --------------
            // --numbat_germline_vcfs: a CSV with columns patient_id, vcf and an
            // optional sample (the VCF's sample column; required when the VCF
            // holds more than one). patient_id values are --numbat_group_col
            // values. An individual with no row keeps the scRNA-inferred
            // genotype, so a cohort can be covered partially.
            def nf_germ = file("${projectDir}/assets/NO_FILE_germline")
            ch_germline = Channel.empty()
            if (params.numbat_germline_vcfs) {
                ch_germ_rows = Channel
                    .fromPath(params.numbat_germline_vcfs, checkIfExists: true)
                    .splitCsv(header: true, strip: true)
                    .map { r ->
                        if (!r.patient_id || !r.vcf)
                            error "--numbat_germline_vcfs needs columns 'patient_id' and 'vcf' " +
                                  "(and optionally 'sample'); got: ${r.keySet().join(', ')}"
                        tuple(r.patient_id.toString(), file(r.vcf, checkIfExists: true),
                              (r.sample ?: '').toString())
                    }
                NUMBAT_GERMLINE(ch_germ_rows, ch_snpvcf)
                ch_germline = NUMBAT_GERMLINE.out.targets.join(NUMBAT_GERMLINE.out.germline)

                // Say which individuals got a germline genotype and which did
                // not. A patient_id typo would otherwise mean "silently the old
                // behaviour" for that patient.
                ch_cellsnp_in.map { it[0].toString() }.unique().map { ['data', it] }
                    .mix(ch_germ_rows.map { ['germ', it[0]] })
                    .collect(flat: false)
                    .subscribe { rows ->
                        def data = rows.findAll { it[0] == 'data' }.collect { it[1] } as Set
                        def germ = rows.findAll { it[0] == 'germ' }.collect { it[1] } as Set
                        println "Numbat: germline WGS genotype for ${(data.intersect(germ)).size()} " +
                            "of ${data.size()} individual(s); the rest infer it from the scRNA."
                        def unknown = germ - data
                        if (unknown)
                            println "WARNING: --numbat_germline_vcfs names individual(s) not in " +
                                "this run's --numbat_group_col (${params.numbat_group_col}): " +
                                "${unknown.sort().join(', ')}. Their VCFs are prepared but used " +
                                "by nothing -- check the spelling against the metadata."
                        def without = data - germ
                        if (without && germ)
                            println "Numbat: no germline VCF for: ${without.sort().join(', ')}."
                    }
            }
            // One map, group -> [targets, germline dir], as a value: every
            // pileup row looks its individual up in it. toList() of an empty
            // channel is [], so without germline VCFs this is an empty map and
            // every row falls back to the 1000G SNP VCF.
            ch_germ_map = ch_germline.toList().map { rows ->
                rows.collectEntries { r -> [ (r[0].toString()): [ r[1], r[2] ] ] }
            }

            ch_cellsnp_rows = ch_cellsnp_in
                .combine(ch_germ_map)
                .combine(ch_snpvcf)
                .map { grp, lib, bam, bai, bc, gmapping, snp ->
                    def g = gmapping[grp.toString()]
                    tuple(grp, lib, bam, bai, bc, g ? g[0] : snp)
                }

            NUMBAT_CELLSNP(ch_cellsnp_rows)

            // How many libraries each individual sent to the pileup. Known as
            // soon as the inputs are, so groupKey can release an individual to
            // NUMBAT_PHASE the moment ITS last pileup lands, instead of waiting
            // for every pileup in the run.
            ch_n_libs = ch_cellsnp_in
                .map { grp, lib, bam, bai, bc -> tuple(grp, lib) }
                .groupTuple()
                .map { grp, libs -> tuple(grp, libs.size()) }

            // Rows are SORTED BY LIBRARY before the lists are built: groupTuple
            // emits in arrival order, which is not deterministic, and the lists
            // correspond to --samples positionally. Sorting makes the
            // correspondence stable and the task hash reproducible.
            //
            // remainder: true, because NUMBAT_CELLSNP gives up on a library
            // after its retries (conf/base.config). Without it, an individual
            // missing one library would be dropped whole -- the same thing as
            // its every library failing. With it, the individual is phased on
            // the libraries that did pile up, and the warning says how many.
            ch_phase_in = NUMBAT_CELLSNP.out.pileup
                .combine(ch_n_libs, by: 0)
                .map { grp, lib, dir, n -> tuple(groupKey(grp, n), lib, dir) }
                .groupTuple(remainder: true)
                .map { key, libs, dirs ->
                    def grp = key.getGroupTarget()
                    if (libs.size() < key.getGroupSize())
                        println "WARNING: Numbat: individual ${grp} is phased on " +
                            "${libs.size()} of ${key.getGroupSize()} libraries -- the " +
                            "pileup failed for the rest. Their cells stay in the count " +
                            "matrix with no allele data, so Numbat calls them from " +
                            "expression alone, as inferCNV would."
                    // .toList() first: an IntRange is immutable and Groovy's
                    // sort() is in-place, so sorting the range directly throws
                    // UnsupportedOperationException at runtime -- which surfaces
                    // only as a bare java.lang.UnsupportedOperationException with
                    // no line number.
                    def order = (0..<libs.size()).toList().sort { libs[it].toString() }
                    tuple(grp, order.collect { libs[it] }, order.collect { dirs[it] })
                }
                .combine(ch_germ_map)
                .map { grp, libs, dirs, gmapping ->
                    def g = gmapping[grp.toString()]
                    tuple(grp, libs, dirs, g ? g[1] : nf_germ)
                }

            NUMBAT_PHASE(ch_phase_in, ch_gmap, ch_snpvcf, ch_paneldir)

            ch_run_in = NUMBAT_PHASE.out.alleles
                .join(ch_counts,  remainder: false)
                .join(ch_cellmap, remainder: false)
                .join(ch_refcells, remainder: false)

            NUMBAT_RUN(ch_run_in)
            ch_calls = NUMBAT_RUN.out.calls
            ch_segments = NUMBAT_RUN.out.segments

            // The reviewer's number, in the run's own console output. It is
            // already in .command.log and in a published table; repeating it
            // here is deliberate, because the one place everybody looks at is
            // the terminal, and "the PBMCs look a bit strange" was a reading of
            // a table nobody had reason to doubt.
            //
            // Both percentages, always. The raw one cannot be near zero for any
            // compartment -- mark_tumor_lineage() always labels its heaviest
            // clade tumour -- so printing it alone would reproduce exactly the
            // misreading this is meant to prevent.
            NUMBAT_RUN.out.by_sample_type
                .splitCsv(header: true, sep: '\t')
                .collect(flat: false)
                .subscribe { rows ->
                    rows.each { r ->
                        def gated = r.containsKey('pct_malignant_gated')
                        println "Numbat [${r.sample_type}]: ${r.n_cells} cells over " +
                            "${r.n_libraries} librar(ies), ${r.pct_malignant_raw}% malignant raw" +
                            (gated ? ", ${r.pct_malignant_gated}% after the lineage gate" : '') + '.'
                    }
                    if (rows.size() > 1 && rows[0].containsKey('pct_malignant_gated')) {
                        def sorted = rows.sort { (it.pct_malignant_gated as Double) }
                        println "Numbat: gated malignant fraction runs from " +
                            "${sorted.first().pct_malignant_gated}% (${sorted.first().sample_type}) to " +
                            "${sorted.last().pct_malignant_gated}% (${sorted.last().sample_type}). " +
                            "A compartment expected to be tumour-free belongs at the low end."
                    }
                }

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
