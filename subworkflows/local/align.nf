/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ALIGN — FASTQ to count matrices
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    Three input modes, unchanged from SCRATCH-QC-multi-mode:

      multi = true    hashtag-multiplexed; cellranger multi demultiplexes natively
      demux = true    pooled FASTQs; barcodes auto-detected from the I1 index
      neither         pre-demultiplexed FASTQs, one row per sample

    Changed on merge: the standalone version reassigned a single
    `ch_cellrange_outs` variable inside `if (modality =~ ...)` blocks, so with
    modality 'GEX+TCR' the TCR branch overwrote the GEX handle and the VDJ
    contigs were never reachable downstream. Emissions are now split into three
    named channels so the TCR and microbiome stages can consume their own
    inputs independently of the GEX branch.
----------------------------------------------------------------------------------------
*/

include { SAMPLESHEET_CHECK } from '../../modules/local/helper/validate/main.nf'
include { AUTO_DEMUX        } from '../../modules/local/auto_demux/main.nf'
include { CELLRANGER_COUNT  } from '../../modules/local/cellranger/count/main.nf'
include { CELLRANGER_VDJ    } from '../../modules/local/cellranger/vdj/main.nf'
include { CELLRANGER_MULTI  } from '../../modules/local/cellranger/multi/main.nf'

// Absolute paths and URIs pass through; a relative path is resolved against the
// SAMPLESHEET's own directory, so a sheet works regardless of where the run was
// launched from.
def resolveSheetPath(p, sheetDir) {
    if (!p) return p
    def t = p.toString()
    return (t.startsWith('/') || t =~ /^[a-zA-Z0-9]+:\/\//) ? t : sheetDir.resolve(t).toString()
}

workflow ALIGN {

    take:
        ch_samplesheet   // channel: path to samplesheet CSV
        modality         // val: GEX | TCR | GEX+TCR
        genome           // val: GRCh38 | GRCm39
        demux            // val: boolean
        multi            // val: boolean

    main:

        ch_versions = Channel.empty()

        // Check the reference BEFORE anything runs. `igenomes_base` defaults to
        // ${projectDir}/references, which is not in the repository and is not
        // something a managed environment such as Cirro will have created — so
        // on a fresh checkout this is missing, and without this check the first
        // sign of it is cellranger exiting inside a container after the queue has
        // already granted a 16-core machine.
        def refs = params.genomes?.get(genome)
        if (!refs)
            error "Unknown --genome '${genome}'. Configured: " +
                  "${params.genomes ? params.genomes.keySet().join(', ') : 'none'}."
        def need = (modality =~ /GEX/) ? ['gex'] : []
        if (modality =~ /TCR/) need << 'vdj'
        def missing = need.findAll { k -> !refs[k] || !file(refs[k]).exists() }
        if (missing)
            error "cellranger reference not found for --genome ${genome}:\n" +
                  missing.collect { k -> "  ${k}: ${refs[k]}" }.join('\n') +
                  "\n\nThese are not shipped with the pipeline (~15 GB for GEX).\n" +
                  "  Download from 10x Genomics, then point --igenomes_base at the\n" +
                  "  directory holding refdata/, or set --genomes.${genome}.gex directly.\n" +
                  "  On S3 an s3:// path works.\n" +
                  "  To skip alignment entirely, start from counts:\n" +
                  "    --from qc --input_gex_matrices_path '<path>/*/outs/*'"

        ch_validated_csv = SAMPLESHEET_CHECK(ch_samplesheet).csv

        // Named, independently-consumable outputs.
        ch_gex_outs = Channel.empty()
        ch_vdj_outs = Channel.empty()
        ch_bam      = Channel.empty()

        if (multi) {

            // -------------------------------------------------------------
            // MULTI: sample_id, multi_config
            // cellranger multi emits per-sample matrices under
            //   <sample_id>/outs/per_sample_outs/<sample>/count/
            // which are expanded into one channel element per demultiplexed
            // sample so downstream QC sees them exactly as in standard mode.
            // -------------------------------------------------------------
            ch_multi_runs = ch_validated_csv
                .splitCsv(header: true, sep: ',')
                .map { row -> tuple(row.sample_id, file(row.multi_config)) }

            ch_multi_result = CELLRANGER_MULTI(ch_multi_runs)

            ch_gex_outs = ch_multi_result.per_sample_outs
                .flatMap { run_id, per_sample_dir ->
                    per_sample_dir.listFiles()
                        .findAll { it.isDirectory() }
                        .collect { sample_dir -> tuple(sample_dir.getName(), sample_dir.listFiles().toList()) }
                }

        } else if (demux) {

            // -------------------------------------------------------------
            // DEMUX: run_id, r1, r2, i1, i2, modality
            // -------------------------------------------------------------
            ch_demux_runs = ch_validated_csv
                .splitCsv(header: true, sep: ',')
                .map { row ->
                    tuple(row.run_id, file(row.r1), file(row.r2), file(row.i1), row.i2 ?: '', row.modality)
                }
                // Group every lane of a run into one AUTO_DEMUX call.
                .groupTuple(by: 0)
                .map { run_id, r1, r2, i1, i2, mod ->
                    def i2_val = i2.find { it?.trim() } ?: ''
                    tuple(run_id, r1, r2, i1, i2_val, mod[0])
                }

            ch_demux_result = AUTO_DEMUX(ch_demux_runs).demux_result

            // detected_samples.csv columns: sample, centroid, modality
            ch_per_sample = ch_demux_result
                .flatMap { demux_dir, manifest ->
                    manifest.readLines()
                        .drop(1)
                        .findAll { it.trim() }
                        .collect { line ->
                            def f  = line.split(',')
                            def s  = f[0]
                            tuple(s,
                                  [ demux_dir.resolve("${s}_S1_L001_R1_001.fastq.gz"),
                                    demux_dir.resolve("${s}_S1_L001_R2_001.fastq.gz") ],
                                  f[2])
                        }
                }

            ch_branches = ch_per_sample.branch {
                gex: it[2] == 'GEX'
                tcr: it[2] == 'TCR'
            }

            if (modality =~ /GEX/) {
                ch_gex_outs = CELLRANGER_COUNT(
                    ch_branches.gex.map { s, reads, mod -> tuple(s, reads) },
                    params.genomes[genome].gex
                ).outs
            }

            if (modality =~ /TCR/) {
                ch_vdj_outs = CELLRANGER_VDJ(
                    ch_branches.tcr.map { s, reads, mod -> tuple(s, reads) },
                    params.genomes[genome].vdj
                ).outs
            }

        } else {

            // -------------------------------------------------------------
            // STANDARD: sample, fastq_1, fastq_2, modality
            // -------------------------------------------------------------
            // Relative FASTQ paths are resolved against the SAMPLESHEET's own
            // directory, not the launch directory. A samplesheet that only works
            // when you happen to launch from the right folder is a portability
            // trap, and absolute paths in a shipped test samplesheet work on
            // exactly one machine. Absolute paths are passed through untouched.
            def sheetDir = file(params.input_samplesheet).parent
            ch_rows = ch_validated_csv
                .splitCsv(header: true, sep: ',')
                // Inlined rather than a `def resolve = { }` closure: Nextflow
                // 26.x's strict parser reports "`resolve` is not defined" for a
                // closure variable referenced later in the same block.
                .map { row ->
                    tuple(row.sample,
                          resolveSheetPath(row.fastq_1, sheetDir),
                          resolveSheetPath(row.fastq_2, sheetDir),
                          row.modality)
                }

            ch_branches = ch_rows.branch {
                gex: it[3] == 'GEX'
                tcr: it[3] == 'TCR'
            }

            if (modality =~ /GEX/) {
                ch_gex_outs = CELLRANGER_COUNT(
                    ch_branches.gex
                        .map { r -> tuple(r[0], r[1], r[2]) }
                        .groupTuple(by: [0])
                        .map { r -> tuple(r[0], r[1..2].flatten()) },
                    params.genomes[genome].gex
                ).outs
            }

            if (modality =~ /TCR/) {
                ch_vdj_outs = CELLRANGER_VDJ(
                    ch_branches.tcr
                        .map { r -> tuple(r[0], r[1], r[2]) }
                        .groupTuple(by: [0])
                        .map { r -> tuple(r[0], r[1..2].flatten()) },
                    params.genomes[genome].vdj
                ).outs
            }
        }

        // The microbiome and Numbat stages consume the aligned BAM, which lives
        // inside the same cellranger `outs` directory as the count matrices.
        //
        // The INDEX travels with it. NUMBAT_PILEUP declares
        // `tuple val(sample_id), path(bam), path(bai), path(barcodes)`, and
        // main.nf's --input_bam_path branch builds that 3-tuple; emitting only
        // (sample, bam) here made the ALIGN-fed path one element short, so
        // `--from align --run_numbat true` — the only route that never needs a
        // pre-existing BAM — could not reach the pileup at all. cellranger writes
        // possorted_genome_bam.bam.bai beside the BAM, so there is nothing to
        // rebuild; a missing index falls back to the same placeholder main.nf uses.
        ch_bam = ch_gex_outs
            .flatMap { it instanceof List && it.size() == 2 ? [it] : [] }
            .map    { sample, files ->
                def fl  = (files instanceof List ? files : [files])
                def bam = fl.find { it.toString().endsWith('possorted_genome_bam.bam') }
                def bai = fl.find { it.toString().endsWith('possorted_genome_bam.bam.bai') }
                bam ? tuple(sample, bam,
                            bai ?: file("${projectDir}/assets/NO_FILE_bai")) : null
            }
            .filter { it != null }

    emit:
        gex_outs = ch_gex_outs
        vdj_outs = ch_vdj_outs
        bam      = ch_bam
        versions = ch_versions
}
