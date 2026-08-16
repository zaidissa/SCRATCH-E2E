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

        // The microbiome stage consumes the aligned BAM, which lives inside the
        // same cellranger `outs` directory as the count matrices.
        ch_bam = ch_gex_outs
            .flatMap { it instanceof List && it.size() == 2 ? [it] : [] }
            .map    { sample, files ->
                def bam = (files instanceof List ? files : [files]).find {
                    it.toString().endsWith('possorted_genome_bam.bam')
                }
                bam ? tuple(sample, bam) : null
            }
            .filter { it != null }

    emit:
        gex_outs = ch_gex_outs
        vdj_outs = ch_vdj_outs
        bam      = ch_bam
        versions = ch_versions
}
