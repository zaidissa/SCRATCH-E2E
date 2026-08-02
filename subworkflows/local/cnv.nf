/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    CNV — inferCNV, SCEVAN, CopyKAT
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    The three callers are independent and run in parallel. Each is individually
    skippable; per-tool `ext.when` in conf/modules.config does the gating so a
    disabled tool leaves an empty channel rather than an error.
----------------------------------------------------------------------------------------
*/

include { INFERCNV        } from '../../modules/local/infercnv/main.nf'
include { SCEVAN          } from '../../modules/local/scevan/main.nf'
include { COPYKAT_PROCESS } from '../../modules/local/copykat/main.nf'

workflow CNV {

    take:
        ch_seurat_object
        ch_page_config

    main:
        ch_infercnv = Channel.empty()
        ch_scevan   = Channel.empty()
        ch_copykat  = Channel.empty()

        // Per-cell calls, kept separate from the raw caller output directories.
        // These are what the stratification step consumes.
        ch_infercnv_meta = Channel.empty()
        ch_scevan_meta   = Channel.empty()

        if (!params.skip_infercnv) {
            INFERCNV(ch_seurat_object,
                     Channel.fromPath(params.notebook_infercnv, checkIfExists: true),
                     ch_page_config)
            ch_infercnv      = INFERCNV.out.results
            ch_infercnv_meta = INFERCNV.out.meta_object
        }

        if (!params.skip_scevan) {
            SCEVAN(ch_seurat_object,
                   Channel.fromPath(params.notebook_scevan, checkIfExists: true),
                   ch_page_config)
            ch_scevan      = SCEVAN.out.results
            ch_scevan_meta = SCEVAN.out.meta_object
        }

        if (!params.skip_copykat) {
            COPYKAT_PROCESS(ch_seurat_object,
                            Channel.fromPath(params.notebook_copykat, checkIfExists: true),
                            ch_page_config)
            ch_copykat = COPYKAT_PROCESS.out.results
        }

    emit:
        infercnv      = ch_infercnv
        scevan        = ch_scevan
        copykat       = ch_copykat
        infercnv_meta = ch_infercnv_meta
        scevan_meta   = ch_scevan_meta
}
