/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    CELL COMMUNICATION — LIANA, CellChat, NicheNet
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    LIANA and CellChat are independent and run concurrently. NicheNet consumes
    both of their outputs.

    Changed on merge: the standalone pipeline declared NicheNet's inputs as
    hardcoded absolute paths built inside the `params` block —

        liana_csv    = "${params.outdir}/${params.project_name}/cellcommunication/liana/data/liana_results.csv"
        cellchat_rds = "${params.outdir}/.../cellchat_object.rds"

    which pointed at files that only exist after the upstream tasks publish, and
    additionally self-referenced `params.outdir` during params evaluation. It
    worked only because a separate `_done` tick channel imposed ordering. Those
    two params are gone; NicheNet now takes the real output channels, so the
    dependency is expressed in the DAG and `-resume` behaves correctly.
----------------------------------------------------------------------------------------
*/

include { asBool } from '../../lib/booleans.nf'
include { CELLCOMM_LIANA    } from '../../modules/local/Cell_Communication/main.nf'
include { CELLCOMM_CELLCHAT } from '../../modules/local/Cell_Communication/main.nf'
include { CELLCOMM_NICHENET  } from '../../modules/local/Cell_Communication/main.nf'
include { NICHENET_FETCH_REFS } from '../../modules/local/Cell_Communication/main.nf'

// Validate NicheNet's reference RDS files. A Git LFS pointer is a small text
// file beginning "version https://git-lfs..."; a real RDS is binary.
// True when the directory holds unresolved LFS pointers rather than real data.
def nichenetNeedsFetch(dir) {
    def d = file(dir, type: 'dir')
    if (!d.exists()) return true
    def rds = d.listFiles().findAll { it.name.endsWith('.rds') }
    if (!rds) return true
    return rds.any { f -> f.size() < 4096 && f.text.startsWith('version https://git-lfs') }
}

def nichenetAssets(dir) {
    def d = file(dir, type: 'dir', checkIfExists: true)
    def rds = d.listFiles().findAll { it.name.endsWith('.rds') }
    if (!rds) error "NicheNet assets directory '${dir}' contains no .rds files."

    def pointers = rds.findAll { f -> f.size() < 4096 && f.text.startsWith('version https://git-lfs') }
    // A -stub run never reads these files, so refusing to start would make the
    // stub suite depend on data it does not use. Same for an explicitly skipped
    // NicheNet: `ext.when` is evaluated per task, AFTER the DAG is built, so
    // this guard would otherwise block a run that never intended to use them.
    if (pointers && (workflow.stubRun || asBool(params.skip_nichenet))) {
        log.warn "NicheNet reference files are unresolved Git LFS pointers " +
                 "(${pointers*.name.join(', ')}). Fine for -stub; a real run needs " +
                 "bin/fetch_nichenet_refs.sh."
    } else if (pointers) {
        error "NicheNet reference files are unresolved Git LFS pointers, not real data:\n  " +
              pointers.collect { it.name + ' (' + it.size() + ' bytes)' }.join('\n  ') +
              "\n\nFetch them with `bin/fetch_nichenet_refs.sh <dest>` and pass " +
              "--nichenet_assets_dir <dest>, or resolve LFS in the source repo with " +
              "`git lfs pull`.\nOtherwise this fails later inside the notebook as " +
              "\"readRDS(): unknown input format\", which does not indicate the cause."
    }
    return d
}

workflow CELL_COMMUNICATION {

    take:
        ch_seurat_object   // value channel
        ch_page_config     // shared quarto template + figure/report library

    main:
        ch_liana_csv    = Channel.empty()
        ch_cellchat_rds = Channel.empty()
        ch_nichenet     = Channel.empty()

        if (!asBool(params.skip_cellcomm)) {

            nb_liana    = Channel.fromPath(params.liana_qmd,    checkIfExists: true)
            nb_cellchat = Channel.fromPath(params.cellchat_qmd, checkIfExists: true)
            nb_nichenet = Channel.fromPath(params.nichenet_qmd, checkIfExists: true)

            CELLCOMM_LIANA(ch_seurat_object.combine(nb_liana), ch_page_config)
            CELLCOMM_CELLCHAT(ch_seurat_object.combine(nb_cellchat), ch_page_config)

            ch_liana_csv    = CELLCOMM_LIANA.out.liana_csv
            ch_cellchat_rds = CELLCOMM_CELLCHAT.out.cellchat_rds

            // Resolve the reference networks BEFORE the DAG needs them.
            //
            // The three .rds files ship as Git LFS pointers. On a laptop that is
            // one `git lfs pull`; in a managed environment such as Cirro the
            // checkout is not yours to fix, and erroring here sank the run at the
            // ninth of eleven stages over a 320 MB download. So: fetch them when
            // they are unresolved and the user has not named a directory of their
            // own, and keep the hard error for the case where they HAVE named one
            // -- that is a wrong path, not a missing convenience.
            def nn_default   = "${projectDir}/assets/nichenet_resources".toString()
            def nn_is_default = params.nichenet_assets_dir.toString() == nn_default

            ch_nichenet_assets =
                // Deliberately NOT excluded from -stub. The process has a stub
                // branch that touches the three filenames, so a stub run proves
                // this wiring; skipping it here would leave the fetch path
                // untested until a real cloud run, which is how the SCDBLFINDER
                // and CellBender h5 bugs survived to production.
                (nichenetNeedsFetch(params.nichenet_assets_dir) && nn_is_default &&
                 params.nichenet_fetch_refs)
                    ? NICHENET_FETCH_REFS().assets
                    : Channel.value(nichenetAssets(params.nichenet_assets_dir))

            // Real data dependency: NicheNet waits on the actual artefacts.
            CELLCOMM_NICHENET(
                ch_seurat_object
                    .combine(nb_nichenet)
                    .combine(ch_liana_csv)
                    .combine(ch_cellchat_rds)
                    .map { seurat, nb, liana_csv, cellchat_rds -> tuple(seurat, nb, liana_csv, cellchat_rds) },
                // `checkIfExists` only proves the directory is present, not
                // that its contents are usable. The shipped .rds files are
                // 130-byte Git LFS POINTERS that were never resolved, so the
                // check passed and NicheNet died deep inside the notebook with
                // "readRDS(): unknown input format" -- a message that says
                // nothing about the real cause. Validate content here, before a
                // container starts.
                ch_nichenet_assets,
                ch_page_config
            )

            // NicheNet's `summary_csv` emit is commented out in the module, so
            // the generic `data` channel is the available handle.
            ch_nichenet = CELLCOMM_NICHENET.out.data
        }

    emit:
        liana_csv    = ch_liana_csv
        cellchat_rds = ch_cellchat_rds
        nichenet     = ch_nichenet
}
