/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    TUMOUR HETEROGENEITY — Leiden -> NMF preprocessing -> NMF -> metaprograms
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    The NMF fit fans out one task per sample. Two interchangeable engines:

      --metaprog_nmf_engine r        NMF::nmf, method "snmf/r"  (original)
      --metaprog_nmf_engine python   scanpy/scikit-learn        (default)
      --metaprog_nmf_engine both     run both; PostNMF keeps the two program
                                     sets side by side for concordance checking

    Both read the same `<sample>_preprocessed.rds` from METAPROG_NMF_PREP, so
    the comparison is solver-versus-solver rather than pipeline-versus-pipeline.

    Changed on merge: the standalone subworkflow imposed ordering between Leiden
    and preprocessing with an artificial `LEI_TICK` barrier
    (`le_rds.map{1}.take(1)` combined back onto the input). Leiden's output is
    not actually consumed by preprocessing, so that barrier only serialised two
    independent steps. They now run concurrently.
----------------------------------------------------------------------------------------
*/

include { METAPROG_LEIDEN     } from '../../modules/local/MetaProg-DE/main.nf'
include { METAPROG_NMF_PREP   } from '../../modules/local/MetaProg-DE/main.nf'
include { METAPROG_NMF        } from '../../modules/local/MetaProg-DE/main.nf'
include { METAPROG_POST       } from '../../modules/local/MetaProg-DE/main.nf'
include { METAPROG_EXPORT_MTX } from '../../modules/local/MetaProg-DE/nmf_python.nf'
include { METAPROG_NMF_PY     } from '../../modules/local/MetaProg-DE/nmf_python.nf'

workflow METAPROG {

    take:
        ch_seurat_object   // value channel

    main:

        def engine = params.metaprog_nmf_engine?.toString()?.toLowerCase() ?: 'python'
        if (!(engine in ['r', 'python', 'both'])) {
            error "Unknown --metaprog_nmf_engine '${engine}'. Expected one of: r, python, both."
        }

        nb_leiden = Channel.fromPath(params.leiden_qmd,  checkIfExists: true)
        nb_prep   = Channel.fromPath(params.nmfprep_qmd, checkIfExists: true)
        nb_nmf    = Channel.fromPath(params.nmf_qmd,     checkIfExists: true)

        // Leiden and NMF preprocessing are independent; run them concurrently.
        METAPROG_LEIDEN(ch_seurat_object.combine(nb_leiden))
        METAPROG_NMF_PREP(ch_seurat_object.combine(nb_prep))

        // Fan out: one task per <sample>_preprocessed.rds
        ch_samples = METAPROG_NMF_PREP.out.preprocessed
            .flatten()
            .map { f -> tuple(f.baseName.replaceFirst(/_preprocessed$/, ''), f) }

        ch_fits = Channel.empty()

        if (engine in ['r', 'both']) {
            METAPROG_NMF(ch_samples.combine(nb_nmf))
            ch_fits = ch_fits.mix(METAPROG_NMF.out.nmf_rds)
        }

        if (engine in ['python', 'both']) {
            METAPROG_EXPORT_MTX(ch_samples)
            METAPROG_NMF_PY(METAPROG_EXPORT_MTX.out.mtx)
            ch_fits = ch_fits.mix(METAPROG_NMF_PY.out.basis)
        }

        // PostNMF scans data/per_sample_mat/nmf_fit/, so every fit artefact is
        // staged back into that layout regardless of which engine produced it.
        ch_post_in = ch_fits
            .collect()
            .map { fits -> tuple(file(params.postnmf_qmd), fits) }

        METAPROG_POST(ch_post_in)

    emit:
        leiden_rds = METAPROG_LEIDEN.out.leiden_rds
        fits       = ch_fits
        post_data  = METAPROG_POST.out.data
}
