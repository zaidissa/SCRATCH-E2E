/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Python NMF engine for the metaprogram stage
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    Drop-in alternative to METAPROG_NMF (R, NMF::nmf with "snmf/r").

    Both engines consume the same `<sample>_preprocessed.rds` produced by
    METAPROG_NMF_PREP, and both hand METAPROG_POST a gene x program basis matrix
    per rank. Only the solver differs, which is what makes the two comparable —
    and what makes `--metaprog_nmf_engine both` a meaningful concordance check
    rather than an apples-to-oranges run.

    Selected with:  --metaprog_nmf_engine python | r | both
----------------------------------------------------------------------------------------
*/

process METAPROG_EXPORT_MTX {

    tag "NMF export: ${sample_id}"
    label 'process_low'

    // Runs in the R image purely to read the .rds; it performs no computation
    // beyond a format change, so it is cheap despite the per-sample fan-out.

    input:
        tuple val(sample_id), path(preprocessed_rds)

    output:
        tuple val(sample_id), path("export/${sample_id}.mtx"),
                              path("export/${sample_id}_genes.tsv"),
                              path("export/${sample_id}_cells.tsv"), emit: mtx

    when:
        task.ext.when == null || task.ext.when

    script:
        """
        set -euo pipefail
        mkdir -p export
        export_preprocessed_mtx.R ${preprocessed_rds} export/${sample_id}
        """

    stub:
        """
        mkdir -p export
        printf '%%%%MatrixMarket matrix coordinate real general\\n1 1 1\\n1 1 1.0\\n' > export/${sample_id}.mtx
        echo GENE1 > export/${sample_id}_genes.tsv
        echo CELL1 > export/${sample_id}_cells.tsv
        """
}


process METAPROG_NMF_PY {

    tag "NMF (python): ${sample_id}"
    label 'process_medium'

    input:
        tuple val(sample_id), path(mtx), path(genes), path(cells)

    output:
        path "data/per_sample_mat/nmf_fit/${sample_id}_basis_matrix.csv", emit: basis
        path "data/per_sample_mat/nmf_fit/${sample_id}_top_genes.csv",    emit: top_genes
        path "data/per_sample_mat/nmf_fit/${sample_id}_nmf_metrics.csv",  emit: metrics
        path "data/per_sample_mat/nmf_fit/${sample_id}_nmf_meta.json",    emit: meta
        path "figures/metaprog/${sample_id}_nmf_rank.png",                emit: figures, optional: true

    when:
        task.ext.when == null || task.ext.when

    script:
        def extras = task.ext.args ?: ''
        """
        set -euo pipefail
        mkdir -p data/per_sample_mat/nmf_fit figures/metaprog

        # Pin BLAS threads to the task allocation. Without this each of the N
        # concurrent per-sample jobs would try to use every core on the node.
        export NMF_NUM_THREADS=${task.cpus}
        export OMP_NUM_THREADS=${task.cpus}
        export OPENBLAS_NUM_THREADS=${task.cpus}
        export MKL_NUM_THREADS=${task.cpus}

        scratch_nmf.py \\
            --prefix ${sample_id} \\
            --sample ${sample_id} \\
            --outdir data/per_sample_mat/nmf_fit \\
            ${extras}

        # Figures live under figures/ by the repo-wide convention.
        if [ -f "data/per_sample_mat/nmf_fit/${sample_id}_nmf_rank.png" ]; then
            mv "data/per_sample_mat/nmf_fit/${sample_id}_nmf_rank.png" figures/metaprog/
        fi
        """

    stub:
        """
        mkdir -p data/per_sample_mat/nmf_fit figures/metaprog
        echo ",3_0,3_1,3_2" > data/per_sample_mat/nmf_fit/${sample_id}_basis_matrix.csv
        echo "GENE1,0.1,0.2,0.3" >> data/per_sample_mat/nmf_fit/${sample_id}_basis_matrix.csv
        echo "sample,program,rank,position,gene,loading" > data/per_sample_mat/nmf_fit/${sample_id}_top_genes.csv
        echo "rank,run,seed,init,reconstruction_err,n_iter,converged,seconds" > data/per_sample_mat/nmf_fit/${sample_id}_nmf_metrics.csv
        echo '{"sample":"${sample_id}","engine":"python-sklearn"}' > data/per_sample_mat/nmf_fit/${sample_id}_nmf_meta.json
        touch figures/metaprog/${sample_id}_nmf_rank.png
        """
}
