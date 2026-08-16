/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    TME CHARACTERISATION — the non-malignant compartment
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    Runs on the batch-corrected non-malignant object, alongside trajectory and
    cell communication. Answers the compositional questions those two do not:
    what the microenvironment is made of, how that differs between conditions,
    and what functional programmes the infiltrate is running.

    New in the merged pipeline — no equivalent module existed in any of the
    eight source repositories.
----------------------------------------------------------------------------------------
*/

process TME_CHARACTERISATION {

    tag "TME characterisation"
    label 'process_high'

    input:
        path seurat_object
        path notebook
        path config

    output:

        // Figures were written by the notebook but never declared as an

        // output, so publishDir had nothing to copy and figures/ stayed empty.

        path("figures/**"), emit: figure_files, optional: true
        path "data/**",                          emit: data,    optional: true
        path "figures/**",                       emit: figures, optional: true
        path "report/${notebook.baseName}.html", emit: report,  optional: true

    when:
        task.ext.when == null || task.ext.when

    script:
        def extras = task.ext.args ? "-P ${task.ext.args}" : ""
        """
        set -euo pipefail
        mkdir -p report data figures

        export OMP_NUM_THREADS=${task.cpus}
        export OPENBLAS_NUM_THREADS=${task.cpus}
        export MKL_NUM_THREADS=${task.cpus}

        quarto render ${notebook} \\
            -P seurat_object:${seurat_object} \\
            -P project_name:${params.project_name} \\
            -P work_directory:\$PWD \\
            ${extras}

        [ -f "${notebook.baseName}.html" ] && cp "${notebook.baseName}.html" report/ || true
        """

    stub:
        """
        mkdir -p report data figures
        echo "celltype,sample,n,frac" > data/tme_composition.csv
        touch figures/tme_stub.png
        echo "<html><body>stub</body></html>" > report/${notebook.baseName}.html
        """
}
