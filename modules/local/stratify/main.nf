/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    STRATIFY — split the cohort into malignant and non-malignant compartments
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    Sits between CNV calling and the two downstream compartments:

        ANNOTATION -> CNV -> STRATIFY -+-> malignant     -> METAPROG
                                       +-> non-malignant -> BATCHCORRECT -> ...

    Evidence is combined per cell from whichever CNV callers ran, with the
    annotation label as a tie-break and as cover for cells no caller scored.
    Cells whose evidence conflicts are labelled `ambiguous` and, by default,
    are excluded from BOTH compartments rather than contaminating either — they
    remain in the assignment table and in the master object so nothing is
    silently discarded.

    This was previously implicit and duplicated: notebook_batch_correctionR.qmd
    built its own `obj_non_tumor` by excluding the "Epithelial" label, and the
    metaprogram stage separately subset to "Epithelial". Two different notebooks
    deciding what a tumour cell is, from a lineage proxy rather than from the
    CNV evidence, with no record of the decision.
----------------------------------------------------------------------------------------
*/

process STRATIFY_COMPARTMENTS {

    tag "Malignant / non-malignant split (${params.stratify_method})"
    label 'process_high'

    input:
        path seurat_object
        path notebook
        path config
        path cnv_evidence, stageAs: 'cnv_evidence/*'

    output:

        // Figures were written by the notebook but never declared as an

        // output, so publishDir had nothing to copy and figures/ stayed empty.

        path("figures/**"), emit: figure_files, optional: true
        path "data/${params.project_name}_malignant_object.RDS",     emit: tumor_rds,      optional: true
        path "data/${params.project_name}_nonmalignant_object.RDS",  emit: nonmalignant_rds, optional: true
        path "data/${params.project_name}_malignancy_assignment.csv", emit: assignment
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
            -P cnv_evidence_dir:cnv_evidence \\
            -P work_directory:\$PWD \\
            ${extras}

        [ -f "${notebook.baseName}.html" ] && cp "${notebook.baseName}.html" report/ || true
        """

    stub:
        """
        mkdir -p report data figures
        touch data/${params.project_name}_malignant_object.RDS
        touch data/${params.project_name}_nonmalignant_object.RDS
        echo "cell_id,malignant_status,evidence" > data/${params.project_name}_malignancy_assignment.csv
        echo "CELL1,malignant,stub" >> data/${params.project_name}_malignancy_assignment.csv
        touch figures/stratify_stub.png
        echo "<html><body>stub</body></html>" > report/${notebook.baseName}.html
        """
}
