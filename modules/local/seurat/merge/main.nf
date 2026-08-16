process SEURAT_MERGE {

    tag "Merging post-QC samples"
    label 'process_high'


    input:
        path(qc_approved)
        path(notebook_merge)
        path(exp_table)
        path(page_config)

    output:

        // Figures were written by the notebook but never declared as an

        // output, so publishDir had nothing to copy and figures/ stayed empty.

        path("figures/**"), emit: figures, optional: true
        path("data/${params.project_name}_qc_merged_object.RDS"), emit: seurat_rds
        // On-disk BPCells counts store backing the object, so the next stage
        // can stage it at the same relative path.
        path("data/bpcells_counts"), emit: bpcells_store
        path("report/${notebook_merge.baseName}.html")

    when:
        task.ext.when == null || task.ext.when
        
    script:
        def param_file = task.ext.args ? "-P input_qc_approved:\'${qc_approved.join(';')}\' -P input_exp_table:${exp_table} -P ${task.ext.args}" : ""
        """
        quarto render ${notebook_merge} ${param_file}
        """
    stub:
        def param_file = task.ext.args ? "-P input_qc_approved:\'${qc_approved.join(';')}\' -P input_exp_table:${exp_table} -P ${task.ext.args}" : ""
        """
        mkdir -p report data/bpcells_counts figures/merge
        touch data/bpcells_counts/matrix.mtx

        touch data/${params.project_name}_qc_merged_object.RDS
        touch report/${notebook_merge.baseName}.html

        """
}
