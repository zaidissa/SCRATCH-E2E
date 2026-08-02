process SEURAT_MERGE {

    tag "Merging post-QC samples"
    label 'process_high'


    input:
        path(qc_approved)
        path(notebook_merge)
        path(exp_table)
        path(page_config)

    output:
        path("data/${params.project_name}_qc_merged_object.RDS"), emit: seurat_rds
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
        mkdir -p report data figures/merge

        touch data/${params.project_name}_qc_merged_object.RDS
        touch report/${notebook_merge.baseName}.html

        """
}
