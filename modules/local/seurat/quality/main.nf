process SEURAT_QUALITY {

    tag "Assessing quality ${sample_id}"
    label 'process_medium'

    input:
        tuple val(sample_id), path(csv_metrics), path(matrices)
        path(notebook_quality)
        path(page_config)

    output:
        tuple val(sample_id), path("objects/*"), path("log/*.txt"), emit: status
        path("${sample_id}_metrics_upgrade.csv"),                   emit: metrics
        path("report/notebook_${sample_id}.html")

    when:
        task.ext.when == null || task.ext.when

    script:
        def param_file = task.ext.args ? "-P sample_name:${sample_id} -P csv_metrics:${csv_metrics} -P input_gex_matrices:${matrices} -P ${task.ext.args}" : ""
        def notebook_sample = "notebook_${sample_id}"
        """
        mv ${notebook_quality} ${notebook_sample}.qmd
        quarto render ${notebook_sample}.qmd ${param_file}
        """
    stub:
        def param_file = task.ext.args ? "-P sample_name:${sample_id} -P csv_metrics:${csv_metrics} -P input_gex_matrices:${matrices} -P ${task.ext.args}" : ""
        def notebook_sample = "notebook_${sample_id}"
        """
        mkdir -p objects log figures report

        touch report/${notebook_sample}.html
        touch ${sample_id}_metrics_upgrade.csv

        touch log/SUCCESS.txt 
        touch objects/${sample_id}_seurat_object.RDS
        echo "${param_file}" > report/param_file.yml
        """
}