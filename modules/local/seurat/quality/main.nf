process SEURAT_QUALITY {

    tag "Assessing quality ${sample_id}"
    label 'process_medium'

    input:
        // ambient_status: cellbender | skipped_by_request | cellbender_failed |
        // disabled — set per sample in qc.nf so the report can say which counts
        // this sample was analysed from.
        tuple val(sample_id), path(csv_metrics), path(matrices), val(ambient_status)
        path(notebook_quality)
        path(page_config)

    output:
        tuple val(sample_id), path("objects/*"), path("log/*.txt"), emit: status
        path("${sample_id}_metrics_upgrade.csv"),                   emit: metrics
        path("report/notebook_${sample_id}.html")
        // This process never declared a figures output, so per-sample QC
        // figures were computed and then discarded. Optional: a sample that
        // fails QC early can legitimately produce none.
        path("figures/**"), emit: figures, optional: true

    when:
        task.ext.when == null || task.ext.when

    script:
        def param_file = task.ext.args ? "-P sample_name:${sample_id} -P csv_metrics:${csv_metrics} -P input_gex_matrices:${matrices} -P ambient_correction:${ambient_status} -P ${task.ext.args}" : ""
        def notebook_sample = "notebook_${sample_id}"
        """
        mv ${notebook_quality} ${notebook_sample}.qmd
        quarto render ${notebook_sample}.qmd ${param_file}
        """
    stub:
        def param_file = task.ext.args ? "-P sample_name:${sample_id} -P csv_metrics:${csv_metrics} -P input_gex_matrices:${matrices} -P ambient_correction:${ambient_status} -P ${task.ext.args}" : ""
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