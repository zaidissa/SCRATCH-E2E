process HELPER_SUMMARIZE {

    tag "Rendering QC Table"
    label 'process_single'


    input:
        path(project_metrics)
        path(notebook_summarize)
        path(page_config)

    output:
        path("report/${notebook_summarize.baseName}.html")        , emit: html

    when:
        task.ext.when == null || task.ext.when

    script:
        def param_file = task.ext.args ? "-P input_metrics_report:\'${project_metrics.join(';')}\' -P ${task.ext.args}" : ""
        """
        quarto render ${notebook_summarize} ${param_file}
        """
    stub:
        def param_file = task.ext.args ? "-P input_metrics_report:\'${project_metrics.join(';')}\' -P ${task.ext.args}" : ""
        """

        mkdir -p report
        touch report/${notebook_summarize.baseName}.html

        echo "${param_file}" > report/param_file.yml
        """

}