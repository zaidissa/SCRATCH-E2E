process GENEVECTOR_ANNOTATION {

    tag "Performing analysis ${notebook.baseName}"
    label 'process_medium'


    input:
        path(notebook)
        path(config)

        val(project_name)
        val(paramA)

    output:
        path("_freeze/${notebook.baseName}"),   emit: cache

    when:
        task.ext.when == null || task.ext.when

    script:
        def project_name = project_name ? "-P project_name:${project_name}" : ""
        def paramA       = paramA       ? "-P paramA:${paramA}" : ""
        """
        quarto render ${notebook} ${project_name} ${paramA}
        """
    stub:
        """
        mkdir -p _freeze/${notebook.baseName}
        """

}
