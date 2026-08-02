process METATIME_ANNOTATION {

    tag ""
    label 'process_medium'


    input:
        path(notebook)
        path(config)

    output:
        path("_freeze/${notebook.baseName}"),   emit: cache

    when:
        task.ext.when == null || task.ext.when

    script:
        def project_name = project_name ? "-P project_name:${project_name}" : ""
        """
        quarto render ${notebook} ${project_name} ${paramA}
        """
    stub:
        """
        mkdir -p _freeze/${notebook.baseName}
        """

}
