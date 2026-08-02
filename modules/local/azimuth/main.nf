process AZIMUTH_ANNOTATION {

    tag "Running Azimuth annotation"
    label 'process_medium'

    

    input:
        path(notebook)
        path(seurat_object)
        path(reference_object)
        path(config)

    output:
        path("_freeze/${notebook.baseName}")                          , emit: cache
        path("data/${params.project_name}_Azimuth_annotation_object.RDS")   , emit: seurat_rds
        path("report/${notebook.baseName}.html")                      , emit: html
        path ("figures/**")                                              , emit: figures

    when:
        task.ext.when == null || task.ext.when

    script:
        def param_file = task.ext.args ? "-P seurat_object:${seurat_object} -P reference_object:'${reference_object}' -P ${task.ext.args}" : ""
        """
        quarto render ${notebook} ${param_file}
        """
    stub:
        def param_file = task.ext.args ? "-P seurat_object:${seurat_object} -P reference_object:'${reference_object}' -P ${task.ext.args}" : ""
        """
        mkdir -p _freeze/${notebook.baseName}
        
        mkdir -p data
        touch data/${params.project_name}_Azimuth_annotation_object.RDS

        mkdir -p report
        touch report/${notebook.baseName}.html

        echo ${param_file} > _freeze/param_file.yml
        mkdir -p figures
        touch figures/azimuth_stub.png
        """

}

