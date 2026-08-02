process SCEVAN {

    tag "Performing SCEVAN analysis"
    label 'process_high'


    input:
        path(seurat_object)
        // path(reference_table)
        path(notebook)
        path(config)

    output:
        path("data/scevan")                             , emit: results
        path("data/${params.project_name}_scevan_meta_object.RDS"), emit: meta_object, optional: true
        path("report/*")
        path("_freeze/${notebook.baseName}")            , emit: cache

    when:
        task.ext.when == null || task.ext.when

    script:
        def additional_args = task.ext.args ? task.ext.args.split(',').collect { "-P ${it.trim()}" }.join(' ') : ""
        def param_file = "-P seurat_object='${seurat_object.toString()}' ${additional_args}"
        
        """
        quarto render ${notebook} ${param_file}
        """
        

    stub:
        def param_file = task.ext.args ? "-P seurat_object:${seurat_object} -P ${task.ext.args}" : ""
        """
        mkdir -p report data/scevan _freeze/${notebook.baseName}
        echo "<html><body>stub</body></html>" > report/${notebook.baseName}.html
        touch data/scevan/${params.project_name}_scevan_stub.txt
        touch data/${params.project_name}_scevan_meta_object.RDS
        touch _freeze/${notebook.baseName}/${notebook.baseName}.html
        """
}
