process INFERCNV {

    tag "Performing INFERCNV analysis"
    label 'process_high'
    label 'process_high_memory'


    input:
        path(seurat_object)
        // path(reference_table)
        path(notebook)
        path(config)


    output:
        path("data/infercnv")                             , emit: results
        path("data/${params.project_name}_infercnv_meta_object.RDS"), emit: meta_object, optional: true
        path("report/*")
        path("_freeze/${notebook.baseName}")              , emit: cache

    when:
        task.ext.when == null || task.ext.when

    script:
        def param_file = task.ext.args ? "-P seurat_object:${seurat_object}  -P ${task.ext.args}" : ""
        // def param_file = task.ext.args ? "-P seurat_object:${seurat_object} -P reference_table:${reference_table} -P ${task.ext.args}" : ""
        """
        quarto render ${notebook} ${param_file}
        """
    stub:
        def param_file = task.ext.args ? "-P seurat_object:${seurat_object}  -P ${task.ext.args}" : ""
        """
        mkdir -p data/infercnv _freeze/${notebook.baseName}
        mkdir -p _freeze/DUMMY/figure-html

        touch _freeze/DUMMY/figure-html/FILE.png

        touch data/infercnv/${params.project_name}_infercnv_meta_object.RDS
        touch data/${params.project_name}_infercnv_meta_object.RDS
        touch _freeze/${notebook.baseName}/${notebook.baseName}.html

        mkdir -p report
        touch report/${notebook.baseName}.html

        echo ${param_file} > _freeze/${notebook.baseName}/params.yml
        """

}
