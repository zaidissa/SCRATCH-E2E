process SEURAT_CLUSTER {

    tag "Clustering cells"
    label 'process_high'


    input:
        path(seurat_object)
        path(notebook_clustering)
        path(page_config)

    output:
        path("data/${params.project_name}_cluster_object.RDS"), emit: seurat_rds
        path("report/${notebook_clustering.baseName}.html")
        path("_freeze/**/figure-html/*.png"), emit: figures

    when:
        task.ext.when == null || task.ext.when
        
    script:
        def param_file = task.ext.args ? "-P seurat_object:${seurat_object} -P ${task.ext.args}" : ""
        """
        quarto render ${notebook_clustering} ${param_file}
        """
    stub:
        """
        mkdir -p report data figures 
        mkdir -p _freeze/DUMMY/figure-html
        
        touch _freeze/DUMMY/figure-html/FILE.png

        touch data/${params.project_name}_cluster_object.RDS
        touch report/${notebook_clustering.baseName}.html

        """

}