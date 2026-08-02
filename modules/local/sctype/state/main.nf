process SCYTPE_STATE_ANNOTATION {

    tag "Running scType annotation - ${cell_population} subsets/states"
    label 'process_medium'

    

    input:
        path(notebook)
        path(seurat_object)
        path(cell_annotation)
        each(cell_population)
        path(config)

    output:
        path("_freeze/notebook_${cell_population}")                                   , emit: cache
        path("data/${params.project_name}_${cell_population}_annotation_object.RDS")  , emit: seurat_rds
        path("data/${params.project_name}_${cell_population}_annotation.csv")         , emit: annotation
        path("report/notebook_${cell_population}.html")                               , emit: html
        path ("figures/**")                                              , emit: figures
        // path("_freeze/**/figure-html/*")                                          , emit: figures
        
    when:
        task.ext.when == null || task.ext.when

    script:
        def param_file = task.ext.args ? "-P seurat_object:${seurat_object} -P input_cell_markers_db:${cell_annotation} -P input_parent_level:'${cell_population}' -P ${task.ext.args}" : ""
        def notebook_cell_type = "notebook_${cell_population}"
        """
        mv ${notebook} ${notebook_cell_type}.qmd
        quarto render ${notebook_cell_type}.qmd ${param_file}
        """
    stub:
        def param_file = task.ext.args ? "-P seurat_object:${seurat_object} -P input_cell_markers_db:${cell_annotation} -P input_parent_level:'${cell_population}' -P ${task.ext.args}" : ""
        def notebook_cell_type = "notebook_${cell_population}"
        """
        mkdir -p data _freeze/${notebook_cell_type}
        mkdir -p _freeze/DUMMY/figure-html

        touch _freeze/DUMMY/figure-html/FILE.png

        touch data/${params.project_name}_${cell_population}_annotation_object.RDS
        touch data/${params.project_name}_${cell_population}_annotation.csv
        
        mkdir -p report
        touch report/${notebook_cell_type}.html

        echo ${param_file} > _freeze/param_file.yml
        mkdir -p figures
        touch figures/sctype_state_${cell_population}_stub.png
        """

}

