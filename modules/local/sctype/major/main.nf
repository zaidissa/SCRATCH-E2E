process SCYTPE_MAJOR_ANNOTATION {

    tag "Running scType annotation - Major cells"
    label 'process_medium'


    input:
        path(notebook_major)
        path(seurat_object)
        path(cell_annotation)
        path(config)

    output:

        // Figures were written by the notebook but never declared as an

        // output, so publishDir had nothing to copy and figures/ stayed empty.

        path("figures/**"), emit: figure_files, optional: true
        path("_freeze/${notebook_major.baseName}")                       , emit: cache
        path("data/${params.project_name}_major_annotation_object.RDS")  , emit: seurat_rds
        path("data/${params.project_name}_major_annotation.csv")         , emit: annotation
        path("data/${params.project_name}_major_annotation.list.txt")    , emit: major_list
        path("report/${notebook_major.baseName}.html")                   , emit: html
        // path("_freeze/**/figure-html/*.png")                             , emit: figures
        path ("figures/**")                                              , emit: figures

    when:
        task.ext.when == null || task.ext.when

    script:
        def param_file = task.ext.args ? "-P seurat_object:${seurat_object} -P input_cell_markers_db:${cell_annotation} -P ${task.ext.args}" : ""
        """
        quarto render ${notebook_major} ${param_file}
        """
    stub:
        """
        mkdir -p data _freeze/${notebook_major.baseName}
        mkdir -p _freeze/DUMMY/figure-html

        touch _freeze/DUMMY/figure-html/FILE.png

        touch data/${params.project_name}_major_annotation_object.RDS
        touch data/${params.project_name}_major_annotation.csv

        mkdir -p figures
        touch figures/FILE.png

        mkdir -p report
        touch report/${notebook_major.baseName}.html

        echo "B_Plasma_Cells\nEndothelial_Cells\nEpithelial\nFibroblast\nMyeloid\nNK_Cells\nT_Cells" > data/${params.project_name}_major_annotation.list.txt
        """

}

