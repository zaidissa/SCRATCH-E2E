process CELLTYPIST_ANNOTATION {

    tag "Running Celltypist annotation"
    label 'process_medium'

    

    input:
        path(notebook)
        path(anndata_object)
        path(config)

    output:

        // Figures were written by the notebook but never declared as an

        // output, so publishDir had nothing to copy and figures/ stayed empty.

        path("figures/**"), emit: figure_files, optional: true
        path("_freeze/${notebook.baseName}")                                  , emit: cache
        path("data/${params.project_name}_celltypist_annotation_object.h5ad") , emit: ann_object
        path("data/Immune_All")                                               , emit: csv_file
        path("report/${notebook.baseName}.html")                              , emit: html
        // path ("figures/**")                                              , emit: figures
                // Jupyter-engine notebooks with embed-resources inline figures as
        // base64 rather than writing _freeze/**/figure-html/*.png, so this
        // was a REQUIRED output that quarto never produces here. The stub
        // branch touches a DUMMY png, which is why -stub never caught it.
        path("_freeze/**/figure-html/*.png")                                  , emit: figures, optional: true
    when:
        task.ext.when == null || task.ext.when

    script:
        def param_file = task.ext.args ? "-P anndata_object:${anndata_object} -P ${task.ext.args}" : ""
        """
        quarto render ${notebook} ${param_file}
        """
    stub:
        def param_file = task.ext.args ? "-P anndata_object:${anndata_object} -P ${task.ext.args}" : ""
        """
        mkdir -p data _freeze/${notebook.baseName}
        mkdir -p _freeze/DUMMY/figure-html

        touch _freeze/DUMMY/figure-html/FILE.png


        touch data/${params.project_name}_celltypist_annotation_object.h5ad
        mkdir -p _freeze/${notebook.baseName} data/Immune_All
        touch _freeze/${notebook.baseName}/${notebook.baseName}.html

        mkdir -p report
        touch report/${notebook.baseName}.html

        echo ${param_file} > _freeze/${notebook.baseName}/params.yml
        """

}
