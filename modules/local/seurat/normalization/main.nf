process SEURAT_NORMALIZE {

    tag "Running normalization and dimensionality reduction"
    label 'process_high'


    input:
        path(seurat_object)
        // BPCells on-disk counts store; staged at the same relative path the
        // object's lazy layers reference, so they resolve without realizing
        // the matrix in memory.
        path(bpcells_store, stageAs: 'data/bpcells_counts')
        path(notebook_normalize)
        path(page_config)

    output:

        // Figures were written by the notebook but never declared as an

        // output, so publishDir had nothing to copy and figures/ stayed empty.

        path("figures/**"), emit: figure_files, optional: true
        path("data/${params.project_name}_reduction_object.RDS"), emit: seurat_rds
        path("report/${notebook_normalize.baseName}.html")
        // optional: the knit_print hook renders ggplots as plotly WIDGETS, so
        // quarto no longer writes figure-html PNGs for those chunks. A stage
        // whose analysis succeeded must not fail because a figure became
        // interactive.
        path("_freeze/**/figure-html/*.png"), emit: figures, optional: true

    when:
        task.ext.when == null || task.ext.when
        
    script:
        def param_file = task.ext.args ? "-P seurat_object:${seurat_object} -P ${task.ext.args}" : ""
        """
        quarto render ${notebook_normalize} ${param_file}
        """
    stub:
        """
        mkdir -p report data figures 
        mkdir -p _freeze/DUMMY/figure-html
        
        touch _freeze/DUMMY/figure-html/FILE.png

        touch data/${params.project_name}_reduction_object.RDS
        touch report/${notebook_normalize.baseName}.html

        """
}
