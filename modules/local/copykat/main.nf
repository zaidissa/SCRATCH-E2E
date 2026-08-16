process COPYKAT_PROCESS {

    tag "Performing COPYKAT analysis"
    label 'process_medium'
    

    input:
        path seurat_rds
        path notebook_copykat
        path page_config

    output:

        // Figures were written by the notebook but never declared as an

        // output, so publishDir had nothing to copy and figures/ stayed empty.

        path("figures/**"), emit: figures, optional: true
        path "copykat_prediction.csv"
        path "seurat_copykat_annotated.rds"
        path "umap_copykat.png"
        path "copykat_result.rds"

    script:
        def param_file = task.ext.args ? "-P ${task.ext.args}" : ""
        """
        cp $seurat_rds input_seurat.rds
        cp $notebook_copykat notebook_copykat.qmd
        quarto render notebook_copykat.qmd --execute --quiet ${param_file}
        """

    stub:
        """
        echo "stub" > copykat_prediction.csv
        touch copykat_result.rds
        touch seurat_copykat_annotated.rds
        touch umap_copykat.png
        """
}

workflow COPYKAT {
    take:
        seurat_rds
        notebook_copykat
        page_config

    main:
        COPYKAT_PROCESS(
            seurat_rds,
            notebook_copykat,
            page_config
        )

}
