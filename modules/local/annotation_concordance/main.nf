process ANNOTATION_CONCORDANCE {

    tag "Comparing annotators"
    label 'process_medium'

    input:
        path(notebook)
        path(celltypist_obs)
        path(sctype_object)
        path(azimuth_object)
        path(config)

    output:
        path("_freeze/${notebook.baseName}")            , emit: cache
        path("report/${notebook.baseName}.html")        , emit: html
        // Every output below is optional on purpose: with fewer than two
        // annotators enabled the notebook renders a "not enough annotators"
        // finding and writes no data. A required output would turn a
        // legitimate configuration into a pipeline failure.
        path("data/*_annotation_concordance.csv")       , emit: concordance, optional: true
        path("data/*_annotator_agreement.csv")          , emit: agreement  , optional: true
        path("figures/**")                              , emit: figures    , optional: true

    when:
        task.ext.when == null || task.ext.when

    script:
        def args = task.ext.args ? "-P ${task.ext.args}" : ""
        """
        quarto render ${notebook} \\
            -P celltypist_obs:${celltypist_obs} \\
            -P sctype_object:${sctype_object} \\
            -P azimuth_object:${azimuth_object} \\
            ${args}
        """
    stub:
        """
        mkdir -p _freeze/${notebook.baseName} report data figures
        touch report/${notebook.baseName}.html
        touch data/${params.project_name}_annotation_concordance.csv
        touch data/${params.project_name}_annotator_agreement.csv
        touch figures/annotation_concordance_stub.png
        """
}
