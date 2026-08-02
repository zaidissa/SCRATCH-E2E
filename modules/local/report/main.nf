process QUARTO_RENDER_PROJECT {

    tag "Creating final report"
    label 'process_low'


    input:
        path(template)
        path(qmd)
        path(cache), stageAs: '_freeze/*'

    output:
        path("report"), emit: project_folder

    shell:
        """
        quarto render .
        """


    stub:
        """
        touch report
        """
}
