process QUARTO_RENDER_PROJECT {

    tag "Creating final report"
    label 'process_low'


    input:
        path(template)
        path(qmd)
        path(cache), stageAs: '_freeze/*'

    output:

        // Figures were written by the notebook but never declared as an

        // output, so publishDir had nothing to copy and figures/ stayed empty.

        path("figures/**"), emit: figures, optional: true
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
