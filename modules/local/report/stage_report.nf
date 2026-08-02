/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    STAGE_REPORT — the smart report, for every stage
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    Runs after an analysis stage and builds its report from what that stage
    PUBLISHED — its tables, figures and any findings.json — rather than from
    inside its notebook.

    Why this way round: the 40 analysis notebooks are validated science that
    cannot be executed here to check an edit. Rewriting each one to emit a smart
    report would be a large unverifiable change. Reading their published outputs
    is generic, testable, and cannot break the analysis. Stages keep their own
    detailed notebook report; this adds the summary-first layer on top.

    Every stage therefore gets: a computed verdict, KPI tiles, findings ranked
    worst-first, searchable tables, an interactive metric explorer with linked
    filters, a figure browser, the offline agent, and findings.json.
----------------------------------------------------------------------------------------
*/

process STAGE_REPORT {

    tag "Report: ${stage}"
    label 'process_single'

    // A reporting failure must never discard a completed analysis stage.
    errorStrategy 'ignore'

    input:
        tuple val(stage), val(stage_title), path(stage_input, stageAs: 'stage_input/*')
        path notebook
        path config

    output:
        tuple val(stage), path("report/${stage}_report.html"), emit: report,   optional: true
        path "data/findings.json",                             emit: findings, optional: true

    when:
        task.ext.when == null || task.ext.when

    script:
        """
        set -euo pipefail
        mkdir -p report data

        quarto render ${notebook} \\
            -P stage:'${stage}' \\
            -P stage_title:'${stage_title}' \\
            -P project_name:'${params.project_name}' \\
            -P stage_dir:stage_input \\
            -P celltype_col:'${params.tme_celltype_col}' \\
            -P sample_col:'${params.tme_sample_col}' \\
            -P work_directory:\$PWD

        [ -f "${notebook.baseName}.html" ] && mv "${notebook.baseName}.html" "report/${stage}_report.html" || true
        """

    stub:
        """
        mkdir -p report data
        echo '<html><body>stub stage report: ${stage}</body></html>' > report/${stage}_report.html
        echo '[{"id":"stub","severity":"good","headline":"stub","stage":"${stage}"}]' > data/findings.json
        """
}
