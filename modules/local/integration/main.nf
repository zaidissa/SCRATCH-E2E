/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Cross-stage integration layer
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    New in the merged pipeline. None of this was possible when the eight
    pipelines ran separately, because no single process ever saw more than one
    stage's output.

      INTEGRATE_OBJECTS  join every stage's per-cell annotation back onto one
                         master object (CNV calls, metaprograms, pseudotime,
                         clonotypes) and export a tidy per-cell table
      CROSS_ANALYSIS     analyses that need two or more stages at once
      RUN_MANIFEST       parameters, versions, timings, per-stage cell counts
      MASTER_REPORT      one navigable Quarto site over all of the above
----------------------------------------------------------------------------------------
*/

process INTEGRATE_OBJECTS {

    tag "Integrating stage annotations"
    label 'process_high'

    input:
        path seurat_object
        path notebook
        path config
        path cnv_infercnv,  stageAs: 'stage_inputs/cnv_infercnv/*'
        path metaprog,      stageAs: 'stage_inputs/metaprog/*'
        path trajectory,    stageAs: 'stage_inputs/trajectory/*'
        path tcr,           stageAs: 'stage_inputs/tcr/*'
        path compartment,   stageAs: 'stage_inputs/compartment/*'
        path annotation,    stageAs: 'stage_inputs/annotation/*'

    output:

        // Figures were written by the notebook but never declared as an

        // output, so publishDir had nothing to copy and figures/ stayed empty.

        path("figures/**"), emit: figure_files, optional: true
        path "data/${params.project_name}_master_object.RDS", emit: master_rds, optional: true
        path "data/${params.project_name}_cell_table.csv",    emit: cell_table, optional: true
        path "data/**",                                       emit: data,    optional: true
        path "figures/**",                                    emit: figures, optional: true
        path "report/${notebook.baseName}.html",              emit: report,  optional: true

    when:
        task.ext.when == null || task.ext.when

    script:
        def extras = task.ext.args ? "-P ${task.ext.args}" : ""
        """
        set -euo pipefail
        mkdir -p report data figures

        export OMP_NUM_THREADS=${task.cpus}
        export OPENBLAS_NUM_THREADS=${task.cpus}
        export MKL_NUM_THREADS=${task.cpus}

        quarto render ${notebook} \\
            -P seurat_object:${seurat_object} \\
            -P project_name:${params.project_name} \\
            -P stage_inputs:stage_inputs \\
            -P work_directory:\$PWD \\
            ${extras}

        [ -f "${notebook.baseName}.html" ] && cp "${notebook.baseName}.html" report/ || true
        """

    stub:
        """
        mkdir -p report data figures
        touch data/${params.project_name}_master_object.RDS
        echo "cell_id,sample,celltype" > data/${params.project_name}_cell_table.csv
        touch report/${notebook.baseName}.html
        """
}


process CROSS_ANALYSIS {

    tag "Cross-stage analysis"
    label 'process_high'

    input:
        path cell_table
        path notebook
        path config
        path stage_data, stageAs: 'stage_inputs/*'

    output:
        path "data/**",                          emit: data,    optional: true
        path "figures/**",                       emit: figures, optional: true
        path "report/${notebook.baseName}.html", emit: report,  optional: true

    when:
        task.ext.when == null || task.ext.when

    script:
        def extras = task.ext.args ? "-P ${task.ext.args}" : ""
        """
        set -euo pipefail
        mkdir -p report data figures

        export OMP_NUM_THREADS=${task.cpus}
        export OPENBLAS_NUM_THREADS=${task.cpus}

        quarto render ${notebook} \\
            -P cell_table:${cell_table} \\
            -P stage_inputs:stage_inputs \\
            -P work_directory:\$PWD \\
            ${extras}

        [ -f "${notebook.baseName}.html" ] && cp "${notebook.baseName}.html" report/ || true
        """

    stub:
        """
        mkdir -p report data figures
        touch report/${notebook.baseName}.html
        touch data/cross_analysis_stub.csv
        """
}


process RUN_MANIFEST {

    tag "Run manifest"
    label 'process_single'

    // Provenance capture. Reads the trace file Nextflow writes for this run
    // (see the `trace` scope in nextflow.config) plus the resolved parameters,
    // so the master report can show what was actually executed rather than what
    // the defaults say.

    input:
        val  ready          // ordering barrier: everything upstream has finished
        path params_json

    output:
        path "manifest/run_manifest.json", emit: manifest
        path "manifest/*.csv",             emit: tables, optional: true

    script:
        """
        set -euo pipefail
        mkdir -p manifest
        cp ${params_json} manifest/resolved_params.json

        cat > manifest/run_manifest.json <<'MANIFEST_EOF'
{
  "pipeline":       "${workflow.manifest.name}",
  "version":        "${workflow.manifest.version}",
  "project_name":   "${params.project_name}",
  "run_name":       "${workflow.runName}",
  "session_id":     "${workflow.sessionId}",
  "start":          "${workflow.start}",
  "profile":        "${workflow.profile}",
  "container_engine": "${workflow.containerEngine ?: 'none'}",
  "nextflow_version": "${nextflow.version}",
  "command_line":   "${workflow.commandLine.replace('"', '\\\\"')}",
  "revision":       "${workflow.revision ?: 'n/a'}",
  "organism":       "${params.organism}",
  "genome":         "${params.genome}",
  "nmf_engine":     "${params.metaprog_nmf_engine}",
  "from":           "${params.from}",
  "to":             "${params.to}"
}
MANIFEST_EOF

        # Per-task resource usage, if the trace file is readable from here.
        TRACE=\$(ls -t ${params.outdir}/pipeline_info/trace_*.txt 2>/dev/null | head -1 || true)
        if [ -n "\$TRACE" ] && [ -f "\$TRACE" ]; then
            cp "\$TRACE" manifest/trace.csv
        fi
        """

    stub:
        """
        mkdir -p manifest
        echo '{"pipeline":"SCRATCH-E2E","stub":true}' > manifest/run_manifest.json
        """
}


process MASTER_REPORT {

    tag "Master report"
    label 'process_single'

    // Assembles the single navigable site. Never allowed to fail the pipeline:
    // conf/base.config sets errorStrategy 'ignore' for this process, because a
    // formatting problem in the report must not discard a completed analysis.

    input:
        path notebook
        path config
        path manifest,     stageAs: 'inputs/manifest/*'
        path stage_reports, stageAs: 'inputs/reports/*'
        path stage_figures, stageAs: 'inputs/figures/*'
        path cross_data,   stageAs: 'inputs/cross/*'

    output:
        // quarto renders the site into report/; master_report/ holds the
        // assembled figures and per-stage pages. Only the latter was declared,
        // so `master_report/index.html` -- the path the completion banner sends
        // users to -- was never published and the run still reported success.
        path "master_report/**", emit: report, optional: true
        path "report/**",        emit: site,   optional: true

    when:
        task.ext.when == null || task.ext.when

    script:
        def extras = task.ext.args ? "-P ${task.ext.args}" : ""
        """
        set -euo pipefail
        mkdir -p master_report

        quarto render ${notebook} \\
            -P project_name:'${params.report_title_final}' \\
            -P inputs_dir:inputs \\
            -P interactive:${params.report_interactive} \\
            -P work_directory:\$PWD \\
            ${extras}

        # Quarto writes <notebook>.html; the site entry point must be index.html.
        if [ -f "${notebook.baseName}.html" ]; then
            cp "${notebook.baseName}.html" master_report/index.html
        fi

        # Carry the per-stage reports and figures into the site so it is
        # self-contained and can be copied or served anywhere.
        [ -d inputs/reports ] && cp -r inputs/reports master_report/stages || true
        [ -d inputs/figures ] && cp -r inputs/figures master_report/figures || true
        """

    stub:
        """
        mkdir -p master_report
        echo '<html><body>stub master report</body></html>' > master_report/index.html
        """
}
