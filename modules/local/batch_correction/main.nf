process BATCHCORRECTION {

    tag "Batch correction (${params.bc_integration_method})"
    label 'process_high'

    // cpus/memory intentionally NOT pinned here — conf/base.config owns the
    // resource model so it can be scaled per site and escalated on retry. The
    // previous hardcoded `cpus 4 / memory '16 GB'` silently capped this step
    // well below what integration of a large merged object needs.

    input:
    path seurat_object
    path notebook
    path config

    output:

        // Figures were written by the notebook but never declared as an

        // output, so publishDir had nothing to copy and figures/ stayed empty.

        path("figures/**"), emit: figure_files, optional: true
    path "report/${notebook.baseName}.html", emit: report,  optional: true
    path "figures/**",                     emit: figures, optional: true
    path "data/**",                        emit: data,    optional: true

    // Explicit chaining handle. The notebook's autosave block writes this
    // method-independent object; the per-method files (postCCA/postRPCA/
    // postHARMONY/postMNN) stay inside `data` for inspection.
    path "data/${params.project_name}_${params.bc_batch_step}_batch_object.RDS", emit: seurat_rds, optional: true

    when:
    task.ext.when == null || task.ext.when

    /*
    * Build Quarto -P params safely (quote values!).
    */
    script:

    // Quoting applied inline. A `def q = { }` closure referenced further down
    // fails on Nextflow 26.x with "`q` is not defined".

    // build flat "-P key:value" list with safe quoting for things containing commas/semicolons
    def parts = []
    parts << "-P seurat_object:${seurat_object}"
    parts << "-P project_name:${params.project_name}"
    parts << "-P input_integration_method:${params.bc_integration_method}"
    // Quote every delimited value with q(). The helper was defined above but
    // never applied: bc_exclude_labels defaults to '', which emitted a bare
    // "-P exclude_labels:" with nothing after the colon. Quarto reads that as
    // NULL, and the notebook's strsplit(NULL, ";") aborts the whole render with
    // "non-character argument". Quoting makes an empty value arrive as ''.
    parts << "-P input_target_variables:'${params.bc_target_variables.replaceAll(',', ';')}'"
    parts << "-P input_batch_step:${params.bc_batch_step}"
    parts << "-P exclude_labels:'${params.bc_exclude_labels.replaceAll(',', ';')}'"
    parts << "-P label_candidates:'${params.bc_label_candidates.replaceAll(',', ';')}'"
    // Evaluation-only downsample. Previously the notebook's own default (0.30)
    // was the only way to set it, and it silently shrank the PUBLISHED object.
    parts << "-P thr_cell_proportion:${params.bc_eval_cell_proportion}"
    parts << "-P n_hvgs:${params.bc_n_hvgs}"
    parts << "-P n_pcs:${params.bc_n_pcs}"
    parts << "-P n_threads:${task.cpus}"
    parts << "-P n_memory:${task.memory.toGiga()}G"
    // IMPORTANT: do NOT quote $PWD so the shell expands it
    parts << "-P work_directory:\$PWD"

    def param_file = parts.join(' ')

    """
    set -euo pipefail

    # ensure these roots exist so the globs always match
    mkdir -p report figures data

    # render (options first, then params)
    quarto render --execute ${notebook} ${param_file}
    """

    stub:
    """
    mkdir -p report figures data
    echo "<html><body>stub</body></html>" > report/${notebook.baseName}.html
    touch figures/batchcorrection_stub.png
    touch data/${params.project_name}_${params.bc_batch_step}_batch_object.RDS
    """
}

