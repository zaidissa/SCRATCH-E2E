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

    // helper to single-quote values and escape any single quotes inside
    def q = { v -> "'${v.toString().replace("'", "'\\''")}'" }

    // build flat "-P key:value" list with safe quoting for things containing commas/semicolons
    def parts = []
    parts << "-P seurat_object:${seurat_object}"
    parts << "-P project_name:${params.project_name}"
    parts << "-P input_integration_method:${params.bc_integration_method}"
    parts << "-P input_target_variables:${params.bc_target_variables.replaceAll(',', ';')}"
    parts << "-P input_batch_step:${params.bc_batch_step}"
    parts << "-P exclude_labels:${params.bc_exclude_labels.replaceAll(',', ';')}"
    // for label_candidates, include quotes because of semicolons
    parts << "-P label_candidates:'${params.bc_label_candidates.replaceAll(',', ';')}'"
    parts << "-P n_hvgs:${params.bc_n_hvgs}"
    parts << "-P n_pcs:${params.bc_n_pcs}"
    parts << "-P n_threads:${task.cpus}"
    parts << "-P n_memory:${task.memory.toGiga()}G"
    // IMPORTANT: do NOT quote $PWD so the shell expands it
    parts << "-P work_directory:$PWD"

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

