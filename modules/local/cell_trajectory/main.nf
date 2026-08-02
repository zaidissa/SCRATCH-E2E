process CELLTRAJECTORY {

    tag "Trajectory inference (${params.traj_species})"
    label 'process_high'

    /*
     * Rewritten on merge. The vendored version had four defects that prevented
     * it running outside its author's environment:
     *
     *  1. `shell:` was used with `${...}` interpolation. Under `shell:`,
     *     `${var}` is a SHELL variable, not a Nextflow one, so `${seurat_object}`
     *     expanded to the empty string and quarto received no input object.
     *     Switched to `script:`, where `${...}` is Nextflow and `\$` is shell.
     *  2. `cpus task.cpus ?: 16` and
     *     `memory { task.memory.toGiga() ? "${task.memory.toGiga()}G GB" : ... }`
     *     were self-referential (a directive reading its own value) and the
     *     memory branch built the malformed string "250G GB". Resources now come
     *     from conf/base.config.
     *  3. `path "${params.outdir}/**"` declared an output glob on an absolute
     *     path outside the task directory, which Nextflow cannot capture.
     *  4. The notebook filename was hardcoded rather than taken from the staged
     *     `notebook` input.
     *
     * The notebook writes into `<outdir>/data/{CytoTraceOut,SlingshotOut,
     * Monocle3Out,Concordance}`; passing `-P outdir:.` lands those under `data/`
     * so this stage matches the repo-wide report/figures/data layout.
     */

    input:
      path seurat_object
      path notebook
      path config

    output:
      path "report/${notebook.baseName}.html", emit: report,  optional: true
      path "data/**",                          emit: rds,     optional: true
      path "figures/**",                       emit: figures, optional: true

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

      quarto render ${notebook} --execute \\
        -P seurat_object:${seurat_object} \\
        -P project_name:${params.project_name} \\
        -P outdir:. \\
        -P n_threads:${task.cpus} \\
        -P n_memory:${task.memory.toGiga()} \\
        -P ct_species:${params.traj_species} \\
        -P seed:${params.seed} \\
        -P work_directory:\$PWD \\
        ${extras}

      [ -f "${notebook.baseName}.html" ] && cp "${notebook.baseName}.html" report/ || true
      """

    stub:
      """
      mkdir -p report data figures
      touch report/${notebook.baseName}.html
      touch data/trajectory_stub.rds
      """
}
