// modules/local/MetaProg-DE/main.nf
nextflow.enable.dsl = 2

/* =============================================================================
   Helper: staged publish paths
   ========================================================================== */
def STAGE_PATH(stage) {
  return "${params.outdir}/${params.project_name}/metaprog/${stage}"
}

/* =============================================================================
   Process: METAPROG_LEIDEN  (single)
   - Expects the annotated Seurat object (RDS)
   - Metaprog_leiden.qmd
   
   ========================================================================== */
process METAPROG_LEIDEN {

  tag "Metaprog: LEIDEN"
  label 'process_medium'
  errorStrategy 'terminate'
  maxRetries 1


  input:
    tuple path(seurat_object), path(notebook)

  output:

      // Figures were written by the notebook but never declared as an

      // output, so publishDir had nothing to copy and figures/ stayed empty.

      path("figures/**"), emit: figure_files, optional: true
    // path "_freeze/${notebook.baseName}"             , emit: cache,      optional: true
    path "report/${notebook.baseName}.html"         , emit: report,     optional: true
    path "figures/metaprog/**"                      , emit: figures,    optional: true
    path "data/**"                                  , emit: data,       optional: true
    path "*_leiden.rds"                             , emit: leiden_rds

  when:
    task.ext.when == null || task.ext.when

  script:
    def extras = task.ext.args ? "-P ${task.ext.args}" : ""
    """
    set -euo pipefail

    export OMP_NUM_THREADS=\${OMP_NUM_THREADS:-${task.cpus}}
    export OPENBLAS_NUM_THREADS=\${OPENBLAS_NUM_THREADS:-${task.cpus}}
    export MKL_NUM_THREADS=\${MKL_NUM_THREADS:-${task.cpus}}
    export R_PARALLEL_NUM_THREADS=\${R_PARALLEL_NUM_THREADS:-${task.cpus}}

    quarto render ${notebook} \\
      -P seurat_object:${seurat_object} \\
      -P project_name:${params.project_name} \\
      -P work_directory:\$PWD \\
      ${extras}

    mkdir -p report
      [ -f "${notebook.baseName}.html" ] && cp "${notebook.baseName}.html" report/
    """


    stub:
        """
        mkdir -p _freeze data figures/metaprog report
        touch _freeze/${notebook.baseName}
        touch data/stub
        touch figures/metaprog/stub
        echo "<html><body>stub</body></html>" > report/${notebook.baseName}.html
        touch stub_leiden.rds
        """
}


/* =============================================================================
   Process: METAPROG_NMF_PREP  (single → returns many *_preprocessed.rds)
   - Runs Metaprog_NMFprepprocessing.qmd
   - Emits one or more *_preprocessed.rds -> drives the fan-out in NMF step
   ========================================================================== */
process METAPROG_NMF_PREP {

  tag "Metaprog: NMF-preprocessing"
  label 'process_medium'
  errorStrategy 'terminate'
  maxRetries 1



  input:
    tuple path(seurat_object), path(notebook)

  output:
    // path "_freeze/${notebook.baseName}"                       , emit: cache,  optional: true
    path "report/${notebook.baseName}.html"                   , emit: report,  optional: true
    path "figures/metaprog/**"                                , emit: figures, optional: true
    path "data/**"                                            , emit: data,    optional: true
    // path "data/per_sample_mat/*_preprocessed.rds"             , emit: preprocessed
    path "data/per_sample_mat/*_preprocessed.rds"        , emit: preprocessed

  when:
    task.ext.when == null || task.ext.when

  script:
    def extras = task.ext.args ? "-P ${task.ext.args}" : ""
    """
    set -euo pipefail

    # threads inside the container follow Nextflow cpus
    export OMP_NUM_THREADS=\${OMP_NUM_THREADS:-${task.cpus}}
    export OPENBLAS_NUM_THREADS=\${OPENBLAS_NUM_THREADS:-${task.cpus}}
    export MKL_NUM_THREADS=\${MKL_NUM_THREADS:-${task.cpus}}
    export R_PARALLEL_NUM_THREADS=\${R_PARALLEL_NUM_THREADS:-${task.cpus}}


    quarto render ${notebook} \\
      -P seurat_object:${seurat_object} \\
      -P project_name:${params.project_name} \\
      -P work_directory:\$PWD \\
      ${extras}

    mkdir -p report
    [ -f "${notebook.baseName}.html" ] && cp "${notebook.baseName}.html" report/
    """


    stub:
        """
        mkdir -p _freeze data data/per_sample_mat figures/metaprog report
        touch _freeze/${notebook.baseName}
        touch data/per_sample_mat/stub_preprocessed.rds
        touch data/stub
        touch figures/metaprog/stub
        echo "<html><body>stub</body></html>" > report/${notebook.baseName}.html
        """
}

/* =============================================================================
   Process: METAPROG_NMF  (parallel: one job per *_preprocessed.rds)
   - Runs Metaprog_NMF.qmd on each preprocessed RDS
   - Emits *_nmf_fit.rds  (matches your QMD saveRDS naming)
   ========================================================================== */


process METAPROG_NMF {

  tag "Metaprog: NMF"
  label 'process_medium'
  errorStrategy 'terminate'
  maxRetries 1



  input:
    tuple val(sample_id), path(preprocessed_rds), path(notebook)
   

  output:
    path "data/per_sample_mat/nmf_fit/*_nmf_fit.rds" , emit: nmf_rds
    path "report/${notebook.baseName}.html"         , emit: report,     optional: true
    // path "report/*.html"                             , emit: report,    optional: true
    path "figures/**"                                , emit: figures,   optional: true
    path "data/**"                                   , emit: data,      optional: true
  
  when:
    task.ext.when == null || task.ext.when
  script:
    def extras = task.ext.args ? "-P ${task.ext.args}" : ""
    """
    set -euo pipefail

    export OMP_NUM_THREADS=\${OMP_NUM_THREADS:-${task.cpus}}
    export OPENBLAS_NUM_THREADS=\${OPENBLAS_NUM_THREADS:-${task.cpus}}
    export MKL_NUM_THREADS=\${MKL_NUM_THREADS:-${task.cpus}}
    export R_PARALLEL_NUM_THREADS=\${R_PARALLEL_NUM_THREADS:-${task.cpus}}

    sample_id=\$(basename "${preprocessed_rds}")
    sample_id=\${sample_id%_preprocessed.rds}

    # render this ONE sample
    quarto render ${notebook} \\
      -P preprocessed_matrix:"${preprocessed_rds}" \\
      -P project_name:${params.project_name} \\
      -P work_directory:\$PWD \\
      ${extras}

    mkdir -p report
      [ -f "${notebook.baseName}.html" ] && cp "${notebook.baseName}.html" report/
    """


    stub:
        """
        mkdir -p data data/per_sample_mat/nmf_fit figures report
        touch data/per_sample_mat/nmf_fit/stub_nmf_fit.rds
        touch data/stub
        touch figures/stub
        echo "<html><body>stub</body></html>" > report/${notebook.baseName}.html
        echo "<html><body>stub</body></html>" > report/stub.html
        """
}



/* =============================================================================
   Process: METAPROG_POST  (single, after all NMF)
   - Runs Metaprog_PostNMF.qmd
   - Your QMD scans work_directory/data/per_sample_mat/nmf_fit/*_nmf_fit.rds
     so we do NOT pass a list; the subworkflow uses collect() only to impose order
   ========================================================================== */
process METAPROG_POST {

  tag "Metaprog: PostNMF"
  label 'process_medium'
  errorStrategy 'terminate'
  maxRetries 1


  input:
    // NOTE: order matches nb_post.combine(nmf_rds.collect())
    // tuple path(notebook), val(nmf_list)
    // tuple path(notebook), path(nmf_fits)
    tuple path(notebook), path(nmf_fits, stageAs: 'data/per_sample_mat/nmf_fit/*')

  output:
    path "report/${notebook.baseName}.html" , emit: report , optional: true
    path "figures/**"                       , emit: figures   , optional: true
    path "data/**"                          , emit: data   , optional: true

  when:
    task.ext.when == null || task.ext.when

  script:
    // // nmf_list might be a single Path; coerce to a Java List first
    // def nmfFiles = (nmf_list instanceof java.util.Collection) ? nmf_list : [ nmf_list ]
    // // Build a bash array literal with proper quoting (paths may have spaces)
    // def filesArray = nmfFiles.collect { f -> '"' + f.toString().replace('"','\\"') + '"' }.join(' ')
    def extras     = task.ext.args ? "-P ${task.ext.args}" : ""
   
    """
    set -euo pipefail


    printf 'Staged fits:\\n' || true
    ls -1 data/per_sample_mat/nmf_fit || true

    # Render the QMD (NOTEBOOK, not an RDS)
    quarto render ${notebook} \\
      -P project_name:${params.project_name} \\
      -P work_directory:\$PWD \\
      ${extras}

    # Normalize HTML location
    mkdir -p report
    [ -f "${notebook.baseName}.html" ] && cp "${notebook.baseName}.html" report/
  
    """


    stub:
        """
        mkdir -p data figures report
        touch data/stub
        touch figures/stub
        echo "<html><body>stub</body></html>" > report/${notebook.baseName}.html
        """
}

/* =============================================================================
   End of file
   ========================================================================== */

   
