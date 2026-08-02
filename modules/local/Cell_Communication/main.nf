// modules/local/Cell_Communication/main.nf
nextflow.enable.dsl = 2

/* =============================================================================
   Process: CELLCOMM_LIANA
   ========================================================================== */
process CELLCOMM_LIANA {

  tag "CellComm: LIANA"
  label 'process_medium'
  errorStrategy 'terminate'
  maxRetries 1

            //  saveAs: { path ->
            //    def p = path.toString()
            //    if (p.startsWith('_freeze/'))         return null
            //    if (p == 'report/index.html')         return "report/liana_index.html"
            //    if (p.endsWith('.html'))              return "report/${file(p).name}"
            //    if (p.startsWith('figures/'))         return "figure/${p - 'figures/'}"
            //    if (p.startsWith('data/'))            return "data/${p - 'data/'}"
            //    if (p.endsWith('_liana_done.rds'))    return "data/${file(p).name}"
            //    return "misc/${file(p).name}"
            //  }

  input:
    tuple path(seurat_object), path(notebook)

  output:
    path "report/${notebook.baseName}.html"         , emit: report,  optional: true
    path "figures/**"                               , emit: figures, optional: true
    path "data/**"                                  , emit: data,    optional: true
    path "data/liana_results.csv"                   , emit: liana_csv
    path "*_liana_done.rds", emit: done_rds

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
      -P assay:${params.cc_assay} \\
      -P slot:${params.cc_slot} \\
      -P celltype_col:${params.cc_celltype_col} \\
      -P min_cells_per_group:${params.cc_min_cells_per_group} \\
      -P liana_resource:${params.cc_liana_resource} \\
      -P liana_top_n_plot:${params.cc_liana_top_n_plot} \\
      ${extras}

    mkdir -p report
    [ -f "${notebook.baseName}.html" ] && cp "${notebook.baseName}.html" report/
    """


    stub:
        """
        mkdir -p data figures report
        echo "stub" > data/liana_results.csv
        touch data/stub
        touch figures/stub
        echo "<html><body>stub</body></html>" > report/${notebook.baseName}.html
        touch stub_liana_done.rds
        """
}


/* =============================================================================
   Process: CELLCOMM_CELLCHAT
   ========================================================================== */
process CELLCOMM_CELLCHAT {

  tag "CellComm: CellChat"
  label 'process_medium'
  errorStrategy 'terminate'
  maxRetries 1

            //  saveAs: { path ->
            //    def p = path.toString()
            //    if (p.startsWith('_freeze/'))            return null
            //    if (p == 'report/index.html')            return "report/cellchat_index.html"
            //    if (p.endsWith('.html'))                 return "report/${file(p).name}"
            //    if (p.startsWith('figures/'))            return "figure/${p - 'figures/'}"
            //    if (p.startsWith('data/'))               return "data/${p - 'data/'}"
            //    if (p.endsWith('_cellchat_done.rds'))    return "data/${file(p).name}"
            //    return "misc/${file(p).name}"
            //  }

  input:
    tuple path(seurat_object), path(notebook)

  output:
    path "report/${notebook.baseName}.html"         , emit: report,  optional: true
    path "figures/**"                               , emit: figures, optional: true
    path "data/**"                                  , emit: data,    optional: true
    path "data/cellchat_object.rds"                 , emit: cellchat_rds
    path "*_cellchat_done.rds", emit: done_rds

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
      -P organism:${params.organism_lower} \\
      -P assay:${params.cc_assay} \\
      -P slot:${params.cc_slot} \\
      -P celltype_col:${params.cc_celltype_col} \\
      -P min_cells_per_group:${params.cc_min_cells_per_group} \\
      ${extras}

    mkdir -p report
    [ -f "${notebook.baseName}.html" ] && cp "${notebook.baseName}.html" report/
    """


    stub:
        """
        mkdir -p data figures report
        touch data/cellchat_object.rds
        touch data/stub
        touch figures/stub
        echo "<html><body>stub</body></html>" > report/${notebook.baseName}.html
        touch stub_cellchat_done.rds
        """
}


/* =============================================================================
   Process: CELLCOMM_NICHENET
   ========================================================================== */
process CELLCOMM_NICHENET {

  tag "CellComm: NicheNet"
  label 'process_medium'
  errorStrategy 'terminate'
  maxRetries 1



  // input:
  //   tuple path(seurat_object), path(notebook), path(liana_csv), path(cellchat_rds)
  
  input:
    tuple path(seurat_object), path(notebook), path(liana_csv), path(cellchat_rds)
    path nichenet_assets_dir


  output:
    path "report/${notebook.baseName}.html"         , emit: report,  optional: true
    path "figures/**"                               , emit: figures, optional: true
    path "data/**"                                  , emit: data,    optional: true
    // path "data/nichenet_summary.csv"                , emit: summary_csv, optional: true
    // path "*_nichenet_done.rds"                      , emit: done_rds

  when:
    task.ext.when == null || task.ext.when

  script:
    def extras = task.ext.args ? "-P ${task.ext.args}" : ""
    """
    set -euo pipefail

    mkdir -p ./quarto_tmp
    export TMPDIR=\$PWD/quarto_tmp
    export QUARTO_TEMP=\$PWD/quarto_tmp
    
    # 2. Wipe any existing metadata
    rm -rf .quarto _freeze

    export OMP_NUM_THREADS=\${OMP_NUM_THREADS:-${task.cpus}}
    export OPENBLAS_NUM_THREADS=\${OPENBLAS_NUM_THREADS:-${task.cpus}}
    export MKL_NUM_THREADS=\${MKL_NUM_THREADS:-${task.cpus}}
    export R_PARALLEL_NUM_THREADS=\${R_PARALLEL_NUM_THREADS:-${task.cpus}}

    quarto render ${notebook} \\
      --to html \\
      --self-contained \\
      -P seurat_object:${seurat_object} \\
      -P project_name:${params.project_name} \\
      -P work_directory:\$PWD \\
      -P nichenet_assets_dir:${nichenet_assets_dir} \\
      -P liana_csv:${liana_csv} \\
      -P cellchat_rds:${cellchat_rds} \\
      -P assay:${params.cc_assay} \\
      -P slot:${params.cc_slot} \\
      -P celltype_col:${params.cc_celltype_col} \\
      -P condition_col:${params.cc_condition_col} \\
      -P min_cells_per_group:${params.cc_min_cells_per_group} \\
      -P receiver_mode:${params.cc_receiver_mode} \\
      -P sender_mode:${params.cc_sender_mode} \\
      -P auto_k_receivers:${params.cc_auto_k_receivers} \\
      -P auto_m_senders:${params.cc_auto_m_senders} \\
      -P manual_receivers:"${params.cc_manual_receivers}" \\
      -P manual_senders:"${params.cc_manual_senders}" \\
      -P condition_ref:${params.cc_condition_ref} \\
      -P condition_test:${params.cc_condition_test} \\
      -P logfc_min:${params.cc_logfc_min} \\
      -P padj_max:${params.cc_padj_max} \\
      -P max_geneset:${params.cc_max_geneset} \\
      -P top_ligands:${params.cc_top_ligands} \\
      -P heatmap_ligands:${params.cc_heatmap_ligands} \\
      -P heatmap_targets_per_ligand:${params.cc_heatmap_targets_per_ligand} \\
      ${extras}

    mkdir -p report
    [ -f "${notebook.baseName}.html" ] && cp "${notebook.baseName}.html" report/
    """

    stub:
        """
        mkdir -p data figures report
        echo "stub" > data/nichenet_summary.csv
        touch data/stub
        touch figures/stub
        echo "<html><body>stub</body></html>" > report/${notebook.baseName}.html
        touch stub_nichenet_done.rds
        """
}


// /* =============================================================================
//    Process: CELLCOMM_REPORT  (final index.html)
//    - Stages upstream figures/data into figures/ and data/
//    ========================================================================== */
// process CELLCOMM_REPORT {

//   tag "CellComm: REPORT"
//   label 'process_small'
//   errorStrategy 'terminate'
//   maxRetries 1
//   container params.container_cellcomm

//   publishDir "${params.outdir}/${params.project_name}/cellcommunication",
//              mode: 'copy',
//              overwrite: true
//             //  saveAs: { path ->
//             //    def p = path.toString()
//             //    if (p.startsWith('_freeze/'))   return null
//             //    if (p == 'report/index.html')   return "report/index.html"
//             //    if (p.endsWith('.html'))        return "report/${file(p).name}"
//             //    if (p.startsWith('figures/'))   return "figure/${p - 'figures/'}"
//             //    if (p.startsWith('data/'))      return "data/${p - 'data/'}"
//             //    return "misc/${file(p).name}"
//             //  }

//   input:
//     // tuple: notebook + lists of files
//     tuple path(notebook),
//           path(li_figs),
//           path(li_data),
//           path(cc_figs),
//           path(cc_data),
//           path(nn_figs),
//           path(nn_data)

//   output:
//     path "report/index.html"   , emit: report_index
//     path "figures/**"          , emit: figures, optional: true
//     path "data/**"             , emit: data,    optional: true

//   script:
//     """
//     set -euo pipefail

//     # Stage (Nextflow already puts inputs here; we normalize into figures/ + data/)
//     mkdir -p figures data report

//     # Copy all figures/data into unified folders for the report
//     # (safe even if some are empty)
//     rsync -a ${li_figs} figures/ || true
//     rsync -a ${cc_figs} figures/ || true
//     rsync -a ${nn_figs} figures/ || true

//     rsync -a ${li_data} data/ || true
//     rsync -a ${cc_data} data/ || true
//     rsync -a ${nn_data} data/ || true

//     quarto render ${notebook} \\
//       -P project_name:${params.project_name} \\
//       -P work_directory:\$PWD

//     # Normalize to report/index.html
//     [ -f "${notebook.baseName}.html" ] && cp "${notebook.baseName}.html" report/index.html
//     """
// }
