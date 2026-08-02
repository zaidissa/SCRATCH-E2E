/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    QC — ambient RNA removal, per-sample filtering, merge, doublet detection
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    Changed on merge: the standalone version computed SCDBLFINDER but then set
    `ch_final_object = ch_merge_object`, discarding the doublet-filtered result,
    so every downstream stage silently analysed the un-deduplicated object.

    That was not arbitrary: the notebook saves TWO objects that both match the
    module's `data/<project>_qc_*.RDS` output glob —

        <project>_qc_dbl_sample_object.RDS    doublets called per sample
        <project>_qc_dbl_cluster_object.RDS   doublets called per cluster

    so `emit: seurat_rds` is a two-element list and chaining it directly would
    hand two objects to a single downstream task. The choice is now explicit via
    `--qc_doublet_filter sample|cluster` (default: sample) and exactly one
    object propagates.
----------------------------------------------------------------------------------------
*/

include { CELLBENDER       } from '../../modules/local/cellbender/main.nf'
include { SEURAT_QUALITY   } from '../../modules/local/seurat/quality/main.nf'
include { HELPER_SUMMARIZE } from '../../modules/local/helper/summarize/main.nf'
include { SEURAT_MERGE     } from '../../modules/local/seurat/merge/main.nf'
include { SCDBLFINDER      } from '../../modules/local/scdblfinder/main.nf'

workflow QC {

    take:
        ch_gex_matrices  // channel: cellranger output files
        ch_exp_table     // channel: per-sample metadata table
        ch_page_config   // channel: quarto assets

    main:

        ch_notebook_quality     = Channel.fromPath(params.notebook_quality,     checkIfExists: true)
        ch_notebook_summarize   = Channel.fromPath(params.notebook_summarize,   checkIfExists: true)
        ch_notebook_merge       = Channel.fromPath(params.notebook_merge,       checkIfExists: true)
        ch_notebook_scdblfinder = Channel.fromPath(params.notebook_scdblfinder, checkIfExists: true)

        // Group cellranger outputs by sample, then pin file order so the process
        // receives (sample, metrics_summary.csv, filtered_..._matrix.h5).
        //
        // Defensive by necessity: a real cellranger directory holds GEX and VDJ
        // runs side by side, so a natural glob such as `*/outs/*` also matches
        // VDJ output that has neither of these two files. Previously those
        // samples produced a (sample, null, null) tuple and the run died with
        // "Path value cannot be null", naming nothing. Incomplete samples are
        // now dropped with a message that says which and why.
        ch_grouped = ch_gex_matrices
            .map { f -> [ f.parent.parent.name, f ] }
            .groupTuple()
            .map { sample, files ->
                [ sample,
                  files.find { it.toString().endsWith('metrics_summary.csv') },
                  files.find { it.toString().endsWith('filtered_feature_bc_matrix.h5') } ]
            }

        ch_grouped
            .filter { sample, csv, h5 -> csv == null || h5 == null }
            .subscribe { sample, csv, h5 ->
                def missing = []
                if (csv == null) missing << 'metrics_summary.csv'
                if (h5  == null) missing << 'filtered_feature_bc_matrix.h5'
                log.warn "QC: skipping '${sample}' — no ${missing.join(' or ')} " +
                         "(a VDJ or incomplete cellranger run will look like this)"
            }

        ch_cell_matrices = ch_grouped.filter { sample, csv, h5 -> csv != null && h5 != null }

        ch_cell_matrices
            .ifEmpty {
                error "No sample under --input_gex_matrices_path has BOTH " +
                      "metrics_summary.csv and filtered_feature_bc_matrix.h5.\n" +
                      "  Expected layout: <sample>/outs/<those files>\n" +
                      "  Typical glob   : --input_gex_matrices_path '/path/to/cellranger/*/outs/*'"
            }

        if (!params.skip_cellbender) {
            ch_cell_matrices = CELLBENDER(ch_cell_matrices)
        }

        SEURAT_QUALITY(ch_cell_matrices, ch_notebook_quality.collect(), ch_page_config)

        HELPER_SUMMARIZE(
            SEURAT_QUALITY.out.metrics.collect(),
            ch_notebook_summarize,
            ch_page_config
        )

        // Keep only samples whose QC status file is SUCCESS.
        ch_qc_approved = SEURAT_QUALITY.out.status
            .filter { sample, object, status -> status.toString().endsWith('SUCCESS.txt') }
            .map    { sample, object, status -> object }
            .collect()

        ch_qc_approved.ifEmpty { error 'No samples matched QC expectations. Loosen the qc_thr_* parameters or inspect the QC summary report.' }

        SEURAT_MERGE(ch_qc_approved, ch_notebook_merge, ch_exp_table, ch_page_config)
        ch_merged = SEURAT_MERGE.out.seurat_rds

        // Propagate the doublet-filtered object, not the pre-filter merge.
        // Select exactly one of the two saved variants by filename.
        if (!params.skip_scdblfinder) {
            SCDBLFINDER(ch_merged, ch_notebook_scdblfinder, ch_page_config)

            def wanted = params.qc_doublet_filter == 'cluster' ? '_qc_dbl_cluster_object.RDS'
                                                              : '_qc_dbl_sample_object.RDS'
            ch_final = SCDBLFINDER.out.seurat_rds
                .flatten()
                .filter { it.name.endsWith(wanted) }
                .ifEmpty { error "SCDBLFINDER produced no object matching '${wanted}'. Check --qc_doublet_filter (sample|cluster)." }
        } else {
            ch_final = ch_merged
        }

    emit:
        seurat_rds  = ch_final
        merged_rds  = ch_merged
        qc_metrics  = SEURAT_QUALITY.out.metrics
}
