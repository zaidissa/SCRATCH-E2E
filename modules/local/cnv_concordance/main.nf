process CNV_CONCORDANCE {

    tag "Comparing CNV callers"
    label 'process_medium'

    // Read-only by design: it compares the callers and reports where they
    // disagree. It does NOT resolve the disagreement — STRATIFY does that, and
    // keeping the two separate is what makes the comparison auditable.

    input:
        path seurat_object
        path infercnv_meta, stageAs: 'caller_inputs/infercnv/*'
        path numbat_calls,  stageAs: 'caller_inputs/numbat/*'
        path notebook
        path page_config

    output:
        // The annotated object with BOTH callers' per-cell calls written in. This
        // is the CNV stage's contribution to the cumulative backbone object:
        // STRATIFY reads this, so nothing upstream has to be re-derived later.
        // NOT optional, deliberately. This object is the CNV stage's handoff to
        // STRATIFY — objects 2 and 3 are subsets of it. An `optional` output that
        // fails to appear yields an EMPTY channel, STRATIFY never executes, and
        // the run reports success with no compartment split at all: the exact
        // silent-failure shape this pipeline keeps producing. Declaring it
        // required turns a missing handoff into a task failure that names itself.
        path("data/${params.project_name}_cnv_object.RDS"), emit: seurat_rds
        path("data/${params.project_name}_cnv_caller_comparison.csv"), emit: comparison, optional: true
        path("data/${params.project_name}_cnv_caller_summary.csv"),    emit: summary,    optional: true
        path("data/*.csv"),                                            emit: data,       optional: true
        path("figures/**"),                                            emit: figures,    optional: true
        path("report/${notebook.baseName}.html"),                      emit: html,       optional: true
        path("_freeze/**/figure-html/*.png"),                          emit: png,        optional: true

    when:
        task.ext.when == null || task.ext.when

    script:
        def args = task.ext.args ? "-P ${task.ext.args}" : ""
        """
        quarto render ${notebook} -P seurat_object:${seurat_object} ${args}
        """

    stub:
        """
        mkdir -p data figures report _freeze/DUMMY/figure-html
        # The stub must emit the names the notebook actually writes — a stub that
        # invents filenames makes -stub prove nothing about downstream wiring.
        touch data/${params.project_name}_cnv_object.RDS
        touch data/${params.project_name}_cnv_caller_comparison.csv
        touch data/${params.project_name}_cnv_caller_summary.csv
        touch figures/cnv_caller_agreement.png
        touch _freeze/DUMMY/figure-html/FILE.png
        touch report/${notebook.baseName}.html
        """
}
