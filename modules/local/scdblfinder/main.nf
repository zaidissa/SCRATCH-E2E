process SCDBLFINDER {

    tag "Removing doublets"
    label 'process_high'

    // NO `executor 'local'` HERE. It was inherited from the standalone pipeline
    // and pinned doublet detection to whatever machine launched the run. On
    // Cirro that is the 2-CPU head node, so the 12 cpus asked for in
    // conf/base.config aborted the whole run with "Process requirement exceeds
    // available CPUs -- req: 12; avail: 2" (2026-10-16). A process-level
    // executor directive overrides `process.executor = 'awsbatch'` from the
    // cirro profile, so nothing in conf/ could correct it. The notebook opens
    // the BPCells store by its STAGED RELATIVE path (data/bpcells_counts), not
    // by the absolute path recorded in the RDS, so this process has never
    // needed to share a filesystem with SEURAT_MERGE.


    input:
        path(seurat_object)
        // BPCells store from SEURAT_MERGE, staged at the path the object's lazy
        // layers reference. Detection runs per sample over this store.
        path(bpcells_store, stageAs: 'data/bpcells_counts')
        path(notebook_scdblfinder)
        path(page_config)

    output:

        // Figures were written by the notebook but never declared as an

        // output, so publishDir had nothing to copy and figures/ stayed empty.

        path("figures/**"), emit: figure_files, optional: true
        path("data/${params.project_name}_qc_*.RDS"), emit: seurat_rds
        // Rewritten to contain only singlets; downstream stages open THIS store.
        path("data/bpcells_counts"), emit: bpcells_store
        path("report/${notebook_scdblfinder.baseName}.html")
                // Optional for the same reason as celltypist: a notebook that inlines
        // its figures must not fail a stage whose analysis succeeded.
        path("_freeze/**/figure-html/*.png"), emit: figures, optional: true

    when:
        task.ext.when == null || task.ext.when
        
    script:
        def param_file = task.ext.args ? "-P seurat_object:${seurat_object} -P ${task.ext.args}" : ""
        """
        quarto render ${notebook_scdblfinder} ${param_file}
        """
    stub:
        """
        mkdir -p report data/bpcells_counts figures
        touch data/bpcells_counts/matrix.mtx
        mkdir -p _freeze/DUMMY/figure-html
        
        touch _freeze/DUMMY/figure-html/FILE.png

        # Must match the notebook's object_dump chunk exactly — a stub that emits
        # names the notebook does not write makes -stub prove nothing about the
        # channel wiring downstream of here.
        touch data/${params.project_name}_qc_dbl_singlet_object.RDS

        touch report/${notebook_scdblfinder.baseName}.html

        """

}