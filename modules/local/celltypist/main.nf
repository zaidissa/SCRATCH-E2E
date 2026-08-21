process CELLTYPIST_ANNOTATION {

    tag "Running Celltypist annotation"
    label 'process_medium'

    

    input:
        path(notebook)
        path(anndata_object)
        path(config)

    output:

        // Figures were written by the notebook but never declared as an

        // output, so publishDir had nothing to copy and figures/ stayed empty.

        path("figures/**"), emit: figure_files, optional: true
        path("_freeze/${notebook.baseName}")                                  , emit: cache
        path("data/${params.project_name}_celltypist_annotation_object.h5ad") , emit: ann_object
        path("data/Immune_All")                                               , emit: csv_file
        path("report/${notebook.baseName}.html")                              , emit: html
        // path ("figures/**")                                              , emit: figures
                // Jupyter-engine notebooks with embed-resources inline figures as
        // base64 rather than writing _freeze/**/figure-html/*.png, so this
        // was a REQUIRED output that quarto never produces here. The stub
        // branch touches a DUMMY png, which is why -stub never caught it.
        path("_freeze/**/figure-html/*.png")                                  , emit: figures, optional: true
    when:
        task.ext.when == null || task.ext.when

    script:
        def param_file = task.ext.args ? "-P anndata_object:${anndata_object} -P ${task.ext.args}" : ""
        """
        # This is the one notebook Quarto renders through its JUPYTER engine, and
        # nbclient gives a new kernel 60 seconds to come up. The docker profile
        # forces linux/amd64, so on Apple silicon the kernel starts under QEMU; with
        # inferCNV or Azimuth running beside it, ipykernel routinely needs longer and
        # the render dies with "Kernel didn't respond in 60 seconds" before executing
        # a single cell. Retrying does not help — it is emulation cost, not a race.
        #
        # Quarto builds its NotebookClient with no config object
        # (share/jupyter/notebook.py: `NotebookClient(nb, resources = resources)`),
        # so neither a jupyter config file nor an env var can reach startup_timeout,
        # and traitlets ignores a reassigned `.default_value`. Passing it through the
        # constructor is what works, so sitecustomize — imported automatically by
        # `site` for any interpreter whose sys.path includes it — supplies the
        # default. Confined to this task's work directory; nothing in the image changes.
        cat > sitecustomize.py <<'PYEOF'
try:
    from nbclient import NotebookClient
    _orig_init = NotebookClient.__init__
    def _init(self, *args, **kwargs):
        kwargs.setdefault('startup_timeout', ${params.celltypist_kernel_timeout})
        _orig_init(self, *args, **kwargs)
    NotebookClient.__init__ = _init
except Exception:
    pass
PYEOF

        # Give the kernel a writable home. The container runs as `-u <uid>:<gid>`
        # with no matching passwd entry, so HOME resolves to '/' and IPython logs
        # "IPython parent '/' is not a writable location, using a temp directory"
        # — after which jupyter has nowhere stable to put the kernel's connection
        # file and the client waits for a kernel that never announces itself. The
        # symptom is indistinguishable from a slow start, which is what sent the
        # timeout above chasing the wrong cause: raising it to 600s only changed
        # the number in the error.
        export HOME="\$PWD"
        export IPYTHONDIR="\$PWD/.ipython"
        export JUPYTER_RUNTIME_DIR="\$PWD/.jupyter/runtime"
        export JUPYTER_DATA_DIR="\$PWD/.jupyter/data"
        mkdir -p "\$IPYTHONDIR" "\$JUPYTER_RUNTIME_DIR" "\$JUPYTER_DATA_DIR"

        PYTHONPATH="\$PWD\${PYTHONPATH:+:\$PYTHONPATH}" \\
        quarto render ${notebook} ${param_file}
        """
    stub:
        def param_file = task.ext.args ? "-P anndata_object:${anndata_object} -P ${task.ext.args}" : ""
        """
        mkdir -p data _freeze/${notebook.baseName}
        mkdir -p _freeze/DUMMY/figure-html

        touch _freeze/DUMMY/figure-html/FILE.png


        touch data/${params.project_name}_celltypist_annotation_object.h5ad
        mkdir -p _freeze/${notebook.baseName} data/Immune_All
        touch _freeze/${notebook.baseName}/${notebook.baseName}.html

        mkdir -p report
        touch report/${notebook.baseName}.html

        echo ${param_file} > _freeze/${notebook.baseName}/params.yml
        """

}
